using ModelingToolkit
using ModelingToolkit: @component, t_nounits as t
using ModelingToolkitStandardLibrary.Blocks
using IfElse
using PowerDynamics
using NetworkDynamics

# Include base system definitions
include("nordic_test_system.jl")

"""
LTCTransformer: ModelingToolkit component for a Load Tap Changer transformer.
Models the continuous approximation of the secondary voltage deadband regulation.
Regulates the 'src' (low-voltage load) terminal by adjusting the 'dst' (high-voltage) tap ratio m.
"""
@component function LTCTransformer(; name, R=0.0, X=0.1, Vref=1.0, db=0.015, T=30.0, m0=1.0, m_min=1/1.20, m_max=1/0.88, defaults...)
    @named src = Terminal()
    @named dst = Terminal()
    
    @parameters begin
        R_p = R
        X_p = X
        Vref_p = Vref
        db_p = db
        T_p = T
        m_min_p = m_min
        m_max_p = m_max
    end
    
    vars = @variables begin
        m(t) = m0
        dm(t) = 0.0
        Vsrc(t)
    end
    
    # Calculate admittance Y = 1/(R + jX) = G_y + jB_y
    den = R_p^2 + X_p^2
    G_y = R_p / den
    B_y = -X_p / den
    
    # Voltage difference across the transformer branch (taking m ratio on dst terminal)
    dV_r = src.u_r - m * dst.u_r
    dV_i = src.u_i - m * dst.u_i
    
    # Admittance current calculations
    imid_r = dV_r * G_y - dV_i * B_y
    imid_i = dV_r * B_y + dV_i * G_y
    
    eqs = [
        # Regulating the 'src' (low-voltage load bus) terminal voltage magnitude
        Vsrc ~ sqrt(src.u_r^2 + src.u_i^2 + 1e-9)
        
        # Deadband integral control on src voltage, driving m (which is 1/n)
        # Delay LTC action until t >= 32.0s to match the PDF timeline (30s delay after fault at t = 2.0s)
        dm ~ ifelse(t < 32.0, 0.0,
                    ifelse(Vsrc < Vref_p - db_p, (Vref_p - Vsrc) / T_p,
                           ifelse(Vsrc > Vref_p + db_p, (Vref_p - Vsrc) / T_p, 0.0)))
                    
        # Tap limiters (stops tapping when reaching min/max limits)
        Differential(t)(m) ~ ifelse(m >= m_max_p, min(0.0, dm),
                                    ifelse(m <= m_min_p, max(0.0, dm), dm))
                                    
        # Terminal currents
        src.i_r ~ -imid_r
        src.i_i ~ -imid_i
        dst.i_r ~ imid_r * m
        dst.i_i ~ imid_i * m
    ]
    
    sys = System(eqs, t; name, systems=[src, dst])
    set_mtk_defaults!(sys, defaults)
    return sys
end

"""
NordicAVRWithOEL: Custom ModelingToolkit AVR type I controller coupled with a takeover
Overexcitation Limiter (OEL) model as defined in the IEEE PES-TR19 Nordic test system.
"""
@component function NordicAVRWithOEL(; name,
    Vref_init=1.0, ilim_fd=1.8991, f=0.0, r=1.0, L1=-11.0, G=70.0,
    Ta=10.0, Tb=20.0, L2=4.0, Texc=0.1, defaults...)
    
    @named v_mag = RealInput()
    @named i_fd = RealInput()
    @named vf = RealOutput()
    
    @parameters begin
        ilim_fd_p = ilim_fd
        f_p = f
        r_p = r
        L1_p = L1
        G_p = G
        Ta_p = Ta
        Tb_p = Tb
        L2_p = L2
        Texc_p = Texc
        Vref, [guess=Vref_init, description="Reference voltage"]
    end
    
    vars = @variables begin
        x_oel(t) = L1
        x_avr(t) = 0.0
        vf_val(t) = Vref_init
    end
    
    # OEL timer/integrator input:
    # If Ifd < Ilim: integrator integrates -1.0 (resets towards lower limit L1)
    # If Ifd >= Ilim: integrator integrates f + r * (Ifd - Ilim), triggering overload accumulation
    oel_in = ifelse(i_fd.u < ilim_fd_p, -1.0, f_p + r_p * (i_fd.u - ilim_fd_p))
    
    # OEL is active when integrator value x_oel becomes positive (> 0.0)
    oel_active = x_oel > 0.0
    
    err_v = Vref - v_mag.u
    err_i = ilim_fd_p - i_fd.u
    
    # Takeover Min-Gate:
    # If OEL is active, the amplifier receives the minimum of the voltage error and current error
    err = ifelse(oel_active,
                 ifelse(err_v < err_i, err_v, err_i),
                 err_v)
                 
    # Lead-lag block (transient gain reduction)
    Vr = G_p * (x_avr + (Ta_p / Tb_p) * (err - x_avr))
    
    eqs = [
        # OEL integrator state with limits [L1, 1.0]
        Differential(t)(x_oel) ~ ifelse(x_oel >= 1.0, min(0.0, oel_in),
                                        ifelse(x_oel <= L1_p, max(0.0, oel_in), oel_in)),
                                        
        # Lead-lag filter state
        Differential(t)(x_avr) ~ (err - x_avr) / Tb_p,
        
        # Exciter state with limits [0.0, L2]
        Differential(t)(vf_val) ~ (ifelse(vf_val >= L2_p, min(0.0, Vr - vf_val),
                                          ifelse(vf_val <= 0.0, max(0.0, Vr - vf_val), Vr - vf_val))) / Texc_p,
                                          
        # Connect output port
        vf.u ~ vf_val
    ]
    
    sys = System(eqs, t; name, systems=[v_mag, i_fd, vf])
    set_mtk_defaults!(sys, defaults)
    return sys
end

"""
Build the Nordic network with dynamic step-down LTCs and dynamic generators (GENROU/GENSAL + AVR + OEL).
"""
function build_nordic_network_dynamic(; G_fault=20.0, t_fault=2.0, t_clear=2.1, trip_line=true, T_ltc=30.0)
    buses = Any[]
    
    for label in BUS_LABELS
        if label != "g20" && haskey(GEN_OP_A, label)
            Pmw, _, Vset = GEN_OP_A[label]
            components = Any[]
            
            # Nominal MVA base (Snom) from GSU transformer definitions
            Snom = 100.0
            for (tname, gname, busname, _, _, trans_Snom) in STEP_UP_TRANSFORMERS
                if gname == label
                    Snom = trans_Snom
                    break
                end
            end
            
            # Determine machine parameters based on Table 2.12
            is_round_rotor = label in ["g6", "g7", "g14", "g15", "g16", "g17", "g18"]
            pmech_val = Pmw / Snom
            
            gen_model = if is_round_rotor
                H = 6.0
                Library.PSSE_GENROU(
                    name = :gen,
                    pmech_input = false,
                    pmech_set = pmech_val,
                    M_b = Snom,
                    S_b = 100.0,
                    fn = 50.0,
                    H = H,
                    D = 15.0,
                    Xd = 2.20,
                    Xq = 2.00,
                    Xpd = 0.30,
                    Xpq = 0.40,
                    Xppd = 0.20,
                    Xppq = 0.20,
                    Xl = 0.15,
                    Tpd0 = 7.0,
                    Tpq0 = 1.5,
                    Tppd0 = 0.05,
                    Tppq0 = 0.05,
                    S10 = 0.1,
                    S12 = 0.3,
                    R_a = 0.0
                )
            else
                # Salient pole machines (g13 vs others)
                is_g13 = label == "g13"
                H = is_g13 ? 2.0 : 3.0
                Xd_val = is_g13 ? 1.55 : 1.10
                Xq_val = is_g13 ? 1.00 : 0.70
                Xpd_val = is_g13 ? 0.30 : 0.25
                Tpd0_val = is_g13 ? 7.0 : 5.0
                
                Library.PSSE_GENSAL(
                    name = :gen,
                    pmech_input = false,
                    pmech_set = pmech_val,
                    M_b = Snom,
                    S_b = 100.0,
                    fn = 50.0,
                    H = H,
                    D = 15.0,
                    Xd = Xd_val,
                    Xq = Xq_val,
                    Xpd = Xpd_val,
                    Xppd = 0.20,
                    Xppq = 0.20,
                    Xl = 0.15,
                    Tpd0 = Tpd0_val,
                    Tppd0 = 0.05,
                    Tppq0 = 0.10,
                    S10 = 0.1,
                    S12 = 0.3,
                    R_a = 0.0
                )
            end
            
            # Determine AVR & OEL parameters based on Table 2.13
            is_fixed_time = label in ["g6", "g7", "g11", "g12"]
            ilim_fd_val = is_round_rotor ? 3.0618 : (label == "g13" ? 2.9579 : 1.8991)
            
            f_val = is_fixed_time ? 1.0 : 0.0
            r_val = is_fixed_time ? 0.0 : 1.0
            
            L1_val = if is_fixed_time
                -20.0
            elseif label == "g13"
                -17.0
            elseif is_round_rotor
                -18.0
            else
                -11.0
            end
            
            G_val = is_round_rotor ? 120.0 : (label == "g13" ? 50.0 : 70.0)
            Ta_val = label in ["g6", "g7", "g14", "g15", "g16", "g17", "g18"] ? 5.0 : (label == "g13" ? 4.0 : 10.0)
            Tb_val = label in ["g6", "g7", "g14", "g15", "g16", "g17", "g18"] ? 12.5 : 20.0
            L2_val = is_round_rotor ? 5.0 : 4.0
            
            avr_model = NordicAVRWithOEL(
                name = :avr,
                Vref_init = Vset,
                ilim_fd = ilim_fd_val,
                f = f_val,
                r = r_val,
                L1 = L1_val,
                G = G_val,
                Ta = Ta_val,
                Tb = Tb_val,
                L2 = L2_val,
                Texc = 0.1
            )
            
            # Build CompositeInjector
            eqs = [
                connect(avr_model.vf, gen_model.EFD_in),
                connect(gen_model.ETERM_out, avr_model.v_mag),
                connect(gen_model.XADIFD_out, avr_model.i_fd)
            ]
            
            comp_gen = CompositeInjector([gen_model, avr_model], eqs; name=Symbol("comp_gen_", label))
            push!(components, comp_gen)
            
            # Compile bus
            bus = compile_bus(MTKBus(components...); vidx=idx(label), name=Symbol("bus_", label))
            pf_model = pfPV(P=Pmw / SBASE_MVA, V=Vset)
            if haskey(INIT_VOLTAGE, label)
                V, deg = INIT_VOLTAGE[label]
                phasor = V * cis(deg2rad(deg))
                set_voltage!(bus, phasor)
                set_voltage!(pf_model, phasor)
            end
            set_pfmodel!(bus, pf_model)
            push!(buses, bus)
        else
            push!(buses, make_bus(label; G_fault=G_fault, t_fault=t_fault, t_clear=t_clear))
        end
    end
    
    lines = Any[]
    
    # 1. Transmission Lines
    line_groups = Dict{Tuple{String, String, Float64}, Vector{Tuple{String, Float64, Float64, Float64}}}()
    for (name, from, to, Rohm, Xohm, Bhalf_microS) in LINE_DATA
        R, X, Bhalf = line_pu(from, Rohm, Xohm, Bhalf_microS)
        t_open = (trip_line && name == "4032-4044") ? t_clear : Inf
        key = (from, to, t_open)
        push!(get!(line_groups, key, Tuple{String, Float64, Float64, Float64}[]), (name, R, X, Bhalf))
    end

    for ((from, to, t_open), group) in line_groups
        y = sum(1 / (R + im * X) for (_, R, X, _) in group)
        z = 1 / y
        Bhalf = sum(item[4] for item in group)
        name = join(first.(group), "_")
        push!(lines, make_line(name, from, to, real(z), imag(z), Bhalf; t_open))
    end

    # 2. Step-up & Interlevel Transformers (remain static PiLine models)
    transformer_groups = Dict{Tuple{String, String, Float64}, Vector{Tuple{String, Float64}}}()
    for (name, from, to, Xown, n, Snom) in vcat(STEP_UP_TRANSFORMERS, INTERLEVEL_TRANSFORMERS)
        X = transformer_x_pu(Xown, Snom)
        key = (from, to, n)
        push!(get!(transformer_groups, key, Tuple{String, Float64}[]), (name, X))
    end

    for ((from, to, n), group) in transformer_groups
        Xeq = 1 / sum(1 / X for (_, X) in group)
        name = join(first.(group), "_")
        @named xf = Library.PiLine(R=0.0, X=Xeq, B_src=0.0, B_dst=0.0, r_src=1.0, r_dst=1 / n)
        push!(lines, compile_line(MTKLine(xf); src=idx(from), dst=idx(to), name=Symbol("xf_", replace(name, "-" => "_"))))
    end

    # 3. Step-down Transformers (replaced with dynamic LTCTransformer models)
    for (name, from, to, Xown, n, Snom) in STEP_DOWN_TRANSFORMERS
        Vref = 1.0
        if haskey(LOAD_OP_A, from)
            _, _, Vref = LOAD_OP_A[from]
        end
        
        X = transformer_x_pu(Xown, Snom)
        m0 = 1 / n
        
        @named ltc_xf = LTCTransformer(R=0.0, X=X, Vref=Vref, m0=m0, T=T_ltc, db=0.015)
        line = compile_line(MTKLine(ltc_xf); src=idx(from), dst=idx(to), name=Symbol("ltc_", replace(name, "-" => "_")))
        
        @named pf_ltc = Library.PiLine(R=0.0, X=X, B_src=0.0, B_dst=0.0, r_src=1.0, r_dst=m0)
        set_pfmodel!(line, compile_line(MTKLine(pf_ltc); src=idx(from), dst=idx(to), name=Symbol("pf_ltc_", replace(name, "-" => "_"))))
        
        push!(lines, line)
    end

    return Network(buses, lines; warn_order=false)
end

# 4. Custom initialization function
function initialize_nordic_dynamic(nw; verbose=false)
    pfnw = PowerDynamics.powerflow_model(nw)
    pfs0 = NWState(pfnw)
    pfs = PowerDynamics.solve_powerflow(nw; pfnw=pfnw, pfs0=pfs0, verbose=false)
    interf = NetworkDynamics.interface_values(pfs)
    pfconstraints = PowerDynamics.specialize_pfinitconstraints(nw, pfs)
    pfformulas = PowerDynamics.specialize_pfinitformulas(nw, pfs)

    # Construct custom overrides to free vf_val and x_avr from defaults list
    defaults = Dict{Any, Any}()
    guesses = Dict{Any, Any}()

    for (k, v) in interf
        defaults[k] = v
    end

    for label in BUS_LABELS
        if label != "g20" && haskey(GEN_OP_A, label)
            _, _, Vset = GEN_OP_A[label]
            vidx_val = idx(label)
            # Remove vf_val and x_avr from defaults by mapping to nothing
            defaults[NetworkDynamics.VIndex(vidx_val, Symbol("comp_gen_", label, "₊avr₊vf_val"))] = nothing
            defaults[NetworkDynamics.VIndex(vidx_val, Symbol("comp_gen_", label, "₊avr₊x_avr"))] = nothing
            
            # Map them as guess overrides instead
            guesses[NetworkDynamics.VIndex(vidx_val, Symbol("comp_gen_", label, "₊avr₊vf_val"))] = Vset
            guesses[NetworkDynamics.VIndex(vidx_val, Symbol("comp_gen_", label, "₊avr₊x_avr"))] = 0.0
        end
    end

    s0 = NetworkDynamics.initialize_componentwise(
        nw;
        default_overrides=defaults,
        guess_overrides=guesses,
        additional_initconstraint=pfconstraints,
        additional_initformula=pfformulas,
        verbose,
        subverbose=false,
        tol=1e-3,
        nwtol=1e-3,
        t=0.0,
    )
    return s0, pfs
end
