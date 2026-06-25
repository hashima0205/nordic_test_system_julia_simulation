using ModelingToolkit
using ModelingToolkit: @component, t_nounits as t
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
    # Vsrc - m * Vdst
    dV_r = src.u_r - m * dst.u_r
    dV_i = src.u_i - m * dst.u_i
    
    # Admittance current calculations
    imid_r = dV_r * G_y - dV_i * B_y
    imid_i = dV_r * B_y + dV_i * G_y
    
    eqs = [
        # Regulating the 'src' (low-voltage load bus) terminal voltage magnitude
        Vsrc ~ sqrt(src.u_r^2 + src.u_i^2 + 1e-9)
        
        # Deadband integral control on src voltage, driving m (which is 1/n)
        dm ~ ifelse(Vsrc < Vref_p - db_p, (Vref_p - Vsrc) / T_p,
                    ifelse(Vsrc > Vref_p + db_p, (Vref_p - Vsrc) / T_p, 0.0))
                    
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
Build the Nordic network, replacing step-down transformers with dynamic LTCTransformers.
"""
function build_nordic_network_with_ltc(; G_fault=20.0, t_fault=1.0, t_clear=1.1, trip_line=true, T_ltc=30.0)
    buses = [make_bus(label; G_fault, t_fault, t_clear) for label in BUS_LABELS]

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
        # Set Vref based on the nominal initial voltage of the distribution bus
        Vref = 1.0
        if haskey(LOAD_OP_A, from)
            _, _, Vref = LOAD_OP_A[from]
        end
        
        X = transformer_x_pu(Xown, Snom)
        m0 = 1 / n
        
        # Compile LTCTransformer
        # src (1st terminal) connected to 'from' (low-voltage load bus)
        # dst (2nd terminal) connected to 'to' (high-voltage transmission bus)
        # This matches the direction in original make_transformer function
        @named ltc_xf = LTCTransformer(R=0.0, X=X, Vref=Vref, m0=m0, T=T_ltc, db=0.015)
        line = compile_line(MTKLine(ltc_xf); src=idx(from), dst=idx(to), name=Symbol("ltc_", replace(name, "-" => "_")))
        
        # Static Power Flow Model for LTC (functions as a static PiLine with r_dst = m0)
        @named pf_ltc = Library.PiLine(R=0.0, X=X, B_src=0.0, B_dst=0.0, r_src=1.0, r_dst=m0)
        set_pfmodel!(line, compile_line(MTKLine(pf_ltc); src=idx(from), dst=idx(to), name=Symbol("pf_ltc_", replace(name, "-" => "_"))))
        
        push!(lines, line)
    end

    return Network(buses, lines; warn_order=false)
end

"""
Custom initialization function. Fits within a strict tolerance since polarities are now aligned.
"""
function initialize_nordic_ltc(nw; verbose=false)
    pfnw = PowerDynamics.powerflow_model(nw)
    pfs0 = NWState(pfnw)
    pfs = PowerDynamics.solve_powerflow(nw; pfnw=pfnw, pfs0=pfs0, verbose=verbose)
    interf = NetworkDynamics.interface_values(pfs)
    pfconstraints = PowerDynamics.specialize_pfinitconstraints(nw, pfs)
    pfformulas = PowerDynamics.specialize_pfinitformulas(nw, pfs)

    s0 = NetworkDynamics.initialize_componentwise(
        nw;
        default_overrides=interf,
        additional_initconstraint=pfconstraints,
        additional_initformula=pfformulas,
        verbose,
        subverbose=false,
        tol=1e-3,      # Restored strict 1e-3 tolerance
        nwtol=1e-3,
        t=0.0,
    )
    return s0, pfs
end
