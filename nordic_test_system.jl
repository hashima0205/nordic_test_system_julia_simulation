using PowerDynamics
using PowerDynamics: Library
using ModelingToolkit
using ModelingToolkit: @component, t_nounits as t
using OrdinaryDiffEq
using NetworkDynamics

const SBASE_MVA = 100.0
const FBASE_HZ = 50.0
const WBASE = 2π * FBASE_HZ
const GRAPH_DIR = "graph"

# IEEE PES-TR19 Nordic test system, operating point A.
# This is a simulation-oriented first model:
# - full 74-bus topology from PES-TR19 Tables 2.2--2.7,
# - voltage-dependent loads from Eq. (2.1): P ~ V^1, Q ~ V^2,
# - generator buses represented by PowerDynamics.Swing machines,
# - the PES-TR19 line 4032-4044 fault-and-trip contingency.
#
# The PES-TR19 long-term collapse is driven by LTCs and OELs.  Those are left
# as explicit extension points below; the base model is kept compact enough to
# initialize and to serve as a reproducible starting point.

@component function NordicLoad(; name, P0=0.0, Q0=0.0, V0=1.0, α=1.0, β=2.0, defaults...)
    @named terminal = Terminal()
    @parameters begin
        P0_p = P0
        Q0_p = Q0
        V0_p = V0
        α_p = α
        β_p = β
    end
    vars = @variables begin
        Vmag(t)
        P(t)
        Q(t)
    end
    eqs = [
        Vmag ~ sqrt(terminal.u_r^2 + terminal.u_i^2 + 1e-9)
        P ~ P0_p * (Vmag / V0_p)^α_p
        Q ~ Q0_p * (Vmag / V0_p)^β_p
        terminal.i_r ~ -(P * terminal.u_r + Q * terminal.u_i) / Vmag^2
        terminal.i_i ~ -(P * terminal.u_i - Q * terminal.u_r) / Vmag^2
    ]
    sys = System(eqs, t; name, systems=[terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end

@component function FaultShunt(; name, G_fault=0.0, t_fault=0.0, t_clear=0.1, sharpness=500.0, defaults...)
    @named terminal = Terminal()
    @parameters begin
        Gf = G_fault
        tf = t_fault
        tc = t_clear
        k = sharpness
    end
    vars = @variables begin
        g(t)
    end
    eqs = [
        g ~ (Gf / 2) * (tanh(k * (t - tf)) - tanh(k * (t - tc)))
        terminal.i_r ~ -g * terminal.u_r
        terminal.i_i ~ -g * terminal.u_i
    ]
    sys = System(eqs, t; name, systems=[terminal])
    set_mtk_defaults!(sys, defaults)
    return sys
end

@component function TimedPiLine(; name, R=0.0, X=0.1, B_src=0.0, B_dst=0.0,
                               G_src=0.0, G_dst=0.0, r_src=1.0, r_dst=1.0,
                               t_open=Inf, initially_active=1.0, sharpness=500.0,
                               defaults...)
    @named src = Terminal()
    @named dst = Terminal()
    @parameters begin
        R_p = R
        X_p = X
        Bs_p = B_src
        Bd_p = B_dst
        Gs_p = G_src
        Gd_p = G_dst
        rs_p = r_src
        rd_p = r_dst
        topen_p = t_open
        active0_p = initially_active
        k_p = sharpness
    end
    vars = @variables begin
        active(t)
    end
    begin
        Z = R_p + im * X_p
        Ysrc = Gs_p + im * Bs_p
        Ydst = Gd_p + im * Bd_p
        Vsrc = src.u_r + im * src.u_i
        Vdst = dst.u_r + im * dst.u_i
        V1 = rs_p * Vsrc
        V2 = rd_p * Vdst
        i1 = Ysrc * V1
        i2 = Ydst * V2
        imid = (V1 - V2) / Z
        isrc = (-imid - i1) * rs_p
        idst = (imid - i2) * rd_p
    end
    eqs = [
        active ~ active0_p * (1 - (1 / 2) * (1 + tanh(k_p * (t - topen_p))))
        src.i_r ~ active * real(isrc)
        src.i_i ~ active * imag(isrc)
        dst.i_r ~ active * real(idst)
        dst.i_i ~ active * imag(idst)
    ]
    sys = System(eqs, t; name, systems=[src, dst])
    set_mtk_defaults!(sys, defaults)
    return sys
end

const TRANSMISSION_BUSES = [
    "1011","1012","1013","1014","1021","1022","1041","1042","1043","1044","1045",
    "2031","2032","4011","4012","4021","4022","4031","4032","4041","4042","4043",
    "4044","4045","4046","4047","4051","4061","4062","4063","4071","4072"
]

const DISTRIBUTION_BUSES = [
    "1","2","3","4","5","11","12","13","22","31","32","41","42","43","46","47",
    "51","61","62","63","71","72"
]

const GENERATOR_BUSES = ["g$i" for i in 1:20]

const BUS_LABELS = vcat(DISTRIBUTION_BUSES, TRANSMISSION_BUSES, GENERATOR_BUSES)
const BUS_INDEX = Dict(label => i for (i, label) in enumerate(BUS_LABELS))

idx(label) = BUS_INDEX[string(label)]

const BASE_KV = Dict{String, Float64}(
    [b => 20.0 for b in DISTRIBUTION_BUSES]...,
    [b => 130.0 for b in ["1011","1012","1013","1014","1021","1022","1041","1042","1043","1044","1045"]]...,
    [b => 220.0 for b in ["2031","2032"]]...,
    [b => 400.0 for b in ["4011","4012","4021","4022","4031","4032","4041","4042","4043","4044","4045","4046","4047","4051","4061","4062","4063","4071","4072"]]...,
    [b => 15.0 for b in GENERATOR_BUSES]...
)

const GEN_OP_A = Dict(
    "g1"=>(600.0, 58.3, 1.0684), "g2"=>(300.0, 17.2, 1.0565), "g3"=>(550.0, 20.9, 1.0595),
    "g4"=>(400.0, 30.4, 1.0339), "g5"=>(200.0, 60.1, 1.0294), "g6"=>(360.0, 138.6, 1.0084),
    "g7"=>(180.0, 60.4, 1.0141), "g8"=>(750.0, 232.6, 1.0498), "g9"=>(668.5, 201.3, 0.9988),
    "g10"=>(600.0, 255.7, 1.0157), "g11"=>(250.0, 60.7, 1.0211), "g12"=>(310.0, 98.3, 1.0200),
    "g13"=>(0.0, 50.1, 1.0170), "g14"=>(630.0, 295.9, 1.0454), "g15"=>(1080.0, 377.9, 1.0455),
    "g16"=>(600.0, 222.6, 1.0531), "g17"=>(530.0, 48.7, 1.0092), "g18"=>(1060.0, 293.4, 1.0307),
    "g19"=>(300.0, 121.2, 1.0300), "g20"=>(2137.4, 377.4, 1.0185)
)

const LOAD_OP_A = Dict(
    "1"=>(600.0,148.2,0.9988), "2"=>(330.0,71.0,1.0012), "3"=>(260.0,83.8,0.9974),
    "4"=>(840.0,252.0,0.9996), "5"=>(720.0,190.4,0.9961), "11"=>(200.0,68.8,1.0026),
    "12"=>(300.0,83.8,0.9975), "13"=>(100.0,34.4,0.9957), "22"=>(280.0,79.9,0.9952),
    "31"=>(100.0,24.7,1.0042), "32"=>(200.0,39.6,0.9978), "41"=>(540.0,131.4,0.9967),
    "42"=>(400.0,127.4,0.9952), "43"=>(900.0,254.6,1.0013), "46"=>(700.0,211.8,0.9990),
    "47"=>(100.0,44.0,0.9950), "51"=>(800.0,258.2,0.9978), "61"=>(500.0,122.5,0.9949),
    "62"=>(300.0,83.8,1.0002), "63"=>(590.0,264.6,0.9992), "71"=>(300.0,83.8,1.0028),
    "72"=>(2000.0,396.1,0.9974)
)

const SHUNTS_MVAR = Dict(
    "1022"=>50.0, "1041"=>250.0, "1043"=>200.0, "1044"=>200.0, "1045"=>200.0,
    "4012"=>-100.0, "4041"=>200.0, "4043"=>200.0, "4046"=>100.0, "4051"=>100.0,
    "4071"=>-400.0
)

const INIT_VOLTAGE = Dict(
    "g1"=>(1.0684,2.59), "g2"=>(1.0565,5.12), "g3"=>(1.0595,10.27), "g4"=>(1.0339,8.03),
    "g5"=>(1.0294,-12.36), "g6"=>(1.0084,-59.42), "g7"=>(1.0141,-68.95), "g8"=>(1.0498,-16.81),
    "g9"=>(0.9988,-1.63), "g10"=>(1.0157,0.99), "g11"=>(1.0211,-29.04), "g12"=>(1.0200,-31.88),
    "g13"=>(1.0170,-54.30), "g14"=>(1.0454,-49.90), "g15"=>(1.0455,-52.19), "g16"=>(1.0531,-64.10),
    "g17"=>(1.0092,-46.85), "g18"=>(1.0307,-43.32), "g19"=>(1.0300,0.03), "g20"=>(1.0185,0.00),
    "1011"=>(1.0618,-6.65), "1012"=>(1.0634,-3.10), "1013"=>(1.0548,1.26), "1014"=>(1.0611,4.26),
    "1021"=>(1.0311,2.64), "1022"=>(1.0512,-19.05), "1041"=>(1.0124,-81.87), "1042"=>(1.0145,-67.38),
    "1043"=>(1.0274,-76.77), "1044"=>(1.0066,-67.71), "1045"=>(1.0111,-71.66), "2031"=>(1.0279,-36.66),
    "2032"=>(1.0695,-23.92), "4011"=>(1.0224,-7.55), "4012"=>(1.0235,-5.54), "4021"=>(1.0488,-36.08),
    "4022"=>(0.9947,-20.86), "4031"=>(1.0367,-39.46), "4032"=>(1.0487,-44.54), "4041"=>(1.0506,-54.30),
    "4042"=>(1.0428,-57.37), "4043"=>(1.0370,-63.51), "4044"=>(1.0395,-64.23), "4045"=>(1.0533,-68.88),
    "4046"=>(1.0357,-64.11), "4047"=>(1.0590,-59.55), "4051"=>(1.0659,-71.01), "4061"=>(1.0387,-57.93),
    "4062"=>(1.0560,-54.36), "4063"=>(1.0536,-50.68), "4071"=>(1.0484,-4.99), "4072"=>(1.0590,-3.98),
    "1"=>(0.9988,-84.71), "2"=>(1.0012,-70.49), "3"=>(0.9974,-79.97), "4"=>(0.9996,-70.67),
    "5"=>(0.9961,-74.59), "11"=>(1.0026,-9.45), "12"=>(0.9975,-5.93), "13"=>(0.9957,-1.58),
    "22"=>(0.9952,-21.89), "31"=>(1.0042,-39.47), "32"=>(0.9978,-26.77), "41"=>(0.9967,-57.14),
    "42"=>(0.9952,-60.22), "43"=>(1.0013,-66.33), "46"=>(0.9990,-66.93), "47"=>(0.9950,-62.38),
    "51"=>(0.9978,-73.84), "61"=>(0.9949,-60.78), "62"=>(1.0002,-57.18), "63"=>(0.9992,-53.49),
    "71"=>(1.0028,-7.80), "72"=>(0.9974,-6.83)
)

const LINE_DATA = [
    ("1011-1013","1011","1013",1.69,11.83,40.841), ("1011-1013b","1011","1013",1.69,11.83,40.841),
    ("1012-1014","1012","1014",2.37,15.21,53.407), ("1012-1014b","1012","1014",2.37,15.21,53.407),
    ("1013-1014","1013","1014",1.18,8.450,29.845), ("1013-1014b","1013","1014",1.18,8.450,29.845),
    ("1021-1022","1021","1022",5.07,33.80,89.535), ("1021-1022b","1021","1022",5.07,33.80,89.535),
    ("1041-1043","1041","1043",1.69,10.14,36.128), ("1041-1043b","1041","1043",1.69,10.14,36.128),
    ("1041-1045","1041","1045",2.53,20.28,73.827), ("1041-1045b","1041","1045",2.53,20.28,73.827),
    ("1042-1044","1042","1044",6.42,47.32,177.50), ("1042-1044b","1042","1044",6.42,47.32,177.50),
    ("1042-1045","1042","1045",8.45,50.70,177.50), ("1043-1044","1043","1044",1.69,13.52,47.124),
    ("1043-1044b","1043","1044",1.69,13.52,47.124), ("2031-2032","2031","2032",5.81,43.56,15.708),
    ("2031-2032b","2031","2032",5.81,43.56,15.708), ("4011-4012","4011","4012",1.60,12.80,62.832),
    ("4011-4021","4011","4021",9.60,96.00,562.34), ("4011-4022","4011","4022",6.40,64.00,375.42),
    ("4011-4071","4011","4071",8.00,72.00,438.25), ("4012-4022","4012","4022",6.40,56.00,328.30),
    ("4012-4071","4012","4071",8.00,80.00,468.10), ("4021-4032","4021","4032",6.40,64.00,375.42),
    ("4021-4042","4021","4042",16.0,96.00,937.77), ("4022-4031","4022","4031",6.40,64.00,375.42),
    ("4022-4031b","4022","4031",6.40,64.00,375.42), ("4031-4032","4031","4032",1.60,16.00,94.248),
    ("4031-4041","4031","4041",9.60,64.00,749.27), ("4031-4041b","4031","4041",9.60,64.00,749.27),
    ("4032-4042","4032","4042",16.0,64.00,625.18), ("4032-4044","4032","4044",9.60,80.00,749.27),
    ("4041-4044","4041","4044",4.80,48.00,281.17), ("4041-4061","4041","4061",9.60,72.00,406.84),
    ("4042-4043","4042","4043",3.20,24.00,155.51), ("4042-4044","4042","4044",3.20,32.00,186.93),
    ("4043-4044","4043","4044",1.60,16.00,94.248), ("4043-4046","4043","4046",1.60,16.00,94.248),
    ("4043-4047","4043","4047",3.20,32.00,186.93), ("4044-4045","4044","4045",3.20,32.00,186.93),
    ("4044-4045b","4044","4045",3.20,32.00,186.93), ("4045-4051","4045","4051",6.40,64.00,375.42),
    ("4045-4051b","4045","4051",6.40,64.00,375.42), ("4045-4062","4045","4062",17.6,128.00,749.27),
    ("4046-4047","4046","4047",1.60,24.00,155.51), ("4061-4062","4061","4062",3.20,32.00,186.93),
    ("4062-4063","4062","4063",4.80,48.00,281.17), ("4062-4063b","4062","4063",4.80,48.00,281.17),
    ("4071-4072","4071","4072",4.80,48.00,937.77), ("4071-4072b","4071","4072",4.80,48.00,937.77)
]

const STEP_UP_TRANSFORMERS = [
    ("g1","g1","1012",0.15,1.00,800.0), ("g2","g2","1013",0.15,1.00,600.0),
    ("g3","g3","1014",0.15,1.00,700.0), ("g4","g4","1021",0.15,1.00,600.0),
    ("g5","g5","1022",0.15,1.05,250.0), ("g6","g6","1042",0.15,1.05,400.0),
    ("g7","g7","1043",0.15,1.05,200.0), ("g8","g8","2032",0.15,1.05,850.0),
    ("g9","g9","4011",0.15,1.05,1000.0), ("g10","g10","4012",0.15,1.05,800.0),
    ("g11","g11","4021",0.15,1.05,300.0), ("g12","g12","4031",0.15,1.05,350.0),
    ("g13","g13","4041",0.10,1.05,300.0), ("g14","g14","4042",0.15,1.05,700.0),
    ("g15","g15","4047",0.15,1.05,1200.0), ("g16","g16","4051",0.15,1.05,700.0),
    ("g17","g17","4062",0.15,1.05,600.0), ("g18","g18","4063",0.15,1.05,1200.0),
    ("g19","g19","4071",0.15,1.05,500.0), ("g20","g20","4072",0.15,1.05,4500.0)
]

const INTERLEVEL_TRANSFORMERS = [
    ("1011-4011","1011","4011",0.10,0.95,1250.0), ("1012-4012","1012","4012",0.10,0.95,1250.0),
    ("1022-4022","1022","4022",0.10,0.93,833.3), ("2031-4031","2031","4031",0.10,1.00,833.3),
    ("1044-4044","1044","4044",0.10,1.03,1000.0), ("1044-4044b","1044","4044",0.10,1.03,1000.0),
    ("1045-4045","1045","4045",0.10,1.04,1000.0), ("1045-4045b","1045","4045",0.10,1.04,1000.0)
]

const STEP_DOWN_TRANSFORMERS = [
    ("11-1011","11","1011",0.10,1.04,400.0), ("12-1012","12","1012",0.10,1.05,600.0),
    ("13-1013","13","1013",0.10,1.04,200.0), ("22-1022","22","1022",0.10,1.04,560.0),
    ("1-1041","1","1041",0.10,1.00,1200.0), ("2-1042","2","1042",0.10,1.00,600.0),
    ("3-1043","3","1043",0.10,1.01,460.0), ("4-1044","4","1044",0.10,0.99,1600.0),
    ("5-1045","5","1045",0.10,1.00,1400.0), ("31-2031","31","2031",0.10,1.01,200.0),
    ("32-2032","32","2032",0.10,1.06,400.0), ("41-4041","41","4041",0.10,1.04,1080.0),
    ("42-4042","42","4042",0.10,1.03,800.0), ("43-4043","43","4043",0.10,1.02,1800.0),
    ("46-4046","46","4046",0.10,1.02,1400.0), ("47-4047","47","4047",0.10,1.04,200.0),
    ("51-4051","51","4051",0.10,1.05,1600.0), ("61-4061","61","4061",0.10,1.03,1000.0),
    ("62-4062","62","4062",0.10,1.04,600.0), ("63-4063","63","4063",0.10,1.03,1180.0),
    ("71-4071","71","4071",0.10,1.03,600.0), ("72-4072","72","4072",0.10,1.05,4000.0)
]

function line_pu(from_bus, R_ohm, X_ohm, Bhalf_microS)
    zbase = BASE_KV[from_bus]^2 / SBASE_MVA
    return R_ohm / zbase, X_ohm / zbase, Bhalf_microS * 1e-6 * zbase
end

function transformer_x_pu(X_own, Snom)
    X_own * SBASE_MVA / Snom
end

function make_bus(label; fault_bus="4032", G_fault=0.0, t_fault=0.0, t_clear=0.1)
    components = Any[]

    if haskey(LOAD_OP_A, label)
        Pmw, Qmvar, V0 = LOAD_OP_A[label]
        @named load = NordicLoad(P0=Pmw / SBASE_MVA, Q0=Qmvar / SBASE_MVA, V0=V0)
        push!(components, load)
    elseif haskey(GEN_OP_A, label)
        Pmw, _, Vset = GEN_OP_A[label]
        if label == "g20"
            @named gen = Library.VδConstraint(V=Vset, δ=0.0)
        else
            H = label == "g13" ? 2.0 : (label in ["g6","g7","g14","g15","g16","g17","g18"] ? 6.0 : 3.0)
            @named gen = Library.Swing(M=2H, D=1.0, V=Vset, Pm=Pmw / SBASE_MVA)
        end
        push!(components, gen)
    end

    if haskey(SHUNTS_MVAR, label)
        @named shunt = Library.StaticShunt(G=0.0, B=SHUNTS_MVAR[label] / SBASE_MVA)
        push!(components, shunt)
    end

    if label == fault_bus
        @named fault = FaultShunt(G_fault=G_fault, t_fault=t_fault, t_clear=t_clear)
        push!(components, fault)
    end

    bus = compile_bus(MTKBus(components...); vidx=idx(label), name=Symbol("bus_", label))

    pf_model = if haskey(LOAD_OP_A, label)
        Pmw, Qmvar, _ = LOAD_OP_A[label]
        pfPQ(P=-Pmw / SBASE_MVA, Q=-Qmvar / SBASE_MVA)
    elseif haskey(GEN_OP_A, label)
        Pmw, _, Vset = GEN_OP_A[label]
        if label == "g20"
            pfSlack(V=Vset, δ=0.0)
        else
            pfPV(P=Pmw / SBASE_MVA, V=Vset)
        end
    elseif haskey(SHUNTS_MVAR, label)
        pfShunt(G=0.0, B=SHUNTS_MVAR[label] / SBASE_MVA)
    else
        pfPQ(P=0.0, Q=0.0)
    end

    if haskey(INIT_VOLTAGE, label)
        V, deg = INIT_VOLTAGE[label]
        phasor = V * cis(deg2rad(deg))
        set_voltage!(bus, phasor)
        set_voltage!(pf_model, phasor)
    end
    set_pfmodel!(bus, pf_model)

    return bus
end

function make_line(name, from, to, R, X, Bhalf; t_open=Inf)
    @named line_model = TimedPiLine(R=R, X=X, B_src=Bhalf, B_dst=Bhalf, t_open=t_open)
    line = compile_line(MTKLine(line_model); src=idx(from), dst=idx(to), name=Symbol("line_", replace(name, "-" => "_")))
    @named pf_line = Library.PiLine(R=R, X=X, B_src=Bhalf, B_dst=Bhalf)
    set_pfmodel!(line, compile_line(MTKLine(pf_line); src=idx(from), dst=idx(to), name=Symbol("pf_line_", replace(name, "-" => "_"))))
    return line
end

function make_transformer(name, from, to, Xown, n, Snom)
    X = transformer_x_pu(Xown, Snom)
    @named xf = Library.PiLine(R=0.0, X=X, B_src=0.0, B_dst=0.0, r_src=1.0, r_dst=1 / n)
    return compile_line(MTKLine(xf); src=idx(from), dst=idx(to), name=Symbol("xf_", replace(name, "-" => "_")))
end

function build_nordic_network(; G_fault=20.0, t_fault=1.0, t_clear=1.1, trip_line=true)
    buses = [make_bus(label; G_fault, t_fault, t_clear) for label in BUS_LABELS]

    lines = Any[]
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

    transformer_groups = Dict{Tuple{String, String, Float64}, Vector{Tuple{String, Float64}}}()
    for (name, from, to, Xown, n, Snom) in vcat(STEP_UP_TRANSFORMERS, INTERLEVEL_TRANSFORMERS, STEP_DOWN_TRANSFORMERS)
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

    return Network(buses, lines; warn_order=false)
end

function initialize_nordic(nw; verbose=false)
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
        tol=1e-3,
        nwtol=1e-3,
        t=0.0,
    )
    return s0, pfs
end

function graph_path(outfile)
    mkpath(GRAPH_DIR)
    return isempty(dirname(outfile)) || dirname(outfile) == "." ? joinpath(GRAPH_DIR, outfile) : outfile
end

function output_prefix(outprefix, outfile)
    if outfile === nothing
        return graph_path(outprefix)
    end
    path, _ = splitext(graph_path(outfile))
    return path
end

function plot_traces(times, values, monitor, outfile; ylabel, title)
    CairoMakie = Base.require(Base.PkgId(Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0"), "CairoMakie"))
    outpath = graph_path(outfile)
    fig = CairoMakie.Figure(size=(1200, 720))
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="Time [s]",
        ylabel=ylabel,
        title=title,
    )

    for (j, bus) in enumerate(monitor)
        CairoMakie.lines!(ax, times, values[:, j], label=bus)
    end

    CairoMakie.axislegend(ax, "Bus"; position=:rb, framevisible=true)
    CairoMakie.save(outpath, fig)
    return outpath
end

function unwrap_radians!(angles)
    for col in axes(angles, 2)
        offset = 0.0
        previous = angles[1, col]
        for row in 2:size(angles, 1)
            raw = angles[row, col] + offset
            jump = raw - previous
            if jump > π
                offset -= 2π
                raw -= 2π
            elseif jump < -π
                offset += 2π
                raw += 2π
            end
            angles[row, col] = raw
            previous = raw
        end
    end
    return angles
end

function collect_bus_traces(sol, monitor)
    voltages = Matrix{Float64}(undef, length(sol.t), length(monitor))
    angles = similar(voltages)
    active_power = similar(voltages)
    reactive_power = similar(voltages)
    currents = similar(voltages)

    for (row, tt) in enumerate(sol.t)
        for (j, b) in enumerate(monitor)
            bus_idx = idx(b)
            voltages[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊u_mag))
            angles[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊u_arg))
            active_power[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊P))
            reactive_power[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊Q))
            currents[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊i_mag))
        end
    end

    unwrap_radians!(angles)
    return (; voltages, angles, active_power, reactive_power, currents)
end

function simulate_nordic(; tspan=(0.0, 3.0), G_fault=20.0, t_fault=1.0, t_clear=1.1,
                         trip_line=true, saveat=0.005, outprefix="nordic_3s",
                         outfile=nothing, monitor=GENERATOR_BUSES)
    nw = build_nordic_network(; G_fault, t_fault, t_clear, trip_line)
    s0, pfs = initialize_nordic(nw)
    prob = ODEProblem(nw, s0, tspan)
    sol = solve(prob, Rodas5P(); saveat, tstops=[t_fault, t_clear], abstol=1e-7, reltol=1e-7)

    traces = collect_bus_traces(sol, monitor)
    prefix = output_prefix(outprefix, outfile)

    plotfiles = Dict(
        :voltage => plot_traces(
            sol.t, traces.voltages, monitor, "$(prefix)_voltage.png";
            ylabel="Voltage magnitude [pu]",
            title="Nordic test system voltage response",
        ),
        :angle => plot_traces(
            sol.t, traces.angles, monitor, "$(prefix)_angle.png";
            ylabel="Voltage angle [rad]",
            title="Nordic test system voltage angle response",
        ),
        :active_power => plot_traces(
            sol.t, traces.active_power, monitor, "$(prefix)_active_power.png";
            ylabel="Active power [pu]",
            title="Nordic test system active power response",
        ),
        :reactive_power => plot_traces(
            sol.t, traces.reactive_power, monitor, "$(prefix)_reactive_power.png";
            ylabel="Reactive power [pu]",
            title="Nordic test system reactive power response",
        ),
        :current => plot_traces(
            sol.t, traces.currents, monitor, "$(prefix)_current.png";
            ylabel="Current magnitude [pu]",
            title="Nordic test system current response",
        ),
    )

    return (; nw, s0, pfs, sol, monitor, traces, plotfiles)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = simulate_nordic()
    println("retcode = ", result.sol.retcode)
    for key in sort(collect(keys(result.plotfiles)); by=String)
        println("saved $(key) plot to ", result.plotfiles[key])
    end
end
