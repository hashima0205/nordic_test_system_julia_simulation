using PowerDynamics
using OrdinaryDiffEq
using NetworkDynamics
using CairoMakie

# Define size and axes on VIndex so it can be evaluated by SciML integrator callbacks
Base.size(::NetworkDynamics.VIndex) = ()
Base.axes(::NetworkDynamics.VIndex) = ()

# Load dynamic network definitions
include("nordic_dynamic_simulation.jl")

function run_dynamic_scenario1()
    println("Building Nordic network with dynamic generators (GENROU/GENSAL + AVR + OEL)...")
    # Fault at t=2.0s, cleared and line 4032-4044 tripped at t=2.1s
    # LTC controller time constant T_ltc = 30.0s
    nw = build_nordic_network_dynamic(G_fault=20.0, t_fault=2.0, t_clear=2.1, T_ltc=30.0)
    
    println("Initializing network state (power flow + component initialization)...")
    s0, pfs = initialize_nordic_dynamic(nw; verbose=true)
    
    println("Preparing ODE simulation problem (0s to 200s)...")
    prob = ODEProblem(nw, s0, (0.0, 200.0))
    
    println("Solving ODE system (this represents Scenario 1 with dynamic OELs)...")
    # Rodas5P solver is ideal for stiff DAE systems.
    sol = solve(prob, Rodas5P(); reltol=1e-5, abstol=1e-5, saveat=0.2, tstops=[2.0, 2.1])
    
    println("Simulation finished with retcode: ", sol.retcode)
    println("Simulation duration: ", sol.t[end], " seconds")
    
    # Track key transmission voltages to observe the collapse
    monitor = ["1041", "1042", "4012", "4062"]
    
    times = sol.t
    voltages = Matrix{Float64}(undef, length(times), length(monitor))
    for (row, tt) in enumerate(times)
        for (j, b) in enumerate(monitor)
            bus_idx = idx(b)
            voltages[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊u_mag))
        end
    end
    
    println("Plotting time series results...")
    fig = Figure(size=(1000, 620))
    ax = Axis(fig[1, 1], 
              xlabel="Time [s]", 
              ylabel="Voltage Magnitude [p.u.]",
              title="Nordic System - Scenario 1 Long-term Voltage Collapse (Dynamic OELs Active)")
    
    # Plot voltage traces
    colors = [:red, :blue, :green, :orange]
    for (j, bus) in enumerate(monitor)
        lines!(ax, times, voltages[:, j], label="Bus $bus", linewidth=2.0, color=colors[j])
    end
    
    # Annotate fault and trip events
    vlines!(ax, [2.0], color=:grey, linestyle=:dash, linewidth=1.2)
    vlines!(ax, [2.1], color=:black, linestyle=:dash, linewidth=1.2)
    text!(ax, 3.0, 0.55, text="Fault & Line Trip (2.0s - 2.1s)", color=:black, fontsize=10)
    
    # Annotate LTC action start
    vlines!(ax, [32.0], color=:blue, linestyle=:dash, linewidth=1.2)
    text!(ax, 33.0, 0.60, text="LTC Action Start (32.0s)", color=:blue, fontsize=10)
    
    axislegend(ax, "Transmission Buses", position=:lb)
    ax.limits = (0.0, max(200.0, maximum(times)), 0.5, 1.1)
    
    mkpath("graph")
    outpath = "graph/nordic_dynamic_scenario1.png"
    save(outpath, fig)
    println("Successfully generated Scenario 1 plot at: ", outpath)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_dynamic_scenario1()
end
