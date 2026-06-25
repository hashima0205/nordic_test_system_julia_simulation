using PowerDynamics
using OrdinaryDiffEq
using NetworkDynamics
using CairoMakie

# Load LTC network definitions
include("nordic_ltc_simulation.jl")

function run_scenario1()
    println("Building Nordic network with LTC controllers...")
    # Fault at t=1.0s, cleared and line 4032-4044 tripped at t=1.1s
    # LTC controller time constant T_ltc = 30.0s
    nw = build_nordic_network_with_ltc(G_fault=20.0, t_fault=1.0, t_clear=1.1, T_ltc=30.0)
    
    println("Initializing network state...")
    s0, pfs = initialize_nordic_ltc(nw)
    
    println("Preparing ODE simulation problem (0s to 180s)...")
    prob = ODEProblem(nw, s0, (0.0, 180.0))
    
    println("Solving ODE system (this represents Scenario 1: No corrective SIPS)...")
    # Rodas5P solver is ideal for stiff DAE systems.
    # The solver will terminate automatically when voltage collapse occurs (singularity).
    sol = solve(prob, Rodas5P(); reltol=1e-6, abstol=1e-6, saveat=0.2, tstops=[1.0, 1.1])
    
    println("Simulation finished with retcode: ", sol.retcode)
    
    # Track key transmission voltages to observe the collapse
    # Bus 1041 is historically the lowest voltage node in this scenario.
    monitor = ["1041", "1042", "4012", "4062"]
    
    times = sol.t
    voltages = Matrix{Float64}(undef, length(times), length(monitor))
    for (row, tt) in enumerate(times)
        for (j, b) in enumerate(monitor)
            bus_idx = idx(b)
            # Fetch the voltage magnitude state variable (u_mag)
            voltages[row, j] = sol(tt, idxs=VIndex(bus_idx, :busbar₊u_mag))
        end
    end
    
    println("Plotting time series results...")
    fig = Figure(size=(1000, 620))
    ax = Axis(fig[1, 1], 
              xlabel="Time [s]", 
              ylabel="Voltage Magnitude [p.u.]",
              title="Nordic System - Scenario 1 Long-term Voltage Collapse (LTC Enabled)")
    
    # Plot voltage traces
    colors = [:red, :blue, :green, :orange]
    for (j, bus) in enumerate(monitor)
        lines!(ax, times, voltages[:, j], label="Bus $bus", linewidth=2.0, color=colors[j])
    end
    
    # Annotate fault and trip events
    vlines!(ax, [1.0], color=:grey, linestyle=:dash, linewidth=1.2)
    vlines!(ax, [1.1], color=:black, linestyle=:dash, linewidth=1.2)
    text!(ax, 1.5, 0.55, text="Fault & Line Trip (1.0s - 1.1s)", color=:black, fontsize=10)
    
    axislegend(ax, "Transmission Buses", position=:lb)
    ax.limits = (0.0, maximum(times), 0.5, 1.1)
    
    mkpath("graph")
    outpath = "graph/nordic_ltc_scenario1.png"
    save(outpath, fig)
    println("Successfully generated Scenario 1 plot at: ", outpath)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_scenario1()
end
