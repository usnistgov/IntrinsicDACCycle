#Program to take in a sorbent name and directory
#Output the intrinsic DAC cycle
# using Pkg
# Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle/")


using IntrinsicDACCycle
using JSON
using Statistics

function Run_an_Opitimize_IntrinsicDACCycle(name)
    T_start = 250.0 #K
    T_end = 375.0 #K
    P_start = 101325.0 #Pa
    P_end = 0.4 * P_start	

    steps = 500
    #Generate Path 1: Change T, then P
    path1_T_A = LinRange(T_start, T_end, steps)
    path1_T_B = LinRange(T_end, T_end, steps)
    path1_T = vcat(path1_T_A, path1_T_B)

    path1_P_A = LinRange(P_start, P_start, steps)
    path1_P_B = LinRange(P_start, P_end, steps)
    path1_P = vcat(path1_P_A, path1_P_B) 


    name = replace(name, ".json" => "")
    name = replace(name, "/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/" => "")
    
    directory = "/users/asm6/DAC_data/"
    alpha = 400/1000000
    results = IntrinsicDACCycle.Intrinisic_refresh_path(directory, name, 
						 path1_T, path1_P, 
					  	 alpha)

    objectives = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(directory, name, 
								       path1_T, path1_P, 
								       alpha, 100)
    results["Intrinsic_capture_efficiency"] = mean(objectives[1])
    results["Intrinsic_capture_efficiency_std"] = std(objectives[1])
    results["Purity_captured_CO2"] = mean(objectives[2])
    results["Purity_captured_CO2_std"] = std(objectives[2])


    filename = directory*"/SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json"
    open(filename,"w") do f
        JSON.print(f, results, 4)
    end
end


name = ARGS[1]
Run_an_Opitimize_IntrinsicDACCycle(name)