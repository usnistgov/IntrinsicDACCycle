#Program to take in a sorbent name and directory
#Output the optimized intrinsic DAC cycle
# using Pkg
# Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle/")


using IntrinsicDACCycle
using JSON

function Run_an_Opitimize_IntrinsicDACCycle(name)
    name = replace(name, ".json" => "")
    name = replace(name, "/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/" => "")
    
    directory = "/users/asm6/DAC_data/"
    α = 400/1000000
    results = IntrinsicDACCycle.Optimize_Intrinsic_Refresh_w_err(directory, name, α)

    filename = directory*"/Optimized_Intrinsic_Cycle/OptIDC_"*name*".json"
    open(filename,"w") do f
        JSON.print(f, results, 4)
    end
end


name = ARGS[1]
Run_an_Opitimize_IntrinsicDACCycle(name)