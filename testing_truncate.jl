### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ a9fd91ee-2da0-11ef-0469-1d5851fbc0fb
using Pkg

# ╔═╡ 67f28ebc-bc4b-43f0-9d28-26a4f12d28c8
Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")

# ╔═╡ c5724d7d-a31c-44fb-9398-cf2206d923f1
using JSON

# ╔═╡ b6be7fd9-36f8-43ed-a4b4-0f649b7732b8
using IntrinsicDACCycle

# ╔═╡ f69abfac-d417-4cc6-ae55-2ae9cf32c6c1
Pkg.instantiate()

# ╔═╡ 716fd2ad-f4d5-4557-91ff-b84e45fbc11c
cd("/users/asm6/Julia_scripts/IntrinsicDACCycle")

# ╔═╡ f5d3aea7-bebd-4e55-bed4-cc572044cc3b
Base_directory = "/users/asm6/DAC_data"

# ╔═╡ 2fe02528-9613-4a4a-b23a-ae510982cd54
begin
	#get all the material files
	list_of_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"/CSD_FEASST_Materials/Materials/"))
	#strip off the .json tag
    list_of_materials = replace.(list_of_material_files, ".json" => "")
	#filter for _clean matierals
	# list_of_clean_materials = filter(x -> occursin.("_clean", x), list_of_materials)
end


# ╔═╡ bb0b7492-78a6-4331-a37c-9bac6323f251
IntrinsicDACCycle.Intrinisic_refresh(Base_directory, list_of_materials[1])

# ╔═╡ d533c96f-2b69-4f98-a507-b0234f388e94
begin
	directory = Base_directory
	name = list_of_materials[1]
	#Read in all the GCMC results
    material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(directory, name)

    # #Test if the Henry constant at 1 atm is close enough to the direct GCMC at 1 atm
    # close_enough_test = Close_enough(material, Kh_N₂, One_atm_N₂)
    # #If not close enough
    # if close_enough_test == false
    #     return nothing
    # end

    #Start Dictionaries for the results
    Results_Dict = sort(Dict{String, Any}("Name" => name))
    Path_Dict = sort(Dict{String, Any}("Refresh_Path" => "Definition of refresh path and material properties along that path (in per kg of sorbent basis)."))
    E_Balance_Dict = sort(Dict{String, Any}("E_Balance" => "Energy balance along path"))
    Step_1_Dict = sort(Dict{String, Any}("Step_1" => "Adsorption"))
    Step_2_Dict = sort(Dict{String, Any}("Step_2" => "Desorption"))
    Step_3_Dict = sort(Dict{String, Any}("Step_3" => "Waste energy recovery"))
    
    #Define refresh cycle (T, P) path and inlet concentration
    t1 = range(0, 100, 101) #progression of desorption [fake time units]
    t2 = range(0, 100, 101) #progression of desorption [fake time units]
	#Isobarically heat 300 K to 350 K
    T1s = 300.0 .+ (0.5 .* t1) #Temperature [K]  
	P1s = 101325 .+ (0 .* t1) #Pressure [Pa] equal to 1 atmosphere of presure
    #Isothermally pull vaccuum from 1 atm to 0.5 atm.
    T2s = T1s[end] .+ (0.0 .* t2) #Temperature [K] 
	P2s = P1s[end] .+ (-101325/200 .* t2) #Pressure [Pa] equal to 1 atmosphere of presure
	#Concatonate the process steps
    Ts = append!(collect(T1s), collect(T2s))
    Ps = append!(collect(P1s), collect(P2s))
    α = 400/1000000 #400 ppm is the concentration of CO2 in ambient air
    βs = IntrinsicDACCycle.T_to_β.(Ts) #[mol/kJ]

    #Test if the Henry constant at 1 atm is close enough to the direct GCMC at 1 atm
    close_enough_test = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
    #If not close enough
    if close_enough_test == false
        #set Ts and Ps and βs to NaN
        Ts = Ts .* NaN
        Ps = Ps .* NaN
        βs = βs .* NaN
    end

    #Extrapolate Henry constants along the path
    #Extrapolate the CO2 isotherm to the βs
    Henry_CO2, Henry_CO2_err = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_CO₂, material) #[mmol/(kg Pa)]

    #Extrapolate the N2 isotherm to the βs
    Henry_N2, Henry_N2_err = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_N₂, material)  #[mmol/(kg Pa)]
	@show length(Henry_CO2)
    #If extrapolating outside a reasonable range (for that material),
    # the Henry Constants will un-physically increase with increasing temperature
    #Keep only the monotonically decreasing parts of the Henry constants
	if ~isnan(Henry_CO2[1])
	
	    mono_Henry_CO2, CO2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_CO2)
	    mono_Henry_N2, N2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_N2)
	    #Choose whichever set of indices is smaller
	    index_lenghts = [length(CO2_indices), length(N2_indices)]
		choice_of_indices = [CO2_indices, N2_indices]
		indices = choice_of_indices[argmin(index_lenghts)]
	
	    #Apply the index truncation to all relevant variables
	    Henry_CO2 = Henry_CO2[indices]
	    Henry_CO2_err = Henry_CO2_err[indices]
		@show length(Henry_CO2)
	    Henry_N2 = Henry_N2[indices]
	    Henry_N2_err = Henry_N2_err[indices]
	
	    Ts = Ts[indices]
	    Ps = Ps[indices]
	    βs = βs[indices]
		
	    #Compare the CO2 Henry constant to Saturation uptake of CO2
	    #Truncate to a sensible range
	    Ts, Ps, βs, Henry_CO2, Henry_CO2_err, Henry_N2, Henry_N2_err = IntrinsicDACCycle.truncate_to_saturation(directory, name, α,
	                                                                                          Ts, Ps, βs, 
	                                                                                          Henry_CO2, Henry_CO2_err, 
	                                                                                          Henry_N2, Henry_N2_err) 
	end 
	@show length(Henry_CO2)
end

# ╔═╡ 42f4b2e9-f9ce-4e2c-9050-e18fb7a2f674
Henry_N2_err

# ╔═╡ Cell order:
# ╠═a9fd91ee-2da0-11ef-0469-1d5851fbc0fb
# ╠═67f28ebc-bc4b-43f0-9d28-26a4f12d28c8
# ╠═f69abfac-d417-4cc6-ae55-2ae9cf32c6c1
# ╠═c5724d7d-a31c-44fb-9398-cf2206d923f1
# ╠═b6be7fd9-36f8-43ed-a4b4-0f649b7732b8
# ╠═716fd2ad-f4d5-4557-91ff-b84e45fbc11c
# ╠═f5d3aea7-bebd-4e55-bed4-cc572044cc3b
# ╠═2fe02528-9613-4a4a-b23a-ae510982cd54
# ╠═bb0b7492-78a6-4331-a37c-9bac6323f251
# ╠═d533c96f-2b69-4f98-a507-b0234f388e94
# ╠═42f4b2e9-f9ce-4e2c-9050-e18fb7a2f674
