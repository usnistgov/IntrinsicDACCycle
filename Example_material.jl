### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ a6b55dd0-ccf9-11ee-2b3c-6967dde4a5f3
begin
	using Pkg
	Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")

	using Plots
	using Distributed
	using JSON
	using DataFrames
	using QHull

	using Revise

	using IntrinsicDACCycle

	Base_directory = "/users/asm6/DAC_data"
end

# ╔═╡ 3d8b08cc-a08a-4312-ac3d-f697c5898348
# name = "BUSQIQ_clean" 
name = "EVUMUE_clean"

# ╔═╡ 60f4cd07-71f0-4fd7-a1f8-8ff36e9fdeb1
begin
	t1 = range(0, 100, 101) #progression of desorption [fake time units]
    t2 = range(0, 100, 101) #progression of desorption [fake time units]
	
	#Heat then pull vacuum 
	#Isobarically heat 300 K to 350 K
    T1s = 300.0 .+ (0.5 .* t1) #Temperature [K]  
	P1s = 101325 .+ (0 .* t1) #Pressure [Pa] equal to 1 atmosphere of presure
    #Isothermally pull vaccuum from 1 atm to 0.5 atm.
    T2s = T1s[end] .+ (0.0 .* t2) #Temperature [K] 
	P2s = P1s[end] .+ (-101325/200 .* t2) #Pressure [Pa] equal to 1 atmosphere of presure

	path1_Ts = append!(collect(T1s), collect(T2s))
    path1_Ps = append!(collect(P1s), collect(P2s))

	#Heat and pull vacuum at same time
    T1s = 300.0 .+ (0.25 .* t1) #Temperature [K]  
	P1s = 101325 .+ (-101325/400 .* t1) #Pressure [Pa] equal to 1 atmosphere of presure
    
    T2s = T1s[end] .+ (0.25 .* t2) #Temperature [K] 
	P2s = P1s[end] .+ (-101325/400 .* t2) #Pressure [Pa] equal to 1 atmosphere of presure
	
	path2_Ts = append!(collect(T1s), collect(T2s))
    path2_Ps = append!(collect(P1s), collect(P2s))


	#Pull vacuum then heat
	#Isothermally pull vaccuum from 1 atm to 0.5 atm.
    T1s = 300.0 .+ (0.0 .* t1) #Temperature [K]  
	P1s = 101325 .+ (-101325/200 .* t1) #Pressure [Pa] equal to 1 atmosphere of presure
    #Isobarically heat 300 K to 350 K
    T2s = T1s[end] .+ (0.5 .* t2) #Temperature [K] 
	P2s = P1s[end] .+ (0.0 .* t2) #Pressure [Pa] equal to 1 atmosphere of presure

	path3_Ts = append!(collect(T1s), collect(T2s))
    path3_Ps = append!(collect(P1s), collect(P2s))
end

# ╔═╡ 634c97eb-d18e-41c0-9254-5e667ec4af07
begin
	plot(path1_Ts, path1_Ps)
	plot!(path2_Ts, path2_Ps)
	plot!(path3_Ts, path3_Ps)
end

# ╔═╡ 123bc411-b068-4d9b-adf3-eb9c50f0ffa3
begin
	α = 400/1000000

	
	path1_dict = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name, 
                                 path1_Ts, path2_Ps, 
                                 α)

	path2_dict = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name, 
                                 path2_Ts, path2_Ps, 
                                 α)

	path3_dict = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name, 
                                 path3_Ts, path3_Ps, 
                                 α)

end

# ╔═╡ 71794428-7509-4ff1-8061-9f36f46213d3
begin
	path1_purity = path1_dict["Purity_captured_CO2"]
	path2_purity = path2_dict["Purity_captured_CO2"]
	path3_purity = path3_dict["Purity_captured_CO2"]

	path1_efficiency = path1_dict["Intrinsic_capture_efficiency"]
	path2_efficiency = path2_dict["Intrinsic_capture_efficiency"]
	path3_efficiency = path3_dict["Intrinsic_capture_efficiency"]

	@show path1_purity
	@show path2_purity
	@show path3_purity

	@show path1_efficiency
	@show path2_efficiency
	@show path3_efficiency
end

# ╔═╡ 22ad58df-db15-49d1-867b-89afda0cea1c
begin
	plot(path1_dict["Refresh_Path"]["Moles_CO2"])
	plot!(path2_dict["Refresh_Path"]["Moles_CO2"])
	plot!(path3_dict["Refresh_Path"]["Moles_CO2"])

end

# ╔═╡ eb6a6d16-ebac-4270-827d-1ae40b1a81c1
begin
	plot(path1_dict["Refresh_Path"]["Moles_N2"])
	plot!(path2_dict["Refresh_Path"]["Moles_N2"])
	plot!(path3_dict["Refresh_Path"]["Moles_N2"])
end

# ╔═╡ 11b21298-be27-4c16-be81-753913f09136
let
	#Heat and pull vacuum at same time
    T1s = 300.0 .+ (0.25 .* t1) #Temperature [K]  
	P1s = 101325 .+ (-101325/400 .* t1) #Pressure [Pa] equal to 1 atmosphere of presure
    
    T2s = T1s[end] .+ (0.25 .* t2) #Temperature [K] 
	P2s = P1s[end] .+ (-101325/400 .* t2) #Pressure [Pa] equal to 1 atmosphere of presure

	plot(T1s, P1s)
	plot!(T2s, P2s)
end

# ╔═╡ 71a38613-eb6d-4b35-aa54-6edcee2e78f8
path2_Ts .- circshift(path2_Ts,1)

# ╔═╡ 1de218f4-8537-426f-a07d-762cad9c9094
path3_Ts[end]

# ╔═╡ b4b10341-d8eb-4bb1-aff3-ae4f947603ec
begin
	 material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(Base_directory, name)

	close_enough_test = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
end

# ╔═╡ Cell order:
# ╠═a6b55dd0-ccf9-11ee-2b3c-6967dde4a5f3
# ╠═3d8b08cc-a08a-4312-ac3d-f697c5898348
# ╠═60f4cd07-71f0-4fd7-a1f8-8ff36e9fdeb1
# ╠═634c97eb-d18e-41c0-9254-5e667ec4af07
# ╠═123bc411-b068-4d9b-adf3-eb9c50f0ffa3
# ╠═71794428-7509-4ff1-8061-9f36f46213d3
# ╠═22ad58df-db15-49d1-867b-89afda0cea1c
# ╠═eb6a6d16-ebac-4270-827d-1ae40b1a81c1
# ╠═11b21298-be27-4c16-be81-753913f09136
# ╠═71a38613-eb6d-4b35-aa54-6edcee2e78f8
# ╠═1de218f4-8537-426f-a07d-762cad9c9094
# ╠═b4b10341-d8eb-4bb1-aff3-ae4f947603ec
