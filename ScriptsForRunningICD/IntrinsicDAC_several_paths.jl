### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 2c0dd6ac-a2b4-11ef-0362-2bac904fde98
begin 
	using Pkg
	Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle/")
end

# ╔═╡ bd748839-5ec3-4c30-9da6-f57ae5d26dde
begin
	using JSON
	using DataFrames
end
	

# ╔═╡ 1aebb92b-9683-46aa-8a1d-7b213c08bff4
begin
	using IntrinsicDACCycle
end

# ╔═╡ 15f41b25-8e16-436d-9d85-543ed89e2828
begin
	using Statistics
end

# ╔═╡ eb6fa2ec-5010-4b3a-8553-df8ef7f17930
begin
	using Plots
end

# ╔═╡ 2f3afb0a-bd61-4240-ad63-0e431b8a737b
html"""<style>
main {
    max-width: 1900px;
}
"""

# ╔═╡ 230a9b22-6dfb-44bc-b89d-3dee927d68f8
begin
	Pareto_Optimized = JSON.parsefile("/users/asm6/Julia_scripts/Optimized_IDC_pareto.json")
end

# ╔═╡ 2c1c993b-567b-4d97-8f3a-10c8a0054eda
list_of_pareto_names = unique(Pareto_Optimized["Name"])

# ╔═╡ 78e61564-6bb3-4bd9-b772-37ce270da58c
begin
	Base_directory = "/users/asm6/DAC_data/"
	Optimized_subdirectory = "Optimized_Intrinsic_Cycle/"  

	name = list_of_pareto_names[2]

	optimized_dict = JSON.parsefile(Base_directory*Optimized_subdirectory*"OptIDC_"*name*".json")


end

# ╔═╡ 99a845f5-dd1b-4e57-aa05-3bd403572c23
name

# ╔═╡ 43faf0f7-fc2c-443e-89b3-b794e81827e6
length(optimized_dict["Refresh_Path"]["Pressures"][1])

# ╔═╡ d38222ab-70bc-4835-be47-274ec4c3a07f
40311/101325

# ╔═╡ 2b4b95da-67eb-4123-a600-2e2615f5f8c2
begin
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

	#Generate Path 2: Change both T and P
	path2_T = LinRange(T_start, T_end, 2*steps)
	path2_P = LinRange(P_start, P_end, 2*steps)

	#Generate Path 3: Change P, then P
	path3_T_A = LinRange(T_start, T_start, steps)
	path3_T_B = LinRange(T_start, T_end, steps)
	path3_T = vcat(path3_T_A, path3_T_B)
	
	path3_P_A = LinRange(P_start, P_end, steps)
	path3_P_B = LinRange(P_end, P_end, steps)
	path3_P = vcat(path3_P_A, path3_P_B)
end

# ╔═╡ 936f4bb4-9e37-4d1a-9e17-4e8ed819877f
begin
	α = 400/1000000
	path1_dict = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name, 
																	path1_T, path1_P, 
																	α)
	path1_dict_w_err = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(Base_directory, name, 
																	path1_T, path1_P, 
																	α, 100)

	path2_dict = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name, 
																	path2_T, path2_P, 
																	α)
	path2_dict_w_err = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(Base_directory, name, 
																	path2_T, path2_P, 
																	α, 100)

	path3_dict = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name, 
																	path3_T, path3_P, 
																	α)
	path3_dict_w_err = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(Base_directory, name, 
																	path3_T, path3_P, 
																	α, 100)
end

# ╔═╡ 58689022-cfbc-4dbd-9973-8e8a67302371
begin
	plot(path1_T, path1_P, label = "Path 1")
	plot!(path2_T, path2_P, label = "Path 2")
	plot!(path3_T, path3_P, label = "Path 3")
	xlabel!("Temperature (K)")
	ylabel!("Pressure (Pa)")
end

# ╔═╡ 1902901a-8152-40d0-a185-fb40a4e1a758
begin
	#Path 1 Performance
	path1_ξ = path1_dict["Intrinsic_capture_efficiency"]
	path1_meanξ = mean(path1_dict_w_err[1])
	path1_stdξ = std(path1_dict_w_err[1])
	
	path1_α = path1_dict["Purity_captured_CO2"]
	path1_meanα = mean(path1_dict_w_err[2])
	path1_stdα = std(path1_dict_w_err[2])
	
	#Path 2 Performance
	path2_ξ = path2_dict["Intrinsic_capture_efficiency"]
	path2_meanξ = mean(path2_dict_w_err[1])
	path2_stdξ = std(path2_dict_w_err[1])

	path2_α = path2_dict["Purity_captured_CO2"]
	path2_meanα = mean(path2_dict_w_err[2])
	path2_stdα = std(path2_dict_w_err[2])

	#Path 2 Performance
	path3_ξ = path3_dict["Intrinsic_capture_efficiency"]
	path3_meanξ = mean(path3_dict_w_err[1])
	path3_stdξ = std(path3_dict_w_err[1])

	path3_α = path3_dict["Purity_captured_CO2"]
	path3_meanα = mean(path3_dict_w_err[2])
	path3_stdα = std(path3_dict_w_err[2])

end

# ╔═╡ 5414f3c2-fcc6-437f-9005-d7ca3528f8ab
begin
	plot()
	# scatter([path1_α], [path1_ξ], label = "Path 1 nominal", c = 1)
	scatter!([path1_meanα], [path1_meanξ],
			xerror = [path1_stdα],
			yerror = [path1_stdξ], label = "Path 1", c = 1)

	# scatter!([path2_α], [path2_ξ], label = "Path 2 nominal", c = 2)
	scatter!([path2_meanα], [path2_meanξ],
			xerror = [path2_stdα],
			yerror = [path2_stdξ], label = "Path 2", c = 2)

	# scatter!([path3_α], [path3_ξ], label = "Path 3 nominal", c = 3)
	scatter!([path3_meanα], [path3_meanξ],
			xerror = [path3_stdα],
			yerror = [path3_stdξ], label = "Path 3", c = 3)

	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end
	

# ╔═╡ 28b2d8fa-cd37-44ce-a63f-078329607907
begin
	plot1 = plot(path1_T, path1_P, label = "Path 1")
	plot1 = plot!(path2_T, path2_P, label = "Path 2")
	plot1 = plot!(path3_T, path3_P, label = "Path 3")
	plot1 = xlabel!("Temperature (K)")
	plot1 = ylabel!("Pressure (Pa)")
	plot1 = title!("a)", titlelocation = :left)

	plot2 = plot()
	plot2 = scatter!([path1_meanα], [path1_meanξ],
			xerror = [path1_stdα],
			yerror = [path1_stdξ], label = "Path 1", c = 1)

	plot2 = scatter!([path2_meanα], [path2_meanξ],
			xerror = [path2_stdα],
			yerror = [path2_stdξ], label = "Path 2", c = 2)

	plot2 = scatter!([path3_meanα], [path3_meanξ],
			xerror = [path3_stdα],
			yerror = [path3_stdξ], label = "Path 3", c = 3)

	plot2 = xlabel!("Purity Captured CO₂")
	plot2 = ylabel!("Intrinsic Capture Efficiency (mol/J)")
	plot2 = title!("b)", titlelocation = :left)

	plot(plot1, plot2, layout = (1, 2))
	plot!(size=(900,400), left_margin = 15*Plots.mm, bottom_margin = 10*Plots.mm)
	savefig("/users/asm6/Julia_scripts/OneSorb_SeveralPaths.svg")
end

# ╔═╡ 48066678-cdde-4b51-b7e2-58f77e55ad6a
begin
	path2_dict
end

# ╔═╡ 76f1dbcb-5933-4a43-83e0-201c40ff91b1
begin
	incriments = collect(1:length(path2_T))./length(path2_T)
	step2_dict  = path2_dict["E_Balance"]["Step_2"]
	plot3 = plot(incriments, step2_dict["Heat_to_desorb_CO2"], label = "Q_CO2")
	plot3 = plot!(incriments, step2_dict["Heat_to_desorb_N2"], label = "Q_N2")
	plot3 = plot!(incriments, step2_dict["Work_to_desorb_CO2"], label = "W_CO2")
	plot3 = plot!(incriments, step2_dict["Work_to_desorb_N2"], label = "W_N2")
	plot3 = plot!(incriments, step2_dict["E_to_heat_adsorbed_CO2"], label = "E_CO2")
	plot3 = plot!(incriments, step2_dict["E_to_heat_adsorbed_N2"], label = "E_N2")
	plot3 = plot!(incriments, step2_dict["E_to_heat_sorbent"], label = "E_sorb")
	plot3 = plot!(incriments, step2_dict["E_to_change_pressure"], label = "E_P")
	plot3 = title!("a)", titlelocation = :left)
	plot3 = plot!(legend = :right)
	plot3 = xlabel!("Step 2 Progress")
	plot3 = ylabel!("Energy (J/kg)")
	

	adsorb_dict = path2_dict["Refresh_Path"]
	n_CO2 = adsorb_dict["Moles_CO2"]
	n_N2 = adsorb_dict["Moles_N2"]
	#plot CO2 loop
	plot4 = plot([0,0], [n_CO2[1], n_CO2[end]], label = "CO2", c = 1) #Step 1
	plot4 = plot!(incriments, n_CO2, label = :none, c = 1) #Step 2
	plot4 = plot!([1,0], [n_CO2[end], n_CO2[end]], label = :none, c = 1) #Step 1
	#plot N2 loop
	plot4 = plot!([0,0], [n_N2[1], n_N2[end]], label = "N2", c = 2) #Step 1
	plot4 = plot!(incriments, n_N2, label = :none, c = 2) #Step 2
	plot4 = plot!([1,0], [n_N2[end], n_N2[end]], label = :none, c = 2) #Step 1
	#Label axis
	plot4 = xlabel!("Step 2 Progress")
	plot4 = ylabel!("Uptake (mol/kg)")
	plot4 = title!("b)", titlelocation = :left)
	plot4 = plot!(legend = :right)

	plot(plot3, plot4, layout = (1, 2))
	plot!(size=(900,400), left_margin = 15*Plots.mm, bottom_margin = 10*Plots.mm)
	savefig("/users/asm6/Julia_scripts/OneSorb_Progress.svg")
end

# ╔═╡ 7773538f-45bf-4d00-be8c-ceba15b8997a
collect(1:length(path2_T))./length(path2_T)

# ╔═╡ cfad142b-5d27-4cc4-81e5-af1a2f211904
begin
	(path1_dict["Intrinsic_capture_efficiency"], path1_dict["Purity_captured_CO2"])
	
end

# ╔═╡ 8c2786ba-959f-4ddd-b6c1-fd5a8ff0bd21
begin
	(path2_dict["Intrinsic_capture_efficiency"], path2_dict["Purity_captured_CO2"])
	
end

# ╔═╡ 2eb4f8c6-2475-425c-a5d3-e827a5eff25e
begin
	(path3_dict["Intrinsic_capture_efficiency"], path3_dict["Purity_captured_CO2"])
	
end

# ╔═╡ Cell order:
# ╠═2c0dd6ac-a2b4-11ef-0362-2bac904fde98
# ╠═bd748839-5ec3-4c30-9da6-f57ae5d26dde
# ╠═1aebb92b-9683-46aa-8a1d-7b213c08bff4
# ╠═15f41b25-8e16-436d-9d85-543ed89e2828
# ╠═eb6fa2ec-5010-4b3a-8553-df8ef7f17930
# ╠═2f3afb0a-bd61-4240-ad63-0e431b8a737b
# ╠═230a9b22-6dfb-44bc-b89d-3dee927d68f8
# ╠═2c1c993b-567b-4d97-8f3a-10c8a0054eda
# ╠═78e61564-6bb3-4bd9-b772-37ce270da58c
# ╠═99a845f5-dd1b-4e57-aa05-3bd403572c23
# ╠═43faf0f7-fc2c-443e-89b3-b794e81827e6
# ╠═d38222ab-70bc-4835-be47-274ec4c3a07f
# ╠═2b4b95da-67eb-4123-a600-2e2615f5f8c2
# ╠═936f4bb4-9e37-4d1a-9e17-4e8ed819877f
# ╠═58689022-cfbc-4dbd-9973-8e8a67302371
# ╠═1902901a-8152-40d0-a185-fb40a4e1a758
# ╠═5414f3c2-fcc6-437f-9005-d7ca3528f8ab
# ╠═28b2d8fa-cd37-44ce-a63f-078329607907
# ╠═48066678-cdde-4b51-b7e2-58f77e55ad6a
# ╠═76f1dbcb-5933-4a43-83e0-201c40ff91b1
# ╠═7773538f-45bf-4d00-be8c-ceba15b8997a
# ╠═cfad142b-5d27-4cc4-81e5-af1a2f211904
# ╠═8c2786ba-959f-4ddd-b6c1-fd5a8ff0bd21
# ╠═2eb4f8c6-2475-425c-a5d3-e827a5eff25e
