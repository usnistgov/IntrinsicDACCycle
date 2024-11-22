### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 081ad48c-ce77-4aa6-85b5-2c739eb8d4ac
using Pkg

# ╔═╡ e78577bb-3ea2-4172-aa82-e3e6dbc02330
Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")

# ╔═╡ 0d9063f6-2313-4794-a9d5-18ca263177e2
using JSON

# ╔═╡ 16682450-cce6-11ee-2cf4-c31f93bd2a9c
Base_directory = "/users/asm6/DAC_data"

# ╔═╡ f8eacc0d-8869-4af6-b9c8-30e1df3e4ab8
begin
	#get all the material files
	list_of_material_files = filter(x -> occursin.("data.",x), readdir(Base_directory*"/CSD_FEASST_Materials/Materials/"))
	#strip off the .json tag
    # list_of_materials = replace.(list_of_material_files, ".json" => "")
	#filter for _clean matierals
	# list_of_clean_materials = filter(x -> occursin.("_clean", x), list_of_materials)
end

# ╔═╡ 41f5e4ae-25da-4667-b17e-fcd6c11fa46c
begin
	for file in list_of_material_files
		full_file = Base_directory*"/CSD_FEASST_Materials/Materials/"*file
		f = open(full_file, "r")
		s = read(f, String)
	
		num_Al = length(collect(eachmatch(r"Al", s)))
	
		if num_Al > 0
			num_O = length(collect(eachmatch(r"O", s)))	
			O_to_Al = num_O / num_Al
			if O_to_Al > 5
				if O_to_Al < 7
					@show file
				end
			end
			
		end
	end
end

# ╔═╡ 139e1b7d-44a4-406d-87ce-cb27131bb303


# ╔═╡ 78e3092d-477d-4c39-8911-519aef78735d


# ╔═╡ 6e2c4b26-8986-4708-8848-27e97f9172e4
begin
	full_file = Base_directory*"/CSD_FEASST_Materials/Materials/"*list_of_material_files[1]
	f = open(full_file, "r")
	s = read(f, String)
	num_C = length(collect(eachmatch(r"C", s)))	
end

# ╔═╡ Cell order:
# ╠═081ad48c-ce77-4aa6-85b5-2c739eb8d4ac
# ╠═e78577bb-3ea2-4172-aa82-e3e6dbc02330
# ╠═0d9063f6-2313-4794-a9d5-18ca263177e2
# ╠═16682450-cce6-11ee-2cf4-c31f93bd2a9c
# ╠═f8eacc0d-8869-4af6-b9c8-30e1df3e4ab8
# ╠═41f5e4ae-25da-4667-b17e-fcd6c11fa46c
# ╠═139e1b7d-44a4-406d-87ce-cb27131bb303
# ╠═78e3092d-477d-4c39-8911-519aef78735d
# ╠═6e2c4b26-8986-4708-8848-27e97f9172e4
