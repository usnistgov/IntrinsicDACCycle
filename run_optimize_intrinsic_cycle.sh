search_dir=/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/
sorbentfiles=$search_dir'*.json'


parallel -j 10 julia --project=Project.toml run_optimize_intrinsic_dac.jl ::: $sorbentfiles


# julia run_optimize_intrinsic_dac.jl ABAYIO_clean