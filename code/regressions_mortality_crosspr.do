
/* 
Regressions, shuffling the precipitation data source 

*/

clear all
set more off

* CONFIGS
global DO_CROSS_PRS 1 // 1 for running the precipitation cross-checks


local code_dir "~/Documents/projects/hle_iv/code/"
local data_dir "/shared/share_hle/data/reg_inputs/mortality/"
local output_dir "/shared/share_hle/data/reg_outputs/mortality/crossprs/"

local degstr 025deg
local ind_suffix _noIND

if $DO_CROSS_PRS {
	import delimited "`data_dir'/reginputs_mortality_`degstr'.csv", clear

	* I misnamed one of the variables
	gen lr_tavg_adm1_avg = lr_tavg_1981_2010 

	do "`code_dir'/prep_data_multclims.do"

	* bysort year: egen tot_pop = total(population)
	* Since the input dataset has multiple obs for each year, loc, agegroup (multiple
	* obs_mods), these weights can get messed up for the ALLYEARS regression, where the
	* number of obs is not constant across years. This is a workaround. 
	bysort adm1_code year: gen first_obs = (_n == 1)
	bysort year: egen tot_pop = total(population * first_obs)
	gen weight = population / tot_pop
	drop first_obs

	* Factor issue - can't use strings as factors
	encode agegroup, gen(agegroup_) 

	* Turn ERA5-025 into ERA5_025 usw. 
	replace obs_mod = subinstr(obs_mod, "-", "_", .)
	* Get names of obs_mods
	levelsof obs_mod, local(mods)

	// Get joint years
	local total = r(r)
	bysort year (obs_mod): gen tag = (obs_mod != obs_mod[_n-1]) | _n == 1
	bysort year: egen n_present = total(tag)
	keep if n_present == `total'

	* Get label for filenames
	levelsof year, local(yrs)
	local first : word 1 of `yrs'
	local n = wordcount("`yrs'")
	local last : word `n' of `yrs'
	local yrlabel = string(`first') + "-" + string(`last')

	preserve
	    keep iso adm0 adm1_id adm2_id agegroup year obs_mod prcp_poly_*
	    rename prcp_poly_* prcp_poly_*_
	    reshape wide prcp_poly_*_, i(iso adm0 adm1_id adm2_id agegroup year) j(obs_mod) string
	    tempfile pr_wide
	    save `pr_wide'
	restore

	merge m:1 iso adm0 adm1_id adm2_id agegroup year using `pr_wide', nogen

	foreach tas_mod of local mods {
		local tas_mod_str : subinstr local tas_mod "_" "-", all
		foreach pr_mod of local mods {
			local pr_mod_str : subinstr local pr_mod "_" "-", all

			local suffix "_tas-`tas_mod_str'_pr-`pr_mod_str'_`yrlabel'`ind_suffix'"

			* if "`tas_src'" == "`pr_src'" continue // actually, run this as well, to keep them all in the same spot
			
			preserve
				keep if obs_mod == "`tas_mod'"
				if "`ind_suffix'" == "_noIND" drop if iso == "IND"
				gen prcp_poly_1_alt = prcp_poly_1_`pr_mod'
				gen prcp_poly_2_alt = prcp_poly_2_`pr_mod'

				reghdfe deathrate_w99 tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
								[pw = weight] ///
								, absorb(adm0_code##c.prcp_poly_1_alt     adm0_code##c.prcp_poly_2_alt ///
										 i.adm2_code#i.CHN_ts#i.agegroup_  i.adm0_code#i.year#i.agegroup_ ) ///
								cluster(adm1_code) residual(e_hat)
							mat b = e(b)
							mat b = b[1,1..4] // keep main variables' coeffs
							mat V = e(V)
							mat V = V[1..4,1..4] // keep main variables' vcvs

						clear
				svmat b
				outsheet using "`output_dir'/regoutputs_coeffs_mortality_pooled_`degstr'`suffix'.csv", comma replace

				clear
				svmat V
				outsheet using "`output_dir'/regoutputs_vcv_mortality_pooled_`degstr'`suffix'.csv", comma replace
			restore
		}
	}
}




