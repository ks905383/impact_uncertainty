/* 
Regressions across products 

TO-DO: remove India from sample, it wasn't used in the original paper either? 

*/

clear all
set more off

* CONFIGS
global DO_LEAVEONEOUT 0 // 1 for run the leave-one-out tests
global DO_BASIC 1 // 1 for running the single global dmgf
global DO_BASIC_COMMONYEARS 1 // 1 for running the single global dmgf across only years with data for all prods 
* global DO_INTERACTED 0 // 1 for running the fully interacted model
* global DO_INTERACTED_NOTOTAL 0 // 1 for running while dropping the total pop group
global DO_JOINT_BASIC 1 // 1 for running joint (all reanalyses at once) pooled global dmgf
* global DO JOINT_INTERACTED 1 // 1 for running the joint (all reanalyses at once) interacted model


local code_dir "~/Documents/projects/hle_iv/code/"
local data_dir "/shared/share_hle/data/reg_inputs/mortality/"
local output_dir "/shared/share_hle/data/reg_outputs/mortality/"

local degstrs 025deg 1deg

if $DO_BASIC {
	foreach degstr of local degstrs {
        import delimited "`data_dir'/reginputs_mortality_`degstr'.csv", clear

		* I misnamed one of the variables
		gen lr_tavg_adm1_avg = lr_tavg_insamp

		do "`code_dir'/prep_data_multclims.do"

		* 1. set weighting schemes
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

		* Get names of obs_mods
		levelsof obs_mod, local(mods)

		foreach mod of local mods {
			foreach ind_suffix in "" "_noIND" {
				local suffix "_`mod'_ALLYEARS`ind_suffix'"

				preserve
	                keep if obs_mod == "`mod'"
	                if "`ind_suffix'" == "_noIND" drop if iso == "IND"

	                reghdfe deathrate_w99 tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
							[pw = weight] ///
							, absorb(adm0_code##c.prcp_poly_1     adm0_code##c.prcp_poly_2 ///
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
}

if $DO_BASIC_COMMONYEARS {
	foreach degstr of local degstrs {
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

		foreach mod of local mods {
			foreach ind_suffix in "" "_noIND" {
				local suffix "_`mod'_`yrlabel'`ind_suffix'"

				preserve
	                keep if obs_mod == "`mod'"
	                if "`ind_suffix'" == "_noIND" drop if iso == "IND"

	                reghdfe deathrate_w99 tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4 ///
							[pw = weight] ///
							, absorb(adm0_code##c.prcp_poly_1     adm0_code##c.prcp_poly_2 ///
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
}


if $DO_JOINT_BASIC {
    foreach degstr of local degstrs {
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
		encode obs_mod, gen(obs_mod_) 

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

        * Get level idxs of mods (before the preserve blocks)
        foreach mod of local mods {
		    quietly summarize obs_mod_ if obs_mod == "`mod'"
		    local mod_clean = subinstr("`mod'", "-", "_", .)
		    local mod_n_`mod_clean' = r(min)
		}

        foreach ind_suffix in "_noIND" "" {

        	preserve
        		if "`ind_suffix'" == "_noIND" drop if iso == "IND"

		        * ---- Reg 1: joint across obs products, but not interacting with FEs (_JOINTnotFEs)
		        reghdfe deathrate_w99 c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.obs_mod_ ///
		                [pw = weight] ///
		                , absorb(adm0_code##c.prcp_poly_1     adm0_code##c.prcp_poly_2 ///
		                         i.adm2_code#i.CHN_ts#i.agegroup_  i.adm0_code#i.year#i.agegroup_ ) ///
		                cluster(adm1_code) residual(e_hat_notFEs)

		        mat b_notFEs = e(b)
		        mat V_notFEs = e(V)

		        * ---- Reg 2: joint across obs products, including interacted with FEs (_JOINT)
		        reghdfe deathrate_w99 c.(tavg_poly_1 tavg_poly_2 tavg_poly_3 tavg_poly_4)#i.obs_mod_ ///
		                [pw = weight] ///
		                , absorb(i.adm0_code#i.obs_mod_  i.adm0_code#c.prcp_poly_1#i.obs_mod_  i.adm0_code#c.prcp_poly_2#i.obs_mod_ ///
		                 i.adm2_code#i.CHN_ts#i.agegroup_#i.obs_mod_  i.adm0_code#i.year#i.agegroup_#i.obs_mod_ ) ///
		                cluster(adm1_code) residual(e_hat_JOINT)

		        mat b_JOINT = e(b)
		        mat V_JOINT = e(V)

		        * ---- Save coefficients, VCVs
		        foreach mod of local mods {
		            local mod_clean = subinstr("`mod'", "-", "_", .)
    				local mod_n = `mod_n_`mod_clean''

		            local p1 "`mod_n'.obs_mod_#c.tavg_poly_1"
		            local p2 "`mod_n'.obs_mod_#c.tavg_poly_2"
		            local p3 "`mod_n'.obs_mod_#c.tavg_poly_3"
		            local p4 "`mod_n'.obs_mod_#c.tavg_poly_4"

		            foreach variant in notFEs JOINT {
		                if "`variant'" == "notFEs" {
		                    local modsuffix "_`mod'_`yrlabel'_JOINTnotFEs`ind_suffix'"
		                    mat b_use = b_notFEs
		                    mat V_use = V_notFEs
		                }
		                else {
		                    local modsuffix "_`mod'_`yrlabel'_JOINT`ind_suffix'"
		                    mat b_use = b_JOINT
		                    mat V_use = V_JOINT
		                }

		                mat b_mod = b_use[1, "`p1'"], b_use[1, "`p2'"], b_use[1, "`p3'"], b_use[1, "`p4'"]

		                mat V_mod = (V_use["`p1'","`p1'"], V_use["`p1'","`p2'"], V_use["`p1'","`p3'"], V_use["`p1'","`p4'"] \ ///
		                             V_use["`p2'","`p1'"], V_use["`p2'","`p2'"], V_use["`p2'","`p3'"], V_use["`p2'","`p4'"] \ ///
		                             V_use["`p3'","`p1'"], V_use["`p3'","`p2'"], V_use["`p3'","`p3'"], V_use["`p3'","`p4'"] \ ///
		                             V_use["`p4'","`p1'"], V_use["`p4'","`p2'"], V_use["`p4'","`p3'"], V_use["`p4'","`p4'"])

		                *preserve
	                    clear
	                    svmat b_mod
	                    outsheet using "`output_dir'/regoutputs_coeffs_mortality_pooled_`degstr'`modsuffix'.csv", comma replace
		                *restore

		                *preserve
	                    clear
	                    svmat V_mod
	                    outsheet using "`output_dir'/regoutputs_vcv_mortality_pooled_`degstr'`modsuffix'.csv", comma replace
		                *restore
		            }
	            }
       		restore
        }
    }
}