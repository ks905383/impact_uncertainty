* GDPPC regressions, adapted from BHM15 replication code
//		Note: users will need to install parmest command (e.g. "ssc install parmest")

clear all
set more off

global DO_BYMOD 1 // 1 for run the regressions for each mod separately
global DO_BYMOD_COMMONYEARS 1 // 1 for run the regressions only in years with all mods
global DO_JOINT 1 // 1 for run the joint regression, only in years with all mods

local data_dir "/shared/share_hle/data/reg_inputs/gdppc/"
local output_dir "/shared/share_hle/data/reg_outputs/gdppc/"

local degstrs 025deg 1deg


// ---------- ESTIMATE GLOBAL RESPONSE WITH BASELINE REGRESSION SPECIFICATION   
// ---------- THEN WRITE OUT FOR CONSTUCTION OF FIGURE 2, panel A

// ---- Separate Regressionen für jedes obs_mod, alle Jahre
if $DO_BYMOD {
    foreach degstr of local degstrs {
        import delimited "`data_dir'/reginputs_gdppc_`degstr'.csv", clear
        encode iso_id, gen(iso_id_)
        gen temp = temp_popweight
        gen precip = precip_popweight
        gen precip2 = precip_popweight_2

        levelsof obs_mod, local(mods)

        foreach mod of local mods {
            local suffix "_`mod'_ALLYEARS"

            preserve
                keep if obs_mod == "`mod'"

                reg growthwdi c.temp##c.temp precip precip2 i.year _yi_* _y2_* i.iso_id_, cluster(iso_id_)
                    mat b = e(b)
                    mat b = b[1,1..2]
                    mat V = e(V)

                clear
                svmat b
                outsheet using "`output_dir'/regoutputs_coeffs_gdppc_`degstr'`suffix'.csv", comma replace

                mat Vtemp = V[1..2,1..2]
                clear
                svmat Vtemp
                outsheet using "`output_dir'/regoutputs_vcv_gdppc_`degstr'`suffix'.csv", comma replace
            restore
        }
    }
}

// ---- Separate Regressionen für jedes obs_mod, nur gemeinsame Jahre
if $DO_BYMOD_COMMONYEARS {
    foreach degstr of local degstrs {
        import delimited "`data_dir'/reginputs_gdppc_`degstr'.csv", clear
        encode iso_id, gen(iso_id_)
        gen temp = temp_popweight
        gen precip = precip_popweight
        gen precip2 = precip_popweight_2

        // Get joint years
        levelsof obs_mod, local(mods)
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
            local suffix "_`mod'_`yrlabel'"

            preserve
                keep if obs_mod == "`mod'"

                // Regression
                reg growthwdi c.temp##c.temp precip precip2 i.year _yi_* _y2_* i.iso_id_, cluster(iso_id_)
                    mat b = e(b)
                    mat b = b[1,1..2]
                    mat V = e(V)

                clear
                svmat b
                outsheet using "`output_dir'/regoutputs_coeffs_gdppc_`degstr'`suffix'.csv", comma replace

                mat Vtemp = V[1..2,1..2]
                clear
                svmat Vtemp
                outsheet using "`output_dir'/regoutputs_vcv_gdppc_`degstr'`suffix'.csv", comma replace
            restore
        }
	}
}

// ---- Gemeinsame Regression für alle obs_mods, interagiert, nur gemeinsame Jahre
if $DO_JOINT {
    foreach degstr of local degstrs {
        import delimited "`data_dir'/reginputs_gdppc_`degstr'.csv", clear
        encode iso_id, gen(iso_id_)
        encode obs_mod, gen(obs_mod_)
        gen temp = temp_popweight
        gen precip = precip_popweight
        gen precip2 = precip_popweight_2

        // Get joint years
        levelsof obs_mod, local(mods)
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

        // Regression
        reg growthwdi c.(c.temp##c.temp precip precip2 _yi_* _y2_*)#i.obs_mod_ i.year#i.obs_mod_ i.iso_id_#i.obs_mod_, cluster(iso_id_)

        mat b = e(b)
        mat V = e(V)

        levelsof obs_mod, local(mods)
        foreach mod of local mods {
            quietly summarize obs_mod_ if obs_mod == "`mod'"
            local mod_n = r(min)

            local temp_name  "`mod_n'.obs_mod_#c.temp"
            local temp2_name "`mod_n'.obs_mod_#c.temp#c.temp"

            mat b_mod = b[1, "`temp_name'"], b[1, "`temp2_name'"]
            mat V_mod = (V["`temp_name'", "`temp_name'"], V["`temp_name'", "`temp2_name'"] \ ///
                         V["`temp2_name'", "`temp_name'"], V["`temp2_name'", "`temp2_name'"])

            local modsuffix "_`mod'_`yrlabel'"
            preserve
                clear
                svmat b_mod
                outsheet using "`output_dir'/regoutputs_coeffs_gdppc_`degstr'`modsuffix'_JOINT.csv", comma replace
            restore

            preserve
                clear
                svmat V_mod
                outsheet using "`output_dir'/regoutputs_vcv_gdppc_`degstr'`modsuffix'_JOINT.csv", comma replace
            restore
        }
    }
}