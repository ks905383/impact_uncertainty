* GDPPC regressions, adapted from BHM15 replication code
//      Note: users will need to install parmest command (e.g. "ssc install parmest")

clear all
set more off

global DO_CROSS_PRS 1 // 1 for running the precipitation cross-checks

local data_dir "/shared/share_hle/data/reg_inputs/gdppc/"
local output_dir "/shared/share_hle/data/reg_outputs/payroll/crossprs/"

local degstr 025deg


// ---- Separate Regressionen für jedes obs_mod, nur gemeinsame Jahre
if $DO_CROSS_PRS {
    import delimited "`data_dir'/reginputs_gdppc_`degstr'.csv", clear
    encode iso_id, gen(iso_id_)
    gen temp = temp_popweight
    gen precip = precip_popweight
    gen precip2 = precip_popweight_2

    * --- Get joint years
    * Turn ERA5-025 into ERA5_025 usw. 
    * (needed for long-to-wide)
    replace obs_mod = subinstr(obs_mod, "-", "_", .)
    * Get names of obs_mods
    levelsof obs_mod, local(mods)
    local total = r(r)
    bysort year (obs_mod): gen tag = (obs_mod != obs_mod[_n-1]) | _n == 1
	bysort year: egen n_present = total(tag)
    keep if n_present == `total'
    
    * Get label for filenames
    levelsof year if n_present == `total', local(jointyrs)
    local first : word 1 of `jointyrs'
    local n = wordcount("`jointyrs'")
    local last : word `n' of `jointyrs'
    local yrlabel = string(`first') + "-" + string(`last')

    preserve
        keep countryname year obs_mod precip precip2
        rename precip* precip*_
        reshape wide precip*_, i(countryname year) j(obs_mod) string
        tempfile pr_wide
        save `pr_wide'
    restore
    
    merge m:1 countryname year using `pr_wide', nogen

    foreach tas_mod of local mods {
        local tas_mod_str : subinstr local tas_mod "_" "-", all
        foreach pr_mod of local mods {
            local pr_mod_str : subinstr local pr_mod "_" "-", all

            local suffix "_tas-`tas_mod_str'_pr-`pr_mod_str'_`yrlabel'"

            preserve
                keep if obs_mod == "`tas_mod'"

                gen precip_alt = precip_`pr_mod'
                gen precip2_alt = precip2_`pr_mod'

                // Regression
                reg growthwdi c.temp##c.temp precip_alt precip2_alt i.year _yi_* _y2_* i.iso_id_, cluster(iso_id_)
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
