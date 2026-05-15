* GDPPC regressions, adapted from BHM15 replication code
//		Note: users will need to install parmest command (e.g. "ssc install parmest")

clear all
set more off

global DO_BYMOD 1 // 1 for run the regressions for each mod separately
global DO_BYMOD_COMMONYEARS 1 // 1 for run the regressions only in years with all mods
global DO_JOINT 1 // 1 for run the joint regression, only in years with all mods

local data_dir "/shared/share_hle/data/reg_inputs/gdppc/"
local output_dir "/shared/share_hle/data/reg_outputs/gdppc/"

local degstr 025deg


import delimited "`data_dir'/reginputs_gdppc_`degstr'.csv", clear
encode iso_id, gen(iso_id_)
gen temp = temp_popweight
gen precip = precip_popweight
gen precip2 = precip_popweight_2

* --- Get joint years
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


* Get label for countries
levelsof countryname, local(reg_ids)

foreach mod of local mods {
	foreach reg_ of local reg_ids { 
	    local suffix "_`mod'_`yrlabel'_no`reg_'"

	    preserve
	        keep if obs_mod == "`mod'"

	        * Drop if ADM0/ADM1 code matches current leave-one-out
	        drop if countryname == "`reg_'"

	        // Regression
	        reg growthwdi c.temp##c.temp precip precip2 i.year _yi_* _y2_* i.iso_id_, cluster(iso_id_)
	            mat b = e(b)
	            mat b = b[1,1..2]
	            mat V = e(V)

	        clear
	        svmat b
	        outsheet using "`output_dir'/loos/regoutputs_coeffs_gdppc_`degstr'`suffix'.csv", comma replace

	        mat Vtemp = V[1..2,1..2]
	        clear
	        svmat Vtemp
	        outsheet using "`output_dir'/loos/regoutputs_vcv_gdppc_`degstr'`suffix'.csv", comma replace
	    restore
	}
}