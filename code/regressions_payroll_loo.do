clear
clear mata
clear matrix
set more off
cap log off
capture log close
pause on

cap ssc install outreg2
cap ssc install estout
cap ssc install reghdfe
cap ssc install ftools
cap ssc install maptile
cap ssc install spmap
cap net install dm79, from(http://www.stata.com/stb/stb56)
cap net install grc1leg, from(http://www.stata.com/users/vwiggins)
cap ssc install regexm

local code_dir "~/Documents/projects/hle_iv/code/"
local data_dir "/shared/share_hle/data/reg_inputs/payroll/"
local output_dir "/shared/share_hle/data/reg_outputs/payroll/"


local degstr 025deg
local wght pop
local yeartyp jointyrs

import delimited "`data_dir'reginputs_payroll_`degstr'_`wght'wght.csv", clear

local a_dry = " prec_dry below40_dry fnd_40to45_dry fnd_45to50_dry fnd_50to55_dry fnd_55to60_dry fnd_60to65_dry fnd_75to80_dry  fnd_65to70_dry "

* Setup for dry bulb regression
rename (days90plus_dry fnd_85to90_dry fnd_80to85_dry  ) (days90plus fnd_85to90 fnd_80to85 )

* Get names of obs_mods
levelsof obs_mod, local(mods)

* --- Get joint years
local total = r(r)
bysort year (obs_mod): gen tag = (obs_mod != obs_mod[_n-1]) | _n == 1
bysort year: egen n_present = total(tag)

* Get label for filenames
levelsof year if n_present == `total', local(jointyrs)
local first : word 1 of `jointyrs'
local n = wordcount("`jointyrs'")
local last : word `n' of `jointyrs'
local yrlabel = string(`first') + "-" + string(`last')

* --- Pre-compute mod_n for each mod (before any preserve/clear)
* (get "clean" version of the names to avoid STATA syntax issues
* from dashes)
encode obs_mod, gen(obs_mod_)
foreach mod of local mods {
    local mod_clean = subinstr("`mod'", "-", "_", .)
    quietly summarize obs_mod_ if obs_mod == "`mod'"
    local mod_n_`mod_clean' = r(min)
}

* Get label for FIPS
gen fips_int = int(fips)
levelsof fips_int, local(reg_ids)

* ---- Separate regressions per obs_mod
foreach mod of local mods {
	foreach reg_ of local reg_ids { 

	    preserve
		keep if obs_mod == "`mod'"

		* Drop if FIPS code matches current leave-one-out
		drop if fips_int == int(`reg_')

		if "`yeartyp'" == "allyrs" {
		    local suffix "_`mod'_ALLYEARS_no`reg_'"
		}
		if "`yeartyp'" == "jointyrs" {
		    keep if n_present == `total'
		    local suffix "_`mod'_`yrlabel'_no`reg_'"
		}

		reghdfe lny_nonag days90plus fnd_85to90 fnd_80to85 `a_dry', absorb(i.fips i.year i.state#c.year) cluster(i.fips)

		mat b = e(b)
		* Remove _cons, which is a protected variable name and can't be a matrix column header
		local last : word `=colsof(b)' of `:colnames b'
		if "`last'" == "_cons" mat b = b[1, 1..`=colsof(b)-1']

		mat V = e(V)
		* Remove _cons, which is a protected variable name and can't be a matrix column header
		local last : word `=colsof(V)' of `:colnames V'
		if "`last'" == "_cons" {
		    local n = colsof(V) - 1
		    mat V = V[1..`n', 1..`n']
		}

		clear
		svmat b, name(col)
		outsheet using "`output_dir'/loos/regoutputs_coeffs_payroll_`degstr'_`wght'`suffix'.csv", comma replace

		clear
		svmat V, name(col)
		outsheet using "`output_dir'/loos/regoutputs_vcv_payroll_pooled_`degstr'_`wght'`suffix'.csv", comma replace

	    restore
	}
}
