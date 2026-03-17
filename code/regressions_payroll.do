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

global DO_MAIN 1
global DO_MAIN_JOINT 1 

local degstrs 025deg 1deg

if $DO_MAIN {
    foreach wght in "raw" "pop" {
        foreach degstr of local degstrs {
            import delimited "`data_dir'reginputs_payroll_`degstr'_`wght'wght.csv", clear

            local a_dry = " prec_dry below40_dry fnd_40to45_dry fnd_45to50_dry fnd_50to55_dry fnd_55to60_dry fnd_60to65_dry fnd_75to80_dry  fnd_65to70_dry "

            * Setup for dry bulb regression
            rename (days90plus_dry fnd_85to90_dry fnd_80to85_dry  ) (days90plus fnd_85to90 fnd_80to85 )

            * Get names of obs_mods
            levelsof obs_mod, local(mods)

            * Get joint years
            local total = r(r)
            bysort year (obs_mod): gen tag = (obs_mod != obs_mod[_n-1]) | _n == 1
            bysort year: egen n_present = total(tag)
        
            * Get label for filenames
            levelsof year, local(yrs)
            local first : word 1 of `yrs'
            local n = wordcount("`yrs'")
            local last : word `n' of `yrs'
            local yrlabel = string(`first') + "-" + string(`last')

            * Pre-compute mod_n for each mod (before any preserve/clear)
            * (get "clean" version of the names to avoid STATA syntax issues
            * from dashes)
            encode obs_mod, gen(obs_mod_)
            foreach mod of local mods {
                local mod_clean = subinstr("`mod'", "-", "_", .)
                quietly summarize obs_mod_ if obs_mod == "`mod'"
                local mod_n_`mod_clean' = r(min)
            }


            foreach yeartyp in allyrs jointyrs {

                * ---- Separate regressions per obs_mod
                foreach mod of local mods {
                    preserve
                        keep if obs_mod == "`mod'"

                        if "`yeartyp'" == "allyrs" {
                            local suffix "_`mod'_ALLYEARS"
                        }
                        if "`yeartyp'" == "jointyrs" {
                            keep if n_present == `total'
                            local suffix "_`mod'_`yrlabel'"
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
                        outsheet using "`output_dir'/regoutputs_coeffs_payroll_`degstr'_`wght'`suffix'.csv", comma replace

                        clear
                        svmat V, name(col)
                        outsheet using "`output_dir'/regoutputs_vcv_payroll_pooled_`degstr'_`wght'`suffix'.csv", comma replace

                    restore
                }

                * ---- Joint regression, only for jointyrs
                if "`yeartyp'" == "jointyrs" {
                    preserve
                        keep if n_present == `total'

                        reghdfe lny_nonag c.(days90plus fnd_85to90 fnd_80to85 `a_dry')#i.obs_mod_ ///
                                , absorb(i.fips#i.obs_mod_ i.year#i.obs_mod_ i.state#c.year#i.obs_mod_) cluster(i.fips)

                        mat b_JOINT = e(b)
                        local last : word `=colsof(b_JOINT)' of `:colnames b_JOINT'
                        if "`last'" == "_cons" mat b_JOINT = b_JOINT[1, 1..`=colsof(b_JOINT)-1']

                        mat V_JOINT = e(V)
                        local last : word `=colsof(V_JOINT)' of `:colnames V_JOINT'
                        if "`last'" == "_cons" {
                            local n = colsof(V_JOINT) - 1
                            mat V_JOINT = V_JOINT[1..`n', 1..`n']
                        }

                        * ---- Save per-mod coefficients and block-diagonal VCVs
                        foreach mod of local mods {
                            local mod_clean = subinstr("`mod'", "-", "_", .)
                            local mod_n = `mod_n_`mod_clean''
                            local modsuffix "_`mod'_`yrlabel'_JOINT"

                            * Build list of all column names for this mod
                            local modcols ""
                            local allcols : colnames b_JOINT
                            local ncols = colsof(b_JOINT)
                            forvalues k = 1/`ncols' {
                                local col1 : word `k' of `allcols'
                                if regexm("`col1'", "^`mod_n'[b]?\.obs\_mod\_\#") {
                                    local modcols "`modcols' `col1'"
                                }
                            }

                            * Extract b
                            local firstcol : word 1 of `modcols'
                            mat b_mod = b_JOINT[1, "`firstcol'"]
                            local nc = wordcount("`modcols'")
                            forvalues k = 2/`nc' {
                                local col1 : word `k' of `modcols'
                                mat b_mod = b_mod, b_JOINT[1, "`col1'"]
                            }

                            * Extract block-diagonal V
                            local nc = wordcount("`modcols'")
                            mat V_mod = J(`nc', `nc', 0)
                            forvalues i = 1/`nc' {
                                forvalues j = 1/`nc' {
                                    local ci : word `i' of `modcols'
                                    local cj : word `j' of `modcols'
                                    mat V_mod[`i', `j'] = V_JOINT["`ci'", "`cj'"]
                                }
                            }

                            local cols : colnames b_mod
                            local newcols ""
                            local ncols = colsof(b_mod)
                            forvalues k = 1/`ncols' {
                                local col : word `k' of `cols'
                                * Alles nach "#c." extrahieren
                                local col_clean = regexr("`col'", ".*#c\.", "")
                                local newcols "`newcols' `col_clean'"
                            }
                            mat colnames b_mod = `newcols'
                            mat colnames V_mod = `newcols'
                            mat rownames V_mod = `newcols'

                            clear
                            svmat b_mod, name(col)
                            outsheet using "`output_dir'/regoutputs_coeffs_payroll_`degstr'_`wght'`modsuffix'.csv", comma replace

                            clear
                            svmat V_mod, name(col)
                            outsheet using "`output_dir'/regoutputs_vcv_payroll_pooled_`degstr'_`wght'`modsuffix'.csv", comma replace
                        }

                        * ---- Save full joint VCV
                        * Have to change the naems of the rows / columns, since STATA can't put "." or "#" into CSV headers for some reason
                        local cols : colnames b_JOINT
                        local newcols ""
                        local ncols = colsof(b_JOINT)
                        forvalues k = 1/`ncols' {
                            local col : word `k' of `cols'

                            * Replace mod encoded number with mod name
                            foreach mod of local mods {
                                local mod_clean = subinstr("`mod'", "-", "_", .)
                                local mod_n = `mod_n_`mod_clean''
                                if regexm("`col'", "^`mod_n'[b]?\.obs\_mod\_\#c\.(.+)$") {
                                    local varname = regexs(1)
                                    local col = "`mod_clean'_`varname'"
                                }
                            }

                            local col_clean = subinstr("`col'", ".", "_", .)
                            local col_clean = subinstr("`col_clean'", "#", "_", .)
                            local col_clean = subinstr("`col_clean'", " ", "_", .)
                            local newcols "`newcols' `col_clean'"
                        }
                        mat colnames b_JOINT = `newcols'
                        mat rownames V_JOINT = `newcols'
                        mat colnames V_JOINT = `newcols'

                        clear
                        svmat V_JOINT, name(col)
                        outsheet using "`output_dir'/regoutputs_vcv_payroll_pooled_`degstr'_`wght'_ALLREA_`yrlabel'_JOINT.csv", comma replace

                        * ---- Save full joint b
                        clear
                        svmat b_JOINT, name(col)
                        outsheet using "`output_dir'/regoutputs_coeffs_payroll_`degstr'_`wght'_ALLREA_`yrlabel'_JOINT.csv", comma replace

                    restore
                }
            }
        }
    }
}