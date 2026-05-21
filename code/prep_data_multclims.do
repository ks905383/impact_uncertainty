/*


Modified version of prep_data.do to deal with multiple climate data sources. 

Intended to be run within the new regression code (so, without loading its own
dataset)

tavg_poly_1* instead of tavg_poly_GMFD_1*

*/


* a. create adm0 LR averages (one year one vote given imbalanced sample)
preserve 
    bysort iso year: keep if _n == 1
    collapse (mean) loggdppc_adm0_avg = loggdppc_adm0, by(iso) 
    tempfile avgadm0inc
    save `avgadm0inc', replace
restore 

* merge in sample average incomes
merge m:1 iso using "`avgadm0inc'", assert(3) nogen 


drop if year > 2010
drop if agegroup == "0"

* 1. create winsorized deathrate with-in country-agegroup
bysort iso agegroup: egen deathrate_p99 = pctile(deathrate), p(99)
gen deathrate_w99 = deathrate
replace deathrate_w99 = deathrate_p99 if deathrate > deathrate_p99 & !mi(deathrate)
drop deathrate_p99

* 2. set up sample 
gen sample = 0
replace sample = 1 if agegroup != "0" & year < = 2010
replace sample = 0 if mi(deathrate_w99)
replace sample = 0 if mi(tavg_poly_1)
replace sample = 0 if mi(prcp_poly_1)
replace sample = 0 if mi(loggdppc_adm1_avg)
replace sample = 0 if mi(lr_tavg_adm1_avg)

keep if sample == 1

* 3. clean up ids
egen adm0_code 			= group(iso)
egen adm1_code 			= group(iso adm1_id)
replace adm2_id 		= adm1_id if iso == "JPN"
egen adm2_code 			= group(iso adm1_id adm2_id)

egen adm0_agegrp_code 	= group(iso agegroup)
egen adm1_agegrp_code	= group(iso adm1_id agegroup)

* 4. assign EU countries same adm0_code
gen adm0_code2 = 41
replace adm0_code2 = 4 if adm0_code == 4 //BRA
replace adm0_code2 = 6 if adm0_code == 6 //CHL
replace adm0_code2 = 7 if adm0_code == 7 //CHN
replace adm0_code2 = 15 if adm0_code == 15 //FRA
replace adm0_code2 = 23 if adm0_code == 23 //JPN
replace adm0_code2 = 27 if adm0_code == 27 //MEX
replace adm0_code2 = 40 if adm0_code == 40 //USA

gen iso2 = iso
replace iso2 = "EUR" if adm0_code2==41

* Assign separate FEs to each time series segment for CHN DSPs that appear in both.
gen CHN_ts = 1
replace CHN_ts = 2 if year >= 2004 & iso=="CHN"
