//Neil Davies 18/06/20
//This simulates some data and attempts to estimate non-transmitted allele effects. 

cd "/Users/ecnmd/Desktop/Research/Henry Dale"

//With a non-transmitted allele effect:

cap: rm "results/trios_power.dta"
forvalues j=10000(10000)100000{
	forvalues i=0.05(0.05)0.5{
		clear
		set obs `j'
		gen f_g_1=rbinomial(1,`i')
		gen f_g_2=rbinomial(1,`i')

		gen m_g_1=rbinomial(1,`i')
		gen m_g_2=rbinomial(1,`i')

		gen m_rand=rbinomial(1,`i')
		gen f_rand=rbinomial(1,`i')

		gen o_g_1=m_g_1*m_rand+m_g_2*(1-m_rand)
		gen o_g_2=f_g_1*f_rand+f_g_2*(1-f_rand)

		gen o_nt_g_1=m_g_1*(1-m_rand)+m_g_2*(m_rand)
		gen o_nt_g_2=f_g_1*(1-f_rand)+f_g_2*(f_rand)

		gen m_g_as=m_g_1+m_g_2
		gen f_g_as=f_g_1+f_g_2
		gen o_g_as=o_g_1+o_g_2

		gen p=o_g_as-.5*m_g_as-.5*f_g_as+rnormal()
		
		cap{
			reg p o_g_as m_g_as f_g_as,ro
			regsave o_g_as m_g_as f_g_as using "results/trios_power",detail(all) pval addvar(maf,`i')
			}

		reg p o_g_as m_g_as f_g_as,ro
		regsave o_g_as m_g_as f_g_as using "results/trios_power",detail(all) pval append  addvar(maf,`i')
		
		//reg p o_g_as o_nt_*,ro
		}
	}
	

use "results/trios_power",clear
drop if var=="m_g_as"
gen maf=coef[_n+2] if var=="o_g_as"
replace maf=coef[_n+1] if var=="f_g_as"

drop if var=="maf"

gen min_detect_fx=1.96*stderr

drop if inlist(N,10000,30000,50000,70000,90000)

separate min_detect_fx, by(N) veryshortlabel

twoway connected min_detect_fx?  maf if var=="f_g_as", legend(title("Number of trios", size(medium)) rows(1) ring(0) position(2)) ///
	xtitle("Minor allele frequncy", size(small)) ytitle("Minimum detectable effect size in SD at 80% power", size(small)) ///
	xlabel(,format(%9.1f) labsize(small)) ylabel(,format(%9.2f) labsize(small))

//Without a non-transmitted allele effect:
clear
set obs 100000

gen f_g_1=rbinomial(1,0.3)
gen f_g_2=rbinomial(1,0.3)

gen m_g_1=rbinomial(1,0.3)
gen m_g_2=rbinomial(1,0.3)

gen m_rand=rbinomial(1,0.5)
gen f_rand=rbinomial(1,0.5)

gen o_g_1=m_g_1*m_rand+m_g_2*(1-m_rand)
gen o_g_2=f_g_1*f_rand+f_g_2*(1-f_rand)

gen o_nt_g_1=m_g_1*(1-m_rand)+m_g_2*(m_rand)
gen o_nt_g_2=f_g_1*(1-f_rand)+f_g_2*(f_rand)

gen m_g_as=m_g_1+m_g_2
gen f_g_as=f_g_1+f_g_2
gen o_g_as=o_g_1+o_g_2

gen p=o_g_as+rnormal()

reg p o_g_as m_g_as f_g_as,ro

reg p o_g_as o_nt_*,ro


