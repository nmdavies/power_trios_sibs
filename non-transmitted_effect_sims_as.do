//Neil Davies 18/06/20
//This simulates some data and attempts to estimate non-transmitted allele effects using allele scores

cd "/Users/ecnmd/Desktop/Research/Henry Dale"

//With a non-transmitted allele effect:

cap: rm "results/trios_power_as.dta"
forvalues j=20000(20000)100000{
	forvalues i=0.005(0.01)0.105{
		clear
		
		set obs `j'
		forvalues k=1(1)1000{
			gen f_g`k'_1=rbinomial(1,`i')
			gen f_g`k'_2=rbinomial(1,`i')

			gen m_g`k'_1=rbinomial(1,`i')
			gen m_g`k'_2=rbinomial(1,`i')

			gen m_rand_`k'=rbinomial(1,`i')
			gen f_rand_`k'=rbinomial(1,`i')
		
			gen o_g`k'_1=m_g`k'_1*m_rand_`k'+m_g`k'_2*(1-m_rand_`k')
			gen o_g`k'_2=f_g`k'_1*f_rand_`k'+f_g`k'_2*(1-f_rand_`k')

			gen o_nt_g`k'_1=m_g`k'_1*(1-m_rand_`k')+m_g`k'_2*(m_rand_`k')
			gen o_nt_g`k'_2=f_g`k'_1*(1-f_rand_`k')+f_g`k'_2*(f_rand_`k')
			}
		egen m_g_as=rowtotal(m_g*_1 m_g*_2)
		egen f_g_as=rowtotal(f_g*_1 f_g*_2)
		egen o_g_as=rowtotal(o_g*_1 o_g*_2)
			
		
		gen p=`i'*(o_g_as-.5*m_g_as-.5*f_g_as)+(1-`i')*rnormal()
		
		cap{
			reg p o_g_as m_g_as f_g_as,ro
			regsave o_g_as m_g_as f_g_as using "results/trios_power_as",detail(all) pval addvar(maf,`i')
			}

		reg p o_g_as m_g_as f_g_as,ro
		regsave o_g_as m_g_as f_g_as using "results/trios_power_as",detail(all) pval append  addvar(maf,`i')
		
		//reg p o_g_as o_nt_*,ro
		}
	}

	

use "results/trios_power_as",clear
drop if var=="m_g_as"
gen maf=coef[_n+2] if var=="o_g_as"
replace maf=coef[_n+1] if var=="f_g_as"

drop if var=="maf"

gen min_detect_fx=1.96*stderr

drop if inlist(N,10000,30000,50000,70000,90000)

separate min_detect_fx, by(N) veryshortlabel

twoway connected min_detect_fx?  r2 if var=="f_g_as" &r2<0.2, legend(title("Number of trios", size(medium)) rows(1) ring(0) position(2)) ///
	xtitle("R-squared for allele scores", size(small)) ytitle("Minimum detectable effect size in SD at 80% power", size(small)) xlabel(,format(%9.2f) labsize(small)) ///
	ylabel(,format(%9.3f) labsize(small))
