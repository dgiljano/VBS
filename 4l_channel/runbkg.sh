for chan in 2e2mu 4mu 4e
	do
	for cate in  1 #0 #2 
		do
		for highmass in 2 
			do
		    for year in 2016 2017 2018 
		           do
		  	      
				#bsub -C 0 -q cmscaf1nd job_bkg.lsf "$chan" $cate $highmass 0
				#bsub -C 0 -q cmscaf1nd job_bkg.lsf "$chan" $cate $highmass 0
				./job_bkg.lsf "$chan" $cate $highmass 1 $year
			done
			done
			done
                        done

