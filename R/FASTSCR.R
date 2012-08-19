#-----------------------------------------------------------------------#
# Package: SAM       	   		   									    #
# Method: Fast Screening Rule				      	                	#
# Authors: Tuo Zhao, Xingguo Li, Han Liu, Lie Wang and Kathryn Roeder   #
# Date: Aug 19th 2012                                                   #
# Version: 1.0.1                                                        #
#-----------------------------------------------------------------------#

fastscr = function(x,y,nscr = NULL, method = "t.test"){
	n = nrow(x)
	d = ncol(x)
	
	if(is.null(nscr))
		nscr = min(n,d)
	
	idx.a = which(y==1)
	idx.b = which(y==-1)


	out = rep(0,d)
	
	if(method == "t.test")
		for(j in 1:d){
			out[j] = t.test(x[idx.a,j],x[idx.b,j],alternative="two.sided")$p.value
		}
	if(method == "wilcox")
		for(j in 1:d){
			out[j] = wilcox.test(x[idx.a,j],x[idx.b,j],alternative="two.sided")$p.value
		}
	
	idx = order(out)[1:nscr]	
	return(idx)		
}
