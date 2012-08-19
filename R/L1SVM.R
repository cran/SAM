#-----------------------------------------------------------------------#
# Package: SAM												      	   	#
# Method: L1 Norm SVM										            #
# Authors: Tuo Zhao, Xinguo Li, Han Liu, Lie Wang and Kathryn Roeder    #
# Date: Aug 19th 2012                                                   #
# Version: 1.0.1                                                        #
#-----------------------------------------------------------------------#

l1svm = function(x,y,lambda = NULL, rho = 1, thol=1e-4, maxIter = 1e4, rescale = T){	
	gcinfo(FALSE)
	fit = list()
	fit$rescale = rescale
	
	x = as.matrix(x)
	y = as.vector(y)
	
	n = nrow(x)
	d = ncol(x)
		
	if(rescale){	
		x.min = apply(x,2,min)
		x.max = apply(x,2,max)
		x.ran = x.max - x.min
		x.min.rep = matrix(rep(x.min,n),nrow=n,byrow=T)
		x.ran.rep = matrix(rep(x.ran,n),nrow=n,byrow=T)
		x = (x-x.min.rep)/x.ran.rep
	}
	fit$x.min = x.min
	fit$x.ran = x.ran	
	
	if(is.null(lambda))
		lambda = exp(seq(0,log(0.2),length=10))/n/log(n)/log(d)
	
	nlambda = length(lambda)
	
	d = d + 1
	
	X = cbind(diag(y)%*%x,y)
	e = rep(1,n)
	G = t(X)%*%X+1/n*diag(d)
	if(n<d)
		L = max(eigen(G,only.values=T)$values)
	if(n>=d)
		L = max(eigen(X%*%t(X),only.values=T)$values)
	
	fit$lambda = lambda
	fit$w = matrix(0,d-1,nlambda)
	fit$b = rep(0,nlambda)
	fit$df = rep(0,nlambda)
	fit$lab = matrix(0,n,nlambda)
	fit$dec = matrix(0,n,nlambda)
	w0 = rep(0,d)
	u0 = rep(0,n)
	r0 = rep(0,n)
	
	for(counter in 1:nlambda){
		ilambda = lambda[counter]
		
		out=.C("SVMc", G=as.double(G),X=as.double(X),lambda0=as.double(ilambda), rho0 = as.double(rho), L0 = as.double(L),n0 = as.integer(n), d0 = as.integer(d), w0 = as.double(w0), u0 = as.double(u0), r0 = as.double(r0), maxIter0 = as.integer(maxIter),thol0 = as.double(thol),PACKAGE="SAM")
		
		w0 = out$w0
		u0 = out$u0
		r0 = out$r0
		
		fit$w[,counter] = out$w0[1:(d-1)]
		fit$b[counter] = out$w0[d]
		fit$df[counter] = sum(out$w0!=0)
		fit$dec[,counter] = cbind(x,rep(1,n))%*%out$w0
		fit$lab[,counter] = sign(fit$dec[,counter])
	}
	
	rm(x,y,X,e,G,w0,u0,r0)
	gc()
	class(fit) = "l1svm"
	return(fit)		
}

print.l1svm = function(x,...){
	cat("Path length:",length(x$df),"\n")
	cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
	if(x$rescale) cat("rescaling: TRUE \n")
	if(!x$rescale) cat("rescaling: FALSE \n")
}

plot.l1svm = function(x,...){
	par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
	plot(1:ncol(x$w),x$w[1,],ylim=range(x$w),col=1,type="l",xlab="Regularization Parameters",ylab = "Coefficient",cex.lab=2)
	for(i in 2:nrow(x$w)){	
		lines(1:ncol(x$w),x$w[i,],col=(i%%256))
	}
}

predict.l1svm = function(object, newdata,...){
	out = list()
	nt = nrow(newdata)
	d = ncol(newdata)
	if(object$rescale){
		x.min.rep = matrix(rep(object$x.min,nt),nrow=nt,byrow=T)
		x.ran.rep = matrix(rep(object$x.ran,nt),nrow=nt,byrow=T)
		
		newdata = (newdata-x.min.rep)/x.ran.rep
		newdata = pmax(newdata,0)
		newdata = pmin(newdata,1)
	}
	
	out$dec = cbind(newdata,rep(1,nt))%*%rbind(object$w,object$b)
	out$lab = sign(out$dec) 
	return(out)
}

