#include "R.h"
#include "math.h"

void gSVMc(double *G, double *X, double *lambda0, double *rho0, double *L0, int *n0, int *d0, int *k0, int *g0, double *w0, double *u0, double *r0, double *gnorm, int *maxIter0, double *thol0){
    
    int n,d,g,k,maxIter;
    int iter,i,j,t,m;
    int col_j,group_m;
    double w_diff,w_all;
    double u_diff,u_all;
    
    double lambda,rho,L,thol,gap,gap1,gap2,shrink;
    
    lambda = lambda0[0];
    rho = rho0[0];
    L = L0[0];
    n = n0[0];
    d = d0[0];
    g = g0[0];
    k = k0[0];
    maxIter = maxIter0[0];
    thol = thol0[0];
	
    iter = 0;
    gap = 1;
    
    double *w1 = (double*) malloc(d*sizeof(double));
    double *u1 = (double*) malloc(n*sizeof(double));
    double *r1 = (double*) malloc(n*sizeof(double));
    double *norm2 = (double*) malloc(g*sizeof(double));
    
    double *v = (double*) malloc(d*sizeof(double));
    double *z = (double*) malloc(n*sizeof(double));
    
    int *idx = (int*) malloc(g*sizeof(int)); //sizes of active sets
	
    for(m=0;m<g;m++){
        if(w0[m*k]!=0)
            idx[m] = 1;
        else
            idx[m] = 0;
    }
    
    //for(m=0;m<g;m++) printf("%d\n",idx[j]);
    
    while(iter<maxIter && gap>thol){
        
        for(i=0;i<n;i++){
            z[i] = r0[i] - u0[i] - 1;
        }
	//for(i=0;i<n;i++) printf("%f\n",z[i]);
    
        //Gradient Descent
        //printf("%d\n",k);
	//printf("%d\n",g);
        for(j=0;j<d;j++){
            v[j] = 0;
            col_j = j*d;
            for(m=0;m<g;m++){
		//printf("idx_m:%d\n",idx[m]);
                if(idx[m]==1){
                    group_m = m*k;
		    //printf("%d\n",m);
                    for(t=0;t<k;t++){
                        v[j] = v[j] + G[col_j+group_m+t]*w0[group_m+t];
			//printf("%f\n",v[j]);
                    }
                }
            }
            v[j] = v[j] + G[col_j+(d-1)]*w0[d-1];
	    //printf("%f\n",v[j]);
        
            col_j = j*n;
            for(i=0;i<n;i++){
                v[j] = v[j] + X[col_j+i]*z[i];
            }
            w1[j] = w0[j] - v[j]/L;
        }
	
	//for(j=0;j<d;j++) printf("%f\n",v[j]);
        //for(j=0;j<d;j++) printf("%f\n",w1[j]);
        
	//Group Soft-thresholding
        
        w_diff = 0;
        w_all = 0;
        for(m=0;m<g;m++){
            group_m = m*k;
            gnorm[m] = 0;
            for(t=0;t<k;t++){
                gnorm[m] = gnorm[m] + pow(w1[group_m+t],2);
            }
            gnorm[m] = sqrt(gnorm[m]);
            
            shrink = (1-lambda/rho/gnorm[m]);
            
            if(shrink<=0){
                for(t=0;t<k;t++){
                    w1[group_m+t] = 0;
                    w_diff = w_diff + fabs(w0[group_m+t]);
                }
                gnorm[m] = 0;
                idx[m] = 0;
            }
            else{
                gnorm[m] = gnorm[m]*shrink;
                for(t=0;t<k;t++){
                    w1[group_m+t] = shrink*w1[group_m+t];
                    w_diff = w_diff + fabs(w1[group_m+t]-w0[group_m+t]);
                    w_all = w_all + fabs(w1[group_m+t]);
                }
                idx[m] = 1;
            }
	    
	    //printf("%f\n",gnorm[m]);
        }
	w_diff = w_diff + fabs(w1[d-1]-w0[d-1]);
	w_all = w_all + fabs(w1[d-1]);
	//for(j=0;j<d;j++) printf("%f\n",w1[j]);
        //printf("%f\n",w_diff);
	//printf("%f\n",w_all);
        u_diff = 0;
        u_all = 0;
        for(i=0;i<n;i++){
            r1[i] = u0[i] + 1;
            for(m=0;m<g;m++){
                if(idx[m]==1){
                    group_m = m*k;
                    for(t=0;t<k;t++){
                        r1[i] = r1[i] - X[(group_m+t)*n+i]*w1[group_m+t];
                    }
                }
            }
            r1[i] = r1[i] - X[(d-1)*n+i]*w1[d-1];
            //printf("%f\n",r1[i]);
            if(r1[i]>1/rho){
                r1[i] = r1[i] - 1/rho;
		u1[i] = 1/rho;
            }
            else if(r1[i]>0){
                u1[i] = r1[i];
                r1[i] = 0;
            }
            else{
                u1[i] = 0;
            }
	    u_diff = u_diff + pow((u1[i] - u0[i]),2);
	    u_all = u_all + pow(u1[i],2);
	    //printf("%f\n",r1[i]);
	    //printf("%f\n",u1[i]);
        }
	//printf("%f\n",u_all);
	//printf("%f\n",u_diff);
	
        
        for(j=0;j<d;j++){
            w0[j] = w1[j];
        }
        for(i=0;i<n;i++){
            r0[i] = r1[i];
            u0[i] = u1[i];
        }
        
	gap1 = sqrt(u_diff)/(sqrt(u_all)+(1e-8));
	gap2 = w_diff/(w_all+(1e-8));
	
	if(gap1>gap2) gap = gap1;
	else gap = gap2;
        
        iter = iter + 1;
    }
    //printf("%d\n",iter);
    free(w1);
    free(r1);
    free(u1);
    free(idx);
    free(v);
    free(z);
}
