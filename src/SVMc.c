#include "R.h"
#include "math.h"

void SVMc(double *G, double *X, double *lambda0, double *rho0, double *L0, int *n0, int *d0, double *w0, double *u0, double *r0, int *maxIter0, double *thol0){
    
    int n,d,maxIter;
    int iter,i,j,k;
    int col_j;
    double w_diff,w_all;
    double u_diff,u_all;
    double tmp;
    
    double lambda,rho,L,thol,gap,gap1,gap2,shrink;
    
    lambda = lambda0[0];
    rho = rho0[0];
    L = L0[0];
    n = n0[0];
    d = d0[0];
    maxIter = maxIter0[0];
    thol = thol0[0];
	
    iter = 0;
    gap = 1;
    gap1 = 1;
    gap2 = 1;
    
    double *w1 = (double*) malloc(d*sizeof(double));
    double *u1 = (double*) malloc(n*sizeof(double));
    double *r1 = (double*) malloc(n*sizeof(double));
    
    double *v = (double*) malloc(d*sizeof(double));
    double *z = (double*) malloc(n*sizeof(double));
    
    int *idx = (int*) malloc((d-1)*sizeof(int)); //sizes of active sets
	
    for(j=0;j<(d-1);j++){
        if(w0[j]!=0)
            idx[j] = 1;
        else
            idx[j] = 0;
    }
        
    while(iter<maxIter && gap>thol){
        
        for(i=0;i<n;i++){
            z[i] = r0[i] - u0[i] - 1;
        }
        
        //Gradient Descent
        for(j=0;j<d;j++){
            v[j] = 0;
            col_j = j*d;
            for(k=0;k<d;k++){
                if(idx[k]==1){
                    v[j] = v[j] + G[col_j+k]*w0[k];
                }
            }
            v[j] = v[j] + G[col_j+(d-1)]*w0[d-1];
           
            
            col_j = j*n;
            for(i=0;i<n;i++){
                v[j] = v[j] + X[col_j+i]*z[i];
            }
            w1[j] = w0[j] - v[j]/L;
        }
	
	        
        //Soft-thresholding
        
        w_diff = 0;
        w_all = 0;
        for(j=0;j<(d-1);j++){
            if(w1[j] >= lambda/rho){
                w1[j] = w1[j] - lambda/rho;
                w_diff = w_diff + fabs(w1[j]-w0[j]);
                w_all = w_all + w1[j];
                idx[j] = 1;
            }
                
            else if(w1[j] <= -lambda/rho){
                w1[j] = w1[j] + lambda/rho;
                w_diff = w_diff + fabs(w1[j]-w0[j]);
                w_all = w_all - w1[j];
                idx[j] = 1;
            }

            else{
                w1[j] = 0;
                w_diff = w_diff + fabs(w0[j]);
                idx[j] = 0;
            }
        }
        
        w_diff = w_diff + fabs(w1[d-1]-w0[d-1]);
        w_all = w_all + fabs(w1[d-1]);
	
        //printf("%f\n",w_all);
        //printf("%f\n",w_diff);
        
        u_diff = 0;
        u_all = 0;
        for(i=0;i<n;i++){
            r1[i] = u0[i] + 1;
            for(j=0;j<(d-1);j++){
                if(idx[j]==1){
                    r1[i] = r1[i] - X[j*n+i]*w1[j];
                }
            }
            r1[i] = r1[i] - X[(d-1)*n+i]*w1[d-1];
            
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
    //printf("%f\n",gap1);
    //printf("%f\n",gap2);
    free(w1);
    free(r1);
    free(u1);
    free(idx);
    free(v);
    free(z);
}
