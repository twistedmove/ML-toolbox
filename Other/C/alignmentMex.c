#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <stdio.h>

double alignment(double *a, double *b,int aLength,int bLength,int width);
double norm(double *a,double *b,int dim);
double min(double a, double b);
double max(double a, double b);

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* more C/C++ code ... */
    mxArray *a_in,*b_in,*c_out;
    double *a,*b,*c;
    const *dims,*dimsB;
    int i,j,res=0;
    int dimy, dimx,dimBy,dimBx;
    
    /* associate inputs*/
    
    a_in = mxDuplicateArray(prhs[0]);
    b_in = mxDuplicateArray(prhs[1]);
    
    
    dims = mxGetDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    
    dimsB = mxGetDimensions(prhs[1]);
    dimBy = (int)dimsB[0]; dimBx = (int)dimsB[1];
    
    c_out = plhs[0] = mxCreateDoubleScalar(0);
    
    /*associate pointers*/
    a = mxGetData(a_in);
    b = mxGetPr(b_in);
    c = mxGetPr(c_out);
    
      
    c[0]=alignment(a,b,dimx,dimBx,dimy);
    
    
}

double alignment(double *a, double *b,int aLength,int bLength,int width){
    mxArray *M_in,*null_in, *var_in,*var_in2;
    double *varA,*varB,*M,*null;
    int i,j,k,kk;
    
    M_in=mxCreateDoubleMatrix(aLength+1,bLength+1,mxREAL);
    M=mxGetPr(M_in);
    
    M[0]=0;
    
    null_in=mxCreateDoubleMatrix(width,1,mxREAL);
    null=mxGetPr(null_in);
    var_in=mxCreateDoubleMatrix(width,1,mxREAL);
    var_in2=mxCreateDoubleMatrix(width,1,mxREAL);
    varA=mxGetPr(var_in);
    varB=mxGetPr(var_in2);
    
    for (i=0;i<width;i++){
        null[i]=0;
    }
   
    
    
    for (j=1;j<=bLength;j++){
        
        for (k=0;k<width;k++){
            
            varA[k]=b[(j-1)*width+k];
            
        }
        
       
        M[j]=M[(j-1)]+norm(varA,null,width);
    }
    
    
    
    
    
    
    
    for (i=1;i<aLength+1;i++){
        
        for (k=0;k<width;k++){
            varA[k]=a[(i-1)*width+k];
            
        }
        
        M[i*(bLength+1)]=M[(i-1)*(bLength+1)]+norm(varA,null,width);
    }
    
     
    for (i=1;i<=aLength;i++){
        
        for (k=0;k<width;k++){
            varA[k]=a[(i-1)*width+k];
        }
        
        for (j=1;j<=bLength;j++){
            for (kk=0;kk<width;kk++){
                varB[kk]=b[(j-1)*width+kk];
            }
           
            M[i*(bLength+1)+j]=max(max(M[(i-1)*(bLength+1)+j-1]+
                    norm(varA,varB,width),M[(i-1)*(bLength+1)+j]+norm(varA,null,width)),
                    M[(i)*(bLength+1)+j-1]+norm(varB,null,width));
            
        }
    }
    
    
          /* print out the matrix
    for (i=0;i<aLength+1;i++){
        for (j=0;j<bLength+1;j++){
            printf("   %f  ",M[i*(bLength+1)+j]);
        }
        printf("\n");
    }
    */
   
    return M[(aLength)*(bLength+1)+bLength];
}

double norm(double *a,double *b,int dim){
    
    int i,j;
    double res=0;
    
    for(j=0;j<dim;j++){
        res+=pow(a[j]-b[j],2);
    }
    
    res=-sqrt(res);
    return (res);
    
    
}

double min(double a, double b){
    
    if (a>b){
        return b;
    }else{
        return a;
    }
    
}

double max(double a, double b){
    
    if (a<b){
        return b;
    }else{
        return a;
    }
    
}

