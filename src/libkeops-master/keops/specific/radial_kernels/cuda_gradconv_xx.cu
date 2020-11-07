/*
*	This cuda routine allows one to compute the derivative wrt the point cloud 'x' of the derivative
*	wrt 'x' of the expression
*		K(x_i,y_j) @ b_j =  sum_j f( |x_i-y_j|^2 ) b_j
*
*
*	We're looking for the gradient with respect to x of
*
*	< e, K(s,a,x,y,b) >  =  \sum_{i,j} f_s'( |x_i-y_j|^2 ) * < a_i, b_j > * 2 < e_i, x_i-y_j>,
*
*	which is an N-by-D array g_i (i from 1 to N), where each line is equal to
*
*	g_i  =  2* \sum_j < a_i, b_j > * [                       f_s'(  |x_i-y_j|^2 ) * e_i
*                                    + 2* < x_i-y_j, e_i > * f_s''( |x_i-y_j|^2 ) * (x_i-y_j) ]
*
*	We will compute this sum over the index 'j' on the GPU, with 'one thread' = 'one index i'.
*	Data will be stored as follow:
*	  - e_i in the thread memory
* 	  - a_i in the thread memory
*	  - x_i in the thread memory
*	  - y_j in the SharedData
*	  - b_j in the SharedData (beta_j, really)
*
*
* Author : Jean Feydy, heavily based on the work of Joan Glaunès and Benjamin Charlier.
*
*/

#include <stdio.h>
#include <assert.h>
#include <cuda.h>

#include "specific/radial_kernels/radial_kernels.h"
#include "specific/radial_kernels/cuda_gradconv_xx.cx"


//////////////////////////////////////////////////////
/////////// CPU -> GPU -> CPU routines ///////////////
//////////////////////////////////////////////////////


template < typename TYPE, KernelFun KernelFp , KernelFun KernelFpp >
int KernelGpuGradConvXX(TYPE ooSigma2,               // 1 / sigma^2
                        TYPE* e_h,                     // N-by-D array (same as x)
                        TYPE* alpha_h, TYPE* x_h,     // N-by-E, N-by-D arrays
                        TYPE* y_h,     TYPE* beta_h,  // M-by-D, M-by-E arrays
                        TYPE* gamma_h,                 // Output: N-by-D (same as x)
                        int dimPoint, int dimVect, int nx, int ny) { // D, E, N, M

    // Data on the device.
    TYPE* e_d;
    TYPE* alpha_d;
    TYPE* x_d;
    TYPE* y_d;
    TYPE* beta_d;
    TYPE* gamma_d;

    // Allocate arrays on device.
    cudaMalloc((void**)&e_d,     sizeof(TYPE)*(nx*dimPoint));
    cudaMalloc((void**)&alpha_d, sizeof(TYPE)*(nx*dimVect ));
    cudaMalloc((void**)&x_d,     sizeof(TYPE)*(nx*dimPoint));
    cudaMalloc((void**)&y_d,     sizeof(TYPE)*(ny*dimPoint));
    cudaMalloc((void**)&beta_d,  sizeof(TYPE)*(ny*dimVect ));
    cudaMalloc((void**)&gamma_d, sizeof(TYPE)*(nx*dimPoint)); // Output: N-by-D (same as x)

    // Send data from host to device.
    cudaMemcpy(e_d,     e_h,     sizeof(TYPE)*(nx*dimPoint), cudaMemcpyHostToDevice);
    cudaMemcpy(alpha_d, alpha_h, sizeof(TYPE)*(nx*dimVect ), cudaMemcpyHostToDevice);
    cudaMemcpy(x_d,     x_h,     sizeof(TYPE)*(nx*dimPoint), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d,     y_h,     sizeof(TYPE)*(ny*dimPoint), cudaMemcpyHostToDevice);
    cudaMemcpy(beta_d,  beta_h,  sizeof(TYPE)*(ny*dimVect ), cudaMemcpyHostToDevice);

    // compute on device.
    dim3 blockSize;
    blockSize.x = CUDA_BLOCK_SIZE; // number of threads in each block
    dim3 gridSize;
    gridSize.x =  nx / blockSize.x + (nx%blockSize.x==0 ? 0 : 1);

    // Copy-paste templating, allowing us to pass the DIMPOINT and DIMVECT at compilation time :
    if(     dimPoint==1 && dimVect==1)
        KernelGpuGradConvXXOnDevice<TYPE,1,1,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else if(dimPoint==2 && dimVect==1)
        KernelGpuGradConvXXOnDevice<TYPE,2,1,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else if(dimPoint==3 && dimVect==1)
        KernelGpuGradConvXXOnDevice<TYPE,3,1,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else if(dimPoint==4 && dimVect==1)
        KernelGpuGradConvXXOnDevice<TYPE,4,1,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else if(dimPoint==2 && dimVect==2)
        KernelGpuGradConvXXOnDevice<TYPE,2,2,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else if(dimPoint==3 && dimVect==3)
        KernelGpuGradConvXXOnDevice<TYPE,3,3,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else if(dimPoint==4 && dimVect==4)
        KernelGpuGradConvXXOnDevice<TYPE,4,4,KernelFp,KernelFpp><<<gridSize,blockSize,blockSize.x*(dimPoint+dimVect)*sizeof(TYPE)>>>
        (ooSigma2, e_d, alpha_d, x_d, y_d, beta_d, gamma_d, nx, ny);
    else {
        printf("GaussGpuGradConvXX error: dimensions of Gauss kernel not implemented in cuda\nYou probably just need a copy-paste in the conda_gradconv_xx.cu file !");
        cudaFree(e_d);
        cudaFree(alpha_d);
        cudaFree(x_d);
        cudaFree(y_d);
        cudaFree(beta_d);
        cudaFree(gamma_d);
        return(-1);
    }

    // block until the device has completed
    cudaThreadSynchronize();

    // Send data from device to host.
    cudaMemcpy(gamma_h, gamma_d, sizeof(TYPE)*(nx*dimPoint),cudaMemcpyDeviceToHost); // Output: N-by-D (same as x)

    // Free memory.
    cudaFree(e_d);
    cudaFree(alpha_d);
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(beta_d);
    cudaFree(gamma_d);

    return 0;
}


// Couldn't find a clean way to give a name to an explicit instantiation :-(
extern "C" int GaussGpuGradConvXX(__TYPE__ ooSigma2, __TYPE__* e_h, __TYPE__* alpha_h, __TYPE__* x_h, __TYPE__* y_h, __TYPE__* beta_h, __TYPE__* gamma_h, int dimPoint, int dimVect, int nx, int ny) {
    return KernelGpuGradConvXX<__TYPE__,GaussFp,GaussFpp>(ooSigma2, e_h, alpha_h, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}
extern "C" int CauchyGpuGradConvXX(__TYPE__ ooSigma2, __TYPE__* e_h, __TYPE__* alpha_h, __TYPE__* x_h, __TYPE__* y_h, __TYPE__* beta_h, __TYPE__* gamma_h, int dimPoint, int dimVect, int nx, int ny) {
    return KernelGpuGradConvXX<__TYPE__,CauchyFp,CauchyFpp>(ooSigma2, e_h, alpha_h, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}
extern "C" int LaplaceGpuGradConvXX(__TYPE__ ooSigma2, __TYPE__* e_h, __TYPE__* alpha_h, __TYPE__* x_h, __TYPE__* y_h, __TYPE__* beta_h, __TYPE__* gamma_h, int dimPoint, int dimVect, int nx, int ny) {
    return KernelGpuGradConvXX<__TYPE__,LaplaceFp,LaplaceFpp>(ooSigma2, e_h, alpha_h, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}
extern "C" int InverseMultiquadricGpuGradConvXX(__TYPE__ ooSigma2, __TYPE__* e_h, __TYPE__* alpha_h, __TYPE__* x_h, __TYPE__* y_h, __TYPE__* beta_h, __TYPE__* gamma_h, int dimPoint, int dimVect, int nx, int ny) {
    return KernelGpuGradConvXX<__TYPE__,InverseMultiquadricFp,InverseMultiquadricFpp>(ooSigma2, e_h, alpha_h, x_h, y_h, beta_h, gamma_h, dimPoint, dimVect, nx, ny);
}

void ExitFcn(void) {
    cudaDeviceReset();
}
