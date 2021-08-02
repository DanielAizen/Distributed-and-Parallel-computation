#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"

__device__ int IsConservativeGPU(char c1, char c2);

__device__ char IsSemiConservativeGPU(char s);

__device__ int MyStrlen(const char *str);

__device__ int MyStrchr(const char *s, char c);

__global__ void CreateMutationGPU(char *d_Seq1, char *d_Mutant_Seq2, double *d_Weights, int h_Offset, int h_Sort_Order, int h_Size2)
{ 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < h_Size2)
    {
        int flagCons = 0;
        char cSC = 0; //char semi conservative
        if (h_Sort_Order == 1)
        {
            flagCons = IsConservativeGPU(d_Seq1[i + h_Offset], d_Mutant_Seq2[i]);
            if (flagCons == 0)
            {
                d_Mutant_Seq2[i] = d_Seq1[i + h_Offset];
            }
        }
        else
        {
            flagCons = IsConservativeGPU(d_Seq1[i + h_Offset], d_Mutant_Seq2[i]);
            if (flagCons == 0)
            {
                if (d_Weights[2] > d_Weights[3])
                {
                    cSC = IsSemiConservativeGPU(d_Seq1[i + h_Offset]);
                    
                    if (cSC != '$')
                    {
                        d_Mutant_Seq2[i] = cSC;
                    }
                }
                else{
                    d_Mutant_Seq2[i] = 'Z';
                }
            }
        }
        
    }
}
__device__ int IsConservativeGPU(char c1, char c2)
{
    const char *conservative_Group[CONSERVATIVE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
    int cmp1, cmp2;
    for (int j = 0; j < CONSERVATIVE; j++)
    {
        cmp1 = MyStrchr(conservative_Group[j],c1);
        cmp2 = MyStrchr(conservative_Group[j],c2);
        if(cmp1 == 0 && cmp2 == 0){
            return 1;
        }
    }
    return 0;
}
__device__ char IsSemiConservativeGPU(char s)
{
    const char *semi_Conservative_Group[SEMI_CONSERVATIVE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
    int cmp;
    for (int j = 0; j < SEMI_CONSERVATIVE; j++)
    {
        cmp = MyStrchr(semi_Conservative_Group[j], s);
        if(cmp == 0){
            for(int k =0; k < MyStrlen(semi_Conservative_Group[j]); k++){
                if(semi_Conservative_Group[j][k] != s && IsConservativeGPU(s , semi_Conservative_Group[j][k]) == 0)
					return semi_Conservative_Group[j][k];
            }
        }
    }
    return '$';
}
__device__ int MyStrlen(const char *str)
{
    int len = 0;
    while (*(str++))
        (len)++;
    return len;
}
__device__ int MyStrchr(const char *s, char c) {
    while (*s != c) {
        if (!*s++) {
            return 1;
        }
    }
    return 0;
}

int computeOnGPU(char *h_Seq1, char *h_Mutant_Seq2, int h_Offset, int h_Sort_Order, double *h_Weights, int h_Len_Seq1, int h_Len_Seq2)
{
    cudaError_t err = cudaSuccess;
    size_t size1 = h_Len_Seq1 * sizeof(char);
    size_t size2 = h_Len_Seq2 * sizeof(char);

    char *d_Seq1;
    char *d_Mutant_Seq2;
    err = cudaMalloc((void **)&d_Seq1, size1);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory for seq1 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&d_Mutant_Seq2, size2);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory for seq2 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    double *d_Weights;
    size_t sizeW = 4 * sizeof(double);
    err = cudaMalloc((void **)&d_Weights, sizeW);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory for seq2 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_Seq1, h_Seq1, size1, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy data from host to device for seq1 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(d_Mutant_Seq2, h_Mutant_Seq2, size2, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy data from host to device for seq2 - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMemcpy(d_Weights, h_Weights, sizeW, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy data from host to device for weights - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
   
    // Launch the Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (h_Len_Seq2 + threadsPerBlock - 1) / threadsPerBlock;
    CreateMutationGPU<<<blocksPerGrid, threadsPerBlock>>>(d_Seq1, d_Mutant_Seq2, d_Weights, h_Offset, h_Sort_Order, h_Len_Seq2);
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch CreateMutationGPU kernel -  %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
	
    // Copy the result from GPU to the host memory.
    err = cudaMemcpy(h_Mutant_Seq2, d_Mutant_Seq2, size2, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    // Free allocated memory on GPU
    if (cudaFree(d_Seq1) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    if (cudaFree(d_Mutant_Seq2) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    if (cudaFree(d_Weights) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    return 0;
}
