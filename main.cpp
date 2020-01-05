// /usr/local/cuda-10.0/bin/nvcc -x cu -I. main.cu
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "vec3.h"
#include "utils.h"

using std::cout;
using std::endl;

#define checkCudaErrors(val) check_cuda( (val), #val, __FILE__, __LINE__ )

void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line){
    if (result) {
        std::cerr << "CUDA error = " <<static_cast<unsigned int>(result) << " at " << file << ":" << line << " '" <<func << "' \n";
        //Make sure we call CUDA Device Resetbefore exiting
        cudaDeviceReset();
        exit(99);
    }
}


__global__ void render(vec3 *fb, int max_x, int max_y){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if ((i >= max_x) || (j>=max_y))
        return;
    int pixel_index = j * max_x + i;
            
    fb[pixel_index] = vec3(float(i) / max_x,
                           float(j) / max_y,
                           0.2);

}

int main (){
    int nx = 200;
    int ny = 100;
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels * sizeof(vec3);

    //allocate FB
    vec3 *fb;
    checkCudaErrors(cudaMallocManaged((void **)&fb, fb_size));

    int tx = 8;
    int ty = 8;

    clock_t start, stop;
    start = clock();

    // Render
    dim3 blocks(nx/tx+1, ny/ty+1);
    dim3 threads(tx,ty);
    render <<<blocks, threads>>>(fb,nx,ny);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);

    checkCudaErrors(cudaFree(fb));
    return 0;
}
