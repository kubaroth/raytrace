// /usr/local/cuda-10.0/bin/nvcc -x cu -I. main.cpp
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>
#include <memory>

#include "vec3.h"
#include "utils.h"

#include <float.h>  // FLT_MAX

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

__device__ vec3 color(const ray& r, hitable **world, curandState *local_rand_state) {
    ray cur_ray = r;
    vec3 cur_attenuation = vec3(1.0,1.0,1.0);
    for(int i = 0; i < 50; i++) {
        hit_record rec;

        if ((*world)->hit(cur_ray, 0.001, MAXFLOAT, rec)) {

            ray scattered;
            vec3 attenuation;
            if (rec.mat_ptr->scatter(r, rec, attenuation,scattered, local_rand_state)) {
                cur_attenuation = cur_attenuation * attenuation;
                cur_ray = scattered;
            }
            else {
                return vec3(0.0,0.0,0.0);
            }

        }
        else {
            vec3 unit_direction = unit_vector(cur_ray.direction());
            float t = 0.5f*(unit_direction.y() + 1.0f);
            vec3 c = (1.0f-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
            return cur_attenuation * c;
        }
    }
    return vec3(0.0,0.0,0.0); // exceeded recursion
}


__global__ void render_init(int max_x, int max_y, curandState *rand_state){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    // printf("Dim: %i %i %i \n", blockDim.x, blockDim.y, blockDim.z);
    if ((i >= max_x) || ( j >=max_y))
        return;
    int pixel_index = j * max_x + i;
    //Each thread gets same seed, a different sequnce number, no offset
    //                vvvvvvvvvvv   regresion on higher resoution 1000
    // curand_init(0, pixel_index, 0 , &rand_state[pixel_index]);  //
    curand_init(pixel_index, 0, 0 , &rand_state[pixel_index]);  // this works as expected
}

__global__ void render(vec3 *fb, int max_x, int max_y, int ns,
                       camera **cam, hitable **world, curandState *rand_state
    ){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if ((i >= max_x) || (j>=max_y))
        return;

    int pixel_index = j * max_x + i;

    curandState local_rand_state = rand_state[pixel_index];
    vec3 col(0,0,0);
    for (int s=0; s<ns; s++) {  // sampling

        float u = float(i + curand_uniform(&local_rand_state)) / float(max_x);
        float v = float(j + curand_uniform(&local_rand_state)) / float(max_y);
        ray r = (*cam)->get_ray(u,v, &local_rand_state);
        col += color(r, world, &local_rand_state);
    }

    fb[pixel_index] = col/float(ns);
}


__global__ void create_world(hitable** d_list,
                  hitable **d_world,  /* ** - set refernce pointer to d_world */
                  camera **d_camera,
                  int nx, int ny){
    if (threadIdx.x == 0 && blockIdx.x ==0) {

        // ground
        int obj_index = 0;
        d_list[obj_index++] = new sphere(vec3(0,-100.5, -1), 100, new lambertian(vec3(0.8,0.8,0.8)));
        d_list[obj_index++] = new sphere(vec3(1,0, -1), 0.5, new lambertian(vec3(0.8,0.3,0.3)));
        d_list[obj_index++] = new sphere(vec3(0,0, -1), 0.5, new dielectric(1.5));
        d_list[obj_index++] = new sphere(vec3(-1,0, -1), 0.5, new metal(vec3(0.5,0.5,0.5), 0.8));

        vec3 look_from(3,3,2);
        vec3 look_at(0,0,-1);
        *d_camera = new camera(look_from,
                               look_at,
                               vec3(0,1,0),
                               20, float(nx)/float(ny),
                               0.5, /* 2.0 big aperture */
                               (look_from-look_at).length()/* distance to focus */
            );
        *d_world = new hitable_list(d_list,4);
    }

}

__global__ void free_world(hitable **d_list, hitable **d_world, camera **d_camera){
    delete *(d_list);
    delete *(d_list+1);
    delete *(d_list+2);
    delete *(d_list+3);
    delete *d_world;
    delete *d_camera;

}

int main (){
    int ns = 10; // number of samples
    int nx = 600;
    int ny = 300;
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels * sizeof(vec3);

    //allocate FB
    vec3 *fb;
    checkCudaErrors(cudaMallocManaged((void **)&fb, fb_size));

    int tx = 8;
    int ty = 8;

    int total_objects = 4;

    // make workd of hittable objects;
    hitable **d_list;
    checkCudaErrors(cudaMalloc((void **)&d_list, total_objects*sizeof(hitable *)));
    hitable **d_world;
    checkCudaErrors(cudaMalloc((void **)&d_world, sizeof(hitable *)));
    camera **d_camera;
    checkCudaErrors(cudaMalloc((void **)&d_camera, sizeof(camera *)));

    create_world<<<1,1>>>(d_list, d_world, d_camera, nx, ny );
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    //create a d_rand_state object for every pixel in our main routine.
    curandState *d_rand_state;
    checkCudaErrors(cudaMalloc((void **)&d_rand_state, num_pixels*sizeof(curandState)));


    clock_t start, stop;
    start = clock();

    dim3 blocks(nx/tx+1, ny/ty+1);
    dim3 threads(tx,ty);

    render_init<<<blocks, threads>>>(nx, ny, d_rand_state);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    render <<<blocks, threads>>>(fb,nx,ny, ns, d_camera, d_world, d_rand_state);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());  // no longer neccessary on Pascal or later

    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);

    // clean up
    checkCudaErrors(cudaDeviceSynchronize());
    free_world<<<1,1>>>(d_list,d_world, d_camera);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaFree(d_list));
    checkCudaErrors(cudaFree(d_world));


    checkCudaErrors(cudaFree(fb));

    return 0;
}
