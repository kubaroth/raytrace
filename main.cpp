// g++ -I. main.cpp
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>
#include <memory>

#include "vec3.h"
#include "utils.h"

#include <float.h>  // FLT_MAX
#include <cassert>

using std::cout;
using std::endl;


vec3 color(const ray& r, hitable **world, int depth){
    hit_record rec;
    if ((*world)->hit(r, 0.001, FLT_MAX, rec)) {  //ignore hits very near zero
        ray scattered;
        vec3 attenuation;  // atenuation gets updated in the mat_ptr->scatterer
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation,scattered)){
            return attenuation * color(scattered, world, depth+1);
        }
        else{
            return vec3(0,0,0);
        }
    }
    else{  // background
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5f*(unit_direction.y() +1.0f);
        return (1.0f-t)*vec3(1.0,1.0,1.0) + t*vec3(0.5,0.7,1.0);
    }
}


void render(vec3 *fb, int max_x, int max_y, int ns,
            camera **camera,
            hitable **world){
    for (int j = max_y-1; j >=0; j--){
        for (int i =0; i < max_x; i++){
            vec3 col(0,0,0);
            int pixel_index = j*max_x + i;
            for (int s=0; s<ns; s++) {  // sampling
                float u = float(i + drand48()) / float(max_x);
                float v = float(j + drand48()) / float(max_y);
                ray r = (*camera)->get_ray(u, v);
                col += color(r, world, 0 /*depth*/);
            }
            fb[pixel_index] = col/float(ns);

        }
    }
}

void create_world(hitable** d_list,
                  hitable **d_world,  /* ** - set refernce pointer to d_world */
                  int* total_objects,
                  camera **d_camera,
                  int nx, int ny){

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
        // At this point we can only verify if number of objects added to d_list matches initial value
        assert (*total_objects == obj_index);
        *d_world = new hitable_list(d_list, obj_index);

}


void free_world(hitable **d_list, hitable **d_world, int* total_objects, camera **d_camera){
    for (int i = 0; i < *total_objects; ++i){
        delete *(d_list + i);
    }
    delete *d_world;
    delete *d_camera;

}

int main (){
    int ns = 10; // number of samples
    int nx = 600;
    int ny = 300;
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels * sizeof(vec3);

    int total_objects = 4;

    //allocate FB
    vec3 *fb;
    fb = new vec3[fb_size];

    hitable **d_list;
    hitable **d_world;
    camera **d_camera;

    d_list = new hitable*[total_objects];
    d_world = new hitable*[1];
    d_camera = new camera*[1];

    create_world(d_list, d_world, &total_objects, d_camera, nx, ny);

    clock_t start, stop;
    start = clock();

    // Render
    render(fb,nx,ny, ns,
           d_camera,
           d_world
        );

    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);
    save_ppm(fb, nx, ny);

    free_world(d_list, d_world, &total_objects, d_camera);;

    delete fb;
    return 0;
}
