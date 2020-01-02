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


vec3 color(const ray& r, hitable *world, int depth){
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) {  //ignore hits very near zero
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
            camera *camera,
            hitable *world){
    for (int j = max_y-1; j >=0; j--){
        for (int i =0; i < max_x; i++){
            vec3 col(0,0,0);
            int pixel_index = j*max_x + i;
            for (int s=0; s<ns; s++) {  // sampling
                float u = float(i + drand48()) / float(max_x);
                float v = float(j + drand48()) / float(max_y);
                ray r = camera->get_ray(u, v);
                col += color(r, world, 0 /*depth*/);
            }
            fb[pixel_index] = col/float(ns);

        }
    }
}


void create_world(vector<hitable*> &d_list,
                  hitable **d_world,  /* ** - set refernce pointer to d_world */
                  camera **d_camera){

    // ground
    d_list.emplace_back(new sphere(vec3(0,-100.5, -1), 100, make_unique<lambertian>(vec3(0.8,0.8,0.0))));
    d_list.emplace_back(new sphere(vec3(0,0,-1), 0.5, make_unique<lambertian>(vec3(0.8,0.3,0.3))));
    d_list.emplace_back(new sphere(vec3(-1,0,-1), 0.5, make_unique<metal>(vec3(0.5,0.5,0.5), 0.8 /*fuzzy*/)));
    d_list.emplace_back(new sphere(vec3(1,0,-1), 0.5, make_unique<dielectric>(1.5)));
        
    *d_world = new hitable_list(d_list);
    *d_camera = new camera();
}


void free_world(vector<hitable*> &d_list,
                hitable **d_world,
                camera **d_camera){

    // for (auto &elem : d_list) delete elem;
    for (int i=0; i < d_list.size(); i++){
        hitable *aa = d_list[i];
        delete aa;
    }
    delete *d_world;
    delete *d_camera;

}

int main (){
    int ns = 1; // number of samples
    int nx = 20;
    int ny = 10;
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels * sizeof(vec3);

    //allocate FB
    vec3 *fb;
    fb = new vec3[fb_size];

    // make world of hittable objects;
    vector<hitable*> d_list;
    hitable *d_world;
    camera *d_camera;
    create_world(d_list, &d_world, &d_camera);
    
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

    free_world(d_list, &d_world, &d_camera);;
    delete d_camera;
    delete fb;
    return 0;
}
