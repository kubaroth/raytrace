#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>

#include "vec3.h"
#include "utils.h"

#include <float.h>  // FLT_MAX

using std::cout;
using std::endl;


bool hit_sphere(const vec3& center, float radius, const ray& r){
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b= 2.0f  * dot(oc, r.direction());
    float c= dot(oc,oc) - radius*radius;
    float discriminant = b*b - 4.0f*a*c;
    return (discriminant > 0.0f);

}
vec3 color(const ray& r, hitable *world){
    hit_record rec;
    bool hit_anything = world->hit(r,0.0, FLT_MAX, rec);
    if(hit_anything){
        return 0.5f*vec3(rec.normal.x() +1.0f, rec.normal.y() +1.0f, rec.normal.z() + 1.0f);
    }
    else{
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
                col += color(r, world);
            }
            fb[pixel_index] = col/float(ns);

        }
    }
}


void create_world(vector<hitable*> &d_list,
                  hitable **d_world,  /* ** - set refernce pointer to d_world */
                  camera **d_camera){
    d_list.emplace_back(new sphere(vec3(0,0,-1), 0.5));
    d_list.emplace_back(new sphere(vec3(0,-100.5, -1), 100));
    *d_world = new hitable_list(d_list);
    *d_camera = new camera();
}


void free_world(hitable **d_list, hitable **d_world){
    delete *(d_list);
    delete *(d_list+1);
    delete d_world;
}

int main (){
    int ns = 10; // number of samples
    int nx = 200;
    int ny = 100;
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

    return 0;
}
