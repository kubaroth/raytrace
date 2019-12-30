#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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


void render(vec3 *fb, int max_x, int max_y,
            vec3 lower_left_corner,
            vec3 horizontal,
            vec3 vertical,
            vec3 origin,
            hitable *world){
    for (int j = max_y-1; j >=0; j--){
        for (int i =0; i < max_x; i++){
            int pixel_index = j*max_x + i;
            float u = float(i)/ float(max_x);
            float v = float(j) / float(max_y);
            ray r(origin, lower_left_corner + u * horizontal + v*vertical);
            fb[pixel_index] = color(r, world);

        }
    }
}


void create_world(hitable **d_list, hitable **d_world){
    *(d_list) = new sphere(vec3(0,0,-1), 0.5);
    *(d_list+1) = new sphere(vec3(0,-100.5, -1), 100);
    *d_world = new hitable_list(d_list,2);
}


void free_world(hitable **d_list, hitable **d_world){
    delete *(d_list);
    delete *(d_list+1);
    delete d_world;
}

int main (){
    int nx = 200;
    int ny = 100;
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels * sizeof(vec3);

    //allocate FB
    vec3 *fb;
    fb = new vec3[fb_size];

    // make workd of hittable objects;
    hitable* d_list[2];
    hitable *d_world;
    create_world(d_list, &d_world);

    clock_t start, stop;
    start = clock();

    // Render
    render(fb,nx,ny,
           vec3(-2.0, -1.0, -1.0),
           vec3(4.0,0.0,0.0),  /* camera decrption*/
           vec3(0.0,2.0, 0.0),
           vec3(0.0, 0.0, 0.0),
           d_world
        );

    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);

    return 0;
}
