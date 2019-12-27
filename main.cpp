#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "vec3.h"
#include "utils.h"

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
vec3 color(const ray& r){
    if (hit_sphere(vec3(0,0,-1), 0.5,r))
        return vec3(1,0,0);
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5f*(unit_direction.y() +1.0f);
    return (1.0f-t)*vec3(1.0,1.0,1.0) + t*vec3(0.5,0.7,1.0);
}


void render(vec3 *fb, int max_x, int max_y,
            vec3 lower_left_corner,
            vec3 horizontal,
            vec3 vertical,
            vec3 origin){
    for (int j = max_y-1; j >=0; j--){
        for (int i =0; i < max_x; i++){
            int pixel_index = j*max_x + i;
            float u = float(i)/ float(max_x);
            float v = float(j) / float(max_y);
            ray r(origin, lower_left_corner + u * horizontal + v*vertical);
            fb[pixel_index] = color(r);

        }
    }
}

int main (){
    int nx = 200;
    int ny = 100;
    int num_pixels = nx*ny;
    size_t fb_size = num_pixels * sizeof(vec3);

    //allocate FB
    vec3 *fb;
    fb = new vec3[fb_size];

    clock_t start, stop;
    start = clock();

    // Render
    render(fb,nx,ny,
           vec3(-2.0, -1.0, -1.0),
           vec3(4.0,0.0,0.0),  /* camera decrption*/
           vec3(0.0,2.0, 0.0),
           vec3(0.0, 0.0, 0.0)
        );

    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);

    return 0;
}
