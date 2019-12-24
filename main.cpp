#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "vec3.h"
#include "utils.h"

using std::cout;
using std::endl;
void render(vec3 *fb, int max_x, int max_y){

    for (int j = max_y-1; j >=0; j--){
        for (int i =0; i < max_x; i++){
            int pixel_index = j*max_x + i;
            fb[pixel_index] = vec3(float(i) / max_x,
                                   float(j) / max_y,
                                   0.2);
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
    render(fb,nx,ny);

    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);

    return 0;
}
