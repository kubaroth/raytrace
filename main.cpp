#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>

#include "vec3.h"
#include "utils.h"

#include <float.h>  // FLT_MAX

#include <embree3/rtcore.h>

#include <tbb/tbb.h>

using std::cout;
using std::endl;
// using namespace tbb;

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
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5*color( ray(rec.p, target-rec.p), world);  // recursion

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

void errorFunction(void* userPtr, enum RTCError error, const char* str)
{
  printf("error %d: %s\n", error, str);
}

class Ray{
private:
        struct RTCRayHit rayhit;

public:
    explicit Ray(float ox, float oy, float oz,
                 float dx, float dy, float dz)
            {
                rayhit.ray.org_x = ox;
                rayhit.ray.org_y = oy;
                rayhit.ray.org_z = oz;
                rayhit.ray.dir_x = dx;
                rayhit.ray.dir_y = dy;
                rayhit.ray.dir_z = dz;
                rayhit.ray.tnear = 0;
                rayhit.ray.tfar = std::numeric_limits<float>::infinity();
                rayhit.ray.mask = -1;
                rayhit.ray.flags = 0;
                rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
                rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            }

    RTCRayHit* operator&() {return &rayhit;};
    RTCRayHit operator*() {return rayhit;};
    RTCRayHit get() {return rayhit;};

};


void waitForKeyPressedUnderWindows()
{
#if defined(_WIN32)
  HANDLE hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);

  CONSOLE_SCREEN_BUFFER_INFO csbi;
  if (!GetConsoleScreenBufferInfo(hStdOutput, &csbi)) {
    printf("GetConsoleScreenBufferInfo failed: %d\n", GetLastError());
    return;
  }

  /* do not pause when running on a shell */
  if (csbi.dwCursorPosition.X != 0 || csbi.dwCursorPosition.Y != 0)
    return;

  /* only pause if running in separate console window. */
  printf("\n\tPress any key to exit...\n");
  int ch = getch();
#endif
}


int main (){
    int ns = 10; // number of samples
    int nx = 500;
    int ny = 250;
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

    ////////// - Embree

    RTCDevice device = rtcNewDevice(NULL);
    RTCScene scene = rtcNewScene(device);

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);


    float* vertices = (float*) rtcSetNewGeometryBuffer(geom,
                                                     RTC_BUFFER_TYPE_VERTEX,
                                                     0,
                                                     RTC_FORMAT_FLOAT3,
                                                     3*sizeof(float),
                                                     3);

    unsigned* indices = (unsigned*) rtcSetNewGeometryBuffer(geom,
                                                          RTC_BUFFER_TYPE_INDEX,
                                                          0,
                                                          RTC_FORMAT_UINT3,
                                                          3*sizeof(unsigned),
                                                          1);

    if (vertices && indices) {
        vertices[0] = 0.f; vertices[1] = 0.f; vertices[2] = 0.f;
        vertices[3] = 1.f; vertices[4] = 0.f; vertices[5] = 0.f;
        vertices[6] = 0.f; vertices[7] = 1.f; vertices[8] = 0.f;

        indices[0] = 0; indices[1] = 1; indices[2] = 2;
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);


    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    auto r1 = Ray(0, 0, -1, 0, 0, 1);
    rtcIntersect1(scene, &context, &r1);


  // printf("%f, %f, %f: ", ox, oy, oz);
  if (r1.get().hit.geomID != RTC_INVALID_GEOMETRY_ID)
  {
    /* Note how geomID and primID identify the geometry we just hit.
     * We could use them here to interpolate geometry information,
     * compute shading, etc.
     * Since there is only a single triangle in this scene, we will
     * get geomID=0 / primID=0 for all hits.
     * There is also instID, used for instancing. See
     * the instancing tutorials for more information */
    printf("Found intersection on geometry %d, primitive %d at tfar=%f\n",
           r1.get().hit.geomID,
           r1.get().hit.primID,
           r1.get().ray.tfar);
  }
  else
    printf("Did not find any intersection.\n");


    // /* This will hit the triangle at t=1. */
    // castRay(scene, 0, 0, -1, 0, 0, 1);

    // /* This will not hit anything. */
    // castRay(scene, 1, 1, -1, 0, 0, 1);


    rtcReleaseScene(scene);
    rtcReleaseDevice(device);


    //////////
    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);

    return 0;
}
