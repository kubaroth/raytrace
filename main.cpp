#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>

#include "vec3.h"
#include "utils.h"

#include <float.h>  // FLT_MAX

#include <embree3/rtcore.h>
// #include <math/vec4.h>

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

void first_test_ray(RTCScene& scene, RTCIntersectContext& context){
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
  else {
    printf("Did not find any intersection.\n");
  }

    // /* This will hit the triangle at t=1. */
    // castRay(scene, 0, 0, -1, 0, 0, 1);

    // /* This will not hit anything. */
    // castRay(scene, 1, 1, -1, 0, 0, 1);


}

void render_embree(vec3 *fb, int max_x, int max_y, int ns,
                   camera *camera,
                   RTCScene& scene, RTCIntersectContext& context) {



    
    for (int j = max_y-1; j >=0; j--){
        for (int i =0; i < max_x; i++){
            vec3 color(0,0,0);
            int pixel_index = j*max_x + i;
            for (int s=0; s<ns; s++) {  // sampling
                float u = float(i + drand48()) / float(max_x);
                float v = float(j + drand48()) / float(max_y);
                // ray r = camera->get_ray(u, v);
                // col += color(r, world);



                /* initialize ray */
                // Ray ray(Vec3fa(camera.xfm.p),
                //         Vec3fa(normalize(x*camera.xfm.l.vx + y*camera.xfm.l.vy + camera.xfm.l.vz)),
                //         0.0f,
                //         inf);

                ray cam_ray = camera->get_ray(u, v);
                auto o = cam_ray.origin();
                auto d = cam_ray.direction();
                auto r1 = Ray(o.x(), o.y(), o.z(), d.x(), d.y(), d.z());
                rtcIntersect1(scene, &context, &r1);
    
                /* intersect ray with scene */
                // rtcIntersect1(data.g_scene,&context,RTCRayHit_(ray));
                // RayStats_addRay(stats);

                // /* shade pixels */
                // Vec3fa color = Vec3fa(0.0f);
                if (r1.get().hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                //         Vec3fa diffuse = data.face_colors[ray.primID];
                //         color = color + diffuse*0.5f;
                //         Vec3fa lightDir = normalize(Vec3fa(-1,-1,-1));

                //         /* initialize shadow ray */
                //         Ray shadow(ray.org + ray.tfar*ray.dir, neg(lightDir), 0.001f, inf, 0.0f);

                //         /* trace shadow ray */
                //         rtcOccluded1(data.g_scene,&context,RTCRay_(shadow));
                //         RayStats_addShadowRay(stats);

                //         /* add light contribution */
                //         if (shadow.tfar >= 0.0f)
                //             color = color + diffuse*clamp(-dot(lightDir,normalize(ray.Ng)),0.0f,1.0f);
                    }

                // /* write color to framebuffer */
                // unsigned int r = (unsigned int) (255.0f * clamp(color.x,0.0f,1.0f));
                // unsigned int g = (unsigned int) (255.0f * clamp(color.y,0.0f,1.0f));
                // unsigned int b = (unsigned int) (255.0f * clamp(color.z,0.0f,1.0f));
                // pixels[y*width+x] = (b << 16) + (g << 8) + r;


            }
            fb[pixel_index] = color/float(ns);

        }
    }
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

    // typedef Vec4<float> Vec4f;
    RTCGeometry spherePoint = rtcNewGeometry(device,   RTC_GEOMETRY_TYPE_SPHERE_POINT);

    int NUM_POINTS = 1;
    float* point_vertices = (float*)rtcSetNewGeometryBuffer(spherePoint,
                                                            RTC_BUFFER_TYPE_VERTEX,
                                                            0,
                                                            RTC_FORMAT_FLOAT4,
                                                            4 * sizeof(float),
                                                            NUM_POINTS);
    point_vertices[0] = 0.f;
    point_vertices[1] = 0.f;
    point_vertices[2] = 0.f;
    point_vertices[3] = 1.f;  // radius?
    
    // TODO set position, st camera
    rtcCommitGeometry(spherePoint);
    rtcAttachGeometry(scene, spherePoint);
    rtcReleaseGeometry(spherePoint);

    
    rtcCommitScene(scene);  // build bvh

    
    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    // first_test_ray(scene, context);
    vec3 *fb_embree;
    fb_embree = new vec3[fb_size];

    render_embree(fb_embree, nx,ny, ns,
                  d_camera,
                  scene, context);
    
    
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);


    //////////
    stop = clock();
    double timer_seconds = ((double)(stop-start)) / CLOCKS_PER_SEC;
    std::cerr << "elapsed time:" << timer_seconds << std::endl;

    save_png(fb, nx, ny);
    save_png(fb_embree, nx, ny, "test_embree.png");

    return 0;
}
