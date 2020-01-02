#pragma once

#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;  // TODO

class vec3{
  public:
    float e[3]; // all public as we have bunch of free functions to access internal v.e[1] etc  

      vec3() {}
      vec3(float e0, float e1, float e2) : e{e0, e1, e2} {}
    //{ e[0] = e0; e[1]=e1; e[2] = e2;}
      inline float x() const {return e[0]; }
      inline float y() const {return e[1]; }
      inline float z() const {return e[2]; }
      inline float r() const {return e[0]; }
      inline float g() const {return e[1]; }
      inline float b() const {return e[2]; }

      inline const vec3& operator+() const {return *this;}
      inline vec3 operator-() const {return vec3(-e[0], -e[1], -e[2]);}
      inline float operator[](int i) const {return e[i];}
      inline float & operator[](int i) {return e[i];}

      inline vec3& operator+=(const vec3 &v2);
      inline vec3& operator-=(const vec3 &v2);
      inline vec3& operator*=(const vec3 &v2);
      inline vec3& operator/=(const vec3 &v2);
      inline vec3& operator*=(const float t);
      inline vec3& operator/=(const float t);

    // not included in the original
      inline vec3 operator+(const float t) { return vec3( x()+t, y()+t, z()+t); } 
      inline vec3 operator-(const float t) { return vec3( x()-t, y()-t, z()-t); } 
      inline vec3 operator*(const float t) { return vec3( x()*t, y()*t, z()*t); } 
      inline vec3 operator/(const float t) { return vec3( x()/t, y()/t, z()/t); } 


      inline float length() const  {return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);}
      inline float squared_length() const {return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];}
      inline void make_unit_vector();

};

  inline vec3& vec3::operator+=(const vec3 &v){
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];
    return *this;
}
  inline vec3& vec3::operator-=(const vec3 &v){
    e[0] -= v.e[0];
    e[1] -= v.e[1];
    e[2] -= v.e[2];
    return *this;
}

  inline vec3& vec3::operator*=(const vec3 &v){
    e[0] *= v.e[0];
    e[1] *= v.e[1];
    e[2] *= v.e[2];
    return *this;
}

  inline vec3& vec3::operator/=(const vec3 &v){
    e[0] /= v.e[0];
    e[1] /= v.e[1];
    e[2] /= v.e[2];
    return *this;
}

  inline vec3& vec3::operator*=(const float t){
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
}

  inline vec3& vec3::operator/=(const float t){
    float k = 1.0/t;
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
    return *this;
}

  inline std::istream& operator>>(std::istream &is, vec3 &t){
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

inline std::ostream& operator<<(std::ostream &os,  vec3 &t){
    os <<t.e[0] << t.e[1] << t.e[2];
    return os;
}

  inline void vec3::make_unit_vector(){
    float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    e[0]*=k; e[1]*=k; e[2]*=2;
}

  inline vec3 operator+(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}

  inline vec3 operator-(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}

  inline vec3 operator*(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
}

  inline vec3 operator/(const vec3 &v1, const vec3 &v2){
    return vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
}

  inline vec3 operator*(float t, const vec3 &v){
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

  inline vec3 operator/(float t, const vec3 &v){
    return vec3(t / v.e[0], t / v.e[1], t / v.e[2]);
}

  inline vec3 operator*(const vec3 &v, float t){
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

  inline float dot(const vec3 &v1, const vec3 &v2){
    return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
}

  inline vec3 cross(const vec3 &v1, const vec3 &v2){
    return vec3 ( (v1.e[1] * v2.e[2] - v1.e[2] * v2.e[1]),
                  (-(v1.e[0] * v2.e[2] - v1.e[2] * v2.e[0])),
                  (v1.e[0] * v2.e[1] - v1.e[1] *v2.e[0]));
}

  inline vec3 unit_vector(vec3 v){  // not a memeber funciton
    return v / v.length();
}


vec3 random_in_unit_sphere() {
    vec3 p;
    do {
        p = 2.0*vec3(drand48(),drand48(),drand48()) - vec3(1,1,1);
    } while (p.squared_length() >= 1.0);
    return p;
}


class ray{
    vec3 A;
    vec3 B;

  public:
     ray(){}
     ray(const vec3& a, const vec3& b) {A = a; B = b; }
     vec3 origin() const { return A; }
     vec3 direction() const { return B; }
     vec3 point_at_parameter(float t) const {return A + t*B;}
};

class material;

struct hit_record{
    float t;
    vec3 p;
    vec3 normal;
    material *mat_ptr;
};

class hitable {
  public:
     virtual bool hit(const ray& r,
                      float t_min,
                      float t_max,
                      hit_record& rec) const = 0;
    virtual ~hitable() {};
};


class hitable_list : public hitable{
    vector<hitable*> list;
    int list_size;
  public:
     hitable_list(){}
     hitable_list(const vector<hitable*> &l) {list = l; list_size=l.size();}
     virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
 
};

 bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {

    hit_record temp_rec;
    bool hit_anything = false;
    float closes_so_far = t_max;
    for (int i = 0; i < list_size; i++){
        if (list[i]->hit(r,t_min, closes_so_far, temp_rec)){
            hit_anything = true;
            closes_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
 }

class sphere: public hitable {
   vec3 center;
   float radius;
   material * mat;
  public:
     sphere() {}
    sphere(vec3 cen, float r, material * m) : center(cen), radius(r), mat(m) {};
     virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;

};
bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const{
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc,oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0){
        float temp = (-b - sqrt(discriminant))/a;
        if (temp < t_max && temp >t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat;
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp<t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat;
            return true;
        }
    }
    return false;
}

class camera{
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
public:
     camera(){
        lower_left_corner = vec3(-2.0, -1.0, -1.0);
        horizontal = vec3(4.0,0.0,0.0);
        vertical= vec3(0.0,2.0,0.0);
        origin = vec3(0.0,0.0,0.0);
    }
     ray get_ray(float u, float v) {
        return ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
    }

};

class material{  /* needs also be added to  hit_racord class */
  public:
    /* when ray hit the surface the hit_record material will be set to point to
       the color() then can look it up
     */
    virtual bool scatter(const ray& r_in,
                         const hit_record& rec, /* to stuff whatever information we want here instead of explicit parms*/
                         vec3 &attenuation,
                         ray & scattered) const = 0;
};

class lambertian : public material {

    vec3 albedo;    

  public:
  lambertian(const vec3& a) : albedo(a){}
    virtual bool scatter(const ray& r_in,
                         const hit_record& rec,
                         vec3 &attenuation,
                         ray & scattered) const {
        vec3 target = rec.p +rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target - rec.p);
        attenuation = albedo;
        return true;
        
    }
};


vec3 reflect(const vec3 &v, const vec3& n){
    return v - 2*dot(v,n) * n;

}


class metal : public material{
    vec3 albedo;    
    float fuzz; // 0: perfect reflection 
  public:
    metal (const vec3& a, float f=0.0) : albedo(a),fuzz(f)  {
        if (f > 1) fuzz = 1;
    }
    virtual bool scatter (const ray& r_in, const hit_record& rec, vec3 &attenuation, ray &scattered) const{
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere()); // larger fuzz more blurry metal
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);

    }
    
};
