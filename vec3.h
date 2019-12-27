#pragma once

#include <math.h>
#include <stdlib.h>
#include <iostream>

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
