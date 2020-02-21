#pragma once

#include <iostream>
#include <fstream>

#include "vec3.h"

// write an uncompressed PNG file from a uint8 RGB buffer
struct Png {
    FILE*f; unsigned int tab[256], crc; ~Png() { fclose(f); }
    Png(const char* fn, int w, int h, const unsigned char* c) {
        crc=0x575e51f5;
        unsigned char d[]={137,80,78,71,13,
                           10,26,10,0,0,0,13,73,72,68,82,73,68,65,84,120,1,0,
                           0,0,73,69,78,68,174,66,96,130};/*chunk headers*/
        for (int i=0;i<256;i++)/*precompute crc tables*/
            for (unsigned int j=8,&v=tab[i]=i;j--;)
                v=(v>>1)^(0xedb88320u&(~(v&1)+1));
        fwrite(d,1,16,f=fopen(fn,"w")); /*write IHDR+IDAT*/
        *this>>w>>h<<8<<2<<0<<0<<0;int a=1,b=0;/*adler-32*/
        *this>>~crc>>(h*(w*3+6)+6);int len=w*3+1;/*nbytes*/
        fwrite(d+16,1,6,f);crc=0x13e5812d;
        for (;h--;) {/*DEFLATE raw block*/
            *this<<!h<<len<<(len>>8)<<~len<<(~len>>8);
            /*filter=0*/*this <<0;/*adler(0)*/b=(a+b)%65521;
            for (int x=w,v;x--;){ /*write & checksum*/
                v=*c++;*this<<v;a=(a+v)%65521;b=(a+b)%65521;
                v=*c++;*this<<v;a=(a+v)%65521;b=(a+b)%65521;
                v=*c++;*this<<v;a=(a+v)%65521;b=(a+b)%65521;
            }
        }
        *this>>((b<<16)|a);*this>>~crc;fwrite(d+22,1,12,f);
    }
    Png& operator<<(int b){crc=(crc>>8)^tab[(crc^fputc(b,f))&255];return *this;}
    Png& operator>>(unsigned int v) {return *this<<(v>>24)<<(v>>16)<<(v>>8)<<v;}
};

void save_ppm(vec3 *fb, int width, int height){
    std::ofstream ofs ("gradient.ppm", std::ofstream::binary);
    // initial bytes
    ofs <<"P6" << std::endl
          << width << " "        << height << std::endl
          << "255" << std::endl;
    // Account for flipping the image
    for (int j = height-1; j>= 0; j--){
          for (int i = width-1; i >= 0; i--){
              size_t pixel_index = j * 1 * width + i;
              float r = fb[pixel_index].r();
              float g = fb[pixel_index].g();
              float b = fb[pixel_index].b();
              int ir = int(255.99*r);
              int ig = int(255.99*g);
              int ib = int(255.99*b);
              // std::cout << ir << " " << ig << " " << ib << '\n';
              ofs << (unsigned char)ir << (unsigned char)ig << (unsigned char)ib;
          }
    }
    ofs.close();
}

void save_png(vec3 *fb, int width, int height){
    int w = width, h = height;
    int total = w * h; // saving Png requires to flip upside down (vertically?)

    unsigned char* image = new unsigned char[3 * w * h];

    for (int j = height-1; j>=0; j--){
        for (int i = 0; i < width; i++){
            size_t pixel_index = j * width + i;

            int ir = int(255.99 * fb[total - pixel_index].r());
            int ig = int(255.99 * fb[total - pixel_index].g());
            int ib = int(255.99 * fb[total - pixel_index].b());

            image[pixel_index*3] = (unsigned char)ir;
            image[pixel_index*3+1] = (unsigned char)ig;
            image[pixel_index*3+2] = (unsigned char)ib;
        }
    }
    Png("test.png", w, h, image);
    delete [] image;

}

