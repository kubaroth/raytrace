Based on "raytracing in one weekend" series (https://www.realtimerendering.com/raytracing/) this version extends and compares 2 approaches of parallel computation. The original CPU version was extended with OpenMP library. The alternative solution adds gpu support (Cuda). Primary goal was to evaluate speedup and estimate maintenance costs having two different implementations in a single code base.

#### cpu version
- branch: master
``g++ -I. main.cpp``

#### cuda version
- branch: master
``/usr/local/cuda-10.0/bin/nvcc -I. main.cu``

#### std:: version and OpenMP support
- tag: openMP
``g++ -std=c++11 -fopenmp -I. main.cpp``
