
#### cpu version
- branch: master
``g++ -I. main.cpp``

#### cuda version
- branch: master
``/usr/local/cuda-10.0/bin/nvcc -I. main.cu``

#### std:: version and OpenMP support
- tag: openMP
``g++ -std=c++11 -fopenmp -I. main.cpp``
