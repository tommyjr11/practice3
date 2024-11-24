#ifndef PARALLEL_H
#define PARALLEL_H
#include <iostream>
#include <omp.h>
#define PARALLEL_FOR _Pragma("omp parallel for if(parallel.choice == 1)")
class Parallel{
    public:
        int choice;
        int thread_num;
        Parallel(int choice, int thread_num);

        
};


#endif