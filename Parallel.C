#include "Parallel.h"
Parallel::Parallel(int choice, int thread_num):choice(choice),thread_num(thread_num){
    if (choice == 1){
        omp_set_num_threads(thread_num);
    }
}