#include <stdio.h>
#include <omp.h>

int main()
{
    int num_threads = 4; // Number of threads to use

    // Set the number of threads to be used for parallelism
    omp_set_num_threads(num_threads);

#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        printf("Hello from thread %d\n", thread_id);
    }

    return 0;
}
