#include <stdio.h>
#include <time.h>
#include <stdlib.h>

// int hello_world()
// {
//     printf("Hello World \n");
//     return 0;
// }

// int small_loop()
// {
//     for (int i = 1; i <= 10000000; i++)
//     {
//         // printf("%d ", i);
//     }

//     printf("\n");
//     return 0;
// }

#define ARRAY_SIZE 10000000 // Adjust the array size as needed

int array_loop()
{
    // Declare a large enough array
    float numbers[ARRAY_SIZE];

    // Seed the random number generator
    srand(time(NULL));

    // Store random floating-point numbers between 0 and 1 in the array
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        numbers[i] = (float)rand() / RAND_MAX; // Generates a random float between 0 and 1
    }

    // Perform the sum of all elements in the array using a loop
    float sum = 0.0;
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        sum += numbers[i];
    }

    // Print the sum
    printf("Sum of random numbers: %f\n", sum);

    return 0;
}

int main()
{
    struct timespec start, end;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);

    // Your code or function calls to measure here
    // printf("Hello, World!\n");
    // hello_world();
    // small_loop();
    array_loop();

    clock_gettime(CLOCK_MONOTONIC, &end);

    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Time taken: %f seconds\n", elapsed);

    return 0;
}
