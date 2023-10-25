// AUTHORS:
// PRITAM SUWAL SHRESTHA (23771397)
// RASPREET KHANUJA (23308425)

#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define LAKE_X_MIN (-100.0)
#define LAKE_X_MAX 100.0
#define LAKE_Y_MIN (-100.0)
#define LAKE_Y_MAX 100.0

#define MIN_SWIM_DISTANCE (-0.1)
#define MAX_SWIM_DISTANCE 0.1

#define MAX_FISH_WEIGHT 2.0
#define MIN_FISH_WEIGHT 0.0

#define NUM_FISHES 5
#define NUM_SIMULATION_STEPS 1

double square(double num) { return num * num; }

typedef struct {
  double x, y;
  double distanceTraveled;
  double weight;
} Fish;

bool isFishOutsideLake(Fish *fish) {
  return (fish->x < LAKE_X_MIN || fish->x > LAKE_X_MAX ||
          fish->y < LAKE_Y_MIN || fish->y > LAKE_Y_MAX);
}

double getRandomNumberInRange(double minValue, double maxValue) {
  double randomDouble = ((double)rand() / RAND_MAX);
  double randomNumber = randomDouble * (maxValue - minValue) + minValue;
  return randomNumber;
}

double calculateDistanceFromOrigin(double x, double y) {
  return sqrt(square(x) + square(y));
}

double calculateObjectiveFunction(Fish *fishes, int numFishes) {
  double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i < numFishes; i++) {
    sum += calculateDistanceFromOrigin((fishes + i)->x, (fishes + i)->y);
  }
  printf("Objective Function: %.2f\n", sum);
  return sum;
}

void swim(Fish *fish) {
  if (!isFishOutsideLake(fish)) {
    fish->x += getRandomNumberInRange(MIN_SWIM_DISTANCE, MAX_SWIM_DISTANCE);
    fish->y += getRandomNumberInRange(MIN_SWIM_DISTANCE, MAX_SWIM_DISTANCE);
  }
}

void updateWeight(Fish *fish, double maxDistanceTraveledInFishSchool) {
  double weightChange =
      fish->distanceTraveled / maxDistanceTraveledInFishSchool;
  fish->weight += weightChange;
  if (fish->weight < MIN_FISH_WEIGHT) {
    fish->weight = MIN_FISH_WEIGHT;
  } else if (fish->weight > MAX_FISH_WEIGHT) {
    fish->weight = MAX_FISH_WEIGHT;
  }
}

void simulationStep(Fish *fishes, int numFishes) {

  double maxDistanceTraveledInFishSchool = 0.0;

#pragma omp parallel for reduction(max : maxDistanceTraveledInFishSchool)
  for (int i = 0; i < numFishes; i++) {
    double prevDistance =
        calculateDistanceFromOrigin((fishes + i)->x, (fishes + i)->y);
    swim(fishes + i);
    double nextDistance =
        calculateDistanceFromOrigin((fishes + i)->x, (fishes + i)->y);
    (fishes + i)->distanceTraveled = nextDistance - prevDistance;
    if (fabs((fishes + i)->distanceTraveled) >
        maxDistanceTraveledInFishSchool) {
      maxDistanceTraveledInFishSchool = fabs((fishes + i)->distanceTraveled);
    }
  }

#pragma omp parallel for
  for (int i = 0; i < numFishes; i++) {
    updateWeight((fishes + i), maxDistanceTraveledInFishSchool);
    //        printf(
    //                "Step Fish %d: x = %.2f, y = %.2f, distanceTraveled =
    //                %.2f, weight = %.2f\n", i, (fishes + i)->x, (fishes +
    //                i)->y, (fishes + i)->distanceTraveled, (fishes +
    //                i)->weight
    //        );
  }
}

void calculateBarycenter(Fish *fishes, int numFishes) {
  double weightSum = 0.0;
  double distanceSum = 0.0;

#pragma omp parallel for reduction(+ : weightSum, distanceSum)
  for (int i = 0; i < numFishes; i++) {
    weightSum += (fishes + i)->weight *
                 calculateDistanceFromOrigin((fishes + i)->x, (fishes + i)->y);
    distanceSum +=
        calculateDistanceFromOrigin((fishes + i)->x, (fishes + i)->y);
  }

  if (distanceSum == 0.0) {
    printf("Distance sum is zero. Cannot calculate barycenter.\n");
    return;
  }

  double barycenter = weightSum / distanceSum;
  printf("Barycenter: %.2f\n", barycenter);
}

void initializeInitialLakeState(Fish *fishes, int numFishes) {
#pragma omp parallel for
  for (int i = 0; i < numFishes; i++) {
    (fishes + i)->x = getRandomNumberInRange(LAKE_X_MIN, LAKE_X_MAX);
    (fishes + i)->y = getRandomNumberInRange(LAKE_Y_MIN, LAKE_Y_MAX);
    (fishes + i)->distanceTraveled = 0.0;
    (fishes + i)->weight = 1.0;
    //        printf(
    //                "Initial Fish %d: x = %.2f, y = %.2f, distanceTraveled =
    //                %.2f, weight = %.2f\n", i, (fishes + i)->x, (fishes +
    //                i)->y, (fishes + i)->distanceTraveled, (fishes +
    //                i)->weight
    //        );
  }
}

int main() {

  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);

  double start_time, end_time, elapsed_time;
  start_time = MPI_Wtime(); // Start recording the elapsed time

  MPI_Datatype MPI_FISH_TYPE;

  // Create MPI type for Fish struct
  int blocklengths[4] = {1, 1, 1, 1}; // One block for each struct element
  MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

  MPI_Aint offsets[4];
  offsets[0] = offsetof(Fish, x);
  offsets[1] = offsetof(Fish, y);
  offsets[2] = offsetof(Fish, distanceTraveled);
  offsets[3] = offsetof(Fish, weight);

  MPI_Type_create_struct(4, blocklengths, offsets, types, &MPI_FISH_TYPE);
  MPI_Type_commit(&MPI_FISH_TYPE);

  if (provided != MPI_THREAD_FUNNELED) {
    fprintf(stderr, "MPI does not support MPI_THREAD_FUNNEL mode.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL) + rank);

  int fishPerProcess = NUM_FISHES / size;
  int remainingFishes = NUM_FISHES % size;
  int localFishCount =
      (rank < remainingFishes) ? fishPerProcess + 1 : fishPerProcess;

  Fish *local_fishes = (Fish *)malloc(localFishCount * sizeof(Fish));
  if (!local_fishes) {
    perror("Memory allocation failed for local_fishes");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if (rank == 0) {
    Fish *all_fishes = (Fish *)malloc(NUM_FISHES * sizeof(Fish));
    if (!all_fishes) {
      perror("Memory allocation failed for all_fishes");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    initializeInitialLakeState(all_fishes, NUM_FISHES);

    int offset = 0;
    for (int i = 0; i < size; i++) {
      int fishToSend =
          (i < remainingFishes) ? fishPerProcess + 1 : fishPerProcess;
      if (i == 0) {
        memcpy(local_fishes, all_fishes, fishToSend * sizeof(Fish));
      } else {
        MPI_Send(all_fishes + offset, fishToSend, MPI_FISH_TYPE, i, 0,
                 MPI_COMM_WORLD);
      }
      offset += fishToSend;
    }
    free(all_fishes);
  } else {
    MPI_Recv(local_fishes, localFishCount, MPI_FISH_TYPE, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }

  for (int i = 0; i < NUM_SIMULATION_STEPS; i++) {

    //        printf("Simulation step %d\n", i);
    simulationStep(local_fishes, localFishCount);
    //        calculateObjectiveFunction(local_fishes, localFishCount);
    //        calculateBarycenter(local_fishes, localFishCount);

    // Objective Function for the Entire Lake
    double localSum = 0.0;
    double globalSum = 0.0;

#pragma omp parallel for reduction(+ : localSum)
    for (int i = 0; i < localFishCount; i++) {
      localSum +=
          calculateDistanceFromOrigin(local_fishes[i].x, local_fishes[i].y);
    }

    MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    if (rank == 0) {
      printf("Objective Function for Entire Lake: %.2f\n", globalSum);
    }

    // Barycenter for the Entire Lake
    double localWeightSum = 0.0;
    double localDistanceSum = 0.0;
    double globalWeightSum = 0.0;
    double globalDistanceSum = 0.0;

#pragma omp parallel for reduction(+ : localWeightSum, localDistanceSum)
    for (int i = 0; i < localFishCount; i++) {
      localWeightSum +=
          local_fishes[i].weight *
          calculateDistanceFromOrigin(local_fishes[i].x, local_fishes[i].y);
      localDistanceSum +=
          calculateDistanceFromOrigin(local_fishes[i].x, local_fishes[i].y);
    }

    MPI_Allreduce(&localWeightSum, &globalWeightSum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&localDistanceSum, &globalDistanceSum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    if (globalDistanceSum == 0.0) {
      if (rank == 0) {
        printf("Distance sum is zero. Cannot calculate barycenter.\n");
      }
    } else {
      double barycenter = globalWeightSum / globalDistanceSum;
      if (rank == 0) {
        printf("Barycenter for Entire Lake: %.2f\n", barycenter);
      }
    }
  }

  if (rank != 0) {
    MPI_Send(local_fishes, localFishCount, MPI_FISH_TYPE, 0, 0, MPI_COMM_WORLD);
  } else {
    Fish *all_fishes = (Fish *)malloc(NUM_FISHES * sizeof(Fish));
    if (!all_fishes) {
      perror("Memory allocation failed for all_fishes");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int offset = 0;
    for (int i = 0; i < size; i++) {
      int fishToReceive =
          (i < remainingFishes) ? fishPerProcess + 1 : fishPerProcess;
      if (i == 0) {
        memcpy(all_fishes, local_fishes, fishToReceive * sizeof(Fish));
      } else {
        MPI_Recv(all_fishes + offset, fishToReceive, MPI_FISH_TYPE, i, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      offset += fishToReceive;
    }

    //        printf("\nAll fish after simulation steps:\n");
    //        for (int i = 0; i < NUM_FISHES; i++) {
    //            printf(
    //                    "Final Fish %d: x = %.2f, y = %.2f, distanceTraveled =
    //                    %.2f, weight = %.2f\n", i, (all_fishes + i)->x,
    //                    (all_fishes + i)->y, (all_fishes +
    //                    i)->distanceTraveled, (all_fishes + i)->weight
    //            );
    //        }

    free(all_fishes);
  }

  free(local_fishes);
  MPI_Type_free(&MPI_FISH_TYPE);

  end_time = MPI_Wtime();
  elapsed_time = end_time - start_time;

  double max_elapsed_time;
  MPI_Reduce(&elapsed_time, &max_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    printf("Maximum time taken by any process: %f seconds\n", max_elapsed_time);
  }

  MPI_Finalize();
  return 0;
}
