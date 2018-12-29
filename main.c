#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>

int dimension;
//int numberOfProcesses;
double precision;

void printArray(double * array) {
    for(int row = 0; row < dimension; row++) {
        for(int column = 0; column < dimension; column++) {
            int index = dimension*row+column;
            printf("%f     ", array[index]);
        }
        printf("\n");
    }
}

double generateRandom(double randMin, double randMax) {
    double rangeOfNumbers = randMax - randMin;
    double divide = RAND_MAX/rangeOfNumbers;
    return randMin + (rand()/divide);
}

double xPlusy(double x, double y) {
    return x+y;
}

double * generatePattern(double * array) {
    for(int row = 0; row < dimension; row++) {
        for(int column = 0; column < dimension; column++) {
            int index = dimension*row+column;
            array[index] = xPlusy(row, column);
        }
    }
    printf("The final array is: \n");
    printArray(array);
    for(int row = 1; row < dimension-1; row++) {
        for(int column = 1; column < dimension-1; column++) {
            int index = dimension*row+column;
            array[index] = generateRandom(0, 100);
        }
    }
    printf("The array to relax is: \n");
    printArray(array);
    return array;
}

void relaxSequentially(double *array) {
    double achievedPrecision;
    double precisionReached = 0;

    double *copyOfArray = (double*)malloc(dimension*dimension* sizeof(double));

    while (precisionReached == 0) {
        for (int row = 1; row < dimension - 1; row++) {
            for (int column = 1; column < dimension - 1; column++) {
                copyOfArray[dimension * row + column] =
                        array[dimension * row + column];
                array[dimension * row + column] =
                        (((array[(dimension * row + column) - 1])
                        + (array[(dimension * row + column) + 1])
                        + (array[(dimension * row + column)
                        - dimension]) + (array[(dimension * row
                        + column) + dimension])) / 4);
                achievedPrecision = copyOfArray[dimension * row + column]
                        - array[dimension * row + column];
                if (achievedPrecision < 0) {
                    achievedPrecision = achievedPrecision * (-1);
                }
                if (achievedPrecision < precision) {
                    precisionReached = 1;
                } else {
                    precisionReached = 0;
                }
            }

        }
    }
    //free(&copyOfArray);
}

int * elementsToRelax(int relaxableElements) {
    int *indicesToRelax = (int*)malloc(relaxableElements * sizeof(int));

    for(int i = 0; i < relaxableElements; i++) {
        indicesToRelax[i] = -1;
    }

    int indexToAddTo = 0;
    for (int row = 1; row < dimension - 1; row++) {
        for (int column = 1; column < dimension - 1; column++) {
            int index = dimension * row + column;
            indicesToRelax[indexToAddTo] = index;
            indexToAddTo++;
        }
    }
    return indicesToRelax;
}

int * divideWorkAmongProcesses(int relaxableElements, int processRank,
        int numberOfProcesses, int* startAndStop) {
    int numberOfElementsPerThread = relaxableElements/numberOfProcesses;
    int startingIndex = processRank*numberOfElementsPerThread;
    int stoppingIndex;
    if(processRank != numberOfProcesses-1) {
        stoppingIndex = startingIndex + numberOfElementsPerThread;
    } else {
        stoppingIndex = relaxableElements;
    }
    startAndStop[0] = startingIndex;
    startAndStop[1] = stoppingIndex;
    return startAndStop;
}

void relaxMpi(double * array, int * startAndStop, int * indicesToRelax, int processRank) {

    double achievedPrecision;

    int precisionReached = 0;
    int startingIndex = startAndStop[0];
    int stoppingIndex = startAndStop[1];
    double *copyOfArray = (double*)malloc(dimension*dimension* sizeof(double));

    while(precisionReached == 0) {
        for(int i = startingIndex; i < stoppingIndex; i++) {

            /**
             * indicesToRelax has all the indices to relax;
             * we are just extracting the ones needed for this process using
             * the starting and stopping indices calculated earlier.
            */
            int indexToRelax = indicesToRelax[i];

            /**
             * making a copy to compare for precision later
             */
            copyOfArray[indexToRelax] = array[indexToRelax];

            /**
             * replacing the element with its border elements;
             * explain the maths
             */
            array[indexToRelax] = (((array[(indexToRelax) - 1]) +
                    (array[(indexToRelax) + 1]) +
                    (array[(indexToRelax) - dimension]) +
                    (array[(indexToRelax) + dimension])) / 4);
            achievedPrecision = copyOfArray[indexToRelax]-array[indexToRelax];

            if(achievedPrecision < 0) {
                achievedPrecision = achievedPrecision * (-1);
            }
            if(achievedPrecision < precision) {
                precisionReached = 1;
            } else {
                precisionReached = 0;
            }
        }
        if(precisionReached == 1) {
            MPI_Bcast(&precisionReached, 1, MPI_INT, processRank, MPI_COMM_WORLD);
        }
    }


}

int checkPrecision() {
    int precisionReached = 0;
    /* Check for precision here */
    return precisionReached;
}

void testForCorrectness(double * array) {
    int correct = 1;
    for(int row = 0; row < dimension; row++) {
        for(int column = 0; column < dimension; column++) {
            int index = dimension*row+column;
            int functionValue = (int) xPlusy(row, column);
            int arrayElement = (int) array[index];
            if(arrayElement == xPlusy(row, column)) {
                correct = 0;
            } else {
                correct = 1;
                printf("The element at (%d,%d) should be %d, is %d\n",
                        row, column, functionValue, arrayElement);
            }
        }
    }
    if(correct == 0) {
        printf("The relaxation technique is correct!\n");
    }
}

int main(int argc, char **argv) {

//    if(argc < 3) {
//        printf("Too few arguments.\n");
//        return 1;
//    }

    clock_t startTime, endTime;
    double seqTimeTaken, parallelTimeTaken;

    dimension = 6;//atoi(argv[1]);
    precision = 0.1;//strtod(argv[2], NULL);

    /**
     * How would I check that I am passing in valid number of processes?
     * if(numberOfProcesses < 1 || numberOfProcesses > 32) {
          numberOfProcesses = 32;
          printf("Invalid number of processes. Setting number of processes"
               " to default: %d.\n", numberOfProcesses);
       }
     */

    double *inputArray = (double*)malloc(dimension*dimension* sizeof(double));

    generatePattern(inputArray);

    printf("The array in main now is: \n");

    printArray(inputArray);

    double *relaxSeqArray = (double*)malloc(dimension*dimension* sizeof(double));

    for(int n = 0; n <= dimension*dimension; n++) {
        relaxSeqArray[n] = inputArray[n];
    }

    startTime = clock();

    relaxSequentially(relaxSeqArray);

    endTime = clock();

    seqTimeTaken = ((double)(endTime-startTime));

    printf("The array after sequential relax is: \n");

    printArray(relaxSeqArray);
    testForCorrectness(relaxSeqArray);
    free(relaxSeqArray);

    printf("Time taken to relax sequentially: %f\n", seqTimeTaken);

    printf("The inputArray is still: \n");

    printArray(inputArray);

    double *relaxMpiArray = (double*)malloc(dimension*dimension* sizeof(double));

    for(int n = 0; n <= dimension*dimension; n++) {
        relaxMpiArray[n] = inputArray[n];
    }

    printf("The array to be relaxed with MPI is: \n");
    printArray(relaxMpiArray);

    int relaxableElements = (int) pow((dimension - 2), 2);
    int* indicesToRelax = elementsToRelax(relaxableElements);

    int rc, processRank, numberOfProcesses, namelength;

    char name[MPI_MAX_PROCESSOR_NAME];

    printf("Reached here1.\n");
    rc = MPI_Init(&argc, &argv);
    printf("Reached here2.\n");
    if(rc != MPI_SUCCESS) {
        printf("Error.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    int *startAndStop = (int*)malloc(2* sizeof(int));

    divideWorkAmongProcesses(relaxableElements, processRank, numberOfProcesses,
                             startAndStop);

    relaxMpi(relaxMpiArray, startAndStop, indicesToRelax, processRank);

    printf("The MPI relaxed array is: \n");
    printArray(relaxMpiArray);

    printf("Hello, World!\n");

    MPI_Finalize();


    return 0;
}