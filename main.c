#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

int dimension, numberOfProcesses;
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
            printf("%f     ", array[index]);
        }
        printf("\n");
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
}

void elementsToRelax() {

}

void divideWorkAmongProcesses() {

}

void relaxMPI() {

}

int checkPrecision() {
    int precisionReached = 0;
    /* Check for precision here */
    return precisionReached;
}

int testForCorrectness() {
    int correct = 0;
    while(correct == 0) {
        /* Check for correctness here */
    }
    printf("The relaxation technique is incorrect.\n");
    return correct;
}

int main(int argc, char **argv) {

    if(argc < 4) {
        printf("Too few arguments.\n");
        return 1;
    }

    clock_t startTime, endTime;
    double seqTimeTaken, parallelTimeTaken;

    dimension = atoi(argv[1]);
    numberOfProcesses = atoi(argv[2]);
    precision = strtod(argv[3], NULL);

    if(numberOfProcesses < 1 || numberOfProcesses > 32) {
        numberOfProcesses = 32;
        printf("Invalid number of processes. Setting number of processes"
               " to default: %d.\n", numberOfProcesses);
    }

    double *relaxArray = (double*)malloc(dimension*dimension* sizeof(double));

    relaxArray = generatePattern(relaxArray);

    printf("The array in main now is: \n");

    printArray(relaxArray);

    double *relaxSeqArray = (double*)malloc(dimension*dimension* sizeof(double));

    for(int n = 0; n <= dimension*dimension; n++) {
        relaxSeqArray[n] = relaxArray[n];
    }

    startTime = clock();

    relaxSequentially(relaxSeqArray);

    endTime = clock();

    seqTimeTaken = ((double)(endTime-startTime));

    printf("The array after sequential relax is: \n)");

    printArray(relaxSeqArray);

    printf("Time taken to relax sequentially: %f\n", seqTimeTaken);

    printf("The relaxArray is still: \n");

    printArray(relaxArray);

    printf("Hello, World!\n");
    return 0;
}