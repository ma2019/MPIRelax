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

double * generatePattern(double * array, int dimension) {
    //double *relaxArray = (double*)malloc(dimension*dimension* sizeof(double));

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
            array[index] = generateRandom(0, 10);
        }
    }
    printf("The array to relax is: \n");
    printArray(array);
    return array;
}

void relaxSequentially() {

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
    precision = atoi(argv[3]);

    if(numberOfProcesses < 1 || numberOfProcesses > 32) {
        numberOfProcesses = 32;
        printf("Invalid number of processes. Setting number of processes"
               " to default: %d.\n", numberOfProcesses);
    }

    double *relaxArray = (double*)malloc(dimension*dimension* sizeof(double));

    relaxArray = generatePattern(relaxArray, dimension);

    printf("The array in main now is: \n");

    printArray(relaxArray);


    printf("Hello, World!\n");
    return 0;
}