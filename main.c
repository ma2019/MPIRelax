#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>
#include <stdbool.h>

int dimension;
int numberOfActiveProcesses;
double targetPrecision;

void printArray(int processRank, double ** array) {
    printf("Rank: %d\n", processRank);
    for(int column = 0; column < dimension; column++) {
        for(int row = 0; row < dimension; row++) {
            printf("%f     ", array[column][row]);
        }
        printf("\n");
    }
    printf("Address: %p\n", &array);
}

double generateRandom(double randMin, double randMax) {
    double rangeOfNumbers = randMax - randMin;
    double divide = RAND_MAX/rangeOfNumbers;
    return randMin + (rand()/divide);
}

double xPlusy(double x, double y) {
    return x+y;
}

double ** generatePattern(int processRank, double ** array) {
    for(int row = 0; row < dimension; row++) {
        for(int column = 0; column < dimension; column++) {
            //int index = dimension*row+column;
            array[row][column] = xPlusy(row, column);
        }
    }
    //printf("The final array is: \n");
    //printArray(processRank, array);
    for(int row = 1; row < dimension-1; row++) {
        for(int column = 1; column < dimension-1; column++) {
            //int index = dimension*row+column;
            array[row][column] = generateRandom(0, 100);
        }
    }
    //printf("The array to relax is: \n");
    //printArray(processRank, array);
    return array;
}

int processInChargeOfThis(int column, int * columnToRelaxArray,
                          int length, int processRank) {
    int relaxableColumns = dimension-2;
    int modulo = column%numberOfActiveProcesses;
    for(int i = 0; i < length; i++) {
        if(columnToRelaxArray[i] == column) {
            return processRank;
        }
    }
    if(column <= numberOfActiveProcesses) {
        return column-1; //column 1 will be by rank 0
    }
    if(modulo == 0) {
        return numberOfActiveProcesses-1; //highest process
    } else {
        return modulo-1;
    }
}

int * divideColumnsAmongProcesses(int processRank, int numberOfProcesses) {
    int relaxableColumns = dimension-2;

    if(numberOfProcesses > relaxableColumns) {
        numberOfActiveProcesses = relaxableColumns;
    } else {
        numberOfActiveProcesses = numberOfProcesses;
    }
    /*
     * 3 options:
     * More processes than relaxable columns
     *    Only some processes will relax
     *    (Introduce activeProcesses for precisionReached)
     * Same processes than relaxable columns
     *    All processes will have 1 column (same as next)
     * More columns than processes - modulo=0
     *    Equal number of columns, cols/procs
     * More columns than processes - modulo!=0
     *    Some will have 1 extra column
     *    (modulo number of processes will have extra column)
     */
    int modulo = relaxableColumns%numberOfActiveProcesses;
    int columnsPerProcess = (relaxableColumns/numberOfActiveProcesses);
    int *columnsToRelax = (int*)malloc((columnsPerProcess+1) * sizeof(int));

    //Fill with dummy numbers
    for(int i = 0; i <= columnsPerProcess; i++) {
        columnsToRelax[i] = -1;
    }

    int processNumber = processRank + 1;

    //First fill the normal ones, ignoring index at the end
    for(int i = 0; i < columnsPerProcess; i++) {
        columnsToRelax[i] = (i*numberOfActiveProcesses)+processNumber;
        /*
         * At the final iteration, i is 2 less than the size
         * since array starts at 0 and we are ignoring the last index.
         */
    }

    if(modulo!=0) { //now populate the edge elements
        /*
         * Go over only the first modulo number of processes
         */
        for(int processesWithExtraColumn = 0;
        processesWithExtraColumn < modulo; processesWithExtraColumn++) {
            if(processesWithExtraColumn == processRank) {
                columnsToRelax[columnsPerProcess] =
                        ((columnsPerProcess)*numberOfActiveProcesses)+processNumber;
                /*
                 * Here columnsPerProcess is the last index, but still -1
                 * as array starts at 0.
                 */
            }
        }
    }


    return columnsToRelax;
}

int checkPrecision(double achievedPrecision) {
    int precisionReached = 0;
    if(achievedPrecision < targetPrecision) {
        precisionReached = 1;
    } else {
        precisionReached = 0;
    }
    return precisionReached;
}

double relax(double ** array, int row, int column) {
    double relaxed = (array[row+1][column] + array[row-1][column]+
            array[row][column-1] + array[row][column+1]) / 4;
    return relaxed;
}

bool valueIsInArray(const int * columnToRelaxArray, int size, int column) {
    for(int everyColumn = 0; everyColumn < size; everyColumn++) {
        if(columnToRelaxArray[everyColumn] == column) {
            return true;
        }
    }
    return false;
}

void relaxMpi(double ** array, int * columnToRelaxArray, int lengthOfColumns,
        int processRank, int numberOfActiveProcesses) {

    double achievedPrecision = -1;

    int precisionReachedByRoot = 0;
    int precisionReachedByAll = 0;
    int toProceed = 1;

    double **copyOfArray = (double**)malloc(dimension* sizeof(double));
    for(int i = 0; i < dimension; i++) {
        copyOfArray[i] = (double*)malloc(dimension* sizeof(double));
    }

    double *myArray = (double*)malloc((dimension-2)* sizeof(double));
    double *leftArray = (double*)malloc((dimension-2)* sizeof(double));
    double *rightArray = (double*)malloc((dimension-2)* sizeof(double));

    for(int i = 0; i < dimension-2; i++) {
        myArray[i] = -1.0;
        leftArray[i] = -1.0;
        rightArray[i] = -1.0;
    }

    int iteration = 0;

    while(toProceed == 1) {
        for(int column = 1; column < dimension-1; column++) {
            if(valueIsInArray(columnToRelaxArray, lengthOfColumns, column)) {
                for(int row = 1; row < dimension-1; row++) {
                    copyOfArray[row][column] = array[row][column];
                    array[row][column] = relax(array, row, column);
                    for(int n = 0; n < dimension-2; n++) {
                        if(myArray[n] == -1) {
                            myArray[n] = array[row][column];
                            break;
                        }
                    }

                    achievedPrecision = copyOfArray[row][column] -
                            array[row][column];

                    if(achievedPrecision < 0) {
                        achievedPrecision = achievedPrecision * (-1);
                    }

                    precisionReachedByRoot = checkPrecision(achievedPrecision);
                }
            }
            int dutyProcess = processInChargeOfThis(
                    column, columnToRelaxArray, lengthOfColumns, processRank);
            MPI_Barrier(MPI_COMM_WORLD);
            for(int row = 0; row < dimension-1; row++) {
                MPI_Bcast(&array[row][column], 1, MPI_DOUBLE,
                        dutyProcess, MPI_COMM_WORLD);
            }
        }

        MPI_Allreduce(&precisionReachedByRoot, &precisionReachedByAll,
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(precisionReachedByAll == numberOfActiveProcesses) {
            toProceed = 0;
        }
        iteration++;
    }
}

void testForCorrectness(const double * array) {
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

    int rc, processRank, numberOfProcesses;

    rc = MPI_Init(&argc, &argv);

    if(rc != MPI_SUCCESS) {
        printf("Error.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    dimension = atoi(argv[1]);
    targetPrecision = strtod(argv[2], NULL);

    /**
     * How would I check that I am passing in valid number of processes?
     * if(numberOfProcesses < 1 || numberOfProcesses > 32) {
          numberOfProcesses = 32;
          printf("Invalid number of processes. Setting number of processes"
               " to default: %d.\n", numberOfProcesses);
       }
     */


    double **inputArray = (double**)malloc(dimension* sizeof(double));
    for(int i = 0; i < dimension; i++) {
        inputArray[i] = (double*)malloc(dimension* sizeof(double));
    }

    generatePattern(processRank, inputArray);

    double **relaxMpiArray = (double**)malloc(dimension * sizeof(double));
    for(int i = 0; i < dimension; i++) {
        relaxMpiArray[i] = (double*)malloc(dimension* sizeof(double));
    }

    for(int row = 0; row < dimension; row++) {
        for(int column = 0; column < dimension; column++) {
            relaxMpiArray[row][column] = inputArray[row][column];
        }
    }

    int* columnsToRelax = divideColumnsAmongProcesses(processRank,
            numberOfProcesses);

    int lengthOfColumnsToRelax = (int) (sizeof(columnsToRelax)/
            sizeof(columnsToRelax[0]));

    if(processRank < numberOfActiveProcesses) {
        relaxMpi(relaxMpiArray, columnsToRelax, lengthOfColumnsToRelax,
                 processRank, numberOfActiveProcesses);
    }
    printArray(processRank, relaxMpiArray);

    MPI_Finalize();

    return 0;
}