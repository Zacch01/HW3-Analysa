# Author : Nitzan Tomer & Zaccharie Attias & Sharon Angdo

# Used Libraries
import math

# Global Variable [Only Used To print the iteration number]
MATRIX_COUNT = -2

# Global Variable To Store The Computer Epsilon Machine
EPSILON = 0.00001


def printIntoFile(data, message, isTrue, final=None):
    """
    Printing the data and the message content into a specified file

    :param data: Data is a list representing matrix or vector
    :param message: Message is a String representing the data explanation
    :param isTrue: If True, The Linear Equation is valid, else False
    """

    # Our Global Variable To Count The Iteration Number
    global MATRIX_COUNT

    # In Case We Are Running A New Linear Equation Calculation, It will erase the lase one
    if MATRIX_COUNT == -2:
        file = open('Calculation.txt', 'w')
        file.close()

    # Open the file and save the data
    with open('Calculation.txt', 'a+') as file:

        # In case the Linear Equation is valid
        if isTrue:

            # Saving the Linear Equation input, and the updated one
            if MATRIX_COUNT < 0:
                file.write(str(message) + '\n')
                for i in range(len(data)):
                    for j in range(len(data[0])):
                        objectData = '{: ^22}'.format(data[i][j])
                        file.write(objectData)
                    file.write('\n')
                file.write('\n')

            # In case we are printing new calculation
            if MATRIX_COUNT == 0:
                file.write('========================================================================================\n')
                for i in range(len(data) + 1):
                    if i == 0:
                        objectData = '{: ^22}'.format('Iteration')

                    else:
                        objectData = '{: ^22}'.format(chr(64 + i))

                    file.write(objectData)
                file.write('\n')

            # Saving the calculation of the Linear Equation
            if MATRIX_COUNT > -1:
                if final:
                    objectData = '{: ^22}'.format('Solution')
                else:
                    objectData = '{: ^22}'.format(str(MATRIX_COUNT + 1))
                file.write(objectData)
                for i in range(len(data)):
                    objectData = '{: ^22}'.format(data[i][0])
                    file.write(objectData)
                file.write('\n')

        # In case Linear Equation is not valid
        else:
            file.write('\n' + str(message) + '\n')

        # Increase Our Global Iteration Counter Variable
        MATRIX_COUNT = MATRIX_COUNT + 1


def GaussSeidelMethod():
    """
    Solving Linear Equation in the Gauss Seidel method

    """
    # Initialize the matrix, and vectorB
    originMatrix, originVectorB = initMatrix()

    # Check if the matrix is Quadratic matrix, and check if the vector is in appropriate size
    if len(originMatrix) == len(originMatrix[0]) and len(originVectorB) == len(originMatrix) and len(originVectorB[0]) == 1:

        # In case the matrix has one solution
        if determinantMatrix(originMatrix):

            # Organize the matrix pivots
            originMatrix, originVectorB = organizeMatrix(originMatrix, originVectorB)

            # Store if the Linear Equation is Diagonal Dominant
            isConvergent = isDiagonalDominant(originMatrix)

            # According message in case the Matrix is Not Diagonal Dominant
            if isConvergent is False:
                printIntoFile(None, 'This Is A Not Diagonal Dominant Matrix', False)

            # Our lists for the Prev iteration values, and our Current iteration values
            prevIteration = [[0 for _ in range(1)] for _ in range(len(originMatrix))]
            currentIteration = [[0 for _ in range(1)] for _ in range(len(originMatrix))]

            # The iteration loop to find the Linear Equation solution
            Counter = 0

            while True:
                if isConvergent is False and Counter > 200:
                    printIntoFile(None, 'The Matrix Is Not Convergent', False)
                    break

                # Calculate the next guess
                for i in range(len(originMatrix)):
                    rowSum = 0
                    for j in range(len(originMatrix)):
                        if i != j:
                            rowSum = rowSum + originMatrix[i][j] * currentIteration[j][0]
                    currentIteration[i][0] = (originVectorB[i][0] - rowSum) / originMatrix[i][i]

                flag = True
                for i in range(len(originMatrix)):
                    if abs(currentIteration[i][0] - prevIteration[i][0]) > EPSILON:
                        flag = False

                # Save the current iteration values into the file, and update the current solution to be the prev
                printIntoFile(currentIteration, 'Linear Equation Final Solution', True)
                prevIteration = [[currentIteration[row][0] for _ in range(len(currentIteration[0]))] for row in range(len(currentIteration))]

                # In case we found our solution, Stop the loop
                if flag:
                    break

                # Stop Condition In case of Not Dominant Diagonal
                Counter = Counter + 1

            # Saving the Linear Equation final solution
            currentIteration = [[format(currentIteration[row][0], '.' + str(int(-math.log10(EPSILON))) + 'f') for _ in range(len(currentIteration[0]))] for row in range(len(currentIteration))]
            printIntoFile(currentIteration, None, True, True)

        # According message In case there is more or less than one solution
        else:
            printIntoFile(None, 'This Is A Singular Matrix', False)

    # In case the input Linear Equation isn't meet the demands
    else:
        printIntoFile(None, "The Input Linear Equation Isn't Match", False)
