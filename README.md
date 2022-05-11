# scientificComputing
Some problems in high-precision scientific computing.


l5 - implementation of solution of a system of linear equations Ax = b, for a given matrix of coefficients A, and the vector of right-hand sides b. The matrix A is sparse, i.e. having a large number of zeros, and blocky.

The study provides an effective way of solving the above problem, taking into account specific structure of matrix A, i.e. its rarity, regularity of occurrence of zero elements and non-zero. It allows to reduce the complexity from to O(n) compared to the standard methods of solving such systems of equations.

In folder l5 there are files:
- l5program.jl - to run the program
- complexityCalculations - to test the program

There are also modules:
- matrixgen - to generate sparse matrix
- Fileutils - to read/save date(matrix, vector) in a file
- Blocksys - to solve equation Ax = b in different ways
- Blocksmatrix - structure of sparse matrix, matrix search functions
