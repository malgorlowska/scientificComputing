# Author: Małgorzata Orłowska
push!(LOAD_PATH,".")

module Blocksys

using SparseArrays
using Blocksmatrix
export basicGauss, GaussWithPivot, findLU, findLUWithPivot, solveWithLU, solveWithLUWithPivot,rightSideVector, solveLU, solveLUWithPivot

# Functions solves equation Ax = b
# where A - matrix, b - right sides vector

# Function implements basic Gauss algorithm(without main element selection)
# optimized for a specific matrix.
#Input:
# A - sparse matrix
# b - right sides vector
#Output:
# x - vector, result of Ax = b
function basicGauss(A::BMatrix, b::Vector{Float64})
    #println(A.size)
     n = A.size
     println(A.matrix[1][1])

 #triangular matrix
    for k in 1:n-1
        for i in (k+1):findNoZeroRow(A,k)
            if isapprox(A[k, k], 0.0; atol=eps())
                println("element of matrix A at index (" , k, ", ",k," is close to zero.")
                return nothing
            end
            z = A[i, k]/A[k, k]
            A[i, k] = 0.0
            b[i] = b[i] - (z*b[k])
            #println("zero kolumn ", findNoZeroColumn(A,i))
            for j in (k+1):findNoZeroColumn(A,i)
                A[i, j] = A[i, j] - (z*A[k, j])
            end
            #println(A)
        end
    end
        #return A, b

#solve equation Ax = b
    x = Vector{Float64}(undef,n)
    x[n] = b[n]\A[n, n]
    for i in n-1:-1:1
        x[i] = b[i]
        for j in (i+1):findNoZeroColumn(A,i)
            x[i] = x[i] - (A[i, j]*x[j])
        end
        x[i] = x[i] / A[i, i]
    end
    return x

end


# Function implements Gauss algorithm with main element selection
# optimized for a specific matrix.
#Input:
# A - sparse matrix
# b - right sides vector
#Output:
# x - vector, result of Ax = b
function GaussWithPivot(A::BMatrix, b::Vector{Float64})

     n = A.size
     p = [1:n;]
 #triangular matrix
    for k in 1:n-1
        noZeroRow = findNoZeroRow(A,k)
        #finding max main row
        max = k
        for l in (k+1):noZeroRow
            if abs(A[p[l], k]) >= abs(A[p[max], k])
                max = l
            end
        end
        pom = p[k]
        p[k] = p[max]
        p[max] = pom

        for i in (k+1):noZeroRow
            #println("wiersz ",i," kolumna ", k, " ",findNoZeroRow(A,k))
            if isapprox(A[p[k], k], 0.0; atol=eps())
                println("element of matrix A at index (" , k, ", ",k," is close to zero.")
                return nothing
            end
            z = A[p[i], k]/A[p[k], k]
            #println("z ", z)
            A[p[i], k] = 0.0
            b[p[i]] = b[p[i]] - (z*b[p[k]])
            #zostaje zamieniony wiersz najdalej z niższego blocku, czyli wierz moze się
            # wydłużyć maksymalnie o A.blockSize
            #println("zero kolumn ", findNoZeroColumn(A,i + A.blockSize))
            for j in (k+1):findNoZeroColumn(A,i + A.blockSize)
                A[p[i], j] = A[p[i], j] - (z*A[p[k], j])
            end
            #println(A)
        end
    end
        #return A, b

#solve equation Ax = b
    x = Vector{Float64}(undef,n)
    x[n] = b[p[n]]\A[p[n], n]
    for i in n-1:-1:1
        x[i] = b[p[i]]
        #zostaje zamieniony wiersz najdalej z niższego blocku, czyli wierz moze się
        # wydłużyć maksymalnie o A.blockSize
        for j in (i+1):findNoZeroColumn(A,i + A.blockSize)
            x[i] = x[i] - (A[p[i], j]*x[j])
        end
        x[i] = x[i] / A[p[i], i]
    end
    return x
end

# Function implements LU decomposition
# optimized for a specific matrix.
#Input:
# A - sparse matrix
#Output:
# LU - U upper triangular matrix
# L - lower triangular matrix
function findLU(A::BMatrix)

     n = A.size

 #triangular matrix
    for k in 1:n-1
        for i in (k+1):findNoZeroRow(A,k)
            if isapprox(A[k, k], 0.0; atol=eps())
                println("element of matrix A at index (" , k, ", ",k," is close to zero.")
                return nothing
            end

            z = A[i, k]/A[k, k]
            A[i, k] = z

            for j in (k+1):findNoZeroColumn(A,i)
                A[i, j] = A[i, j] - (z*A[k, j])
            end
        end
    end

    return A
end


# Function implements LU decomposition
# with main element selection
# optimized for a specific matrix.
#Input:
# A - sparse matrix
#Output:
# p - vector of permutation
# LU - U upper triangular matrix
# L - lower triangular matrix
function findLUWithPivot(A::BMatrix)

     n = A.size
     p = [1:n;]
 #triangular matrix
    for k in 1:n-1
        noZeroRow = findNoZeroRow(A,k)
        #finding max main row
        max = k
        for l in (k+1):noZeroRow
            if abs(A[p[l], k]) >= abs(A[p[max], k])
                max = l
            end
        end
        pom = p[k]
        p[k] = p[max]
        p[max] = pom

        for i in (k+1):noZeroRow
            if isapprox(A[p[k], k], 0.0; atol=eps())
                println("element of matrix A at index (" , k, ", ",k," is close to zero.")
                return nothing
            end

            z = A[p[i], k]/A[p[k], k]
            A[p[i], k] = z

            for j in (k+1):findNoZeroColumn(A,i + A.blockSize)
                A[p[i], j] = A[p[i], j] - (z*A[p[k], j])
            end
        end
    end

    return p, A

end

#Function solves equation Ax = b using LU distribution,
#in order to do that solves equations
#Ly = b i Ux = y
#Input:
# LU - distribution of matrix A
# b - right sides vector
#Output:
# x - vector, result of Ax = b
function solveWithLU(LU::BMatrix, b::Vector{Float64})
    n = LU.size

    #solve Ly = b z macierza trójkątną dolną
    for k in 1:(n-1)
        for i in (k+1):findNoZeroRow(LU, k)
            b[i] = b[i] - (b[k]*LU[i, k])
        end
    end
    #solve Ux = y z macierzą trójkątną górną
    x = Vector{Float64}(undef,n)
    x[n] = b[n]/LU[n, n]
    for i in n-1:-1:1
        x[i] = b[i]
        for j in (i+1):findNoZeroColumn(LU,i)
            x[i] = x[i] - (LU[i, j]*x[j])
        end
        x[i] = x[i] / LU[i, i]
    end
    return x

end

#Function solves equation Ax = b using LU distribution,
#Input:
# A - sparse matrix
# b - right sides vector
#Output:
# x - vector, result of Ax = b
function solveLU(A::BMatrix, b::Vector{Float64})
    LU = findLU(A)
    x = solveWithLU(LU, b)
    return x

end

#Function solves equation Ax = b using LU distribution,
# with main element selection,
#in order to do that solves equations
#Ly = b i Ux = y
#Input:
# LU - distribution of matrix A
# b - right sides vector
#Output:
# x - vector, result of Ax = b
function solveWithLUWithPivot(LU::BMatrix, b::Vector{Float64}, p::Vector{Int})
    n = LU.size

    #solve Ly = b z macierza trójkątną dolną, tu optymalizacja??
    for k in 1:(n-1)
        for i in (k+1):findNoZeroRow(LU, k)
            b[p[i]] = b[p[i]] - (b[p[k]]*LU[p[i], k])
        end
    end

    #solve Ux = y z macierzą trójkątną górną
    x = Vector{Float64}(undef,n)
    x[n] = b[p[n]]/LU[p[n], n]
    for i in n-1:-1:1
        x[i] = b[p[i]]
        for j in (i+1):findNoZeroColumn(LU,i + LU.blockSize)
            x[i] = x[i] - (LU[p[i], j]*x[j])
        end
        x[i] = x[i] / LU[p[i], i]
    end
    return x

end

#Function solves equation Ax = b using LU distribution,
# with main element selection,
#Input:
# A - sparse matrix
# b - right sides vector
#Output:
# x - vector, result of Ax = b
function solveLUWithPivot(A::BMatrix, b::Vector{Float64})
    p,LU = findLUWithPivot(A)
    x = solveWithLUWithPivot(LU, b, p)
    return x

end

#Function calculates right side vector
#Input:
# A - sparse matrix
#Output:
# b - right sides vector
function rightSideVector(A::BMatrix)
    n = A.size
    #x = ones(Float64, n)
    b = zeros(n)

    for j in 1:n
        first = findFirstColumn(A, j)
        last = findNoZeroColumn(A, j)
        #println(first, " ", last)
        for i in first:last
            b[j] += A[j, i]
        end
    end

return b

end


end
