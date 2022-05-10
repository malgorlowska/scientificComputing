# Author: Małgorzata Orłowska
push!(LOAD_PATH,".")
module Blocksmatrix

using SparseArrays
export BMatrix, createMatrix, findNoZeroRow, findNoZeroColumn, findFirstColumn

#Structure of sparse block matrix
mutable struct BMatrix
    matrix::SparseMatrixCSC{Float64, Int}
    size::Int
    blockSize::Int
end

function Base.setindex!(A::BMatrix, v::Float64, i::Int, j::Int)
    A.matrix[i, j] = v
end

function Base.getindex(A::BMatrix, i::Int, j::Int)
    return A.matrix[i, j]
end

function createMatrix( matrixSize::Int, bSize::Int, sMatrix::SparseMatrixCSC{Float64, Int})
    matrix = sMatrix
    size = matrixSize
    blockSize = bSize

    return BMatrix(matrix, size, blockSize)
end

function findNoZeroRow(A::BMatrix, column::Int)# czy prościej??

    l = A.blockSize
    #println("reszta ", column%l)
    if mod(column, l) in 1:(l-2)
        return (floor(Int,(column/l))+1)*l
    elseif mod(column,l) == 0
        return column + l
    else
        return min(A.size,((floor(Int,(column/l))+1)*l+l))
    end

end

function findNoZeroColumn(A::BMatrix, row::Int)
    return min(A.size, row + A.blockSize)
end

function findFirstColumn(A::BMatrix, row::Int)
    l = A.blockSize
    #if mod(row, l) in 1:(l-2)
        newv = floor(Int,((row-1)/l))*l-1
        #println("nowa ", newv)
        return max(1,newv)

end

end
