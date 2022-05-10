# Author: Małgorzata Orłowska
push!(LOAD_PATH,".")

module Fileutils

using LinearAlgebra
using Blocksmatrix
using SparseArrays
using DelimitedFiles
export readMatrix, readRightVector, saveToFile, saveToFileWithRelativeError

#Functions reads the data from file and writes it to a sparse matrix.
#Input:
# fileName - file name to read
#Output:
# M - sparse matrix (BMatrix)
function readMatrix(fileName::String)

    open(fileName) do file
        firstLine = split(readline(file))
        size = parse(Int, firstLine[1])
        blockSize = parse(Int, firstLine[2])
        I = Vector{Int64}(undef,0)
        J = Vector{Int64}(undef,0)
        V = Vector{Float64}(undef,0)

        for ln in eachline(file)
            line = split(ln)
            push!(I,parse(Int, line[1]))
            push!(J,parse(Int, line[2]))
            push!(V,parse(Float64, line[3]))
        end
        M = sparse(I, J, V)
        return createMatrix(size, blockSize, M)
    end

end

#Functions reads the data from file and writes it to a vector.
#Input:
# fileName - file name to read
#Output:
# b - right sides vector
function readRightVector(fileName::String)
    open(fileName) do file
        vectorSize = parse(Int, readline(file))
        b = Vector{Float64}(undef, 0)

        for ln in eachline(file)
            push!(b,parse(Float64, ln))
        end

        return vectorSize, b
    end

end

#Function saves a vector in a file.
#Input:
# fileName - file name
# x - result of equation Ax = b
function saveToFile(fileName::String, x::Vector{Float64})
    writedlm(fileName, x)
end

#Function saves a relative error and a vector in a file.
#Input:
# fileName - file name
# x - result of equation Ax = b
function saveToFileWithRelativeError(fileName::String, x::Vector{Float64})
    a = ones(Float64,size(x,1))
    #println("size ", size(x,1))
    open(fileName, "w") do f
            var = (norm(x-a))/(norm(a))
            write(f, "$var")
            println(var)
            write(f, "\n\r")
            for i in 1:size(x,1)
                var = x[i]
                write(f,"$var \n")
        end
    end
end

end
