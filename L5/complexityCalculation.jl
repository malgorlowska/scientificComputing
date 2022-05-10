# Author: Małgorzata Orłowska
push!(LOAD_PATH,".")

using Blocksmatrix
using Blocksys
using Fileutils
using matrixgen
using Plots

matrixSizes = [2000, 4000, 6000, 8000, 10000, 12000, 14000]
blockSize = 10
numberOfFunction = 5
namesAlgorithms = ["pierwsza", "druga", "trzecia", "czwarta", "piąta"]
solutions = [ (namesAlgorithms[i], Int64[], Float64[], Int64[]) for i in 1:numberOfFunction]

function f(A::BMatrix,b::Vector{Float64})
    return Array(A.matrix)\b
end

#for whatever reason re-using the same table resulted in having the
#same exact benchmark for all algorithms...
#3 hours spent on debugging and this seems to be a solution...?
times = Vector{Float64}(undef,size(matrixSizes))
memories = Vector{Int64}(undef,size(matrixSizes))

times2 = Vector{Float64}(undef,size(matrixSizes))
memories2 = Vector{Int64}(undef,size(matrixSizes))

times3 = Vector{Float64}(undef,size(matrixSizes))
memories3 = Vector{Int64}(undef,size(matrixSizes))

times4 = Vector{Float64}(undef,size(matrixSizes))
memories4 = Vector{Int64}(undef,size(matrixSizes))

times5 = Vector{Float64}(undef,size(matrixSizes))
memories5 = Vector{Int64}(undef,size(matrixSizes))

testData = Array{Tuple{BMatrix, Vector{Float64}}}(undef, size(matrixSizes))

#generating matrixes and right sides vectors
for(i, matrixSize) in enumerate(matrixSizes)
    blockmat(matrixSize, blockSize, 10.0, "generatedMatrix.txt")
    A = readMatrix("generatedMatrix.txt")
    b = rightSideVector(A)
    testData[i] = (A, b)
end

#first function
for(i, matrixSize) in enumerate(matrixSizes)
    A = deepcopy(testData[i][1])
    b = deepcopy(testData[i][2])
    statistics1 = @timed f(A,b)
    times[i] = statistics1.time
    memories[i] = statistics1.bytes
end

solutions[1] = ("pierwsza", matrixSizes, times, memories)

#second function
for(i, matrixSize) in enumerate(matrixSizes)
    A = deepcopy(testData[i][1])
    b = deepcopy(testData[i][2])
    statistics2 = @timed basicGauss(A,b)
    times2[i] = statistics2.time
    memories2[i] = statistics2.bytes
end

solutions[2] = ("druga", matrixSizes, times2, memories2)


#third function
for(i, matrixSize) in enumerate(matrixSizes)
    A = deepcopy(testData[i][1])
    b = deepcopy(testData[i][2])
    statistics3 = @timed GaussWithPivot(A,b)
    times3[i] = statistics3.time
    memories3[i] = statistics3.bytes
end

solutions[3] = ("trzecia", matrixSizes, times3, memories3)

#fourth function
for(i, matrixSize) in enumerate(matrixSizes)
    A = testData[i][1]
    b = testData[i][2]
    statistics = @timed solveLU(A,b)
    times4[i] = statistics.time
    memories4[i] = statistics.bytes
end

solutions[4] = ("czwarta", matrixSizes, times4, memories4)

#fifth function
for(i, matrixSize) in enumerate(matrixSizes)
    A = deepcopy(testData[i][1])
    b = deepcopy(testData[i][2])
    statistics = @timed solveLUWithPivot(A,b)
    times5[i] = statistics.time
    memories5[i] = statistics.bytes
end

solutions[5] = ("piąta", matrixSizes, times5, memories5)

println(solutions)
p = plot(matrixSizes, [solutions[1][3] solutions[2][3] solutions[3][3] solutions[4][3] solutions[5][3]], label=["function f" "Gauss" "GaussWithPivot" "LU" "LUWithPivot"], title = "time")
savefig("l5times5.png")
p = plot(matrixSizes, [solutions[1][4] solutions[2][4] solutions[3][4] solutions[4][4] solutions[5][4]], label=["function f" "Gauss" "GaussWithPivot" "LU" "LUWithPivot"], title = "memory")
savefig("l5memoriess5.png")
p = plot(matrixSizes, [ solutions[2][3] solutions[3][3] solutions[4][3] solutions[5][3]], label=["Gauss" "GaussWithPivot" "LU" "LUWithPivot"], title = "time")
savefig("l5times6.png")
p = plot(matrixSizes, [ solutions[2][4] solutions[3][4] solutions[4][4] solutions[5][4]], label=["Gauss" "GaussWithPivot" "LU" "LUWithPivot"], title = "memory")
savefig("l5memoriess4.png")
