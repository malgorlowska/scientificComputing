# Author: Małgorzata Orłowska

push!(LOAD_PATH,".")

using Blocksmatrix
using Blocksys
using Fileutils
using SparseArrays

global m = readMatrix("A.txt")
println("block size ", m.blockSize)
#println(m.matrix)

global vectorSize, b = readRightVector("b.txt")
println("vector size ", vectorSize)
#println(b)

#zad1 bez wyboru elementu głównego
# obliczenie x dla danego A i b
#x = basicGauss(m,b)
#saveToFile("l5wynikA1.txt", x)

#zad1 bez wyboru elementu głównego
# obliczenie b na podstawie A i dla x=[1,..,1]
# i porównanie go z wyznaczonym b
newVector = rightSideVector(m)
x = basicGauss(m,newVector)
saveToFileWithRelativeError("l5wynikA2.txt", x)

#zad1 z częściowym wyborem elementu głównego
# obliczenie x dla danego A i b
#x = GaussWithPivot(m,b)
#saveToFile("l5wynik3.txt", x)

#zad1 z częściowym wyborem elementu głównego
# obliczenie b na podstawie A i dla x=[1,..,1]
# i porównanie go z wyznaczonym b
#newVector = rightSideVector(m)
#x = GaussWithPivot(m,newVector)
#saveToFileWithRelativeError("l5wynikA4.txt", x)

#zad2 bez wyboru elementu głównego
#wyznaczenie rozkładu LU
#LU = findLU(m)

#zad2 z częściowym wyborem elementu głównego
#wyznaczenie rozkładu LU
#LU = findLUWithPivot(m)

#zad3 bez wyboru elementu głównego
#rozwiązanie układu z wyznaczonym wcześniej LU

#newVector = rightSideVector(m)
#LU = findLU(m)
#x = solveWithLU(LU, newVector)
#saveToFileWithRelativeError("l5wynikA5.txt", x)

#zad3 z częściowym wyborem elementu głównego
#rozwiązanie układu z wyznaczonym wcześniej LU
#newVector = rightSideVector(m)
#p,LU = findLUWithPivot(m)
#x = solveWithLUWithPivot(LU, newVector, p)
#saveToFileWithRelativeError("l5wynikA6.txt", x)

#=statistics = @timed solveLUWithPivot(m,b)
#println("time")
println(statistics.time)=#
