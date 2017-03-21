using JuMP
using GLPKMathProgInterface

#
#  @author Paweł Otrębski 188383
#
#  zadanie 1, zestaw 1.
#

#
# function to generate a hilbert matrix of given size
# size - size of the matrix: size x size
#
function Hilbert(size::Int)
  A = zeros(Float64,size,size)
  for i=1:size, j=1:size  A[i,j] = 1/(i+j-1) end
  return Matrix{Float64}(A)
end

#
# helper function to generate the coeffients for the vector c and b
# size - size of the vector
#
function generateCoefficients(size::Int)
  c = Vector{Float64}(size)
  for i = 1:size
    sum = 0
    for j = 1:size
      sum = sum+(1/(i+j-1))
    end
    c[i] = sum
  end
  return c
end

#
# A - linear contraints matrix
# b - solution Vector
# c - constants vector
#
function LinearProg(A::Matrix{Float64},b::Vector{Float64},c::Vector{Float64})
#
# m - number of rows, n number of columns
#
  (m,n) = size(A)
#
#  Define the solver
#
  model = Model(solver = GLPKSolverLP())
#
#  x - decision varibales
#
  @variable(model,x[1:n]>=0)
  @objective(model,:Min,vecdot(c,x))
  @constraint(model,A*x .==b)

  solution = solve(model, suppress_warnings=true)

  if solution==:Optimal
    return solution,getobjectivevalue(model),getvalue(x)
  end
  return solution,nothing,nothing
end

#
# minN - minimal number of N
# maxN - maximum number N should take
#
function runTests(minN::Int,maxN::Int)
  println("n ","condition ", "objectivevalue ", "x-value ", "RelativeError\n")
  for i=minN:maxN
    b = generateCoefficients(i)
    c = b
    A = Hilbert(i)
    condition = cond(A)
    trueAnswer = ones(i)
    s,m,x = LinearProg(A,b,c)
    relativeError = norm(x-trueAnswer)/norm(trueAnswer)
    println(i," ",condition," ",m," ",x," ",relativeError)
  end
end

try
  min = parse(ARGS[1])
  max = parse(ARGS[2])
  runTests(min,max)
catch
  println("ARGUMENT ERROR")
end
