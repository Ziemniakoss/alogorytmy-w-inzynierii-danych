#!/usr/bin/julia

include("./Bfgs.jl")

M = [1 .4
    .4  1]
g(x) = x'M*x + [10,20]'x
x = [420.0, 20.0]
Bfgs(g,x )
bfgsWithPreAllocation(g, [420.0,20.0])
@time result1 = Bfgs(g, [420.0,20.0])
@time result2 = bfgsWithPreAllocation(g, [420.0,20.0])
println(result1)
println(result2)
