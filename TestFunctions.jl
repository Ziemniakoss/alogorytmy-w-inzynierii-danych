a = 1;
b = 100;

rosenbrock(x::Array{Float64, 1})::Float64 = (a - x[1])^2 + b * (x[2] - x[1]^2)^2
