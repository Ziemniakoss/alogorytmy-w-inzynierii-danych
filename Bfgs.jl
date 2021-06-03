using Zygote
include("./LineSearch.jl")
include("./MathWithPreallocations.jl")


"""

# Parametry

## Opcjonalne

- ϵ -> dokładność rozwiązania
- max_iter -> maksymalna liczba iteracji BFGS
- line_search_max_iter -> maksymalna liczba iteracji wybranego linesearcha

"""

function Bfgs(f, x::Vector{Float64}; max_iter::Int64=1000, ϵ::Float64=1e-6, line_search_max_iter::Int64=1000)
	B = rand(length(x), length(x)) # Tutaj lepsze zrobić lepsze wyszukanie
	∇f = f'(x)
	i = 0
	while i < max_iter
		p = B \ -∇f
		s = Armijo(f, x, p; maxIterations=line_search_max_iter) * p
		
		if s's < ϵ
			break
		end

		x .+= s
		∇fn = f'(x)
		y = ∇fn - ∇f
		∇f = ∇fn
		Bs = B * s
		B += y * y' ./ y's - Bs * Bs' ./ s'Bs
		i += 1
	end
	return x
end

# https://scicomp.stackexchange.com/questions/11323/effect-of-initial-guess-b-approximate-hessian-on-bfgs-algorithm
#
function bfgsWithPreAllocation(f, x::Vector{Float64}; max_iter::Int64=1000, ϵ::Float64=1e-10, maxLineSearchIterations::Int64=1000)
	∇f = f'(x)
	B .= ∇f /  #rand(length(x), length(x)) # I # Tutaj lepsze zrobić lepsze wyszukanie
	i = 0
	y = zeros(length(x))
	BS = zeros(length(x), 1)
	while i < max_iter
		p = B \ -∇f
		s = Armijo(f, x, p; maxIterations=maxLineSearchIterations) * p
		
		if s's < ϵ
			break
		end

		x += s
		∇fn = f'(x)
		sub!(∇fn, ∇f, y) # y = ∇fn - ∇f
		∇f = ∇fn
		mul!(BS, B, s)
		B .+= y * y' ./ y's - BS * BS' ./ s'BS
		i += 1
	end
	return x
end

