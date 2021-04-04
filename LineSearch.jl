using Zygote
using LinearAlgebra

"""

# Paramtery

- f -> funkcja dla której szukamy minumum. Jeżeli ma więcej niż jeden parametr, musi je przyjmować jako wektor
- x -> wstępne x
- p -> 

## Opcjonalne

- α ->
- τ ->
- c ->
- max_iter -> maksymalna liczba iteracji, domyślnie 10 000

# Bibliografia

- https://en.wikipedia.org/wiki/Backtracking_line_search
"""
function Armijo(f, x, p; α=0.99, τ=0.95, c=0.3, max_iter=10000)::Float64
	m = f'(x)'p
	t = -c * m
	i = 0
	while i < max_iter
		if f(x) - f(x + α * p) >= α * t
			break
		end
		α *= τ
		i += 1
	end
	return α
end
