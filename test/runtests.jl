using DriftPolynomials
using Test


if isempty(ARGS)
    tests = [ 
        "drifts.jl", 
       # "drift_polys.jl",
    	"draw_drifts.jl"
	]
else
    tests = ARGS
end

for test in tests
    include(test)
end

