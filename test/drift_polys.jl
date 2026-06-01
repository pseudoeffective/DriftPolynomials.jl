@testset "DriftPolys" begin

# --- single-configuration drift class -------------------------------------
# When nothing can drift, the drift polynomial is just the dominant monomial of
# the shape (row i contributes one factor of x_i per box).
R2 = drift_ring(2)
x = gens(R2)
d = partition2drift([2,1])
@test length(collect(drift_class(d))) == 1
@test drift_poly(d; ring=R2) == x[1]^2*x[2]

# --- genuine multi-term drift class ---------------------------------------
# partition [2,1] flagged by [2,3] has a 5-element drift class; the drift
# polynomial is the corresponding 5-term sum (hand-checked).
R3 = drift_ring(3)
x = gens(R3)
d = partition2drift([2,1],[2,3])
@test length(collect(drift_class(d))) == 5
@test drift_poly(d; ring=R3) ==
    x[1]^2*x[2] + x[1]^2*x[3] + x[1]*x[2]^2 + x[1]*x[2]*x[3] + x[2]^2*x[3]

# keyword and positional forms of partition2drift agree
@test partition2drift([2,1],[2,3]) == partition2drift([2,1]; ff=[2,3])

end
