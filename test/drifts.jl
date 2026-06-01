@testset "Drifts" begin

# constructors

mtx = Matrix( Int8[ 0 0 2 ; 0 2 1 ; 2 1 1 ] )
mtxd = Matrix( Int8[ 0 0 8 ; 0 8 1 ; 8 1 1 ] )
b = BPD( mtx )
d = bpd2drift(b)

@test d.mtx == mtxd

mtx2 = Matrix( [ "O" "O" ""; "O" "" "+"; "" "+" "+" ] )
d2 = Drift( mtx2 )

@test d==d2



# iterators

w=[1,2,7,5,3,4,6]
bps=collect(all_bpds(w))
b=bps[57]
d=bpd2drift(b)
dc = collect(drift_class(d))

@test length(dc)==66

la=[3,2,1]
ff=[2,4,5]
d=partition2drift(la,ff)
dc=collect(drift_class(d))

@test length(dc)==44


# is_rkmtx: regression for the nested-loop fix (was `for (i,j) in (1:n-1,1:m-1)`,
# which never iterated the interior grid). A rank matrix from an integrable drift
# is accepted; a matrix violating the interior diamond condition is rejected.
@test is_rkmtx( drift2rkmtx( partition2drift([2,1]) ) )
@test !is_rkmtx( Matrix( [ 0 0 0 ; 0 0 0 ; 0 0 2 ] ) )


end
