# Definition and functions for double polynomial rings
# David Anderson, June 2025


"""
    drift_ring(n, m=0; coeff=ZZ, xname=:x, yname=:y) -> MPolyRing

Polynomial ring with variables x1..xn then y1..ym over `coeff` (default Nemo ZZ).
Use `extract_vars` to recover the x- and y-families.
"""
function drift_ring(n::Int, m::Int=0; coeff=ZZ, xname::Symbol=:x, yname::Symbol=:y)
    names = vcat(["$(xname)$(i)" for i in 1:n], ["$(yname)$(j)" for j in 1:m])
    R, _ = polynomial_ring(coeff, names)
    return R
end

# keep xy_ring for backward compatibility, but deprecate it in favor of drift_ring
@deprecate xy_ring(n::Int, m::Int=0; kwargs...) drift_ring(n, m; kwargs...)


"""
    extract_vars(R::MPolyRing; varname::Symbol) -> Vector

Return the generators of `R` whose printed name begins with `varname`, in the order
`gens(R)` returns them. Relies on the ring having been built with its generators in
index order (x1,x2,...,xn then y1,...,ym). Prefix match means `varname=:x` also catches
any `x`-prefixed family; fine for this package's rings.
"""
function extract_vars(R::MPolyRing; varname::Symbol)
  pre = string(varname)
  z = elem_type(R)[]
  for v in gens(R)
    startswith(string(v), pre) && push!(z, v)
  end
  return z
end
