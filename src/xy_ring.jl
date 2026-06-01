# Definition and functions for double polynomial rings
# David Anderson, June 2025


function xy_ring(xx::Vector{String}, yy::Vector{String}; coeff::Ring=ZZ)

  return polynomial_ring(coeff, vcat(xx, yy))

end


function xy_ring(xx::Vector{String}; kwargs...)
    return xy_ring(xx, String[]; kwargs...)
end


function xy_ring(n::Int, m::Int; varnames::Tuple{Symbol,Symbol}=(:x,:y), coeff::Ring=ZZ)

  (x,y) = varnames
  local xvars = ["$(x)$(i)" for i=1:n]
  local yvars = ["$(y)$(j)" for j=1:m]

  return polynomial_ring(coeff, vcat(xvars, yvars))

end

function xy_ring(n::Int; kwargs...)
  return xy_ring(n, 0; kwargs...)
end


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
