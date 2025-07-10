# Definition and functions for double polynomial rings
# David Anderson, June 2025


function xy_ring(xx::Vector{String},yy::Vector{String}; coef_ring::Ring=ZZ)

  return polynomial_ring(coef_ring,vcat(xx,yy))

end


function xy_ring(xx::Vector{String}; coef_ring::Ring=ZZ)
    return xy_ring(xx,String[];kwargs...)
end


function xy_ring(n::Int,m::Int; varnames::Tuple{Symbol,Symbol}=(:x,:y), coef_ring::Ring=ZZ)

  (x,y) = varnames
  local xvars = ["x$(i)" for i=1:n]
  local yvars = ["y$(i)" for i=1:m]

  return polynomial_ring(coef_ring,vcat(xvars,yvars))

end


function xy_ring(n::Int; kwargs...)
  return xy_ring(n,0;kwargs...)
end

function extract_vars(R::MPolyRing; varname::Symbol)
  z_variables = MPolyRingElem[]
  for v in gens(R)
    if startswith(string(v),"$varname")
      push!(z_variables,v)
    end
  end
  return z_variables
end
