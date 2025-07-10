# Tools for computing drift polynomials
# David Anderson, May 2025.


#####
# Drift polynomials
#####

function drift2bin( d::Drift; double::Bool=false, ring::MPolyRing=xy_ring( size(d)[1], size(d)[2] )[1] )
# product of binomials for d
# default is single polynomial, use version="d[ouble]" for double
  n,m=size(d)
  bin = ring(1)

  x = extract_vars(ring,varname=:x)
  aa = length(x)

  if double
    y = extract_vars(ring,varname=:y)
    bb = length(y)
  else
    bb = 0
  end

  for i=1:n
    for j=1:m

      if d.mtx[i,j]==0 || d.mtx[i,j]>=10
        p=ring(0)
        if i<=aa
          p=p+x[i]
        end
        if j<=bb
          p=p+y[j]
        end
        bin = bin*p
      end

    end
  end

  return bin

end


#function drift_poly( d::Drift; version::String="single", R::MPolyRing=xy_ring( size(d)[1], size(d)[2] )[1]  )
function drift_poly( d::Drift; ring::MPolyRing=xy_ring( size(d)[1], size(d)[2] )[1], kwargs...  )
# compute drift pol by iterator
  dc=drift_class(d)

  pol=ring(0)

  for dd in dc
    pol = pol+drift2bin(dd; ring=ring, kwargs...)
  end

  return(pol)

end


