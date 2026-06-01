# Computing Schubert polynomials via the drift class formula
# David Anderson, May 2025.

export schub_drifts

# TO DO: clarify the logic and clean up tableau_components


# TO DO: make sure this works now that schur_poly is decoupled from SchubertPolynomials

function schub_drifts( w::Vector{Int}; double::Bool=true, ring::MPolyRing=drift_ring( max(length(w)-1,1), max(length(w)-1,1) ) )
# compute schubert pol by drift class formula


  fbpds = flat_bpds(w)

  pol = ring(0)

  for b in fbpds
    b=markconfig(b)
    pol = pol+dc2sd(b; double=double, ring=ring)
  end

  return pol

end


function dc2sd( dc::Drift; double::Bool=true, ring::MPolyRing=drift_ring( size(dc)[1], size(dc)[2] ) )
# drift configuration to s-polynomial
# must take marked configuration as input

  local n,m=size(dc)

  for k=(n+m):-1:2
    for i=maximum([1,k-m]):minimum([n,k-1])
      if is_marked( dc.mtx[i,k-i] ) && collides( dc.mtx[i,k-i] )
        (dc1,dc2)=drift_split( dc, i, k-i )
        return ( dc2sd( dc1;  double=double, ring=ring ) + dc2sd( dc2;  double=double, ring=ring ) )
      end
    end
  end

  sd = ring(1)

  tcomps = tableau_components(dc)

  for tt in tcomps
    sd = sd*schur_poly( tt[1], tt[2]; double=double, ring=ring, mu=tt[3], xoffset=tt[4][1], yoffset=tt[4][2], rowmin=true )
  end

  return sd

end




function tableau_components(dc::Drift)
# return labelled tableaux for a drift configuration dc
# must take marked configuration as input

  local n,m=size(dc)

  if !isflat(dc)
    return( tableau_components( nw_reset(dc) ) )
  end

  local lyds=Vector{Vector}([])

  local corners=Vector{Tuple{Int,Int}}([])

  for i=1:n
    for j=1:m
        if !( (i,j) in corners) && is_marked(dc.mtx[i,j]) && ((i,j)==(1,1) || (i>1 && j>1 && !is_marked( dc.mtx[i-1,j] ) && !is_marked( dc.mtx[i,j-1] )  )) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while i+s<=n && is_marked( dc.mtx[i+s,j] )

          local k=0
          while j+k<=m && is_marked( dc.mtx[i+s,j+k] )  # find SE boxes
            k +=1
          end
          push!(la,k)


          local kk=0
          while j-kk-1>0 && is_marked( dc.mtx[i+s,j-kk-1] )  # find skew boxes
            kk +=1
          end

          mu=mu+fill(kk,length(mu))
          la=la+fill(kk,length(la))
          push!(mu,0)

          if s>0 && i+s>1 && j-kk>1 && ( dc.mtx[i+s-1,j-kk-1]==1 || dc.mtx[i+s-1,j-kk-1]==9 || dc.mtx[i+s-1,j-kk-1]==6  )
            push!(corners,(i+s,j-kk) ) # record new corner
          end
          j=j-kk
          s +=1
        end

        rr=Vector{Vector{Int}}([])
        for el=1:length(la)
          push!(rr, fill(0,mu[el]) )
          for mm=mu[el]+1:la[el]
            push!(rr[end], drift_distance(dc.mtx[i-1+el,j-1+mm])+el )
          end
        end

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end


#=
# not used
# TO DO: bring all drift poly computations from schub_polys and put them here

function tabcomps(bpd)
# return labelled tableaux for a flat bpd

  local n=size(bpd)[1]

  if !isflat(bpd)
    return( tableau_components( flatten(bpd) ) )
  end

  local lyds=Vector{Vector}([])

  local corners=Vector{Tuple{Int,Int}}([])

  for i=1:n
    for j=1:n
      if !( (i,j) in corners) && bpd[i,j]==0 && ((i,j)==(1,1) || (i>1 && j>1 &&bpd[i-1,j-1]==1)) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while i+s<=n && bpd[i+s,j]==0

          local k=0
          while j+k<=n && bpd[i+s,j+k]==0  # find SE boxes
            k +=1
          end
          push!(la,k)

          local el=1
            while i+s+el<=n && j+k-1+el<=n && bpd[i+s+el,j+k-1+el]=="%" || bpd[i+s+el,j+k-1+el]=="|"  # || bpd[i+s+el,j+k-1+el]=="-"
              el +=1
            end
          push!(rr,el-1)


          local kk=0
          while j-kk-1>0 && bpd[i+s,j-kk-1]==0  # find skew boxes
            kk +=1
          end

          mu=mu+fill(kk,length(mu))
          la=la+fill(kk,length(la))
          push!(mu,0)

          if s>0 && i+s>1 && j-kk>1 && bpd[i+s-1,j-kk-1]==1
            push!(corners,(i+s,j-kk) ) # record new corner
          end
          j=j-kk
          s +=1
        end

#=
        for k=length(la):-1:2
          if la[k-1]==la[k]
            rr[k-1]=rr[k]
          end
        end
=#

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end

=#
