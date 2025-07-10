# Tools for computing drift polynomials
# David Anderson, May 2025.

export unmarkconfig, markbox, drift_split
# TO DO: review function drift2bin(d::Drift, R::DoublePolyRing) which acts like bpd2bin, records a product of binomials x[i]+y[j]
# TO DO: review function drift_poly(d::Drift, R::DoublePolyRing) which produces the drift polynomial from the iterator
# TO DO: clarify the logic in markbox
# TO DO: clarify the logic and clean up tableau_components

#####
# Towards drift transition
#####


 
function drift_split( dc::Drift, i::Int, j::Int )
    # do split of dc at (i,j)

    # assumes input dc is marked

    dr=dc.mtx

    # only split at colliding box
  if dr[i,j]<100
    return dc
  end

    # find colliding box at end of column or row
  if dr[i+1,j]>=100
    return drift_split(dc,i+1,j)
  end

  if dr[i,j+1]>=100
    return drift_split(dc,i,j+1)
  end


  dr1=copy(dr)

  k=dr[i,j]-100
  dr1[i+k+1,j+k+1]=Int8(6)
  dr1[i,j]=k+10

  dc1=Drift(dr1)

  dr2=copy(dr)

 # dr2[i+k,j+k]=Int8(7)
  dr2[i+k+1,j+k+1]=Int8(7)
  dr2[i,j]=Int8(8)

# changed 6 to 9
 # dr2[i+k,j+k]=Int8(6)

  b=1
  while dr2[i,j+b]>=10
    if dr2[i+k+1,j+b+k+1]==6
      dr2[i+k+1,j+b+k+1]=10
    else
      dr2[i+k+1,j+b+k+1]=dr[i,j+b][1]-k-1 +10
    end
    dr2[i,j+b]=Int8(8)
    b+=1
  end

  b=1
  while dr2[i+b,j]>=10
    if dr2[i+b+k+1,j+k+1]==6
      dr2[i+b+k+1,j+k+1]=10
    else
      dr2[i+b+k+1,j+k+1]=dr[i+b,j][1]-k-1 +10
    end
    dr2[i+b,j]=Int8(8)
    b+=1
  end

  dc2=Drift(dr2)
  dc2=unmarkconfig(dc2)
  dc2=markconfig(dc2)

  return (dc1,dc2)

end


# better: label boxes with pairs (int,boolean), int=distance box can drift, boolean=collides or not
# probably should adapt the marking functions so they can work on marked diagrams, not just the unmarked drift diagrams
# no, this leads to errors -- for now an inefficient "unmarking"
# done?  check!

# new approach, keep everything in Int
# >=10 : marked 
# 10<= x < 100 : collides=false, blocked = false
# 100 <= x < 1000 : collides=true, blocked = false
# 1000 <= x < 2000 : collides=false, blocked = true
# >=2000 : collides=true, blocked = true
function markbox( dc::Drift, i::Int, j::Int )

  n,m = size(dc)

  dr=dc.mtx

  # quick returns
  #if dr[i,j]==7  # 7 should be obsolete
  #  return 10
  #end

  if dr[i,j]>0 && dr[i,j]<10 && dr[i,j]!=7
    return dr[i,j]
  end

  if i==n || j==m
    return 10
  end

  blocked=false
  if dr[i,j]==7 || dr[i,j]>=1000
    blocked=true
  end
  collides=false
    # this case should never occur for bpds
  if (dr[i+1,j]==8 || dr[i+1,j]==6) && (dr[i,j+1]==8 || dr[i,j+1]==6) && (dr[i+1,j+1]==0 || dr[i+1,j+1]>=10)
      res =blocked ? 2000 : 100
      return res
  end

  if (dr[i,j+1]==0 || dr[i,j+1]==7 || dr[i,j+1]>=10 ) && ( dr[i+1,j]>0 && dr[i+1,j]<10 && dr[i+1,j]!=7)
        k1=markbox( dc, i, j+1 )
        k2=0
   
        if (k1>=10 && k1<100)
            k1 = k1-10
        elseif (k1>=100 && k1<1000)
            k1 = k1-100
            collides=true
        elseif (k1>=1000 && k1<2000)
            k1 = k1-1000
        else
            k1=k1-2000
            collides=true
        end

        #  while i+k2+1<n && j+k2<n && dr[i+k2+1,j+k2+1]==8
        while k2<k1 && i+k2+1<n && j+k2<m
            if dr[i+k2+2,j+k2+1]==0 || dr[i+k2+2,j+k2+1]>=10
                res = blocked ? k2+2000 : k2+100
                return res
            end
            k2+=1
            if dr[i+k2+1,j+k2+1]==6
                res = blocked ? k2+1000 : k2+10
                return res
            end
        end
        res = blocked ? k2+1000 : k2+10
        return res
    end

  if ( dr[i+1,j]==0 || dr[i+1,j]==7 || dr[i+1,j]>=10 ) && dr[i,j+1]>0 && dr[i,j+1]<10 && dr[i+1,j]!=7
        k1=markbox( dc, i+1, j )
        k2=0

        if (k1>=10 && k1<100)
            k1 = k1-10
        elseif (k1>=100 && k1<1000)
            k1 = k1-100
            collides=true
        elseif (k1>=1000 && k1<2000)
            k1=k1-1000
        else
            k1=k1-2000
            collides=true
        end    
        #    while i+k2<n && j+k2+1<n && dr[i+k2+1,j+k2+1]==8
        while k2<k1 && i+k2<n && j+k2+1<m
            if dr[i+k2+1,j+k2+2]==0 || dr[i+k2+1,j+k2+2]>=10
              res = blocked ? k2+2000 : k2+100
              return res
            end
            k2+=1
            if dr[i+k2+1,j+k2+1]==6
                res = blocked ? k2+1000 : k2+10
                return res
            end
        end
        res = blocked ? k2+1000 : k2+10
        return res
  end

  if ( dr[i+1,j]==0 || dr[i+1,j]>=10 ) && ( dr[i,j+1]==0 || dr[i,j+1]>=10 )
        k1 = markbox( dc, i+1, j )
        if (k1>=10 && k1<100)
            k1 = k1-10
        elseif (k1>=100 && k1<1000)
            k1 = k1-100
        elseif (k1>=1000 && k1<2000)
            k1 = k1-1000
        else
            k1 = k1-2000
        end    

        k2 = markbox( dc, i, j+1 )
        if (k2>=10 && k2<100)
            k2 = k2-10
        elseif k2>=100
            k2 = k2-100
        end    

        res = blocked ? min(k1,k2)+1000 : min(k1,k2)+10
        return res
    end


  if can_drift( dc, i, j)
        if i<n-1
            aa = dr[i+2,j+1]
            if aa==0 || aa>=10 || aa==7
                res = blocked ? 2000 : 100
                return res
            end
        end

        if j<m-1
            aa = dr[i+1,j+2]
            if aa==0 || aa>=10 || aa==7
                res = blocked ? 2000 : 100
                return res
            end
        end

        dc1 = drift(dc,i,j)

        kk=markbox( dc1, i+1, j+1 )

        return kk+1
   end

   res = blocked ? 1000 : 10

  return res

end




function markconfig( dc::Drift )
# mark interfering boxes in drift config dc

  n,m=size(dc)
  mm = [Int(0) for i=1:n, j=1:m]

  for k=n+m-1:-1:m
      for i=(k-m+1):n
         mm[i,k-i+1] = markbox(dc, i, k-i+1)
      end
  end

  for k=m-1:-1:1
      for i=1:min(k,n)
         mm[i,k-i+1] = markbox(dc, i, k-i+1)
      end
  end

  return Drift(mm)

end


function markconfig( bpd::BPD )

  return markconfig( bpd2drift(bpd) )

end


function unmarkbox( dc::Drift, i, j)

  if dc.mtx[i,j]>=1000
    return Int8(7)
  elseif dc.mtx[i,j]>=10
    return Int8(0)
  else
    return dc.mtx[i,j]
  end

end

function unmarkconfig(dc::Drift)

  n,m=size(dc)

  mm = [unmarkbox(dc, i, j) for i=1:n, j=1:m]

  return Drift(mm)

end


# Turn this off until schur_poly is decoupled from SchubertPolynomials
#=

function schub_drifts( w::Vector{Int}, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by drift class formula

  fbpds = flat_bpds(w)

  pol = R.ring(0)

  for b in fbpds
    b=markconfig(b)
    pol = pol+dc2sd(b,R)
  end

  return pol

end


function dc2sd( dc::Drift, R::DoublePolyRing=xy_ring( size(dc)[1], size(dc)[2] )[1]  )
# drift configuration to s-polynomial
# must take marked configuration as input

  local n,m=size(dc)

  for k=(n+m):-1:2
    for i=maximum([1,k-m]):minimum([n,k-1])
      if isa( dc.mtx[i,k-i], Tuple ) && dc.mtx[i,k-i][2]
        (dc1,dc2)=drift_split( dc, i, k-i )
        return ( dc2sd( dc1, R ) + dc2sd( dc2, R ) )
      end
    end
  end

  sd = R.ring(1)

  tcomps = tableau_components(dc)

  for tt in tcomps
    sd = sd*schur_poly( tt[1], tt[2], R; mu=tt[3], xoffset=tt[4][1], yoffset=tt[4][2], rowmin=true )
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
        if !( (i,j) in corners) && isa(dc.mtx[i,j],Tuple) && ((i,j)==(1,1) || (i>1 && j>1 && !isa( dc.mtx[i-1,j], Tuple) && !isa( dc.mtx[i,j-1],Tuple )  )) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while i+s<=n && isa( dc.mtx[i+s,j], Tuple )

          local k=0
          while j+k<=m && isa( dc.mtx[i+s,j+k], Tuple )  # find SE boxes
            k +=1
          end          
          push!(la,k)


          local kk=0
          while j-kk-1>0 && isa( dc.mtx[i+s,j-kk-1], Tuple )  # find skew boxes
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
            push!(rr[end], dc.mtx[i-1+el,j-1+mm][1]+el )
          end
        end

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end

=#


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