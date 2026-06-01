# Tools for computing drift polynomials
# David Anderson, May 2025.

export unmarkconfig, markbox, drift_split
# TO DO: review function drift2bin(d::Drift, R::DoublePolyRing) which acts like bpd2bin, records a product of binomials x[i]+y[j]
# TO DO: review function drift_poly(d::Drift, R::DoublePolyRing) which produces the drift polynomial from the iterator
# TO DO: clarify the logic in markbox

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
      dr2[i+k+1,j+b+k+1]=drift_distance(dr[i,j+b])-k-1 +10
    end
    dr2[i,j+b]=Int8(8)
    b+=1
  end

  b=1
  while dr2[i+b,j]>=10
    if dr2[i+b+k+1,j+k+1]==6
      dr2[i+b+k+1,j+k+1]=10
    else
      dr2[i+b+k+1,j+k+1]=drift_distance(dr[i+b,j])-k-1 +10
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




"""
    markconfig(dc::Drift) -> Drift
    markconfig(bpd::BPD)   -> Drift

Mark the interfering boxes of a drift configuration `dc` (or the drift configuration of a
BPD), returning a new `Drift` whose entries use the integer marking encoding described at
the top of `drifts.jl` (values `>=10` record each box's drift distance, collision and
blocked status). Marking proceeds anti-diagonal by anti-diagonal via `markbox`.
"""
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
