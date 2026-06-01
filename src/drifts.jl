# Tools for working with drift configurations in Julia
# David Anderson, May 2025.


# TO DO: improve Base.show overload, better cross symbol
# export can_drift, can_undrift, drift, undrift
#########

struct Drift
    mtx::Matrix{<:Integer}
end

# Integer encoding of drift-configuration entries
# ------------------------------------------------
# 0..9 are raw tiles (0 empty box, 1 droplet, 2..6 pipe glyphs, 7 blocked droplet,
#       8 blank, 9 marked-blocked droplet).
# Values >= 10 are *marked* boxes produced by `markconfig`/`markbox`, encoded as a
# base symbol plus a range that records the box's drift distance and status:
#     10   <= x < 100   : marked, collides=false, blocked=false
#     100  <= x < 1000  : marked, collides=true,  blocked=false
#     1000 <= x < 2000  : marked, collides=false, blocked=true
#            x >= 2000   : marked, collides=true,  blocked=true



# Symbol to integer mapping
# to do: pare this down
const DRIFT_TO_INT = Dict(
    "O" => Int8(0),
    "+" => Int8(1),
    "." => Int8(6),
    "*" => Int8(7),
    "" => Int8(8)
)

function Drift(matrix::Matrix{String})
    int_matrix = map(x -> DRIFT_TO_INT[x], matrix)
    return Drift(int_matrix)
end



# convert integers back to symbols for display
function int_to_symbol_drift(i::Integer)
    symbols = ["\u25A1",  #0 "□ "
        "\u002B",    #1 \u271A "✚ "
        "\u256D\u2500",    #2 "╭─"
        "\u256F",         #3 "╯ "
        "\u2502",         #4 "│ "
        "\u2500\u2500",    #5 "──"
        'o', #   #6 \u2022 "• "
        "\u25A0",  #7 "■ " # '*', 
        "\u00B7",   #8 '·'
        "\u2715",  #9 '✕'
        "\u2B27"  #10 '⬧'  

    ]
    if i<10            # raw tile 0..9
      return symbols[i+1]
    elseif i<100       # marked, collides=false  (10..99)   -> □
      return symbols[1]
    elseif i<1000      # marked, collides=true   (100..999) -> ⬧
      return symbols[11]
    else               # marked, blocked         (>=1000)   -> ■
      return symbols[8]
    end
end

function int_to_symbol_drift(t::Tuple)
    return t[1]
end




# add method to display Drift
function Base.show(io::IO, dc::Drift)
    println(io)
    for i in 1:size(dc.mtx, 1)
        for j in 1:size(dc.mtx, 2)
            print(io, int_to_symbol_drift(dc.mtx[i, j]), " ")
        end
        println(io)
    end
end

# overload identity for Drift type
Base.:(==)(d1::Drift, d2::Drift) = d1.mtx == d2.mtx

# add method to Base.size for Drift
function Base.size(d::Drift)
  size(d.mtx)
end

# add method to Base.copy for Drift type
function Base.copy(dc::Drift)
  dcm = copy(dc.mtx)
  return Drift(dcm)
end



function can_drift(dc::Drift,i::Int,j::Int)
  if i==size(dc)[1] || j==size(dc)[2]
    return false
  end

 # check corners
    if dc.mtx[i,j] > 0 && dc.mtx[i,j]<10 && dc.mtx[i,j]!=7
      return(false)
    end

    if dc.mtx[i+1,j] == 0 || dc.mtx[i+1,j]>=10 || dc.mtx[i+1,j] == 1 || dc.mtx[i+1,j] == 7
      return(false)
    end

    if dc.mtx[i,j+1] == 0 || dc.mtx[i,j+1]>=10 || dc.mtx[i,j+1] == 1 || dc.mtx[i,j+1] == 7
      return(false)
    end

    if dc.mtx[i+1,j+1] == 0 || dc.mtx[i+1,j+1]>=10 || dc.mtx[i+1,j+1] == 1 || dc.mtx[i+1,j+1]==6 || dc.mtx[i+1,j+1]==7
      return(false)
    end

  return(true)

end

function drift( dc::Drift, i::Int, j::Int)
# perform drift move of dc at i,j if possible
# assume unmarked

  if can_drift(dc,i,j)

    dc2=copy(dc.mtx)

    dc2[i+1,j+1]=Int8(0)

    if dc2[i,j]==7
      dc2[i,j]=Int8(9)
    else
      dc2[i,j]=Int8(8)
    end

    return Drift(dc2)

  end

end 


function can_undrift(dc::Drift,i::Int,j::Int)
  # assume unmarked
  if i==1 || j==1
    return false
  end

 # check corners
    if dc.mtx[i,j] != 0
      return(false)
    end

    a = dc.mtx[i-1,j]
    if a == 0 || a == 1 || a==6 || a==7
      return(false)
    end

    a = dc.mtx[i,j-1]
    if a == 0 || a == 1 || a==7
      return(false)
    end

    a = dc.mtx[i-1,j-1]
    if a == 0 || a == 1 || a==7 || a==6
      return(false)
    end

  return(true)

end

function undrift( dc::Drift, i::Int, j::Int)
# perform undrift move of dc at i,j if possible
# assume unmarked

  if can_undrift(dc,i,j)

    dc2=copy(dc.mtx)

    if dc2[i-1,j-1]!=9
     dc2[i-1,j-1]=Int8(0)
    else
      dc2[i-1,j-1]=Int8(7)
    end
    
    dc2[i,j]=Int8(8)
  
    return Drift(dc2)

  end

end 


function step_drifts(dc::Drift)
# produce all one-step drifts of dc
   local n,m=size(dc)

   local dfts = []

   for i=1:n-1
     for j=1:m-1
            if can_drift(dc,i,j)
              local dc2=drift(dc,i,j)
              push!(dfts,dc2)
       end
     end
   end

   return(dfts)
end



struct DriftIterator
  stack::Vector{Any}
  seen::Set{Matrix}
end


function DriftIterator(dc::Drift)
  # Initialize
  seen = Set([dc.mtx])
  dfts = step_drifts(dc)
  stack = [(dc, dfts)]
  return DriftIterator(stack,seen)
end


Base.eltype(::Type{DriftIterator}) = Drift


Base.IteratorSize(::Type{<:DriftIterator}) = Base.SizeUnknown()


function Base.iterate(iter::DriftIterator, state=nothing)

    while !isempty(iter.stack)
        current, drfts = pop!(iter.stack)

        unseen_drfts = filter( b -> !(b.mtx in iter.seen), drfts )

        for b in unseen_drfts
          push!(iter.seen, b.mtx)  # mark new drift as seen
          push!( iter.stack, (b, step_drifts(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end


function drift_class(dc::Drift)
# iterator of all diagrams in drift class of dc

  dc2 = nw_reset(dc)

  iter = DriftIterator(dc2)

#  empty!(hash_all_drifts_SE)  # clear the lookup table

  return iter

end


function isflat(dc::Drift)
  dcm = dc.mtx
  local n=size(dcm)[1]

  for i=2:n-1
    for j=2:n-1
      if dcm[i,j]==0 && dcm[i-1,j]!=0 && dcm[i,j-1]!=0 && dcm[i-1,j-1]==2
         return(false)
      end
    end
  end
  return(true)
end


function nw_reset(dc::Drift)
# returns flat diagram in drift class of dc
   n,m=size(dc)

   for i=2:n
     for j=2:m
       if can_undrift(dc,i,j)
         dc2=undrift(dc,i,j)
         return nw_reset(dc2)
       end
     end
   end

   return(dc)   

end   


function se_reset(dc::Drift)
# returns sharp diagram in drift class of dc
   n,m=size(dc)

   for i=1:n-1
     for j=1:m-1

       if can_drift(dc,i,j)
         dc2=drift(dc,i,j)
         return se_reset(dc2)
       end
     end
   end

   return(dc)   

end

#=
function can_cancel(dc::Drift, i::Int,j::Int)

  if i>=size(dc)[1] || j>=size(dc)[2] return false end
 
  if dc.mtx[i,j]!=0 return false end

  if can_drift(dc,i,j) return can_cancel(drift(dc,i,j),i+1,j+1) end
 
  if dc.mtx[i+1,j+1]!=1 return false end
  if (dc.mtx[i+1,j] in [0,1]) || (dc.mtx[i,j+1] in [0,1]) return false end

  return true

end


function cancel_drift!(dc::Drift, i::Int, j::Int)

  if !can_cancel(dc,i,j) return dc end

  dc.mtx[i,j]=Int8(8)
  for s=1:min(size(dc)[1]i,size(dc)[2]-j)
    if dc.mtx[i+s,j+s]==1
      dc.mtx[i+s,j+s]=Int8(8)
      return dc
    end
  end

end

function cancel_all_drift(dc::Drift)

   n,m=size(dc)
   if n!=m throw(error("Drift configuration must be square")) end
  
  if !all(x -> x == 0 || x == 1 || x==8, dc.mtx)
    return dc
  end

  dc2=deepcopy(dc)

  while countboxes(dc2)>0
    ctbx=countboxes(dc2)
    for s=-n+1:0
      for t=0:s+n-1
        cancel_drift!(dc2,n-t,n+s-t)
      end
    end
    for s=1:n-1
      for t=0:n-s-1
        cancel_drift!(dc2,n-s-t,n-t)
      end
    end
    if ctbx==countboxes(dc2) return dc2 end
  end

  return dc2
end

function iscancelable(dc::Drift)
  n,m=size(dc)
  if n!=m return false end

  return (cancel_all_drift(dc)==empty_drift(n))
end
=#

function drift2rkmtx(dc::Drift)

  n,m=size(dc)
  rkmtx = fill(0,n,m)
  #nw corner
  if !(dc.mtx[1,1] in [0,1])
    rkmtx[1,1]=1
  elseif dc.mtx[1,1]==1
    rkmtx[1,1]=2
  end

  #first row
  for j=2:m
    if !(dc.mtx[1,j] in [0,1])
      dc.mtx[1,j-1]==0 ? rkmtx[1,j]=1 : rkmtx[1,j]=rkmtx[1,j-1]
    elseif dc.mtx[1,j]==1
      dc.mtx[1,j-1]==0 ? rkmtx[1,j]=2 : rkmtx[1,j]=rkmtx[1,j-1]+1
    end
  end
  #first column
  for i=2:n
    if !(dc.mtx[i,1] in [0,1])
      dc.mtx[i-1,1]==0 ? rkmtx[i,1]=1 : rkmtx[i,1]=rkmtx[i-1,1]
    elseif dc.mtx[i,1]==1
      dc.mtx[i-1,1]==0 ? rkmtx[i,1]=2 : rkmtx[i,1]=rkmtx[i-1,1]+1
    end
  end
  #rest
  for i=2:n
    for j=2:m
      if dc.mtx[i,j]==0
        rkmtx[i,j]=rkmtx[i-1,j-1]
      elseif dc.mtx[i,j]==1
        rkmtx[i,j]=rkmtx[i-1,j-1]+2
      else
        rkmtx[i,j]=rkmtx[i-1,j-1]+1
      end
    end
  end

  return rkmtx
end

function rkmtx2asm(mtx::Matrix{<:Integer})
  n,m=size(mtx)
  aa = fill(0,n,m)
  if n==0 || m==0 return mtx end
  #first row
  aa[1,1]=mtx[1,1]
  for j=2:m
    aa[1,j]=mtx[1,j]-mtx[1,j-1]
  end
  #first column
  for i=2:n
    aa[i,1]=mtx[i,1]-mtx[i-1,1]
  end
  #rest
  for i=2:n
    for j=2:m
      aa[i,j]=mtx[i-1,j-1]+mtx[i,j]-mtx[i-1,j]-mtx[i,j-1]    
    end
  end
  return aa
end

function is_rkmtx(mtx::Matrix{<:Integer})
  n,m=size(mtx)
  for i in 1:n-1, j in 1:m-1
    #check monotone
    if !(mtx[i+1,j]-mtx[i,j] in [0,1]) || !(mtx[i,j+1]-mtx[i,j] in [0,1])
      return false
    end
    #check diamond condition
    if !(mtx[i,j]+mtx[i+1,j+1]-mtx[i+1,j]-mtx[i,j+1] in [-1,0,1])
      return false
    end
  end
  # check monotone on last row and column
  for j=1:m-1
    if !(mtx[n,j+1]-mtx[n,j] in [0,1])
      return false
    end
  end
  for i=1:n-1
    if !(mtx[i+1,m]-mtx[i,m] in [0,1])
      return false
    end
  end
  return true
end

"""
    isintegrable(dc::Drift)

Check if drift configuration `dc` is integrable, i.e., comes from part of a BPD.
"""
function isintegrable(dc::Drift)
  return is_rkmtx( drift2rkmtx(dc) )
end




#####
# Generating drift configurations
#####

function empty_drift( n::Int, m::Int=n )

  mtx = fill( Int8(8),n,m )

  return Drift(mtx)

end

# random drift config of size n,m
function random_drift( n::Int, m::Int=n ; extended::Bool=false)
  if extended 
    possible_entries = Int8[0, 1, 8, 6, 7]
  else
    possible_entries = Int8[0, 1, 8]
  end

  return Drift(rand( possible_entries, n,m ))

end


"""
    partition2drift(lambda::Vector{Int}; ff=..., n=...)
    partition2drift(lambda::Vector{Int}, ff::Vector{Int}, n::Int=...)

Build the drift configuration with partition `lambda` placed in the NW corner, whose
boxes drift down to the rows bounded by the flag `ff`. `n` is the side length of the
(square) configuration. By default `ff` flags every part to the partition length and
`n` is large enough to contain the diagram.

The keyword form keeps `lambda` as the only positional argument (preferred); the
positional form `partition2drift(lambda, ff)` is retained for back-compatibility.
"""
function partition2drift( lambda::Vector{Int};
                          ff::Vector{Int}=fill(length(lambda),length(lambda)),
                          n::Int=maximum( [ maximum(ff), length(lambda), lambda[1] ] )  )
# drift config with lambda in NW corner

  local mtx = fill( Int8(8), n, n)

  for i=1:length(lambda)
    for j=1:lambda[i]
      mtx[i,j]=Int8(0)
    end
    if i<=length(ff) && ff[i]+1 <= n && lambda[i]+ff[i]+1-i <= n
      mtx[ff[i]+1,lambda[i]+ff[i]-i+1]=Int8(1)
    end
  end

  return Drift(mtx)

end

# positional forwarder for back-compatibility (e.g. `partition2drift(la, ff)`)
function partition2drift( lambda::Vector{Int}, ff::Vector{Int},
                          n::Int=maximum( [ maximum(ff), length(lambda), lambda[1] ] )  )
  return partition2drift( lambda; ff=ff, n=n )
end



"""
    bpd2drift(bpd::BPD) -> Drift

Convert a bumpless pipe dream `bpd` to its drift configuration: empty boxes (`0`) and
droplets (`1`) are preserved, every other tile becomes the blank drift entry `8`.
"""
function bpd2drift( bpd::BPD )
# generate drift config from BPD

  local n = size(bpd)

  local bpd2=copy(bpd.mtx)


  for i=1:n
    for j=1:n
      if !isa(bpd2[i,j],Tuple) && bpd2[i,j]!=0 && bpd2[i,j]!=1
        bpd2[i,j]=Int8(8)
      end
    end
  end


  return Drift(bpd2)

end

function countboxes( dc::Drift )

  ct = count(==(0),dc.mtx)
  return ct

end

