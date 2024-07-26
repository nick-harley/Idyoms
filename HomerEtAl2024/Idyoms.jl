module Idyoms

export Bounded, Unbounded
export View
export A, B, C, D, X
export Interpolated, Backoff

export Distribution, Prediction
export infcontent, entropy

using Chakra


# ALPHABETS

export get_alphabet

function get_alphabet(v::Viewpoint{T},s::Vector{Chakra.Constituent})::Set{T} where T

    # RETURN THE VIEWPOINT ALPHABET OF A SEQUENCE
    
    return Set{T}(filter(e->e!=none,vp_map(v,s)))

end

function get_alphabet(v::Viewpoint{T},ss::Vector{Vector{Chakra.Constituent}})::Set{T} where T

    # RETURN THE VIEWPOINT ALPHABET OF A VECTOR OF SEQUENCES
    
    return union([get_alphabet(v,s) for s in ss]...)

end



### ORDER BOUNDS

abstract type OrderBound end

struct Bounded{h} <: OrderBound
    Bounded(h::Int) = h < 0 ? error("Can't have bound less than 0.") : new{h}()
end

struct Unbounded <: OrderBound end


# TYPE OF SEQUENCES

export Seq

Seq{S} = SubArray{S, 1, Vector{S}, Tuple{UnitRange{Int64}}, true}


# TRIM SEQUENCES

function trim(x::Seq{X},l::Int) where X    
    length(x) <= l && return x
    return @views x[end-l+1:end]
end

trim(p::Pair{Seq{X},Y},l::Int) where {X,Y} = trim(p[1],l)=>p[2]

function itrim(z::Seq{Z},l::Int) where Z
    length(z) <= l && return z
    return @views z[1:l]
end

itrim(p::Pair{Y,Seq{Z}},l::Int) where {Y,Z} = p[1]=>itrim(p[2],l)



# VIEWS


struct View{X,Y,Z}

    X::Vector{X}
    Y::Vector{Y}
    Z::Vector{Z}
    xn::Vector{Int}
    zn::Vector{Int}
    
    function View(elems::Vector{Tuple{Option{X},
                                      Option{Y},
                                      Option{Z}}}) where {X,Y,Z}

        # CONSTRUCT VIEW FROM SOURCE AND TARGET SEQUENCES
        
        x = [e[1] for e in elems] # sequences including none
        y = [e[2] for e in elems] 
        z = [e[3] for e in elems]
        
        xind = findall(e->e != none,x) # indices of attribtues 
        yind = findall(e->e != none,y)
        zind = findall(e->e != none,z)

        xns = [length(findall(x->x<i,xind)) for i in yind]
        zns = [length(findall(z->z<=i,zind))+1 for i in yind]
        
        x = x[xind] # sequences
        y = y[yind] 
        z = z[zind]
        
        return new{X,Y,Z}(x,y,z,xns,zns)
    end

    function View(seq::Vector,
                  vpx::Viewpoint{X},
                  vpy::Viewpoint{Y},
                  vpz::Viewpoint{Z}) where {X,Y,Z}

        # CONSTRUCT VIEW FROM SEQUENCE AND TWO VIEWPOINTS
        
        return View(Tuple{Option{X},
                          Option{Y},
                          Option{Z}}[zip(vp_map(vpx,seq),
                                         vp_map(vpy,seq),
                                         vp_map(vpz,seq))...])
    end
end

Base.length(v::View) = length(v.Y)

export getX, getY, getZ

getX(v::View,i::Int) = @views v.X[1:v.xn[i]]
getY(v::View,i::Int) = v.Y[i]
getZ(v::View,i::Int) = @views v.Z[v.zn[i]:end]


### NGRAMS

NMG{X,Y,Z} = Tuple{Seq{X},Y,Seq{Z}}

export nmg, nmgs, NMgs

nmg(v::View,i::Int) = getX(v,i) , getY(v,i) , getZ(v,i)
nmg(v::View,i::Int,n::Int,m::Int) = begin
    x , y , z = nmg(v,i)
    trim(x,n), y, itrim(z,m)
end

nmg(v::View,i::Int,::Bounded{n},::Bounded{m}) where {n,m} = nmg(v,i,n,m)
    
nmgs(v::View,n::Int,m::Int) = [nmg(v,i,n,m) for i in n+1:length(v)-m]

NMgs(v::View,N::Int,M::Int) = begin
    gs = vcat([nmgs(v,n,m) for n in 0:N for m in 1:M]...)
    l = length(v)
    edge_gs = [nmg(v,l,n,0) for n in 0:N]
    return vcat(gs,edge_gs)
end
    
NMgs(v::View) = NMgs(v,length(v)-1,length(v)-1)
NMgs(v::View,::Bounded{n},::Bounded{m}) where {n,m} = NMgs(v,n,m)
NMgs(v::View,::Unbounded,::Bounded{m}) where m = NMgs(v,length(v)-1,m)
NMgs(v::View,::Bounded{n},::Unbounded) where n = NMgs(v,n,length(v)-1)
NMgs(v::View,::Unbounded,::Unbounded) = NMgs(v)

NMgs(vs::Vector{View{X,Y,Z}},n,m) where {X,Y,Z} = vcat([NMgs(v,n,m) for v in vs]...)

Base.length(p::Pair{Seq{X},Y}) where {X,Y} = length(p[1])


# TODO: Which of these functions is correct?
# if a and b are 0, a/b is nan so return 0
# if a/b is inf should it return 0 or 1?

#non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) ? 0 : x)
#non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) || isinf(x) ? 0 : x)
non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) ? 0 : isinf(x) ? 1 : x)




using DataStructures

# TYPE OF COUNTERS

Counter{T} = DefaultDict{T,Int,Int}    

emptycounter(T::DataType) = Counter{T}(0)

function incrementcounter(c::Counter{T},nxt::T) where T

    # INCREMENT COUNTER
    
    c[nxt] += 1
end

function count(c::Counter{T},nxt::T) where {T}

    # OCCURRENCE COUNT OF NXT IN A COUNTER
    
    return Base.get(c,nxt,0)
end

function count(c::Counter{T},nxt::T,k::Int) where {T}

    # OCCURRENCE COUNT OF NXT WITH INITIAL COUNT K
    
    x = count(c,nxt)
    x == 0 && return 0
    return x + k
end

function sumcount(c::Counter{T},k::Int,ex::Set{T}) where T

    # TOTAL OCCURRENCE COUNT NOT INCLUDING EX
    
    return sum([count(c,e,k) for e in keys(c) if !(e in ex)])
end

function symset(c::Counter{T}) where T

    # SET OF SYMBOLS COUNTED
    
    return Set{T}(keys(c))
end

function symcount(c::Counter{T}) where T

    # NUMBER OF SYMBOLS COUNTED
    
    return length(c)
end

function nsymset(n::Int, c::Counter{T}) where T

    # SET OF SYMBOLS COUNTED N TIMES
    
    return Set{T}([e for e in keys(c) if count(c,e) == n])
end

function nsymcount(n::Int, c::Counter{T}) where T

    # NUMBER OF SYMBOLS COUNTED N TIMES
    
    return length(nsymset(n,c))
end



### ESCAPE METHODS

abstract type Escape end

struct A <: Escape end
struct B <: Escape end
struct C <: Escape end
struct D <: Escape end
struct X <: Escape end

ppmk(::A) = 0
ppmk(::B) = -1
ppmk(::C) = 0
ppmk(::D) = -0.5
ppmk(::X) = 0

function typecount(::A, c::Counter{T}) where T
    return 1
end

function typecount(::B, c::Counter{T}) where T
    return symcount(c)
end

function typecount(::C, c::Counter{T}) where T
    return symcount(c)
end

function typecount(::D, c::Counter{T}) where T
    return symcount(c)/2
end

function typecount(::X, c::Counter{T}) where T
    return nsymcount(1,c)+1
end


### TALLIES

export TALLY

struct TALLY{X,Y,Z}
    
    X::Dict{Seq{X},Int}
    Y::Dict{Y,Int}
    Z::Dict{Seq{Z},Int}
    
    counts::Array{Int,3}
    
    function TALLY(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

        aX::Vector{Seq{X}} = unique([g[1] for g in gs])
        aY::Vector{Y} = unique([g[2] for g in gs])
        aZ::Vector{Seq{Z}} = unique([g[3] for g in gs])
    
        iX = Dict([x=>i for (i,x) in enumerate(aX)])
        iY = Dict([y=>i for (i,y) in enumerate(aY)])
        iZ = Dict([z=>i for (i,z) in enumerate(aZ)])
    
        counts = zeros(Int,length(aX),length(aY),length(aZ))
        
        for g in gs
            x, y, z = g
            counts[iX[x],iY[y],iZ[z]] += 1
        end

        new{X,Y,Z}(iX,iY,iZ,counts)
    end
end

export getCount, counter

function getCount(t::TALLY{X,Y,Z},g::NMG{X,Y,Z}) where {X,Y,Z}
    x,y,z = g
    !haskey(t.X,x) && return 0
    !haskey(t.Y,y) && return 0
    !haskey(t.Z,z) && return 0
    t.counts[t.X[x],t.Y[y],t.Z[z]]
end

function getCount(t::TALLY{X,Y,Z},x::Seq{X}) where {X,Y,Z}
    !haskey(t.X,x) && return 0
    sum(t.counts[t.X[x],:,:])
end

function getCount(t::TALLY{X,Y,Z},y::Y) where {X,Y,Z}
    !haskey(t.Y,y) && return 0
    sum(t.counts[:,t.Y[y],:])
end

function getCount(t::TALLY{X,Y,Z},z::Seq{Z}) where {X,Y,Z}
    !haskey(t.Z,z) && return 0
    sum(t.counts[:,:,t.Z[z]])
end

function getCount(t::TALLY{X,Y,Z},p::Pair{Seq{X},Y}) where {X,Y,Z}
    x,y = p
    !haskey(t.X,x) && return 0
    !haskey(t.Y,y) && return 0
    sum(t.counts[t.X[x],t.Y[y],:])
end

function getCount(t::TALLY{X,Y,Z},p::Pair{Seq{X},Seq{Z}}) where {X,Y,Z}
    x,z = p
    !haskey(t.X,x) && return 0
    !haskey(t.Z,z) && return 0
    sum(t.counts[t.X[x],:,t.Z[z]])
end

function getCount(t::TALLY{X,Y,Z},p::Pair{Y,Seq{Z}}) where {X,Y,Z}
    y,z = p
    !haskey(t.Y,y) && return 0
    !haskey(t.Z,z) && return 0
    sum(t.counts[:,t.Y[y],t.Z[z]])
end

function counter(t::TALLY{X,Y,Z},p::Pair{Seq{X},Y})::Counter{Seq{Z}} where {X,Y,Z}
    c = emptycounter(Seq{Z})
    [c[z] = getCount(t,p...,z) for z in keys(t.Z)]
    return c
end

function counter(t::TALLY{X,Y,Z},p::Seq{X},y::Y)::Counter{Seq{Z}} where {X,Y,Z}
    c = emptycounter(Seq{Z})
    [c[z] = getCount(t,p...,z) for z in keys(t.Z)]
    return c
end



Tally{S,T} = Dict{S,Counter{T}}

emptytally(S,T) = Tally{S,T}()

function incrementtally!(tally::Tally{S,T},s::S,t::T) where {S,T}
    
    Base.get!(tally,s,Counter{T}(0))[t] += 1
    
end

function maketally(gs::Vector{Pair{S,T}}) where {S,T}
    
    tally = emptytally(S,T)
    
    [incrementtally!(tally,p...) for p in gs]
    
    return tally
end

function tally_xY(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

    return maketally([g[1]=>g[2] for g in gs])
        
end

function tally_xyZ(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

    return maketally([(g[1]=>g[2])=>g[3] for g in gs])
    
end

function tally_yZ(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

    return @views maketally([[g[2]][1:end]=>g[3] for g in gs])

end

function tally_xZ(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

    return maketally([g[1]=>g[3] for g in gs])
    
end

function tally_zY(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

    return maketally([g[3]=>g[2] for g in gs])
    
end

function tally_xzY(gs::Vector{NMG{X,Y,Z}}) where {X,Y,Z}

    return maketally([(g[1]=>g[3])=>g[2] for g in gs])

end

# TODO: WHAT TO DO WITH UPDATE TALLY???

function updatetally(tally::Tally{S,T},s::S,t::T) where {S,T}
    
    for n in 0:length(s)
        incrementtally!(tally,s,t)
    end
end



### SMOOTHING METHODS

abstract type Smoothing end

struct Interpolated <: Smoothing end

struct Backoff <: Smoothing end



### ESTIMATE PROBABILITY

function estimate(s::S,
                  t::T,
                  tally::Tally{S,T},
                  seen::Set{T},
                  A::Set{T},
                  B::Backoff,
                  E::Escape,
                  U::Bool,
                  ex::Set{T} = Set{T}()) where {S,T}

    # ESTIMATE PROBABILITY WITH BACKOFF SMOOTHING

    n = length(s)
    
    k = ppmk(E)

    counter = Base.get(tally,s,Counter{T}(0))

    trans_count = count(counter,t,k)

    state_count = sumcount(counter,k,ex)

    type_count = typecount(E,counter)

    p = non_nan( trans_count,
                 state_count )
    
    w = non_nan( state_count,
                 state_count + type_count)

    p > 0 && return w * p, n
    
    e = 1 - w
    
    n == 0 && return e/(length(A)+1-length(seen)), -1         

    U && union!(ex,symset(counter))
    
    p, o = estimate(trim(s,n-1),t,tally,seen,A,B,E,U,ex)

    return e * p, o
end

function estimate(s::S,
                  t::T,
                  tally::Tally{S,T},
                  seen::Set{T},
                  A::Set{T},
                  B::Interpolated,
                  E::Escape,
                  U::Bool,
                  ex::Set{T}=Set{T}()) where {S,T}

    # ESTIMATE PROBABILITY WITH INTERPOLATED SMOOTHING
    
    n = length(s)
    
    k = ppmk(E)

    counter = Base.get(tally,s,Counter{T}(0))

    trans_count = count(counter,t,k)

    state_count = sumcount(counter,k,ex)

    type_count = typecount(E,counter)

    p = non_nan( trans_count,
                 state_count )
    
    w = non_nan( state_count,
                 state_count + type_count)

    e = 1 - w

    e == 0 && error("SOMETHING WRONG")
    
    if n==0

        p == 0 && return e / (length(A)+1-length(seen)), -1

        return w*p + (e/(length(A)+1-length(seen))), 0
    end
 
    U && union!(ex,symset(counter))
    
    _p, o = estimate(trim(s,n-1),t,tally,seen,A,B,E,U,ex)

    p == 0 && return e*_p, o
    
    return w*p + e*_p, n
end


### TYPE OF PROBABILITY DISTRIBUTIONS

Distribution{T} = Dict{T,Float64}

# GET ALPHABET OF DISTRIBUTION
domain(d::Distribution) = keys(d)

# NORMALISE A DISTRIBUTION
normalise!(d::Distribution) = (total = sum(values(d)); map!(x->x/total,values(d)))

# LOOKUP PROBABILITY OF ELEMENT
(d::Distribution{T})(e::T) where T = Base.get(d,e,0)

function mult(d1,d2::Distribution{T}) where T

    return Distribution{T}([t=>d1(t)*d2(t) for t in keys(d1)]...)
    
end

function logdiv(d1,d2::Distribution{T}) where T

    return Distribution{T}([t=>log(2,d1(t)/d2(t)) for t in keys(d2)]...)
    
end



function estimate_dist(s::S,
                       tally::Tally{S,T},
                       seen::Set{T},
                       A::Set{T},
                       B::Smoothing,
                       E::Escape,
                       U::Bool) where {S,T}

    # ESTIMATE THE PROBABILITY DISTRIBUTION OVER ALPHABET A
    
    dist = Distribution{T}()
    ords = Dict{T,Int}()

    for t in A
        p, o = estimate(s,t,tally,seen,A,B,E,U)
        dist[t] = p
        ords[t] = o
    end

    normalise!(dist)

    return dist, ords
end




using Random

function sample(dist::Distribution{T}) where T

    # RETURN A RANDOM ELEMENT OF DISTRIBUTION
    
    ps = shuffle(collect(dist))
    a = first.(ps)
    p = cumsum(last.(ps))
    n = rand(Float64)
    q = findfirst(q->n<=q,p)
    return first(ps[q])
end








function select_order(ctx::Vector{S},
                      tally::Tally{S,T},
                      seen::Set{T},
                      O::Bounded{h}) where {S,T,h}

    # RETURN BOUNDED CONTEXT
    
    length(ctx) <= h && return @views ctx[1:end]
    return @views ctx[end-h+1:end]
end

function select_order(ctx::Vector{S},
                      tally::Tally{S,T},
                      seen::Set{T},
                      O::Unbounded) where {S,T}

    # RETURN UNBOUNDED CONTEXT
    
    length(ctx) == 0 && return @views ctx[1:end]

    for l in 0:length(ctx)

        @views _ctx = ctx[end-l+1:end]

        c = Base.get(tally,_ctx,Counter{T}(0))

        tc = symcount(c)

        # RETURN SHORTEST DETERMINISTIC CONTEXT
        tc == 1 && return _ctx 

        # OR, RETURN THE LONGEST MATCHING CONTEXT
        tc == 0 && return @views _ctx[2:end]

    end
        
    return ctx

end

# TODO: Types here???

function select_order(v::View{X,Y,Z},
                      i::Int,
                      tally::Tally{Seq{X},Y},
                      seen::Set,
                      O::Bounded{h}) where {X,Y,Z,h}
    
    # RETURN BOUNDED NGRAM i FROM VIEW v
    
    return getngram(v,i,h)
end


# TODO: What does PPM* mean in contrast information

function select_order(v::View{X,Y,Z},
                      i::Int,
                      tally::Tally{Seq{X},Y},
                      seen::Set{Y},
                      O::Unbounded) where {X,Y,Z}

    # RETURN UNBOUNDED NGRAM AT i FROM VIEW v
    
    x, y, z = getngram(v,i)

    length(x) == 0 && return x, y, z

    for l in 0:length(x)

        @views _x = x[end-l+1:end]

        c = Base.get(tally,_x,Counter{Y}(0))

        tc = symcount(c)

        # RETURN SHORTEST DETERMINISTIC CONTEXT
        tc == 1 && return _x, y, z 

        # OR, RETURN LONGEST MATCHING CONTEXT
        tc == 0 && return @views _x[2:end], y, z

    end
        
    return x, y, z

end






### TYPE OF PREDICTIONS 

struct Prediction{T}
    symbol::T
    estimate::Float64
    order::Int
    distribution::Distribution{T}
end



# INFORMATION CONTENT OF AN ELEMENT FROM A DISTRIBUTION
infcontent(d::Distribution{T},e::T) where T = - log(2,d(e))
infcontent(x::Prediction) = infcontent(x.distribution,x.symbol)

# ENTROPY OF A DISTRIBUTION
entropy(d::Distribution) = sum([d(e)*infcontent(d,e) for e in domain(d)])
entropy(x::Prediction) = entropy(x.distribution)

# MAX ENTROPY OF A DISTRIBUTION
max_entropy(d::Distribution) = log(2,length(domain(d)))
max_entropy(x::Prediction) = max_entropy(x.distribution)

# RELATIVE ENTROPY OF A DISTRIBUTION
relative_entropy(d::Distribution) = (hm=max_entropy(d); hm>0 ? entropy(d) / hm : 1 )
relative_entropy(x::Prediction) = relative_entropy(x.distribution)

# CALCULATE WEIGHT OF A DISTRIBUTION
weight(d::Distribution,b::Int) = relative_entropy(d) ^ (-b)
weight(x::Prediction,b::Int) = weight(x.distribution,b)


function combine_dist(ds::Vector{Distribution{T}},
                      b::Int=0) where T

    # RETURN WEIGHTED COMBINATION OF DISTRIBUTIONS
    
    A = union([domain(d) for d in ds]...)

    ws = [weight(d,b) for d in ds]

    sum_weights = sum(ws)

    estimates = [(e=>sum([ws[m] * ds[m](e) for m in 1:length(ds)]) / sum_weights) for e in A]

    dist = Distribution{T}(estimates)

    normalise!(dist)
    
    return dist
    
end

function combine_dist(ps::Vector{Prediction{T}},
                      b::Int=0) where T

    # RETURN WEIGHTED COMBINATION OF PREDICTION DISTRIBUTIONS
    
    sym = ps[1].symbol

    ds = [p.distribution for p in ps]

    os = [p.order for p in ps]

    new_dist = combine_dist(ds,b)

    p = new_dist(sym)

    o = max(os...)

    Prediction(sym,p,o,new_dist)
end


function combine_predictions(ps::Vector{Vector{Prediction{T}}},b::Int) where T

    # COMBINE SEQUENCE PREDICTIONS
    
    return Prediction{T}[Idyoms.combine_dist([ep...],b) for ep in zip(ps...)]

end

function combine_predictions(pss::Vector{Vector{Vector{Prediction{T}}}},b::Int) where T

    # COMBINE A VECTOR OF SEQUENCE PREDICTIONS
    
    return Vector{Prediction{T}}[combine_predictions([ps...],b) for ps in zip(pss...)]
    
end



function mean_infcontent(ps::Vector{Prediction{T}}) where T

    # RETURN MEAN INFORMATION CONTENT OF A SEQUENCE OF PREDICTIONS
    
    return sum([infcontent(p) for p in ps])/length(ps)

end

function mean_infcontent(pss::Vector{Vector{Prediction{T}}}) where T

    # RETURN MEAN INFORMATION CONTENT OF A SET OF PREDICTION SEQUENCES
    
    return sum([mean_infcontent(ps) for ps in pss])/length(pss)

end


# PREDICT

function ppm(s::S,
             t::T,
             tally::Tally{S,T},
             seen::Set{T},
             A::Set{T},
             B::Smoothing,
             E::Escape,
             U::Bool) where {S,T}

    # RETURN PPM PREDICTION FOR AN ELEMENT IN A CONTEXT
    
    dist, ords = estimate_dist(s,tally,seen,A,B,E,U)
    
    return Prediction(nxt,dist[nxt],ords[nxt],dist)
    
end


end

