module PPM

export A, B, C, D, X
export Interpolated, Backoff
export Bounded, Unbounded

export generate_ngrams, generate_hgrams, hgram_sequence

export Distribution, Prediction
export infcontent, entropy

export ppm, ppm_seq, ppm_seq_inc


using Chakra, Viewpoints


### ORDER BOUNDS

abstract type OrderBound end

struct Bounded{h} <: OrderBound
    Bounded(h::Int) = h < 0 ? error("Can't have bound less than 0.") : new{h}()
end

struct Unbounded <: OrderBound end


# TYPE OF CONTEXTS

Context{S} = SubArray{S, 1, Vector{S}, Tuple{UnitRange{Int64}}, true}

function trim(ctx::Context{S},l::Int) where {S,T}

    # TRIM CONTEXT TO MAX LENGTH L
    
    length(ctx) <= l && return ctx
    return @views ctx[end-l+1:end]
end

struct View{S,T}

    # VIEW OF A SEQUENCE
    
    target::Vector{Option{T}}
    source::Vector{Option{S}}
    targetindex::Vector{Int}
    sourceindex::Vector{Int}
    targetelements::Vector{T}
    sourceelements::Vector{S}

    function View(elems::Vector{Tuple{Option{S},
                                      Option{T}}}) where {S,T}

        # CONSTRUCT VIEW FROM SOURCE AND TARGET SEQUENCES
        
        t = last.(elems) # target sequence
        s = first.(elems) # source sequence
        tind = findall(x->x != none,t) # indices of target attribtues 
        sind = findall(x->x != none,s) # indices of source attribtues
        telems = t[tind] # target attribtues
        selems = s[sind] # source attribtues

        return new{S,T}(t,s,tind,sind,telems,selems)
    end

    function View(seq::Vector,
                  src::Viewpoint{S},
                  trg::Viewpoint{T}) where {S,T}

        # CONSTRUCT VIEW FROM SEQUENCE AND TWO VIEWPOINTS
        
        return View(Tuple{Option{S},Option{T}}[zip(vp_map(src,seq),vp_map(trg,seq))...])
    end
end

function Base.length(v::View)

    # LENGHT OF A VIEW
    
    length(v.targetindex)

end
    
function getnext(v::View{S,T},i::Int) where {S,T}

    # ITH TARGET ATTRIBTUE FROM VIEW
    
    return v.targetelements[i]
end

function getcontext(v::View{S,T},i::Int) where {S,T}

    # ITH CONTEXT FROM VIEW
    
    return @views v.sourceelements[1:findfirst(s->s==v.targetindex[i],v.sourceindex)-1]
end

function getngram(v::View{S,T},i::Int) where {S,T}

    # ITH TARGET-CONTEXT PAIR FROM VIEW 
    
    return getnext(v,i) => getcontext(v,i)
end

function getngram(v::View{S,T},i::Int,n::Int) where {S,T}

    # ITH N-GRAM FROM VIEW
    
    return getnext(v,i) => trim(getcontext(v,i),n-1)
end

function generate_ngrams(v::View{S,T},n::Int) where {S,T}

    # ALL N-GRAMS FROM VIEW
    
    return Pair{T,Context{S}}[getngram(v,i,n) for i in n:length(v)]
end

function generate_hgrams(v::View{S,T},h::Int)::Vector{Pair{T,Context{S}}} where {S,T}

    # ALL 1:H-GRAMS FROM VIEW
    
    vcat([generate_ngrams(v,n) for n in 1:h]...)
end

generate_hgrams(v::View,::Bounded{b}) where b = generate_hgrams(v,b+1)

generate_hgrams(v::View) = generate_hgrams(v,length(v)) 

generate_hgrams(v::View,::Unbounded) = generate_hgrams(v)





non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) ? 0 : x)





using DataStructures

# TYPE OF COUNTERS

Counter{T} = DefaultDict{T,Int,Int}    

function emptycounter(T::DataType)

    # EMPTY COUNTER

    Counter{T}(0)
end

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
    
    return Set(keys(c))
end

function symcount(c::Counter{T}) where T

    # NUMBER OF SYMBOLS COUNTED
    
    return length(c)
end

function nsymset(n::Int, c::Counter{T}) where T

    # SET OF SYMBOLS COUNTED N TIMES
    
    return Set([e for e in keys(c) if count(c,e) == n])
end

function nsymcount(n::Int, c::Counter{T}) where {S,T}

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
    return nsymcount(1,ctx,tally,seen)+1
end



# The TYPE OF TALLYS

Tally{S,T} = Dict{Context{S},Counter{T}}

emptytally(S,T) = Tally{S,T}()

function incrementtally(tally::Tally{S,T},nxt::T,ctx::Context{S}) where {S,T}
    Base.get!(tally,ctx,Counter{T}(0))[nxt] += 1
end

function maketally(gs::Vector{Pair{T,Context{S}}}) where {S,T}
    tally = emptytally(S,T)
    for p in gs
        incrementtally(tally,p...)
    end
    return tally
end

function updatetally(tally::Tally{S,T},nxt::T,ctx::Context{S}) where {S,T}
    for n in 0:length(ctx)
        incrementtally(tally,nxt,trim(ctx,n))
    end
end








### SMOOTHING METHODS

abstract type Smoothing end

struct Interpolated <: Smoothing end

struct Backoff <: Smoothing end


function estimate(nxt::T,
                  ctx::Context{S},
                  tally::Tally{S,T},
                  seen::Set{T},
                  A::Set{T},
                  B::Backoff,
                  E::Escape,
                  U::Bool,
                  ex::Set{T} = Set{T}()) where {S,T}

    # ESTIMATE PROBABILITY WITH BACKOFF SMOOTHING

    n = length(ctx)
    
    k = ppmk(E)

    counter = Base.get(tally,ctx,Counter{T}(0))

    trans_count = count(counter,nxt,k)

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
    
    @views p, o = estimate(nxt,ctx[2:end],tally,seen,A,B,E,U,ex)

    return e * p, o
end

# PPM INTERPOLATED
function estimate(nxt::T,
                  ctx::Context{S},
                  tally::Tally{S,T},
                  seen::Set{T},
                  A::Set{T},
                  B::Interpolated,
                  E::Escape,
                  U::Bool,
                  ex::Set{T}=Set{T}()) where {S,T}

    n = length(ctx)
    
    k = ppmk(E)

    counter = Base.get(tally,ctx,Counter{T}(0))

    trans_count = count(counter,nxt,k)

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
    
    @views _p, o = estimate(nxt,ctx[2:end],tally,seen,A,B,E,U,ex)

    p == 0 && return e*_p, o
    
    return w*p + e*_p, n
end


function select_order(v::View{S,T},
                      i::Int,
                      tally::Tally{S,T},
                      seen::Set{T},
                      O::Bounded{h}) where {S,T,h}
    return getngram(v,i,h)
end


function select_order(v::View{S,T},
                      i::Int,
                      tally::Tally{S,T},
                      seen::Set{T},
                      O::Unbounded) where {S,T}

    nxt, ctx = getngram(v,i)

    length(ctx) == 0 && return nxt, ctx

    for l in 0:length(ctx)

        @views _ctx = ctx[end-l+1:end]

        c = Base.get(tally,_ctx,Counter{T}(0))

        tc = symcount(c)

        tc == 1 && return nxt,_ctx 

        tc == 0 && return @views nxt, _ctx[2:end]

    end
        
    return nxt, ctx

end




Distribution{T} = Dict{T,Float64}

struct Prediction{T}
    symbol::T
    estimate::Float64
    order::Int
    distribution::Distribution{T}
end

domain(d::Distribution) = keys(d)

normalise!(d::Distribution) = (total = sum(values(d)); map!(x->x/total,values(d)))

(d::Distribution{T})(e::T) where T = Base.get(d,e,0)

infcontent(d::Distribution{T},e::T) where T = - log(2,d(e))
infcontent(x::Prediction) = infcontent(x.distribution,x.symbol)


entropy(d::Distribution) = sum([d(e)*infcontent(d,e) for e in domain(d)])
entropy(x::Prediction) = entropy(x.distribution)


max_entropy(d::Distribution) = log(2,length(domain(d)))
max_entropy(x::Prediction) = max_entropy(x.distribution)


relative_entropy(d::Distribution) = (hm=max_entropy(d); hm>0 ? entropy(d) / hm : 1 )
relative_entropy(x::Prediction) = relative_entropy(x.distribution)


weight(d::Distribution,b::Int) = relative_entropy(d) ^ (-b)
weight(x::Prediction,b::Int) = weight(x.distribution,b)


function combine(ds::Vector{Distribution{T}},
                 b::Int=0) where T

    A = union([domain(d) for d in ds]...)

    ws = [weight(d,b) for d in ds]

    sum_weights = sum(ws)

    estimates = [(e=>sum([ws[m] * ds[m](e) for m in 1:length(ds)]) / sum_weights) for e in A]

    dist = Distribution{T}(estimates)
    normalise!(dist)
    
    return dist
    
end

function combine(ps::Vector{Prediction{T}},
                 b::Int=0) where T
    
    sym = ps[1].symbol
    ds = [p.distribution for p in ps]
    os = [p.order for p in ps]
    new_dist = combine(ds,b)
    p = new_dist(sym)
    o = max(os...)
    Prediction(sym,p,o,new_dist)
end


function mean_infcontent(ps::Vector{Prediction{T}}) where T

    return sum([infcontent(p) for p in ps])/length(ps)

end

function mean_infcontent(pss::Vector{Vector{Prediction{T}}}) where T

    return sum([mean_infcontent(ps) for ps in pss])/length(pss)

end


# PREDICT

function ppm(nxt::T,
             ctx::Context{S},
             tally::Tally{S,T},
             seen::Set{T},
             A::Set{T},
             B::Smoothing,
             E::Escape,
             U::Bool) where {S,T}

    dist = Dict{T,Float64}()
    ords = Dict{T,Int}()
    for e in A
        p, o = estimate(e,ctx,tally,seen,A,B,E,U)
        dist[e] = p
        ords[e] = o
    end

    normalise!(dist)
    
    return Prediction(nxt,dist[nxt],ords[nxt],dist)
    
end

function ppm_seq(v::View{S,T},
                 tally::Tally{S,T},
                 seen::Set{T},
                 A::Set{T},
                 B::Smoothing,
                 E::Escape,
                 U::Bool,
                 O::OrderBound) where {S,T}

    # MODEL A SEQUENCE
    
    predictions = Prediction{T}[]
    
    for i in 1:length(v)
        nxt,ctx = select_order(v,i,tally,seen,O)
        pred = ppm(nxt,ctx,tally,seen,A,B,E,U)
        push!(predictions,pred)
        #push!(predictions,ppm(select_order(v,i,tally,seen,O)...,tally,seen,A,B,E,U))
    end
    
    return predictions
end


function ppm_seq_inc(v::View{S,T},
                     A::Set{T},
                     B::Smoothing,
                     E::Escape,
                     U::Bool,
                     O::OrderBound;
                     tally = emptytally(S,T),
                     seen = Set{T}()) where {S,T,h}

    # MODEL A SEQUENCE INCREMENTALLY
    
    predictions = Prediction{T}[]
    
    for i in 1:length(v)
        nxt, ctx = select_order(v,i,tally,seen,O)
        push!(predictions,ppm(nxt,ctx,tally,seen,A,B,E,U))
        updatetally(tally,getngram(v,i)...)
        push!(seen,nxt)
    end
    
    return predictions
    
end



# DATASET MANAGEMENT



function folddataset(data::Vector{View{S,T}},
                     nfolds::Int) where {S,T}

    # PARTITION DATASET INTO N FOLDS

    size = length(data)
    
    if !(0 < nfolds <= size)
        error("The number of folds must be less than the size of the data set size.")
    end

    # ASSIGN EACH VIEW TO A FOLD
    folds = Int[Int(round(i/nfolds % 1 *nfolds)) for i in 1:size]
    
    folds[findall(x->x==0,folds)] .= nfolds

    folds = sort(folds)

    # CREATE A TRAINING SET FOR EACH FOLD
    training_sets = Vector{View{S,T}}[data[findall(x->x != i,folds)] for i in 1:nfolds]
    
    return folds, training_sets

end

function train(data::Vector{View{S,T}},
               O::OrderBound) where {S,T}

    # CREATE A TALLY FROM A SET OF VIEWS
    
    gs = vcat([generate_hgrams(ex,O) for ex in data]...)

    seen = Set([g.first for g in gs])

    tally = maketally(gs)

    return tally, seen
end



# BUILD MODELS


function ppm_stm(data::Vector{View{S,T}},
                 A::Set{T},
                 B::Smoothing,
                 E::Escape,
                 U::Bool,
                 O::OrderBound) where {S,T}

    # STM
    
    return Vector{Prediction{T}}[ppm_seq_inc(v,A,B,E,U,O) for v in data]

end


function ppm_ltm(data::Vector{View{S,T}},
                 A::Set{T},
                 B::Smoothing,
                 E::Escape,
                 U::Bool,
                 O::OrderBound,
                 nfolds::Int=10)::Vector{Vector{Prediction{T}}} where {S,T}

    # LTM

    folds, training = folddataset(data,nfolds)
    
    db = [train(t,O) for t in training]

    tally = first.(db)
    seen = last.(db)

    return Vector{Prediction{T}}[ppm_seq(data[i],tally[folds[i]],seen[folds[i]],A,B,E,U,O) for i in 1:length(data)]

end


function ppm_ltm_plus(data::Vector{View{S,T}},
                      A::Set{T},
                      B::Smoothing,
                      E::Escape,
                      U::Bool,
                      O::OrderBound,
                      nfolds::Int=10) where {S,T}

    # FIX: LTM+

    folds, training = folddataset(data,nfolds)

    db = [train(t,O) for t in training]
    # Should fresh tallies be used for each test set? 
    return Vector{Prediction{T}}[ppm_seq_inc(data[i],A,B,E,U,O; tally = db[folds[i]][1], seen = db[folds[i]][2]) for i in 1:length(data)]
end



function ppm_both(data::Vector{View{S,T}},
                  A::Set{T},
                  B::Smoothing,
                  E::Escape,
                  U::Bool,
                  O::OrderBound,
                  nfolds::Int=10,
                  b::Int=0) where {S,T}

    # FIX: LTM-STM 
    
    stm = ppm_stm(data,A,B,E,U,O)
    ltm = ppm_ltm(data,A,B,E,U,O,nfolds)

    return [[combine(Prediction{T}[p1,p2],b) for (p1,p2) in zip(s,l)] for (s,l) in zip(stm,ltm)]
end

function ppm_both_plus(data::Vector{View{S,T}},
                       A::Set{T},
                       B::Smoothing,
                       E::Escape,
                       U::Bool,
                       O::OrderBound,
                       nfolds::Int=10,
                       b::Int=0) where {S,T}

    # FIX: LTM+-STM
    
    stm = ppm_stm(data,A,B,E,U,O)
    ltm = ppm_ltm_plus(data,A,B,E,U,O,nfolds)

    return [[combine(Prediction{T}[p1,p2],b) for (p1,p2) in zip(s,l)] for (s,l) in zip(stm,ltm)]
end





using DataFrames

function todataframe(id::Int,ps::Vector{Prediction{T}}) where {T}
    seqid = repeat([id],length(ps))
    eventid = [1:length(ps)...]
    sym = [p.symbol for p in ps]
    ord = [p.order for p in ps]
    prob = [p.distribution(p.symbol) for p in ps]
    ic = infcontent.(ps)
    en = entropy.(ps)

    return DataFrame(SeqID = seqid,
                     EventID = eventid,
                     Symbol = sym,
                     Order = ord,
                     Prob = prob,
                     IC = ic,
                     H = en)

end

function todataframe(ps::Vector{Vector{Prediction{T}}}) where {T}
    vcat([todataframe(i,p) for (i,p) in enumerate(ps)]...)
end


# function getalphabet(data::Vector{View{S,T}})::Set{T} where {S,T}

struct IdyomParameters
    E::Escape
    O::OrderBound
end


function getviews(seq::Vector,
                  trg::Viewpoint{T},
                  srcs::Vector{Viewpoint{S} where S}) where T

    return map(s->View(seq,s,trg),srcs)   

end


function idyom(seq::Vector,
               trg::Viewpoint{T},
               srcs::Vector{Viewpoint{S} where S}) where T

    # construct views

    views = getviews(seq,trg,srcs)
    
    # construct models

    # combine models
    

end
    
# end of module
end
