module Idyoms

export Bounded, Unbounded
export View
export A, B, C, D, X
export Interpolated, Backoff

export generate_ngrams, generate_hgrams, hgram_sequence

export Distribution, Prediction
export infcontent, entropy

export ppm, ppm_seq, ppm_seq_inc


using Chakra


### ORDER BOUNDS

abstract type OrderBound end

struct Bounded{h} <: OrderBound
    Bounded(h::Int) = h < 0 ? error("Can't have bound less than 0.") : new{h}()
end

struct Unbounded <: OrderBound end


# TYPE OF CONTEXTS

Context{S} = SubArray{S, 1, Vector{S}, Tuple{UnitRange{Int64}}, true}

function trim(ctx::Context{S},l::Int) where S

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
        
        t = last.(elems) # target sequence including none
        s = first.(elems) # source sequence including none
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

function get_alphabet(v::Viewpoint{T},s::Vector{Chakra.Constituent})::Set{T} where T

    # RETURN THE VIEWPOINT ALPHABET OF A SEQUENCE
    
    return Set{T}(filter(e->e!=none,vp_map(v,s)))

end

function get_alphabet(v::Viewpoint{T},ss::Vector{Vector{Chakra.Constituent}})::Set{T} where T

    # RETURN THE VIEWPOINT ALPHABET OF A VECTOR OF SEQUENCES
    
    return union([get_alphabet(v,s) for s in ss]...)

end
    
function getnext(v::View{S,T},i::Int)::T where {S,T}

    # iTH TARGET ATTRIBUTE FROM VIEW
    
    return v.target[v.targetindex[i]]
end

function getcontext(v::View{S,T},i::Int)::Context{S} where {S,T}

    # iTH CONTEXT FROM VIEW

    return @views v.sourceelements[1:length(findall(j->j<v.targetindex[i],v.sourceindex))]
end

function getngram(v::View{S,T},i::Int) where {S,T}

    # iTH TARGET-CONTEXT PAIR FROM VIEW 
    
    return getnext(v,i) => getcontext(v,i)
end

function getngram(v::View{S,T},i::Int,n::Int) where {S,T}

    # iTH N-GRAM FROM VIEW
    
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



# TODO: Which of these functions is correct?
# if a and b are 0, a/b is nan so return 0
# if a/b is inf should it return 0 or 1?

#non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) ? 0 : x)
#non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) || isinf(x) ? 0 : x)
non_nan(a::Int,b::Int)::Float64 = (x = a/b; isnan(x) ? 0 : isinf(x) ? 1 : x)





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


### NGRAM TALLIES

# THE TYPE OF NGRAM TALLIES
Tally{S,T} = Dict{Context{S},Counter{T}}


# CONSTRUCT AN EMPTY TALLY
emptytally(S,T) = Tally{S,T}()

function incrementtally!(tally::Tally{S,T},nxt::T,ctx::Context{S}) where {S,T}

    # INCREMENT A TALLY WITH NGRAM
    
    Base.get!(tally,ctx,Counter{T}(0))[nxt] += 1
end

function maketally(gs::Vector{Pair{T,Context{S}}}) where {S,T}

    # RERURN A NEW TALLY CONSTRUCTED FROM A VECTOR OF NGRAMS
    
    tally = emptytally(S,T)
    
    for p in gs
        incrementtally!(tally,p...)
    end
    
    return tally
end

function updatetally(tally::Tally{S,T},nxt::T,ctx::Context{S}) where {S,T}

    # INCREMENT A TALLY WITH ALL NGRAMS
    
    for n in 0:length(ctx)
        incrementtally!(tally,nxt,trim(ctx,n))
    end
end



### SMOOTHING METHODS

abstract type Smoothing end

struct Interpolated <: Smoothing end

struct Backoff <: Smoothing end



### ESTIMATE PROBABILITY

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

function estimate(nxt::T,
                  ctx::Context{S},
                  tally::Tally{S,T},
                  seen::Set{T},
                  A::Set{T},
                  B::Interpolated,
                  E::Escape,
                  U::Bool,
                  ex::Set{T}=Set{T}()) where {S,T}

    # ESTIMATE PROBABILITY WITH INTERPOLATED SMOOTHING
    
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


function select_order(v::View{S,T},
                      i::Int,
                      tally::Tally{S,T},
                      seen::Set{T},
                      O::Bounded{h}) where {S,T,h}
    
    # RETURN BOUNDED NGRAM i FROM VIEW v
    
    return getngram(v,i,h+1)
end


function select_order(v::View{S,T},
                      i::Int,
                      tally::Tally{S,T},
                      seen::Set{T},
                      O::Unbounded) where {S,T}

    # RETURN UNBOUNDED NGRAM AT i FROM VIEW v
    
    nxt, ctx = getngram(v,i)

    length(ctx) == 0 && return nxt, ctx

    for l in 0:length(ctx)

        @views _ctx = ctx[end-l+1:end]

        c = Base.get(tally,_ctx,Counter{T}(0))

        tc = symcount(c)

        # RETURN SHORTEST DETERMINISTIC CONTEXT
        tc == 1 && return nxt,_ctx 

        # OR, RETURN LONGEST MATCHING CONTEXT
        tc == 0 && return @views nxt, _ctx[2:end]

    end
        
    return nxt, ctx

end



### TYPE OF PROBABILITY DISTRIBUTIONS

Distribution{T} = Dict{T,Float64}

# GET ALPHABET OF DISTRIBUTION
domain(d::Distribution) = keys(d)

# NORMALISE A DISTRIBUTION
normalise!(d::Distribution) = (total = sum(values(d)); map!(x->x/total,values(d)))

# LOOKUP PROBABILITY OF ELEMENT
(d::Distribution{T})(e::T) where T = Base.get(d,e,0)


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

function estimate_dist(ctx::Context{S},
                       tally::Tally{S,T},
                       seen::Set{T},
                       A::Set{T},
                       B::Smoothing,
                       E::Escape,
                       U::Bool) where {S,T}

    # ESTIMATE THE PROBABILITY DISTRIBUTION OVER ALPHABET A
    
    dist = Dict{T,Float64}()
    ords = Dict{T,Int}()

    for x in A
        p, o = estimate(x,ctx,tally,seen,A,B,E,U)
        dist[x] = p
        ords[x] = o
    end

    normalise!(dist)

    return dist, ords
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


Model{T} = Vector{Vector{Prediction{T}}}


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

function ppm(nxt::T,
             ctx::Context{S},
             tally::Tally{S,T},
             seen::Set{T},
             A::Set{T},
             B::Smoothing,
             E::Escape,
             U::Bool) where {S,T}

    # RETURN PPM PREDICTION FOR AN ELEMENT IN A CONTEXT
    
    dist, ords = estimate_dist(ctx,tally,seen,A,B,E,U)
    
    return Prediction(nxt,dist[nxt],ords[nxt],dist)
    
end

struct ViewModel{S,T}
    view::View{S,T}
    predictions::Vector{Prediction{T}}
end

function ppm_seq(v::View{S,T},
                 tally::Tally{S,T},
                 seen::Set{T},
                 A::Set{T},
                 B::Smoothing,
                 E::Escape,
                 U::Bool,
                 O::OrderBound) where {S,T}

    # RETURN PPM PREDICTION SEQUENCE FOR THE ELEMENTS OF A VIEW
    
    predictions = Prediction{T}[]
    
    for i in 1:length(v)

        # SELECT iTH NGRAM CONTEXT LENGTH
        nxt,ctx = select_order(v,i,tally,seen,O)

        # PREDICT iTH ELEMENT
        push!(predictions,ppm(nxt,ctx,tally,seen,A,B,E,U))
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
                     seen = Set{T}()) where {S,T}

    # RETURN PPM PREDICTION SEQUENCE FOR A VIEW WITH INCREMENTAL TALLY UPDATE
    
    predictions = Prediction{T}[]
    
    for i in 1:length(v)

        # SELECT iTH NGRAM CONTEXT LENGTH
        nxt, ctx = select_order(v,i,tally,seen,O)

        # PREDICT iTH ELEMENT
        push!(predictions,ppm(nxt,ctx,tally,seen,A,B,E,U))

        # UPDATE NGRAM TALLY
        updatetally(tally,getngram(v,i)...)

        # UPDATE SEEN ELEMENT SET
        push!(seen,nxt)
    end
    
    return predictions
    
end



### DATASET MANAGEMENT

function folddataset(data::Vector{View{S,T}},
                     nfolds::Int) where {S,T}

    # RETURN FOLD INDEXES AND TRAINING SETS FOR EACH FOLD
    
    # PARTITION DATASET INTO N FOLDS:
    size = length(data)
    
    if !(0 < nfolds <= size)
        error("The number of folds must be less than the size of the data set size.")
    end

    # ASSIGN EACH SEQUENCE VIEW TO A FOLD:
    folds = Int[Int(round(i/nfolds % 1 *nfolds)) for i in 1:size]
    
    folds[findall(x->x==0,folds)] .= nfolds

    folds = sort(folds)

    # CREATE A TRAINING SET FOR EACH FOLD:
    training_sets = Vector{View{S,T}}[data[findall(x->x != i,folds)] for i in 1:nfolds]
    
    return folds, training_sets

end

function train(data::Vector{View{S,T}},
               O::OrderBound) where {S,T}

    # RETURN AN NGRAM TALLY AND A SET OF SEEN ELEMENTS FOR A SET OF VIEWS
    
    # GENERATE ALL NGRAMS: 
    gs = vcat([generate_hgrams(ex,O) for ex in data]...)

    # GENERATE SET OF SEEN ELEMENTS:    
    seen = Set{T}([g.first for g in gs])

    # TALLY NGRAMS:  
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

    # RETURN SHORT-TERM MODEL (STM) PREDICTIONS FOR A DATASET:
    
    return Vector{Prediction{T}}[ppm_seq_inc(v,A,B,E,U,O) for v in data]

end


function ppm_ltm(data::Vector{View{S,T}},
                 A::Set{T},
                 B::Smoothing,
                 E::Escape,
                 U::Bool,
                 O::OrderBound,
                 nfolds::Int=10)::Vector{Vector{Prediction{T}}} where {S,T}

    # RETURN LONG-TERM MODEL (LTM) PREDICTIONS FOR A DATASET:

    # FOLD DATASET:
    folds, training = folddataset(data,nfolds)

    # CREATE NGRAM TALLIES AND SEEN ELEMENTS FOR EACH FOLD:
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

    # RETURN LONG-TERM MODEL WITH INCREMENTAL UPDATE (LTM+) PREDICTIONS FOR A DATASET

    # FOLD DATASET:
    
    folds, training = folddataset(data,nfolds)

    # CREATE NGRAM TALLIES AND SEEN ELEMENTS FOR EACH FOLD:
    
    db = [train(t,O) for t in training]

    # TODO: SHOULD FRESH TALLIES BE USED FOR EACH TEST SET?

    return Vector{Prediction{T}}[ppm_seq_inc(data[i],A,B,E,U,O;tally=db[folds[i]][1],seen=db[folds[i]][2]) for i in 1:length(data)]
end



function ppm_both(data::Vector{View{S,T}},
                  A::Set{T},
                  B::Smoothing,
                  E::Escape,
                  U::Bool,
                  O::OrderBound,
                  nfolds::Int=10,
                  b::Int=0) where {S,T}

    # RETURN COMBINED LTM AND STM PREDICTIONS FOR A DATASET
    
    stm = ppm_stm(data,A,B,E,U,O)
    ltm = ppm_ltm(data,A,B,E,U,O,nfolds)

    return [[combine_dist(Prediction{T}[p1,p2],b) for (p1,p2) in zip(s,l)] for (s,l) in zip(stm,ltm)]
end

function ppm_both_plus(data::Vector{View{S,T}},
                       A::Set{T},
                       B::Smoothing,
                       E::Escape,
                       U::Bool,
                       O::OrderBound,
                       nfolds::Int=10,
                       b::Int=0) where {S,T}

    # RETURN COMBINED LTM+ AND STM PREDICTIONS FOR A DATASET
    
    stm = ppm_stm(data,A,B,E,U,O)
    ltm = ppm_ltm_plus(data,A,B,E,U,O,nfolds)

    return [[combine_dist(Prediction{T}[p1,p2],b) for (p1,p2) in zip(s,l)] for (s,l) in zip(stm,ltm)]
end







### GENERATION

function gen(ctx::Context{S},
             tally::Tally{S,T},
             seen::Set{T},
             A::Set{T},
             B::Smoothing,
             E::Escape,
             U::Bool) where {S,T}
    
    # RETURN A RANDOM ELEMENT FROM TEH DISTRIBUTION GENERATED BY THE CONTEXT
    
    dist, ords = estimate_dist(ctx,tally,seen,A,B,E,U)
    
    nxt = sample(dist)
    
    return Prediction(nxt,dist[nxt],ords[nxt],dist)
end

function gen_seq(len::Int,
                 tally::Tally{T,T},
                 seen::Set{T},
                 A::Set{T},
                 B::Smoothing,
                 E::Escape,
                 U::Bool,
                 O::OrderBound) where T

    # RETURN A SEQUENCE OF LENGTH len OF RANDOM ELEMENTS
    
    preds = Prediction{T}[]
    seq = T[]
    
    for i in 1:len
        ctx = select_order(seq,tally,seen,O)
        p = gen(ctx,tally,seen,A,B,E,U,)
        push!(preds,p)
        push!(seq,p.symbol)
    end
    
    return preds, seq
end

function gen_seq_inc(len::Int,
                     A::Set{T},
                     B::Smoothing,
                     E::Escape,
                     U::Bool,
                     O::OrderBound,
                     tally::Tally{T,T} = emptytally(T,T),
                     seen::Set{T} = Set{T}()) where T

    # RETURN A SEQUENCE OF LENGTH len OF RANDOM ELEMENTS WITH INCREMENTAL UPDATE
    
    preds = Prediction{T}[]
    seq = T[]
    
    for i in 1:len

        # SELECT PREDICTION ORDER
        ctx = select_order(seq,tally,seen,O)

        # GENERATE NEXT ELEMENT
        p = gen(ctx,tally,seen,A,B,E,U,)

        # UPDATE SEQUENCE
        push!(preds,p)
        push!(seq,p.symbol)

        # UPDATE NGRAM TALLY
        updatetally(tally,getngram(v,i)...)

        # UPDATE SEEN ELEMENT SET
        push!(seen,nxt)
    end
    
    return preds, seq
end

### DATAFRAMES

using DataFrames

function todataframe(id::Int,ps::Vector{Prediction{T}}) where {T}

    # RETURN A DATAFRAME GENERATED FROM A SEQUENCE OF PREDICTIONS
    
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

    # RETURN A DATAFRAME GENERATED FROM A LIST OF PREDICTION SEQUENCES
    
    vcat([todataframe(i,p) for (i,p) in enumerate(ps)]...)
end

include("AutomaticViewpointSelection.jl")
    
# end of module
end
