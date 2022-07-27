module Viewpoints 

using Chakra

export Viewpoint
export vp, link, compose, delay, thread, vp_map


abstract type Viewpoint{T} end

returntype(v::Viewpoint{T}) where T = T

struct AtomicViewpoint{T} <: Viewpoint{T}
    attribute::Symbol
    returntypes::List{DataType}
    AtomicViewpoint(a::Symbol) = new{Chakra.__typ__(a)}(a,[Chakra.__typ__(a)])
end

function (v::AtomicViewpoint{T})(s::List)::Option{T} where T 
    obind(lpeek(s), o->Chakra.geta(v.attribute,o))
end

struct LinkedViewpoint{T} <: Viewpoint{T}
    components::List{Viewpoint}
    returntypes::List{DataType}
    LinkedViewpoint(v1::Viewpoint,v2::Viewpoint,vs::Viewpoint...) = begin
        components = [v1,v2,vs...]
        returntypes = [returntype(c) for c in components]
        new{Tuple{returntypes...}}(components,returntypes)
    end
end

function (v::LinkedViewpoint{T})(s::List)::Option{T} where T
    res = []
    for c in v.components
        val = c(s)
        if val == none
            return none
        end
        push!(res,val)
    end
    #print(res)
    return Tuple(res)
end

struct DerivedViewpoint{T} <: Viewpoint{T}
    base::Viewpoint
    modifier::Function
    returntypes::Vector{DataType}
    DerivedViewpoint(v::Viewpoint,f::Function) = begin
        t = Base._return_type(f,Tuple(v.returntypes))
        if t == Union{}
            error("Type mismatch: The function $f is not composable with the viewpoint $v")
        end
        new{t}(v,f)
    end
end

function (v::DerivedViewpoint{T})(s::Vector)::Option{T} where T
    n = length(v.base.returntypes)
    if n == 1
        obind(v.base(s),(x)->v.modifier(x))
    end
    obind(v.base(s),(x)->v.modifier(x...))
end

struct DelayedViewpoint{T} <: Viewpoint{T}
    base::Viewpoint{T}
    lag::Int64
end

function (v::DelayedViewpoint{T})(s::Vector)::Option{T} where T 
    v.base(lpopn(s,v.lag))
end

struct ThreadedViewpoint{T} <: Viewpoint{T}
    base::Viewpoint{T}
    test::Viewpoint{Bool}
end

function (v::ThreadedViewpoint{T})(s::Vector)::Option{T} where T
    v.test(s) ? v.base(s) : none
end


# Viewpoint constructor interface

vp(x::Symbol) = AtomicViewpoint(x)

link(v1::Viewpoint,v2::Viewpoint,vs::Viewpoint...) = LinkedViewpoint(v1,v2,vs...)

compose(v::Viewpoint,f::Function) = DerivedViewpoint(v,f)

delay(v::Viewpoint,l::Int) = DelayedViewpoint(v,l)

thread(b::Viewpoint,t::Viewpoint) = ThreadedViewpoint(b,t)



# Additional viewpoint constructors

diff(v::Viewpoint{T}) where T = compose(link(v,delay(v,1)),(x,y)->x-y)







function vp_map(v::Viewpoint{T},s::Vector)::List{Option{T}} where T

    # Maps a viewpoint of a sequence

    return [v(s[1:n]) for n in 1:length(s)]
end



# end of module
end
