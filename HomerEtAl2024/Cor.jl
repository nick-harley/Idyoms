using Statistics, StatsBase, JLD2, FileIO

function main(seqs,v,n,m,b,e,u)

    res = load("correlations.jld2")["res"]

    dsname = string(v,n,m,b,e,u)

    haskey(res,dsname) && return res[dsname]
    
    views = [View(seq,v,v,v) for seq in seqs];

    gs = NMgs(views,n,m);

    txyZ = Idyoms.tally_xyZ(gs);
    txzY = Idyoms.tally_xzY(gs);
    txZ = Idyoms.tally_xZ(gs);
    tyZ = Idyoms.tally_yZ(gs);
    txY = Idyoms.tally_xY(gs);

    dX = Set([g[1] for g in gs]);
    dY = Set([g[2] for g in gs]);
    dZ = Set([g[3] for g in gs]);

    dist(x,t,s,a,b,e,u) = Idyoms.estimate_dist(x,t,s,a,b,e,u)[1];

    pZxy(x::Seq{X},y::Y) where {X,Y} = dist((x=>y),txyZ,dZ,dZ,b,e,u);
    pZy(y::Y) where Y = @views dist([y][1:end],tyZ,dZ,dZ,b,e,u);
    pZx(x::Seq{X}) where X = dist(x,txZ,dZ,dZ,b,e,u);
    pYx(x::Seq{X}) where X = dist(x,txY,dY,dY,b,e,u);
    pYxz(x::Seq{X},z::Z) where {X,Z} = dist((x=>z),txzY,dY,dY,b,e,u);

    function FCI(x::Seq{X},y::Y,z::Seq{Z}) where {X,Y,Z}

        pzxy = pZxy(x,y)
        pzx = pZx(x)
        pzy = pZy(y)
        pyxz = pYxz(x,z)
        pyx = pYx(x)

        fpc = sum([pzxy(z)*log(2,pzxy(z)/pzx(z)) for z in dZ])
        fcc = sum([pzxy(z)*log(2,pzxy(z)/pzy(z)) for z in dZ])
        frc = sum([pyxz(y)*log(2,pyxz(y)/pyx(y)) for y in dY])
        ic = Idyoms.infcontent(pyx,y)
        h = Idyoms.entropy(pyx)

        efpc = 0
        for y in dY
            p = pZxy(x,y)
            ci = sum([p(z)*log(2,p(z)/pzx(z)) for z in dZ])
            efpc += pyx(y)*ci
        end
        
        return fpc, fcc, frc, pyx(y), pzxy(z), pzx(z), pzy(z), pyxz(y), ic, h, efpc
        
    end

    function FCI(v::View{X,Y,Z}) where {X,Y,Z}
        
        fci = [FCI(nmg(v,i,n,m)...) for i in 1:length(v)]
        fpc = [t[1] for t in fci];
        fcc = [t[2] for t in fci];
        frc = [t[3] for t in fci];
        pyx = [t[4] for t in fci];
        pzxy = [t[5] for t in fci];
        pzx = [t[6] for t in fci];
        pzy = [t[7] for t in fci];
        pyxz = [t[8] for t in fci];
        ic = [t[9] for t in fci];
        h = [t[10] for t in fci];
        efpc = [t[11] for t in fci];
        
        return fpc, fcc, frc, pyx, pzxy, pzx, pzy, pyxz, ic, h, efpc 

    end

    fci = [FCI(view) for view in views]

    fpc = vcat([t[1] for t in fci]...);
    fcc = vcat([t[2] for t in fci]...);
    frc = vcat([t[3] for t in fci]...);
    pyx = vcat([t[4] for t in fci]...);
    pzxy = vcat([t[5] for t in fci]...);
    pzx = vcat([t[6] for t in fci]...);
    pzy= vcat([t[7] for t in fci]...);
    pyxz= vcat([t[8] for t in fci]...);
    ic = vcat([t[9] for t in fci]...);
    h = vcat([t[10] for t in fci]...);
    efpc = vcat([t[11] for t in fci]...);
        
    cors = Dict(
        1=>Dict(
            1=>Dict(
                1=> cor(fpc,ic),
                2=> mean([cor(d[1],d[9]) for d in fci])),         
            2=>Dict(
                1=> corspearman(fpc,ic),
                2=> mean([corspearman(d[1],d[9]) for d in fci])),
            3=>Dict(
                1=> corkendall(fpc,ic),
                2=> mean([corkendall(d[1],d[9]) for d in fci]))),
        2=>Dict(
            1=>Dict(
                1=> cor(efpc,h),
                2=> mean([cor(d[11],d[10]) for d in fci])),
            2=>Dict(
                1=> corspearman(efpc,h),
                2=> mean([corspearman(d[11],d[10]) for d in fci])),
            3=>Dict(
                1=> corkendall(efpc,h),
                2=> mean([corkendall(d[11],d[10]) for d in fci]))))

    res[dsname] = cors
    save("correlations.jld2","res",res)

    return cors

end

