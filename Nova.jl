module Nova

using MIDI, DataFrames
using Chakra
using Viewpoints

export data

function loadmidi(filepaths)

    df = DataFrame(id=[],particles=[],pitch=[],onset=[],duration=[],velocity=[])
    
    for fp in filepaths

        file = readMIDIFile(fp)

        tracks = file.tracks
        trackids = []
        
        for (t,track) in enumerate(tracks)
            trackid = string(fp,"/track",t)
            push!(trackids,trackid)

            notes = getnotes(track,file.tpq)
            noteids =[]
            
            for (n,note) in enumerate(notes)

                noteid = string(trackid,"/note",n)
                push!(noteids,noteid)
                
                pitch = Int(note.pitch)
                onset = Int(note.position)
                duration = Int(note.duration)
                velocity = Int(note.velocity)
                push!(df,Dict(:id=>noteid,
                              :particles=>[],
                              :pitch=>pitch,
                              :onset=>onset,
                              :duration=>duration,
                              :velocity=>velocity))
            end

            push!(df,Dict(:id=>trackid,
                          :particles=>noteids,
                          :pitch=>missing,
                          :onset=>missing,
                          :duration=>missing,
                          :velocity=>missing))
            
        end

        fileid = fp
        push!(df,Dict(:id=>fileid,
                      :particles=>trackids,
                      :pitch=>missing,
                      :onset=>missing,
                      :duration=>missing,
                      :velocity=>missing))
        
    end
    return df
    
end

export ID, C, H

ID = String
C = DataFrameRow{DataFrame,DataFrames.Index}
H = DataFrame

Chakra.@Attribute :pitch Int
Chakra.@Attribute :onset Int
Chakra.@Attribute :duration Int
Chakra.@Attribute :velocity Int

function Chakra.fnd(x::ID,h::H)::Option{C}
    c = h[h.id .== x,:]
    isempty(c) && return none
    return c[1,:]
end

function Chakra.geta(::Att{a,T},c::C)::Option{T} where {a,T}
    val = c[a]
    val isa Missing && return none
    return val
end

function Chakra.dom(h::H)::List{ID}
    h.id
end

function Chakra.pts(c::C)::List{ID}
    c[:particles]
end

pitchvp = vp(:pitch)
onsetvp = vp(:onset)
ioivp = compose(link(onsetvp,delay(onsetvp,1)),(x,y)->x-y)

filenames = readdir("nova")
paths = map(fn->string("nova/",fn),filenames)
data = loadmidi(paths)

end
