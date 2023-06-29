module AutomaticViewpointSelection

using Chakra
using Charm
using ..Idyoms
using Combinatorics
using LinearAlgebra
using Random
using StatsBase
using Statistics
using Distributions
using ProgressMeter


# Define a struct to represent a pair of a viewpoint and its name
struct viewpoint_pair
    viewpoint::Any
    name::String
end

# Given a list of viewpoints and a list of names, return a list of viewpoint pairs
function initial_viewpoint_pair(viewpoints,
    viewpoints_name)

    n = length(viewpoints_name)
    viewpoint_pairs = []
    for i in 1:n
        push!(viewpoint_pairs, viewpoint_pair(viewpoints[i], viewpoints_name[i]))
    end
    return viewpoint_pairs
end

function initial_viewpoint_pair_2(viewpoints,
    viewpoints_name)

    n = length(viewpoints_name)
    viewpoint_pairs = []
    for i in 1:n
        push!(viewpoint_pairs, (viewpoints[i], viewpoints_name[i]))
    end
    return viewpoint_pairs
end

# Given a starting list of viewpoints, a list of other viewpoints, and a maximum length,
# return a list of all possible mixtures of the starting list and the other viewpoints
function all_possible_mixtures(other_viewpoints,
    max_length)

    mixtures = []
    for i in 2:max_length
        combs = combinations(other_viewpoints, i)
        for comb in combs
            new_mixture = []
            append!(new_mixture, comb)
            push!(mixtures, new_mixture)
        end
    end
    mixtures_result = []
    for i in mixtures
        new_mixture = []
        name = ""
        for j in i
            push!(new_mixture, j[1])
            name = (name * j[2] * "_X_")
        end
        name = name[1:end-3] 
        push!(mixtures_result, (Chakra.link(new_mixture...), name))
    end
    return mixtures_result
end


function main_fitness(src,
    trg,
    sequences,
    a,
    b,
    e,
    u,
    o,
    weight,
    nfold)

    score = 0
    views = [[View(seq,s[1],trg[1]) for seq in sequences] for s in src];
    models = [Idyoms.ppm_both(v,a,b,e,u,o,nfold) for v in views];
    m = Idyoms.combine_predictions(models,weight);
    MIC = Idyoms.mean_infcontent(m)
    score = 1/MIC
    return score
end

function best_first_search_main(trg,
    trg_name,
    src,
    src_name,
    sequences,
    max_linked_size,
    a,
    b,
    e,
    u,
    o,
    weight,
    nfold)
    
    # Primary initialization
    src_pairs = initial_viewpoint_pair_2(src, src_name)
    trg_pair = (trg, trg_name)
    Best_Viewpoint = [src_pairs[1]]
    Best_fitness = main_fitness(Best_Viewpoint,trg_pair,sequences,a,b,e,u,o,weight,nfold)
    Final_mixture = []
    VP_to_avoid = []
    println("=========================================")
    
    
    # Now doing the process for no linked viewpoint
    for i in 1:length(src_pairs)
        available_vps = setdiff(src_pairs, VP_to_avoid)
        push!(Final_mixture, available_vps[1])
        if i == 1
            Best_fitness_for_halt_check = 0
        else
            Best_fitness_for_halt_check = Best_fitness
        end
        for v in available_vps
            println("----------------------------")
            println("mixture length: ", length(Best_Viewpoint))
            saved_vp = Final_mixture[end]
            Final_mixture[end] = v
            fitness = main_fitness(Final_mixture,trg_pair,sequences,a,b,e,u,o,weight,nfold)
            if fitness > Best_fitness
                println("Fitness improvment: ", Best_fitness, " -> ", fitness)
                Best_fitness = fitness
                Best_Viewpoint = Final_mixture
            else
                Final_mixture[end] = saved_vp
            end
            println("----------------------------")
        end
        println("=========================================")
        if Best_fitness_for_halt_check == Best_fitness
            pop!(Best_Viewpoint)
            break
        else
            push!(VP_to_avoid, Best_Viewpoint[end])
        end
    end
    
    # Get the possible mixtures for linking
    if max_linked_size == 1
        linked_possible_mixtures = []
    else
        linked_possible_mixtures = all_possible_mixtures(src_pairs, max_linked_size)
    end
    
    # linked
    for i in 1:length(linked_possible_mixtures)
        available_vps = setdiff(linked_possible_mixtures, VP_to_avoid)
        push!(Final_mixture, available_vps[1])
        Best_fitness_for_halt_check = Best_fitness
        for v in available_vps
            println("----------------------------")
            println("mixture length: ", length(Best_Viewpoint))
            saved_vp = Final_mixture[end]
            Final_mixture[end] = v
            fitness = main_fitness(Final_mixture,trg_pair,sequences,a,b,e,u,o,weight,nfold)
            if fitness > Best_fitness
                println("Fitness improvment: ", Best_fitness, " -> ", fitness)
                Best_fitness = fitness
                Best_Viewpoint = Final_mixture
            else
                Final_mixture[end] = saved_vp
            end
            println("----------------------------")
        end
        println("=========================================")
        if Best_fitness_for_halt_check == Best_fitness
            pop!(Best_Viewpoint)
            break
        else
            push!(VP_to_avoid, Best_Viewpoint[end])
        end
    end
    
    return [s[2] for s in Final_mixture], Best_fitness
    
end
end
