# using Polyhedra, CDDLib
# using GLPK;
# using SparseArrays, LinearAlgebra, Random
# solver = GLPK.Optimizer;


function generate_adjacency_list(P::Polyhedron)
    println("called")
    ϵ = 1e-7
    ver_rep = removevredundancy(vrep(P), solver)
    vr = MixedMatVRep(ver_rep)
    h_rep = removehredundancy(hrep(P), solver)
    hr = MixedMatHRep(h_rep)
    A, b = hr.A, hr.b
    nonzero_indices = [i for i in 1:size(A, 1)] |> arr -> filter(i -> A[i, :]' * A[i, :] > 1e-6, arr)
    A = A[nonzero_indices, :]
    b = b[nonzero_indices]
    no_of_vertics_and_dirs = size(vr.V, 2) + size(vr.R, 2)

    brepearted = repeat(b, outer=(1, size(vr.V, 2)))
    zeros_repeadted_R_dim = zeros(eltype(b), (size(b, 1), size(vr.R, 2)))
    rows_V, rows_R = [vr.V, vr.R] .|> size .|> first
    if rows_V != rows_R
        return nothing
    end
    bt = sparse(abs.(A * hcat(vr.V, vr.R) - hcat(brepearted, zeros_repeadted_R_dim)) .< ϵ)

    adjdatabit = sparse(bt' * bt .>= fulldim(h_rep) - 1)
    adjdatabit = adjdatabit - diagm(diag(adjdatabit))
    adjdatalist = Dict{Int,Vector{Int}}()
    adjdatalist = Dict(j => findall(i -> i == 1, adjdatabit[j, :]) for j in 1:no_of_vertics_and_dirs)

    vinx = [pi for pi in eachindex(points(vr))]
    rinx = [pi for pi in eachindex(rays(vr))]
    vr_inx = vec(hcat(vinx, rinx))
    d_inx = Dict(vr_inx[j] => vr_inx[findall(i -> i == 1, adjdatabit[j, :])] for j in 1:no_of_vertics_and_dirs)
    return Dict("bit" => adjdatabit, "list" => adjdatalist, "indices" => d_inx)
end

function generate_adjacency_list(h::HRepresentation)
    P = polyhedron(h, CDDLib.Library())
    generate_adjacency_list(P)
end

function generate_adjacency_list(v::VRepresentation)
    P = polyhedron(v, CDDLib.Library())
    generate_adjacency_list(P)
end


# function chooseDirection(vi::Int, P::Polyhedron, dim::Int)
#     adj = generate_adjacency_list(P)
#     adj_list = get(adj["list"], vi, [])
#     ee = ones(dim)
#     if (length(adj_list) < dim)
#         return ee
#     end
#     adj_list = shuffle(adj_list)[1:dim]
#     v = vrep(P) |> MixedMatVRep
#     extreme = hcat(v.V, v.R)
#     A = extreme[:, adj_list]
#     d = A \ ee
#     d
# end


function fixed_direction(v; first_elm=1)
    d = ones(eltype(v), length(v))
    d[1] = convert(eltype(v), first_elm)
    d
end

function fixed_point(phat)
    function direction(v)
        d = phat - v
        d
    end
    direction
end

function find_verex_index(v, A)
    _, n = size(A)
    findfirst([A[:, col] == v for col in 1:n])
end

function adj_direction(vi::Int, P::Polyhedron, dim::Int)
    ee = ones(dim)
    adj = generate_adjacency_list(P)
    if isnothing(adj)
        println("🅰")
        return ee / norm(ee)
    end
    adj_list = get(adj["list"], vi, [])
    if (length(adj_list) < dim)
        println("🅱")
        return ee / norm(ee)
    end
    adj_list = shuffle(adj_list)[1:dim]
    v = vrep(P) |> MixedMatVRep
    extreme = hcat(v.V, v.R)
    A = extreme[:, adj_list]
    d = A \ ee
    d
end