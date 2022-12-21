

rng = MersenneTwister(1425)


function generate_adjacency_list(P::Polyhedron)
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
    bt = sparse(abs.(A * hcat(vr.V, vr.R) - hcat(brepearted, zeros_repeadted_R_dim)) .< ϵ)

    adjdatabit = sparse(bt' * bt .>= fulldim(h_rep) - 1)
    adjdatabit = adjdatabit - diagm(diag(adjdatabit))
    adjdatalist = Dict{Int,Vector{Int}}()
    adjdatalist = Dict(j => findall(i -> i == 1, adjdatabit[j, :]) for j in 1:no_of_vertics_and_dirs)

    vinx = [pi for pi in eachindex(points(v))]
    rinx = [pi for pi in eachindex(rays(v))]
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


function chooseDirection(vi::Int, P::Polyhedron, dim::Int)
    adj = generate_adjacency_list(P)
    adj_list = get(adj["list"], vi, [])
    ee = ones(dim)
    if (length(adj_list) < dim)
        return ee
    end
    adj_list = shuffle(adj_list)[1:dim]
    v = vrep(P) |> MixedMatVRep
    extreme = hcat(v.V, v.R)
    A = extreme[:, adj_list]
    d = A \ ee
    d
end


# # H₁ = {x ∈  𝐑² | x₁ + x₂ ≥ 2 }
# H1 = HalfSpace([-1, -1], -2)

# # H₂ = {x ∈  𝐑² | x₁  ≥ 0 }
# H2 = HalfSpace([-1, 0], 0)

# # H₃ = {x ∈  𝐑² | x₂  ≥ 0 }
# H3 = HalfSpace([0, -1], 0)

# h = hrep([H1, H2, H3]);
# v_of_h = polyhedron(h) |> vrep;
# h_of_v = polyhedron(v_of_h) |> hrep;
# ad = generate_adjacency_list(h)
# add = generate_adjacency_list(h_of_v)
# dd = generate_adjacency_list(v_of_h)

v = convexhull([2, 0], [0, 2]) + conichull([1, 0], [0, 1])
# nH = HalfSpace([-1.0, -0.0], 0.0) ∩ HalfSpace([-0.0, -0.0], 1.0) ∩ HalfSpace([-0.0, -1.0], 0.0) ∩ HalfSpace([-1.0, -1.0], -2.0)
nP = polyhedron(v, CDDLib.Library())
# h2 = hrep(nP)
# removehredundancy(nH,solver)
# vMat = MixedMatVRep(v)
# extreme = hcat(vMat.V, vMat.R)
# nhMat.A
# adj_v = generate_adjacency_list(nP)
d = chooseDirection(1, nP, 2)