function isVertixInVerices(v::Vector{<:Real}, V::Matrix{<:Real})
    v = convert(Vector{Float64}, v)
    V = convert(Matrix{Float64}, V)
    m, _ = size(V)
    for i in 1:m
        ov = V[i, 1:end]
        if v == ov
            return true
        end
    end
    return false

end


function getUnusedVertices(Vk, Vused)
    if isempty(Vk)
        return nothing
    end
    m, _ = size(Vk)
    Vunused = filter(i -> !isVertixInVerices(Vk[i, :], Vused), 1:m) |> inx -> Vk[inx, :]
    Vunused
end


function algorith1(f, ws, ps, dps, W;
    chooseDirection=fixed_direction,
    chooseVertex=chooseVertex,
    newDirection=nothing,
    ϵ=0.05)
    SCCounter = 0
    dim = length(W[1])
    X = W .|> x -> ws(x)
    SCCounter += 1
    Vused = Matrix{Float64}(undef, 1, dim)

    Xbar0 = Matrix{Float64}(undef, 1, dim)
    for i in 1:dim
        Xbar0 = vcat(Xbar0, X[i]')
    end
    Xbar0 = Xbar0[2:end, :]


    #=
    P₀ = ∪ {y ∈ Rᵈ |  -(wⁱ)ᵀy≤ -(wⁱ)ᵀf(xⁱ)} 
    =#
    hr = zip(W, X) .|> wx -> HalfSpace(-first(wx), -first(wx)' * f(last(wx)))
    hr = hrep(hr)
    P0 = polyhedron(hr, CDDLib.Library())

    # V0 = vrep(P0) |> MixedMatVRep |> v -> v.V
    V0 = doubledescription(hr) |> MixedMatVRep |> v -> v.V
    Vunused = getUnusedVertices(V0, Vused)


    while !isempty(Vunused)
        v = chooseVertex(Vunused)
        d = chooseDirection(v)
        if ~isnothing(newDirection)
            adj_d = newDirection(2, P0, dim)
            println(adj_d)
        end
        ##
        # define a function to add a list of adjaceny vertices 
        ###
        (xv, zv) = ps(v, d)
        SCCounter += 1
        wv = dps(v, d)
        SCCounter += 1
        Vused = vcat(Vused, v')
        Xbar0 = vcat(Xbar0, xv')
        if (zv > ϵ)
            H = HalfSpace(-wv, -wv' * v - zv)
            hr = hr ∩ H
            ddvrep = doubledescription(hr)
            P0 = polyhedron(ddvrep, CDDLib.Library())
            # Vtemp = vrep(P0) |> MixedMatVRep |> s -> s.V
            Vtemp = ddvrep |> MixedMatVRep |> s -> s.V
            V0 = Vtemp
            Vunused = getUnusedVertices(V0, Vused)
        else
            break
        end

    end
    Xbar0, SCCounter, V0, Vunused, P0, hr, Vused
end