# using Random
# rng = MersenneTwister(1234);


function random_vertex(V::Matrix{Float64})
    m, _ = size(V)
    ind = rand(rng, 1:m)
    Vector(V[ind, :])
end

function next_vertex(V::Matrix{Float64})
    Vector(V[1, :])
end

