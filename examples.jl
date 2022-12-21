
# using LinearAlgebra

# Example 1
function example1(; dim=2, explicit=false)
    f(x) = x
    # W = [[1, 0], [0, 1]]
    W = map(i -> [i == j for j in 1:dim], 1:dim)


    # export ws
    function ws(w)
        dim = length(w)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, x[i=1:dim] >= 0)
        @NLconstraint(m, sum((x[i] - 1)^2 for i in 1:dim) <= 1)
        @objective(m, Min, w' * f(x))
        optimize!(m)
        value.(x)
    end

    function ps(v, d)
        dim = length(v)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, z)
        @variable(m, x[i=1:dim] >= 0)
        @NLconstraint(m, sum((x[i] - 1)^2 for i in 1:dim) <= 1)
        @NLconstraint(m, [i = 1:dim], f(x)[i] <= v[i] + z * d[i])
        @objective(m, Min, z)
        optimize!(m)
        map(x -> value.(x), (x, z))
    end

    # 
    # @return w^*
    # 
    function dps(v, d)
        dim = length(v)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, w[i=1:dim] >= 0)
        # ---
        @variable(m, t)
        @variable(m, x[i=1:dim] >= 0)
        @NLconstraint(m, sum((x[i] - 1)^2 for i in 1:dim) <= 1)
        @NLconstraint(m, t <= sum(w[i] * f(x)[i] for i in 1:dim))
        # ---
        @NLconstraint(m, sum(w[i] * d[i] for i in 1:dim) == 1)
        @objective(m, Max, t - w' * v)
        optimize!(m)
        value.(w)
    end
    # 
    # @return w^*
    # explicitly...
    # 
    function explicit_dps(v, d)
        dim = length(v)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, w[i=1:dim] >= 0)
        @NLconstraint(m, sum(w[i] * d[i] for i in 1:dim) == 1)
        @NLobjective(m, Max, (sum(w[i] * (1 + w[i] / (sqrt(sum(w[j]^2 for j in 1:dim)))) for i in 1:dim)) - sum(w[i] * v[i] for i in 1:dim))
        optimize!(m)
        value.(w)
    end

    return explicit ? (f, ws, ps, explicit_dps, W) : (f, ws, ps, dps, W)
end



# Example 3
# a ∈ {5, 7, 10, 20}
function example3(; a=5)
    f(x) = x
    # W = [[1, 0], [0, 1]]
    dim = 3
    W = map(i -> [i == j for j in 1:dim], 1:dim)


    # export ws
    function ws(w)
        dim = 3
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, x[i=1:dim] >= 0)
        @NLconstraint(m, (x[1] - 1)^2 + ((x[2] - 1) / a)^2 + ((x[3] - 1) / 5)^2 <= 1)
        @objective(m, Min, w' * f(x))
        optimize!(m)
        value.(x)
    end

    function ps(v, d)
        dim = length(v)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, z)
        @variable(m, x[i=1:dim] >= 0)
        @NLconstraint(m, (x[1] - 1)^2 + ((x[2] - 1) / a)^2 + ((x[3] - 1) / 5)^2 <= 1)
        @NLconstraint(m, [i = 1:dim], f(x)[i] <= v[i] + z * d[i])
        @objective(m, Min, z)
        optimize!(m)
        map(x -> value.(x), (x, z))
    end

    # 
    # @return w^*
    # 
    function dps(v, d)
        dim = length(v)
        m = Model(Ipopt.Optimizer)
        set_silent(m)
        @variable(m, w[i=1:dim] >= 0)
        # ---
        @variable(m, t)
        @variable(m, x[i=1:dim] >= 0)
        @NLconstraint(m, (x[1] - 1)^2 + ((x[2] - 1) / a)^2 + ((x[3] - 1) / 5)^2 <= 1)
        @NLconstraint(m, t <= sum(w[i] * f(x)[i] for i in 1:dim))
        # ---
        @NLconstraint(m, sum(w[i] * d[i] for i in 1:dim) == 1)
        @objective(m, Max, t - w' * v)
        optimize!(m)
        value.(w)
    end


    f, ws, ps, dps, W
end
