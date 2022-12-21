function plot_xbar(Xbar0)
    _, n = size(Xbar0)
    x = Xbar0[:, 1]
    y = Xbar0[:, 2]
    ar = (x, y)
    if n > 2
        z = Xbar0[:, 3]
        ar = (ar..., z)
    end
    scatter(ar...)
end

function example_1_lhs(x)
    return x .|> (v -> (v - 1)^2) |> sum
end

function check_example_1(Xbar)

    m, _ = size(Xbar)
    check = map(i -> ~(example_1_lhs(Xbar[i, :]) ≈ 1), 1:m) |> sum

    if (check == 0)
        println("All $m solutions are feasible.")
    else
        println("$check solutions are not feasible.")

    end

end

function example_3_lhs(x, a)
    return (x[1] - 1)^2 + ((x[2] - 1) / a)^2 + ((x[3] - 1) / 5)^2
end

function check_example_3(Xbar, a)

    m, _ = size(Xbar)
    check = map(i -> ~(example_3_lhs(Xbar[i, :], a) ≈ 1), 1:m) |> sum

    if (check == 0)
        println("All $m solutions are feasible.")
    else
        println("$check solutions are not feasible.")

    end

end