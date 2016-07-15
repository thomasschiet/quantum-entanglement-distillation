facts("Rains") do
  context("2->1 distillation") do
    n = 2

    ρ = sortAB(copies(full(ronaldState(1)), n), 2, n)
    (problem, F) = Rains(ρ, 4, 4, 2, verbose = false)
    @fact problem.optval --> roughly(1, TOL)

    ρ = sortAB(copies(full(ronaldState(0)), n), 2, n)
    (problem, F) = Rains(ρ, 4, 4, 2, verbose = false)
    @fact problem.optval --> roughly(0.5, TOL)

    ρ = sortAB(copies(full(ronaldState(0.7)), n), 2, n)
    (problem, F) = Rains(ρ, 4, 4, 2, verbose = false)
    @fact problem.optval --> greater_than(0.5)
    @fact problem.optval --> less_than(1)
    oldoptval = problem.optval
  end

  context("3->1 distillation") do
    n = 2
    ρ = sortAB(copies(full(ronaldState(0.7)), n), 2, n)
    (problem, F) = Rains(ρ, 4, 4, 2, verbose = false)
    oldoptval = problem.optval

    n = 3
    ρ = sortAB(copies(full(ronaldState(0.7)), n), 2, n)
    (problem, F) = Rains(ρ, 8, 8, 2, verbose = false)
    @fact problem.optval --> greater_than_or_equal(oldoptval)
  end
end
