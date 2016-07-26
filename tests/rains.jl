facts("RainsProb") do
  context("2->1 distillation") do
    n = 2
    δ = 0.5
    ρ = sortAB(copies(full(ronaldState(1)), n), 2, n)
    (problem, F, p_succ) = RainsProb(ρ, 4, 4, 2, 0, δ, verbose = false)
    @fact F --> roughly(1, TOL) "max entangled state should output 1"
    @fact p_succ --> roughly(δ, TOL)

    δ = 0.5
    ρ = sortAB(copies(full(ronaldState(0.3)), n), 2, n)
    (problem, F, p_succ) = RainsProb(ρ, 4, 4, 2, 0, δ, verbose = false)
    @fact F --> roughly(0.785907642566565, TOL)
    @fact p_succ --> roughly(δ, TOL)

    δ = 0.5
    ρ = sortAB(copies(full(ronaldState(0)), n), 2, n)
    (problem, F, p_succ) = RainsProb(ρ, 4, 4, 2, 0, δ, verbose = false)
    @fact F --> roughly(0.5, TOL) " |11⟩⟨11| should give 0.5"
    @fact p_succ --> roughly(δ, TOL)

    ρ = sortAB(copies(full(ronaldState(0.7)), n), 2, n)
    (problem, F, p_succ) = RainsProb(ρ, 4, 4, 2, 0, 0.2, verbose = false)
    @fact F --> greater_than(0.5)
    @fact F --> less_than(1)
    @fact p_succ --> roughly(0.2, TOL)
    oldoptval = F
  end

  context("3->1 distillation") do
    n = 2
    ρ = sortAB(copies(full(ronaldState(0.5)), n), 2, n)
    (problem, F, p_succ) = RainsProb(ρ, 4, 4, 2, 0, 0.2, verbose = false)
    @fact p_succ --> roughly(0.2, TOL)
    oldoptval = F

    n = 3
    ρ = sortAB(copies(full(ronaldState(0.5)), n), 2, n)
    (problem, F, p_succ) = RainsProb(ρ, 8, 8, 2, 0, 0.2, verbose = false)
    @fact p_succ --> roughly(0.2, TOL)
    @fact F --> greater_than_or_equal(oldoptval)
  end
end
