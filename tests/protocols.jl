facts("Protocols") do
  (p1, p2) = sort([rand(), rand()])

    state1 = ronaldState(p1)
    state2 = ronaldState(p2)
    @fact deutsch(state1)[1] --> less_than_or_equal(deutsch(state2)[1])
    @fact bennett(state1)[1] --> less_than_or_equal(bennett(state2)[1])

    state1 = ronald2State(p1)
    state2 = ronald2State(p2)
    @fact deutsch(state1)[1] --> less_than_or_equal(deutsch(state2)[1])
    @fact bennett(state1)[1] --> less_than_or_equal(bennett(state2)[1])

    state1 = wernerState(p1)
    state2 = wernerState(p2)
    @fact deutsch(state1)[1] --> less_than_or_equal(deutsch(state2)[1])
    @fact bennett(state1)[1] --> less_than_or_equal(bennett(state2)[1])
end
