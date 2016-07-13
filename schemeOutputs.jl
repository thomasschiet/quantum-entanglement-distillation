function bennett(; state = "Werner", p = false, F = false)
  if state == "Ronald"
    if p == false
      error("Ronald state given, expected parameter");
    end
    if F ≠ false
      error("Ronald state given, unexpected input fidelity")
    end

    F = (1+p)/2
    p = false
    # p = (2p + 1)/3
  end

  if p ≠ false
    F = (3p+1)/4
  end

  if F ≠ false
    p_succ = (F^2+2F*(1-F)/3 + 5((1-F)/3)^2)
    F_out = (F^2+((1-F)/3)^2)/(F^2+2F*(1-F)/3 + 5((1-F)/3)^2)
    return (p_succ, F_out)
  end
end

function deutsch(; state = "Werner", p = false, F = false)
  if state == "Werner"
    λ_00 = (3p + 1)/4
    λ_01 = λ_10 = λ_11 = (1 - p)/4
  end
  if state == "Ronald"
    λ_00 = (p + 1)/2
    λ_01 = (1 - p)/2
    λ_10 = λ_11 = 0
  end
  p_succ = N = (λ_00 + λ_11)^2 + (λ_01 + λ_10)^2
  λ_prime_11 = (λ_00^2+λ_11^2)/N

  return (p_succ, λ_prime_11)
end

state = "Ronald"
println("Computing output for ", state, " state:")
println("Bennett");
for p = collect(0.1:0.1:1)
  println("p=",p)
  println(bennett(state = state, p = p))
end

println("Deutsch");
for p = collect(0.1:0.1:1)
  println("p=",p)
  println(deutsch(state = state, p = p))
end

state = "Werner"
for p = collect(0.1:0.1:1)
  println("\\begin{tikzpicture}
      \\begin{axis}[
          title={Werner state \\(p=", p , "\\)},
          xlabel=Chance of success \\(p_{succ}\\),
          ylabel=Output fidelity \\(F'\\),
          xmin=0.5,
          xmax=1,
          ymin=0.4,
          ymax=1,
          minor xtick={0.1,0.3,...,0.9},
          minor ytick={0.45,0.55,...,0.95},
          grid=both,
          legend pos=south east]


          \\legend{Upper bound, Bennett, Deutsch}
          \\addplot +[mark=none] table[x=p_succ, y=fidelity] {data/werner/p0_", int(p*10) , ".dat};
          \\addplot +[only marks] coordinates {", bennett(state = state, p = p), "};
          \\addplot +[only marks] coordinates {", deutsch(state = state, p = p), "};
      \\end{axis}
  \\end{tikzpicture}
  \\clearpage")
end
