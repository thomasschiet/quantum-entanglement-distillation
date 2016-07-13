# Calculate the distance to Q-space as
# described http://link.springer.com/chapter/10.1007%2F978-3-642-31594-7_67
# doesn't actually work yet though; it gives negative values.

function q(ρ, Q_base)
  M = length(Q_base)
  result = Vector(M)

  for i = 1:M
    result[i] = trace(Q_base[i]*ρ)
  end

  return result
end

# algorithm to find distance
function dis(p, ϵ, d, Q_base, w)
  M = length(Q_base)

  γ = ϵ/(8M*w)
  T = int(ceil(log(d) / γ^2))

  print(T)
  W = eye(d)
  ρ = Vector(T)
  N = Vector(T)

  for t = 1:T
    # part a
    ρ[t] = W / trace(W)
    q_vec = q(ρ[t], Q_base)
    print(p)
    print(q_vec)
    s_vec = p - q_vec

    z = Vector(length(p))
    y = Vector(length(p))
    i = 1
    for c in s_vec
      z[i] = exp(-im*angle(c))
      y[i] = conj(z[i])
      i = i + 1
    end

    # part b
    Q = zeros(size(Q_base[1]))
    for i = 1:M
      Q += y[i] * Q_base[i]
    end

    N[t] = real(p⋅z) * eye(d) - 0.5 * (Q + conj(Q)) + 2M*w*eye(d)

    # part c
    W_exponent = zeros(size(W))
    for τ = 1:t
      W_exponent = W_exponent - γ * N[τ]
    end
    W = expm(W_exponent)
  end

  result = 0
  for t = 1:T
    result = result + trace(conj(ρ[t]) * N[t] - 2M*w*eye(d))
  end

  return result
end

id = [1 0; 0 1]
p_1 = [0 1; 1 0]
p_2 = [0 -im; im 0]
p_3 = [1 0; 0 -1]

Pauli_base = (id, p_1, p_2, p_3)

p = [1, 2, 5, 3]
for c in p
  print(c)
end


Q = zeros(size(Pauli_base[1]))

dis(p, 2, 2, Pauli_base, 2)
