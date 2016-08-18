v0 = zeros(8)
v = Any[]
for i = 1:8
  x = copy(v0)
  x[i] = 1
  push!(v, x)
end

u = Any[]
for i = 1:8
  x = kron( v[i] , v[i] )
  push!(u, x)
end

for j = 1:7
  for i = 1:8-j
    println(j*(17-j)/2+i)
    # u[j*(17-j)/2+i] = ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) + kron( v[i+j] , v[j] ) )
    push!(u, ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) + kron( v[i+j] , v[j] ) ))
  end
end

for j= 1: 7
  for i = 1:8-j
    push!(u, ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) - kron( v[i+j] , v[j] )))
  end
end

Pcb = spzeros(512, 512)
for i = 1:64
  for j = 1:8
    Pcb[:, 64 * (j - 1) + i] = kron(v[j], u[i])
  end
end

Pcb
