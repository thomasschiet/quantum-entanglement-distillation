#
# Compute the partial transpose of a bipartite state on the system index by j
#
export partialtranspose
export PT

function partialtranspose(ρ, nA::Int64, nB::Int64, j::Int64)
	n = ρ.size[1]

	# make ρ a n^2×1 vector.
	ρ_vec = reshape(ρ, n^2, 1)
	ρ_vec = eye(n^2, n^2) * ρ_vec
	ρ = reshape(ρ_vec, n, n)
	return ρ;

	# Check the dimensions of ρ
	totalDim = nA * nB;
	(da,db) = size(ρ);

	@assert da == db "The state must be square matrix."
	@assert da == totalDim "Dimensions of the input state do not match."

	# Shuffle the entries around
	outρ = ρ;

	print(outρ)

	for r = 0:(nA-1)
	for s = 0:(nA-1)
		offsetRow = r * nA;
		offsetCol = s * nA;
		for k = 1:nB
		for l = 1:nB
			p1 = offsetRow + k;
			p2 = offsetCol + l;
			q1 = offsetRow + l;
			q2 = offsetCol + k;
			outρ[q1,q2] == ρ[p1,p2];
			outρ[p1,p2] == ρ[q1,q2];
		end
		end
	end
	end

	return outρ;
end

function addSwap(matrix, coord1, coord2)
  col = round(Int, (coord1[1] - 1)*sqrt(size(matrix)[1]) + coord1[2])
  row = round(Int, (coord1[2] - 1)*sqrt(size(matrix)[1]) + coord2[2])

  matrix[row, col] = 1
  matrix[col, row] = 1

  return matrix
end

function PTmatrix(shape, blocksize)
  println(shape)
  println(blocksize)
  if (shape[1] % blocksize[1] ≠ 0 || shape[2] % blocksize[2] ≠ 0)
    error("incorrect blocksizeensions/shape")
  end

  if (shape[1] ≠ shape[2])
    error("currently only square matrices are supported")
  end

  n = shape[1]
  m = shape[2]

  matrix = spzeros(n*m, n*m)

  for i = 1:n
    for j = 1:m
      block = (round(Int, floor((i-1)/blocksize[1])), round(Int, floor((j-1)/blocksize[2])))
      B = blocksize[1]
      coord = (i, j)
      newcoords = (j - B*block[2] + B * block[1], i - B*block[1] + B * block[2])
      # println(coord, " in block ", block)
      # println("swapping ", coord, " with ", newcoords)
      matrix = addSwap(matrix, coord, newcoords)
    end
  end

  return matrix
end

function PT(matrix, blocksize = false)
  # return matrix
  # if (blocksize ≡ false)
    # blocksize = (sqrt(size(matrix)[1]), sqrt(size(matrix)[1]))
  # end

  M_size = size(matrix)
  M_vec = reshape(matrix, prod(size(matrix)), 1)
  # println("Make permutation matrix")
  PTM = PTmatrix(size(matrix), blocksize)
  # println("Do calculation")
  M_vec =  PTM * M_vec
  # println("Finished calculation")
  return reshape(M_vec, M_size[1], M_size[2])
end
