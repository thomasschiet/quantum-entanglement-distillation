#
# input: real 2n x 2n matrix transformed like in
# http://arxiv.org/pdf/1007.2905v2.pdf
#
# output: complex nxn matrix
function realToComplex(A)
  # check if the matrix is formatted properly!
  # that is, the northeast and southwest blocks should be the same
  # whereas the northwest and southeast block have opposing signs
  # furthermore, the dimensions must be divisable by 2

  dimensions = size(A)

  if dimensions[1] != dimensions[2]
    error("expected square matrix")
  end

  if dimensions[1] % 2 != 0
    error("dimension not divisable by 2")
  end



  realPart = A[1:dimensions[1]/2, 1:dimensions[1]/2]
  imagPart = A[dimensions[1]/2+1:dimensions[1], 1:dimensions[1]/2]

  return realPart + imagPart * im

end
