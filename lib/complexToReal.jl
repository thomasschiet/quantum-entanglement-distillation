#
# Makes a real matrix from a complex one
# http://arxiv.org/pdf/1007.2905v2.pdf
#
# input: Complex matrix
#
# output: Real matrix in the form
# [Re -Im; Im Re]
#
function complexToReal(A)
  realA = real(A + conj(A))/2
  imagA = real(A - conj(A))/2
  return [realA -imagA; imagA realA]
end
