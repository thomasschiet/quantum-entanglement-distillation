#
# Outputs a ronald state of the form p EPR + (1-p) max mixed
#
# Inputs: p

function ronald2StateCorrPhase(p)

	integrand(φ) = ( 1/(2*pi) ) * ronald2State(p, φ) ⊗ ronald2State(p, φ) ;

	eps = 10.0^(-4)

	out = real( quadgk( integrand, 0, 2π; reltol=sqrt(eps), abstol=0, maxevals=10^7, order=7, norm=vecnorm)[1])

  return out

end
