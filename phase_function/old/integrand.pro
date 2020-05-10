function integrand,i

;Calculate the integrand used in the dust latitute distribution
;h(beta). This assumes that the inclination distribution g(i)
;has the form g(i)=(PI/2)*sin(i)*[a0+a1*exp(-0.5*(i/sigma1)^2/2)].
;Set a1=0 to isolate the isotropic contribution to h(beta), and
;a0=0 to isolate the gaussian contribution.

common share, b, a_iso, a_gss, s_gss

g=a_iso+a_gss*exp(-0.5d0*(i/s_gss)^2)
g=0.6366197723676d0*g*sin(i)
integrand=g/sqrt(sin(i)^2-sin(b)^2)

return,integrand
end

