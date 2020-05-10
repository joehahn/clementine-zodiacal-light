function integrand,i

;Calculate the integrand used in the dust latitute distribution
;h(beta). This assumes that the inclination distribution g(i)
;has the form g(i)=(PI/2)*sin(i)*exp(-0.5*(i/sigma)^2/2).
;Set the common variables, latitude beta and standard deviation
;sigma, before calling qromo (which in turn calls this funciton).

common share, beta, sigma

intrgnd=1.0
if (sigma ne 0.0) then intrgnd=exp(-0.5*(i/sigma)^2)
intrgnd=0.636620*intrgnd*sin(i)/sqrt(sin(i)^2-sin(beta)^2)

return,intrgnd
end

