function latitude_integrand,i

;Calculate the integrand used in the dust latitute distribution
;h(beta). This assumes that the inclination distribution g(i)
;has the form g(i)=(PI/2)*sin(i)*exp(-0.5*(i/sigma)^2/2).
;Set the common variables, latitude beta and standard deviation
;sigma, before calling qromo (which in turn calls this function).

common share, beta, sigma

intgrnd=0.0
if (sigma le 0.0) then intgrnd=1.0
if ((sigma gt 0.0) and (i lt 7.0*sigma)) then $
  intgrnd=exp(-0.5*(i/sigma)^2)
factor=sin(i)^2-sin(beta)^2
if (factor gt 0.0) then begin $
  factor=sin(i)/sqrt(factor) & $
endif else begin $
  factor=1.0 & $
endelse
intgrnd=intgrnd*0.636620*factor

return,intgrnd
end

