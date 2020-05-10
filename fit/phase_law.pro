function phase_law, theta, nu

;compute the phase function Phi used by Hong(1985)
;as a function of scattering angle theta for the dust radial
;power-law index nu. This function is normalized so that
;Phi(theta=pi)=1.0 at the phase-angle pi-theta=0

;parameters for Henyey-Greenstein functions
g = [  0.700, -0.200, -0.810 ]
w = [  0.665,  0.330,  0.005 ]

;Phi1 &  Z for nu=1
Phi1=0.0
Z=0.0
p4=4.0*!pi
if (theta gt 1.0e-6) then begin $
  for k=0,2 do begin $
    if (theta lt !pi-1.0e-6) then begin $
      f1=g(k)*(sin(theta)^2)/(1.0-g(k)) & $
      f2=(1.0+g(k))/sqrt(1.0+g(k)^2-2.0*g(k)*cos(theta))-1.0 & $
      Z=Z+w(k)*f2/f1/p4 & $
    endif else begin $
      Z=Z+w(k)*0.5*(1.0-g(k))/((1.0+g(k))^2)/p4
    endelse & $
    Phi1=Phi1+w(k)*(1.0-g(k)^2)/(p4)/ $
      (1.0+g(k)^2-2.0*g(k)*cos(theta))^1.5 & $
  endfor & $
endif

;convert to nu<>1
Phi=(Phi1-(nu-1.0)*cos(theta)*Z)>0.0

;normalize such that Phi(theta=pi)=1.0
Phi1=0.0
Z=0.0
for k=0,2 do begin $
  Z=Z+w(k)*0.5*(1.0-g(k))/((1.0+g(k))^2)/p4 & $
  Phi1=Phi1+ w(k)*(1.0-g(k))/(p4)/(1.0+g(k))^2 & $
endfor
Phi_pi=(Phi1+(nu-1.0)*Z)>0.0
Phi=Phi/Phi_pi

return,Phi
end
