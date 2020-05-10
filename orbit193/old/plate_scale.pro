pro plate_scale, xc, yc, x, y, RA, DEC, plate_scl, errors

;Input the optical center coordinates (xc,yc), an array
;containing the stars' x-coordinates, y-coordinate,
;right ascension RA, and declination DEC. There must be data
;for at least 2 stars. The outputs are the plate_scale
;and the rms errors (in units of plate_scale, or pixels).
;All angles have units of radians.

;solve for the plate scale using these trial solutions ps
Nps=5001
ps0=0.00137
range=0.2
ps_try=ps0*(1.0+range*(2.0*findgen(Nps)/(Nps-1)-1.0))
errors=fltarr(Nps)

Nstars=n_elements(x)
for i=0,Nps-1 do begin $
  dx=ps_try(i)*(x-xc) & $
  dy=ps_try(i)*(y-yc) & $
  theta=asin(dy) & $		;for spherical lens
  phi=asin(dx/cos(theta)) & $	;for spherical lens
  for star0=0,Nstars-1 do begin $
    for star1=star0+1,Nstars-1 do begin $
      angle=acos(cos(DEC(star0))*cos(DEC(star1))* $
        cos(RA(star0)-RA(star1)) + $
        sin(DEC(star0))*sin(DEC(star1))) & $
      angle=acos(cos(DEC(star0))*cos(DEC(star1))* $
        cos(RA(star0)-RA(star1)) + $
        sin(DEC(star0))*sin(DEC(star1))) & $
      angle_try=acos(cos(theta(star0))*cos(theta(star1))* $
        cos(phi(star0)-phi(star1)) + $
        sin(theta(star0))*sin(theta(star1))) & $
      errors(i)=errors(i)+total((angle-angle_try)^2) & $
    endfor & $
  endfor & $
endfor
j=where(errors eq min(errors))
plate_scl=ps_try(j(0))
errors=sqrt(errors)/plate_scl

end

