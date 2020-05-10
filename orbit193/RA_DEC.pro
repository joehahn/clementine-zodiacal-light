pro RA_DEC, plate_scl, x_opt, y_opt, x_star, y_star, RA_star, $
  DEC_star, twist_angle, RA_opt, DEC_opt, Nx, Ny, RA_field, $
  DEC_field

;Input the plate-scale (plate_scl), the x-y coordinates for the optical
;axis (x_opt,y_opt), the x-y coordinates for a list of stars
;(x_stars,y_stars), their right ascension & declination
;(RA_star,DEC_star), and initial guesses at the twist_angle, the
;optical center's ra & dec (RA_opt,DEC_opt), and the array size
;(Nx,Ny). The outputs are more accurate estimates of 
;(twist_angle, RA_opt, DEC_opt) and the fields ra & dec
;(RA_field,DEC_field).

delta_angle=0.45 		;radians
twist_angle0=twist_angle
RA_opt0=RA_opt
DEC_opt0=DEC_opt
theta_star=asin(plate_scl*(y_star-y_opt))
phi_star=asin(-plate_scl*(x_star-x_opt)/cos(theta_star))
N=21
errors=fltarr(N,N,N)
Nstars=n_elements(x_star)
r2d=180./!pi

;scan array of possible solutions for 
; RA_opt, DEC_opt, and twist_angle
while (delta_angle gt 0.015*plate_scl) do begin $
  twist_angle=twist_angle0+delta_angle*2.0*(findgen(N)/(N-1)-0.5) & $
  RA_opt=RA_opt0+delta_angle*2.0*(findgen(N)/(N-1)-0.5) & $
  DEC_opt=DEC_opt0+delta_angle*2.0*(findgen(N)/(N-1)-0.5) & $
  for i=0,N-1 do begin $		;twist_angle index
    for j=0,N-1 do begin $		;RA_opt index
      for k=0,N-1 do begin $		;DEC_opt index

        f=cos(twist_angle(i))*sin(phi_star)*cos(theta_star)- $
           sin(twist_angle(i))*sin(theta_star) & $
        g= sin(twist_angle(i))*sin(phi_star)*cos(theta_star)+ $
           cos(twist_angle(i))*sin(theta_star) & $
        h= cos(phi_star)*cos(theta_star) & $
        kk=-g*sin(DEC_opt(k))+h*cos(DEC_opt(k)) & $

        DEC_star_try=g*cos(DEC_opt(k))+h*sin(DEC_opt(k)) & $
        DEC_star_try=asin(DEC_star_try) & $

        y= f*cos(RA_opt(j))+kk*sin(RA_opt(j)) & $
        x=-f*sin(RA_opt(j))+kk*cos(RA_opt(j)) & $
        RA_star_try=atan(y,x) & $
        l=where(RA_star_try lt 0.0) & $
        if (l(0) ne -1) then $
          RA_star_try(l)=RA_star_try(l)+2.0*!pi & $

        errors(i,j,k)=sqrt(total((RA_star-RA_star_try)^2+ $
          (DEC_star-DEC_star_try)^2))/plate_scl & $

      endfor & $
    endfor & $
  endfor & $
  index=where(errors eq min(errors)) & $
  k=index(0)/(N*N) & $
  j=index(0)/N-N*k & $
  i=index(0)-N*j-N*N*k & $
  if ((i eq 0) or (i eq N-1)) then print,'RA_DEC: did not converge to solution.'
  if ((j eq 0) or (j eq N-1)) then print,'RA_DEC: did not converge to solution.'
  if ((k eq 0) or (k eq N-1)) then print,'RA_DEC: did not converge to solution.'
  twist_angle0=twist_angle(i) & $
  RA_opt0=RA_opt(j) & $
  DEC_opt0=DEC_opt(k) & $
  ;print,delta_angle/plate_scl,twist_angle0*r2d,RA_opt0*r2d, $
  ;  DEC_opt0*r2d,errors(i,j,k) & $
  delta_angle=3.0*delta_angle/N & $
endwhile
twist_angle=twist_angle0
RA_opt=RA_opt0
DEC_opt=DEC_opt0

;print errors in solution
for l=0,Nstars-1 do begin $
  dra=RA_star(l)-RA_star_try(l) & $
  ddec=DEC_star(l)-DEC_star_try(l) & $
  print,'star ',l,' registration error = ', $
    sqrt(dra^2+ddec^2)/plate_scl,' pixels' & $
endfor

;get field RA & DEC
xx=findgen(Nx)
x=fltarr(Nx,Ny)
for i=0,Ny-1 do x(0,i)=xx
yy=transpose(findgen(Ny))
y=fltarr(Nx,Ny)
for i=0,Nx-1 do y(i,0)=yy
theta=asin(plate_scl*(y-y_opt))
phi=asin(-plate_scl*(x-x_opt)/cos(theta))
f= cos(twist_angle)*sin(phi)*cos(theta)- $
   sin(twist_angle)*sin(theta)
g= sin(twist_angle)*sin(phi)*cos(theta)+ $
   cos(twist_angle)*sin(theta)
h= cos(phi)*cos(theta)
k=-g*sin(DEC_opt)+h*cos(DEC_opt)
DEC_field=g*cos(DEC_opt)+h*sin(DEC_opt)
DEC_field=asin(DEC_field)
y= f*cos(RA_opt)+k*sin(RA_opt)
x=-f*sin(RA_opt)+k*cos(RA_opt)
RA_field=atan(y,x)
l=where(RA_field lt 0.0)
if (l(0) ne -1) then RA_field(l)=RA_field(l)+2.0*!pi

end
