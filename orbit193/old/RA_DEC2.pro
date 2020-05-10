pro RA_DEC2, plate_scl, x_opt, y_opt, x_star, y_star, RA_star, $
  DEC_star, twist_angle, RA_opt, DEC_opt, Nx, Ny, RA_field, $
  DEC_field

theta_star=asin(plate_scl*(y_star-y_opt))
phi_star=asin(plate_scl*(x_star-x_opt)/cos(theta_star))

f=-cos(twist_angle)*sin(phi_star)*cos(theta_star)- $
   sin(twist_angle)*sin(theta_star)
g= sin(twist_angle)*sin(phi_star)*cos(theta_star)- $
   cos(twist_angle)*sin(theta_star)
h= cos(phi_star)*cos(theta_star)
k=-g*sin(DEC_opt)+h*cos(DEC_opt)

DEC_star_try=g*cos(DEC_opt)+h*sin(DEC_opt)
DEC_star_try=-asin(DEC_star_try)		;sign problem?

y= f*cos(RA_opt)+k*sin(RA_opt)
x=-f*sin(RA_opt)+k*cos(RA_opt)
RA_star_try=atan(y,x)
l=where(RA_star_try lt 0.0)
if (l(0) ne -1) then RA_star_try(l)=RA_star_try(l)+2.0*!pi
r2d=180./!pi
for l=0,n_elements(x_star)-1 do $
  print,l,ra_star(l)*r2d,ra_star_try(l)*r2d, $
    dec_star(l)*r2d,dec_star_try(l)*r2d

end
