;density.pro

pflag=0			;plotting flag
sigma_1=4.6e-21		;cm^2/cm^3
r_max=3.0		;AU
z_max=3.0		;AU

;sigma's for three different dust populations, radians
sigma_pop = [0.000, 0.122, 0.576]

;ecliptic contributions by the three different dust populations
;OCCs, asteroids, JFCs
contrib   = [0.050, 0.450, 0.500]

;each population's radial power-law
nu_pop    = [2.000, 1.000, 1.45]

print,'Calculating the latitude distribution h(beta)...'
r2d=180.0/!pi
Po2=!pi/2.0
Nbeta=601
beta_tab=Po2*findgen(Nbeta)/(Nbeta-1)
common share, beta, sigma
.run latitude_integrand.pro
Npop=n_elements(sigma_pop)
h_tab=fltarr(Nbeta,Npop)
for j=0,Npop-1 do begin $
  sigma=sigma_pop(j) & $
  for i=0,Nbeta-1 do begin $
    beta=beta_tab(i) & $
    h_tab(i,j)=qromo('latitude_integrand',beta,Po2,double=0, $
      eps=1.0e-5) & $
  endfor & $
  h_tab(0,j)=h_tab(*,j)/h_tab(0,j) & $    ;normalize to h(0)=1
  print,float(j+1)/Npop & $
endfor
;some fixes
h_tab(Nbeta-1,0)=h_tab(Nbeta-2,*)
j=where(sigma_pop eq 0.0)
h_tab(*,j)=1.0
h_total=fltarr(Nbeta)
for j=0,Npop-1 do begin $
  h_tab(0,j)=contrib(j)*h_tab(*,j) & $
  h_total=h_total+h_tab(*,j) & $
endfor
window,xs=550,ys=550
plot,beta_tab*r2d,h_total,thick=1,xstyle=1,yrange=[0,1], $
  xtitle='latitude !7b!3    (degrees)',ytitle='h(!7b!3)', $
  charsize=1.3
for j=0,Npop-1 do oplot,beta_tab*r2d,h_tab(*,j)

;generate map of dust density(x,z)
Ngrid=101
x=r_max*findgen(Ngrid)/(Ngrid-1)
z=z_max*findgen(Ngrid)/(Ngrid-1)
density=fltarr(Ngrid,Ngrid)
beta=fltarr(Ngrid,Ngrid)
r=fltarr(Ngrid,Ngrid)
for i=0,Ngrid-1 do begin $
  for j=0,Ngrid-1 do begin $
    beta(i,j)=atan(z(j),x(i)) & $
    r(i,j)=sqrt(x(i)^2+z(j)^2) & $
  endfor & $
endfor 
for p=0,Npop-1 do begin $
  for i=0,Ngrid-1 do begin $
    for j=0,Ngrid-1 do begin $
      if ((i ne 0) OR (j ne 0)) then begin $
        db=abs(beta_tab-beta(i,j)) & $
        b=where(db eq min(db)) & $
        d=h_tab(b,p)/r(i,j)^nu_pop(p) & $
        density(i,j)=density(i,j)+d & $
      endif & $
    endfor & $
  endfor & $
  print,p & $
endfor
density(0,0)=3.0*max(density)
d2=fltarr(2*Ngrid-1,2*Ngrid-1)
d2(0,0)=rotate(density,2)
d2(Ngrid-1,Ngrid-1)=density
d2(Ngrid-1,0)=transpose(rotate(density,1))
d2(0,Ngrid-1)=transpose(rotate(density,3))
d2=d2*sigma_1
x2=[-rotate(x(1:*),2),x]
z2=[-rotate(z(1:*),2),z]
thck=2*pflag+1
lvls=sigma_1*1.5^(findgen(21)-10)
c_style=lvls*0.0
j=where(lvls eq sigma_1)
c_style(j)=1
if (pflag eq 0) then window,xs=550,ys=0.5*550 
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=16,ysize=9.0,xoffset=3,yoffset=15 & $
endif
contour,d2,x2,z2,xstyle=1,ystyle=1,charsize=1.3, $
  xtitle='x    (AU)',ytitle='z    (AU)',c_linestyle=c_style, $
  xrange=3*[-1,1],yrange=1.5*[-1,1],levels=lvls, $
  thick=thck,xthick=thck,ythick=thck,charthick=thck
if (pflag eq 1) then output_plot,'density.ps'
print,lvls/sigma_1

stop

;compare to fan model
nu=1.25
k1=1.5		;ok agreement is acheived for 1<k<2
k2=1.0
r(0)=r(1)
arg=-k1*(abs(sin(beta))^k2)
d_approx=sigma_1*exp(arg)/r^nu
contour,density*sigma_1,x,z,xstyle=1,ystyle=1,charsize=1.3, $
  xtitle='x    (AU)',ytitle='z    (AU)',c_linestyle=c_style, $
  xrange=2*[0,1],yrange=1.0*[0,1],levels=lvls, $
  thick=thck,xthick=thck,ythick=thck,charthick=thck
contour,d_approx,x,z,xstyle=5,ystyle=5,charsize=1.3, $
  xtitle='x    (AU)',ytitle='z    (AU)',c_linestyle=c_style, $
  xrange=2*[0,1],yrange=1.0*[0,1],levels=lvls,color=128, $
  thick=2,noerase=1
