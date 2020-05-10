;latitude_check.pro

;check the computed latitude distributions h(beta)

;Set pflag=1 to generate Postscript output
pflag=0

;fractional error allowed in numerical integration
err=1.0e-5

;sigma's for three different dust distributions
sigma_is=0.000			;isotropic dust
sigma_lo=0.122			;sigma=7.0 degrees
sigma_hi=0.576			;sigma=33 degrees

;fractional contribution by each population
is=0.05				;isotropic dust
lo=0.45				;sigma=7.0 degrees
hi=0.50				;sigma=33 degrees

;Get h(beta)
common share, beta, sigma
.run latitude_integrand.pro
Po2=!pi/2.0
r2d=180./!pi
N=201
b=Po2*findgen(N)/(N-1)
b(N-1)=0.999*Po2
h_is=fltarr(N)
h_lo=fltarr(N)
h_hi=fltarr(N)
print,'Computing the latitude distribution'
for j=0,N-1 do begin $
  beta=b(j) & $
  sigma=sigma_is & $
  h_is(j)=qromo('latitude_integrand',beta,Po2,double=0,eps=err)&$
  sigma=sigma_lo & $
  h_lo(j)=qromo('latitude_integrand',beta,Po2,double=0,eps=err)&$
  sigma=sigma_hi & $
  h_hi(j)=qromo('latitude_integrand',beta,Po2,double=0,eps=err)&$
  print,float(j+1)/N & $
endfor

;be sure to normalize h(beta=0)=1 !!!
c_is=1.0/h_is(0)
c_lo=1.0/h_lo(0)
c_hi=1.0/h_hi(0)
h_is=c_is*h_is
h_lo=c_lo*h_lo
h_hi=c_hi*h_hi
print,'c_is, c_lo, c_hi = ',c_is,c_lo,c_hi

;smoooth out slight jags in h_lo
h_lo=smooth(h_lo,5)

;plot h(beta) for each population
window,xs=550,ys=550
plot,b*r2d,h_is,xrange=[0,90],yrange=[0,1.1], $
  ystyle=1,xtitle='latitude !7b!3    (degrees)',charsize=1.3, $
  ytitle='latitude distribution h!lj!n(!7b!3)',xstyle=1
oplot,b*r2d,h_lo
oplot,b*r2d,h_hi

;plot sum over all contributions
thck=2*pflag+1
h_is=is*(fltarr(N)+1.0)
h_lo=lo*h_lo
h_hi=hi*h_hi
h=h_is+h_lo+h_hi
if (pflag eq 1) then setplot
plot,b*r2d,h,xrange=[0,90],yrange=[0,1.1],thick=2*thck, $
  ystyle=1,xtitle='latitude !7b!3    (degrees)',charsize=1.3, $
  ytitle='latitude distribution f!lj!nh!lj!n(!7b!3)',xstyle=1, $
  xthick=thck,ythick=thck,charthick=thck
oplot,b*r2d,h_is,thick=thck
oplot,b*r2d,h_lo,thick=thck
oplot,b*r2d,h_hi,thick=thck
xyouts,13.4,0.61,'total',charthick=thck,charsize=1.3
xyouts,25.5,0.275,'high',charthick=thck,charsize=1.3
xyouts,35.1,0.068,'isotropic',charthick=thck,charsize=1.3
xyouts,10.5,0.17,'low',charthick=thck,charsize=1.3
if (pflag eq 1) then output_plot,'latitude.ps'
f_comet=(h_is+h_hi)/h
k=where(f_comet gt 0.9)
print,'comets contribute > 90% at latitude beta (degrees) > ', $
  b(k(0))*r2d

;fractional integrated dust contributions
gamma_lo=int_tabulated(b,h_lo*cos(b))/lo
gamma_hi=int_tabulated(b,h_hi*cos(b))/hi
gamma_is=int_tabulated(b,h_is*cos(b))/is
print,'fractional integrated contributions' & $
print,'gamma_lo, gamma_hi, gamma_is = ', $
  gamma_lo, gamma_hi, gamma_is

;plot inclination distribution
g_is=(2.0/!pi)*fltarr(N)+is
g_lo=(2.0/!pi)*c_lo*lo*exp(-0.5*(b/sigma_lo)^2)
g_hi=(2.0/!pi)*c_hi*hi*exp(-0.5*(b/sigma_hi)^2)
g_total=g_is+g_lo+g_hi
plot,b*r2d,g_total*sin(b),yrange=[0,0.40],xstyle=1,thick=3, $
  xtitle='inclination i    (degrees)',charsize=1.3, $
  ytitle='inclination distribution g(i)    (radians!u-1!n)'
oplot,b*r2d,g_is*sin(b)
oplot,b*r2d,g_lo*sin(b)
oplot,b*r2d,g_hi*sin(b)
if (pflag eq 1) then setplot
plot,b*r2d,g_total,yrange=[0,4],xstyle=1,thick=2*thck, $
  charsize=1.3,xtitle='inclination i    (degrees)',ystyle=1, $
  ytitle='inclination distribution f!lj!ng!lj!n(i)/sin(i)    (radians!u-1!n)', $
  xthick=thck,ythick=thck,charthick=thck
oplot,b*r2d,g_is,thick=thck
oplot,b*r2d,g_lo,thick=thck
oplot,b*r2d,g_hi,thick=thck
xyouts,7.41,2.96,'total',charthick=thck,charsize=1.15
xyouts,2.20,0.75,'high',charthick=thck,charsize=1.15
xyouts,2.30,1.69,'low',charthick=thck,charsize=1.0
xyouts,30.0,0.11,'isotropic',charthick=thck,charsize=1.15
if (pflag eq 1) then output_plot,'inclination_dist.ps'

;compare to Figure 6 in Divine (1993)
plot,b*r2d,g_total*Po2,yrange=[0,6],xstyle=1,thick=2*thck, $
  charsize=1.3,xtitle='inclination i    (degrees)',ystyle=1, $
  ytitle='inclination distribution g(i)/sin(i)    (radians!u-1!n)', $
  xthick=thck,ythick=thck,charthick=thck
oplot,b*r2d,g_is*Po2,thick=thck
oplot,b*r2d,g_lo*Po2,thick=thck
oplot,b*r2d,g_hi*Po2,thick=thck


