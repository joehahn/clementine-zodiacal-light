;inc_distribution3.pro

;set pflag=1 to generate postscript figures
pflag=0
if (pflag eq 0) then window,xs=501,ys=501

;constants
Po2=1.570796326795d0
r2d=180.0/!pi

;set parameters
k=1.15			;isophot ellipticity from 3D_phase_integral.pro
sigma_lo=10.0/r2d	;1-sigma inclination for low-i population
sigma_hi=30.0/r2d	;1-sigma inclination for  hi-i population

;set up integrand function
common share, b, a_iso, a_gss, s_gss
.run integrand.pro

;compute h(beta) for isotropic, low-i, and high-i populations
N=101			
beta=Po2*dindgen(N)/(N-1)
h_iso=dblarr(N)
h_lo=dblarr(N)
h_hi=dblarr(N)
for i=0,N-2 do begin $

  b=beta(i) & $

  ;get isotropic contribution
  a_iso=1.d0 & $
  a_gss=0.d0 & $
  s_gss=2.0d0 & $
  h_iso(i)=qromo('integrand',b,Po2,double=0) & $

  ;get low-i gaussian contribution
  a_iso=0.d0 & $
  a_gss=1.d0 & $
  s_gss=sigma_lo & $
  h_lo(i)=qromo('integrand',b,Po2,double=0) & $

  ;get hi-i gaussian contribution
  a_iso=0.d0 & $
  a_gss=1.d0 & $
  s_gss=sigma_hi & $
  h_hi(i)=qromo('integrand',b,Po2,double=0) & $

  print,i,b*r2d,h_iso(i),h_lo(i),h_hi(i) & $
endfor
h_iso(N-1)=h_iso(N-2)
h_lo(N-1)=h_lo(N-2)
h_hi(N-1)=h_hi(N-2)

;plot fit to h_exp
a_is=0.305		;isotropic coefficient
a_lo=1.73		;low-inclination coefficient
a_hi=1.0		; hi-inclination coefficient
h_i=a_is*h_iso		;isotropic contribution
h_l=a_lo*h_lo		;low-inclination contribution
h_h=a_hi*h_hi		; hi-inclination contribution
h_fit=h_i+h_l+h_h	;total fitted h
h_exp=exp(-k*sin(beta))	;expected h
thck=1
if (pflag eq 1) then setplot
plot,beta*r2d,h_exp,xstyle=1,yrange=[0.0,1.1], $
  ystyle=1,xtitle='latitude !7b!3    (degrees)', $
  xthick=thck,ythick=thck,charthick=thck, $
  ytitle='latitude distribution h(!7b!3)',charsize=1.3,nodata=1
oplot,beta*r2d,h_exp,thick=10+thck,color=128
oplot,beta*r2d,h_fit,thick=3+thck
oplot,beta*r2d,h_i,linestyle=0,thick=thck
oplot,beta*r2d,h_l,linestyle=0,thick=thck
oplot,beta*r2d,h_h,linestyle=0,thick=thck
xyouts,19.8,0.74,'total',charsize=1.3,charthick=thck
xyouts,35.2,0.32,'isotropic',charsize=1.3,charthick=thck
xyouts,55.0,0.11,'high-i',charsize=1.3,charthick=thck
xyouts,13.3,0.11,'low-i',charsize=1.3,charthick=thck
if (pflag eq 1) then output_plot,'latitude_dist.ps'

;plot inclination distributiong_is=a_is*sin(beta)/Po2
g_is=a_is*sin(beta)/Po2+beta*0.0
g_lo=a_lo*sin(beta)*exp(-0.5*(beta/sigma_lo)^2)/Po2
g_hi=a_hi*sin(beta)*exp(-0.5*(beta/sigma_hi)^2)/Po2
g_tot=g_is+g_lo+g_hi
if (pflag eq 1) then setplot
plot,beta*r2d,g_tot,xstyle=1,charsize=1.3, $
  xtitle='inclination !7i!3    (degrees)',thick=3, $
  xthick=thck,ythick=thck,charthick=thck, $
  ytitle='inclination distribution g(i)',yrange=[0,0.4],ystyle=1
oplot,beta*r2d,g_is,thick=thck
oplot,beta*r2d,g_lo,thick=thck
oplot,beta*r2d,g_hi,thick=thck
xyouts,30.0,0.305,'total',charsize=1.3,charthick=thck
xyouts,73.0,0.17,'isotropic',charsize=1.3,charthick=thck
xyouts,24.0,0.205,'high-i',charsize=1.3,charthick=thck
xyouts,25.1,0.03,'low-i',charsize=1.3,charthick=thck
if (pflag eq 1) then output_plot,'inclination_dist.ps'

;get the relative abundances in the ecliptic
h0_total=h_i(0)+h_l(0)+h_h(0)
print,'isotropic contribution at beta=0 ',h_i(0)/h0_total & $
print,'low-incln contribution at beta=0 ',h_l(0)/h0_total & $
print,' hi-incln contribution at beta=0 ',h_h(0)/h0_total

;get the total abundances
N_tot=int_tabulated(beta,cos(beta)*h_fit)
N_iso=int_tabulated(beta,cos(beta)*h_i)
N_low=int_tabulated(beta,cos(beta)*h_l)
N_hi=int_tabulated(beta,cos(beta)*h_h)
print,'integrated isotropic contribution = ',N_iso/N_tot & $
print,'integrated low-incln contribution = ',N_low/N_tot & $
print,'integrated  hi-incln contribution = ',N_hi/N_tot

stop

;plot inclination distribution/sin(i)
g_is=a_is/Po2+beta*0.0
g_lo=a_lo*exp(-0.5*(beta/sigma_lo)^2)/Po2
g_hi=a_hi*exp(-0.5*(beta/sigma_hi)^2)/Po2
g_tot=g_is+g_lo+g_hi
plot,beta*r2d,g_tot,xstyle=1,charsize=1.3, $
  xtitle='inclination !7i!3    (degrees)', $
  ytitle='g(i)/sin(i)',thick=3
oplot,beta*r2d,g_is
oplot,beta*r2d,g_lo
oplot,beta*r2d,g_hi


