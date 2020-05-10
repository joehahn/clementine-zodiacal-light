;inc_distribution3.pro

Po2=1.570796326795d0
r2d=180.0/!pi
k=1.15
sigma_lo=10.0/r2d
sigma_hi=30.0/r2d

common share, b, a_iso, a_gss, s_gss
.run integrand.pro
window,xs=501,ys=501

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
a_is=0.305		;isotropic contribution
a_lo=1.73		;low-inclination contribution
a_hi=1.0		; hi-inclination contribution
h_fit=a_is*h_iso+a_lo*h_lo+a_hi*h_hi	;total fitted h
h_exp=exp(-k*sin(beta))	;expected h
plot,beta*r2d,h_exp,xstyle=1,yrange=[0.0,1.1], $
  ystyle=1,xtitle='latitude !7b!3    (degrees)', $
  ytitle='h(!7b!3)',charsize=1.3,nodata=1
oplot,beta*r2d,h_exp,thick=8,color=128
oplot,beta*r2d,h_fit,thick=4
oplot,beta*r2d,a_is*h_iso,linestyle=0,color=196
oplot,beta*r2d,a_lo*h_lo,linestyle=0,color=196
oplot,beta*r2d,a_hi*h_hi,linestyle=0,color=196
xyouts,19.8,0.74,'total',charsize=1.3
xyouts,35.2,0.32,'isotropic',charsize=1.3
xyouts,55.0,0.11,'high-i',charsize=1.3
xyouts,13.3,0.11,'low-i',charsize=1.3

;plot inclination distribution
g_is=a_is*sin(beta)/Po2
g_lo=a_lo*sin(beta)*exp(-0.5*(beta/sigma_lo)^2)/Po2
g_hi=a_hi*sin(beta)*exp(-0.5*(beta/sigma_hi)^2)/Po2
g_tot=g_is+g_lo+g_hi
plot,beta*r2d,g_tot,xstyle=1,charsize=1.3, $
  xtitle='inclination !7i!3    (degrees)', $
  ytitle='g(i)',thick=3,yrange=[0,0.35],ystyle=1
oplot,beta*r2d,g_is
oplot,beta*r2d,g_lo
oplot,beta*r2d,g_hi
xyouts,31.5,0.31,'total',charsize=1.3
xyouts,69.1,0.17,'isotropic',charsize=1.3
xyouts,25.2,0.20,'high-i',charsize=1.3
xyouts,25.4,0.03,'low-i',charsize=1.3

;plot inclination distribution
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


