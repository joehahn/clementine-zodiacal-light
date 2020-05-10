;inclination_dist.pro

;set up integrand function
common share, beta, sigma
.run integrand.pro

window,xs=550,ys=550
Po2=!pi/2.0
r2d=180./!pi
N=21
b=Po2*findgen(N)/N

;isotropic
h0=fltarr(N)
sigma=0.0
for j=0,N-1 do begin $
  beta=b(j) & $
  h0(j)=qromo('integrand',beta,Po2,double=0,eps=1.0e-6) & $
endfor
h0_approx=b*0+1.0

;10 degrees
h10=fltarr(N)
sigma=10.0/r2d
for j=0,N-1 do begin $
  beta=b(j) & $
  h10(j)=qromo('integrand',beta,Po2,double=0,eps=1.0e-6) & $
endfor
h10_approx=0.8*sigma*exp(-0.5*(b/sigma)^2)

;30 degrees
h30=fltarr(N)
sigma=30.0/r2d
for j=0,N-1 do begin $
  beta=b(j) & $
  h30(j)=qromo('integrand',beta,Po2,double=0,eps=1.0e-6) & $
endfor
h30_approx=0.8*sigma*exp(-0.5*(b/sigma)^2)

;60 degrees
h60=fltarr(N)
sigma=60.0/r2d
for j=0,N-1 do begin $
  beta=b(j) & $
  h60(j)=qromo('integrand',beta,Po2,double=0,eps=1.0e-6) & $
endfor
h60_approx=0.8*sigma*exp(-0.5*(b/sigma)^2)

;plots
plot,b*r2d,h10,psym=8,xrange=[0,90],yrange=[0,1], $
  ystyle=1,xtitle='latitude !7b!3    (degrees)',charsize=1.3, $
  ytitle='latitude distribution h!lj!n(!7b!3)',xstyle=1
oplot,b*r2d,h10_approx
oplot,b*r2d,h30,psym=8
oplot,b*r2d,h30_approx
oplot,b*r2d,h60,psym=8
oplot,b*r2d,h60_approx
;oplot,b*r2d,h0,psym=8
;oplot,b*r2d,h0


