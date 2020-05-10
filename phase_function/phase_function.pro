;phase_function.pro

;plot the phase function given in Hong (1985) A&A v146, p67
;versus scattering angle.

;set pflag=1 to generate postscript figure
pflag=0

;parameters for Henyey-Greenstein functions
g = [  0.700, -0.200, -0.810 ]
w = [  0.665,  0.330,  0.005 ]

;scattering angle theta
N=1001
theta=!pi*findgen(N)/(N-1)
theta(0)=0.5*theta(1)
r2d=180.0/!pi
p4=4.0*!pi
t90=where(theta ge !pi/2.0)
t90=t90(0)
theta(t90)=!pi/2.0

;Phi1 & Z
Phi1=fltarr(N)
Z=fltarr(N)
for k=0,2 do begin $
  f1=g(k)*(sin(theta)^2)/(1.0-g(k)) & $
  f2=(1.0+g(k))/sqrt(1.0+g(k)^2-2.0*g(k)*cos(theta))-1.0 & $
  ;get value at theta=!pi
  f1(N-1)=1.0 & $
  f2(N-1)=0.5*(1-g(k))/(1.0+g(k))^2 & $
  Z=Z+w(k)*f2/f1/p4 & $
  Phi1=Phi1 +  w(k)*(1.0-g(k)^2)/(p4)/ $
    (1.0+g(k)^2-2.0*g(k)*cos(theta))^1.5 & $
endfor

;compare Hong's Fig. 2 to the Clementine Z(e)
window,xs=650,ys=650
!p.multi=[0,1,2,0]
plot_io,theta*r2d,Z*200/Z(N-1),xrange=[0,180],xstyle=1, $
  yrange=[100,20000],ystyle=1,title='Fig. 2 in Hong (1985)', $
  xtitle='elongation !7e!3    (degrees)',charsize=1.5, $
  ytitle='ecliptic surface brightness Z(!7e!3)'
;clementine surface brightness normalized to Hong's at eln=15 deg
Z_clem=1.0/abs(sin(theta))^2.45
t15=where(theta*r2d ge 15.0)
t15=t15(0)
Z_clem=Z_clem*Z(t15)/Z_clem(t15)
plot_io,theta*r2d,Z/Z(t15),xrange=[0,30],xstyle=1, $
  yrange=[0.01,100],ystyle=1,charsize=1.5, $
  xtitle='elongation !7e!3    (degrees)', $
  ytitle='ecliptic surface brightness Z(!7e!3)'
oplot,theta*r2d,Z_clem/Z(t15),linestyle=1
!p.multi=0

;generate Hong's Fig. 5
plot,theta*r2d,alog10(Phi1),yrange=[-1.7,0],xstyle=1,ystyle=1, $
  xrange=[0,180]
nu=1.2
Phi=(Phi1-(nu-1.0)*cos(theta)*Z)
Phi=Phi*Phi1(t90)/Phi(t90)
oplot,theta*r2d,alog10(Phi>1.0e-6)
nu=1.3
Phi=(Phi1-(nu-1.0)*cos(theta)*Z)
Phi=Phi*Phi1(t90)/Phi(t90)
oplot,theta*r2d,alog10(Phi>1.0e-6)
nu=1.45
Phi=(Phi1-(nu-1.0)*cos(theta)*Z)
Phi=Phi*Phi1(t90)/Phi(t90)
oplot,theta*r2d,alog10(Phi>1.0e-6)

;Phi for dust nu=1.45
thck=2*pflag+1
nu=1.45
Phi=(Phi1-(nu-1.0)*cos(theta)*Z)>0.0
Phi=Phi/Phi(N-1)
if (pflag eq 1) then setplot
plot,theta*r2d,Phi,xstyle=1,xrange=[0,180],yrange=[0,1], $
  charsize=1.3,xtitle='scattering angle !9P!3    (degrees)', $
  ytitle='phase law !7w(!9P!3)    and    !7w(!9P!3)sin!u!7m!n(!9P!3)', $
  thick=thck,xthick=thck,ythick=thck,charthick=thck
xyouts,142,0.735,'!7w(!9P!3)',charthick=thck,charsize=1.0
j=where(Phi gt 0.0)
oplot,theta(j)*r2d,Phi(j)*abs(sin(theta(j)))^nu,linestyle=2, $
  thick=thck
xyouts,142,0.37,'!7w(!9P!3)sin!u!7m!n(!9P!3)', $
  charthick=thck,charsize=1.0;,color=128
if (pflag eq 1) then output_plot,'phase_function.ps'

;compare to Herb's phase law
Phi_zook=exp(-0.77*(!pi-theta))
;oplot,theta*r2d,Phi_zook,linestyle=1
;oplot,theta*r2d,Phi_zook*abs(sin(theta(j)))^nu,linestyle=2

;
intgrnd=Phi*abs(sin(theta))^nu
intgrl=fltarr(N/2)
chi=0.5*!pi*findgen(N/2)/(N/2-1)
for i=0,N/2-2 do begin $
  j=where((theta ge chi(i)) and (theta le !pi-chi(i))) & $
  if (j(0) ne -1) then $
    intgrl(i)=int_tabulated(theta(j),intgrnd(j)) & $
endfor
intgrl=intgrl/intgrl(0)
plot,chi*r2d,intgrl,xstyle=1,charsize=1.3, $
  xtitle='!7v!3    (degrees)', $
  ytitle='fractional contribution to integral Z(!7v!3)'
oplot,chi*r2d,chi*0+0.9,color=128
j=where(intgrl lt 0.9)
print,'90% contribution to Z occurs at chi (degrees) < ', $
  chi(j(0))*r2d

;compare the surface brightness inferred from Hong's
;phase function to the clementine surface brightness,
;normalized at 15 degrees.
I=fltarr(N)
for j=0,N-2 do begin $
  x=theta(j:*) & $
  y=Phi(j:*)*(abs(sin(theta(j:*)))^nu) & $
  I(j)=int_tabulated(x,y) & $
endfor
Z_fit=I/abs(sin(theta))^(nu+1.0)>0.1
!p.multi=[0,1,2,0]
plot_oo,theta*r2d,Z_fit/Z_fit(t15),xstyle=1, $
  xrange=[1,30],xtitle='elongation !7e!3    (degrees)',$
  ytitle='surface brightness Z(!7e!3)',charsize=1.3
oplot,theta*r2d,Z_clem/Z_clem(t15),psym=8,symsize=0.6
plot,theta*r2d,I,xrange=[0,30],xstyle=1, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3, $
  ytitle='Z(!7e!3)sin!u!7m!n(e!3)'
!p.multi=0
