;lamy.pro

;Generate an approximate form for the forward-scattering
;phase function given in Lamy & Perrin, 1986, A&A, 163, 269.
;This is achieved by starting with the Hong (1985) phase-fn
;with nu=1, and then adding a 1/theta^3 term that is presumably
;due to diffration by large particles

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

;Phi1
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

;generate Hong's Fig. 5
window,xs=550,ys=550
plot,theta*r2d,alog10(Phi1),yrange=[-1.7,0.3],xstyle=1,ystyle=1,$
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
Phi1d45=(Phi1-(nu-1.0)*cos(theta)*Z)
Phi1d45=Phi1d45*Phi1(t90)/Phi1d45(t90)
oplot,theta*r2d,alog10(Phi1d45>1.0e-6),linestyle=1

;compare Hong's nu=1 phase function to the volume scattering
;function of Lamy & Perrin (1986), all normalized at
;theta=103 degrees
theta_lp = [ 5.0,         8.0,    17.0,    29.0,    67.0,   103.0,   148.0, 180.0   ]
Phi_lp1  = [ 1.0e-20, 4.0e-21, 1.0e-21, 4.0e-22, 1.0e-22, 7.0e-23, 1.0e-22, 1.5e-22 ]
theta_lp=theta_lp/r2d
;normalize Hong's nu=1 function to unity at theta=180 degrees
j_hong=where(theta ge 180.0/r2d)
j_hong=j_hong(0)
Phi1=Phi1/Phi1(j_hong)
Phi1d45=Phi1d45/Phi1d45(j_hong)
;adjust lp function to agree with Hong's at theta=103 degrees
j_lp=where(theta_lp ge 103.0/r2d)
j_lp=j_lp(0)
j_hong=where(theta ge 103.0/r2d)
j_hong=j_hong(0)
Phi_lp1=Phi_lp1*Phi1(j_hong)/Phi_lp1(j_lp)
phi_lp=Phi1+0.04/theta^3
plot_oo,theta_lp*r2d,Phi_lp1,psym=8,xstyle=1,ystyle=1, $
  yrange=[0.1,100],xrange=[3,180],charsize=1.3, $
  xtitle='scattering angle !9P!3    (degrees)', $
  ytitle='phase law !7w(!9P!3)'
oplot,theta*r2d,Phi_lp
oplot,theta*r2d,Phi1,linestyle=1
oplot,theta*r2d,Phi1d45,linestyle=2
