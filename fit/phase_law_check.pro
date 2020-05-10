;phase_law_check.pro

;double-check the computation for the phase-law

;power-law index for the dust radial distribution
nu=1.45

N=1001
theta=!pi*findgen(N)/(N-1)
r2d=180.0/!pi

.run phase_law
Phi=fltarr(N)
for i=0,N-1 do Phi(i)=phase_law(theta(i), nu)

window,xs=550,ys=550
r2d=180.0/!pi
plot,theta*r2d,Phi,psym=10,xrange=[0,180],xstyle=1
oplot,theta*r2d,Phi*abs(sin(theta))^nu,linestyle=2
