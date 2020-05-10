;phase_funtion.pro
;Compare's Zook's phase function to Hong's.

N=1001
g = [0.700, -0.200, -0.810]
w = [0.665,  0.330,  0.005]

;Hong's phase function
theta=!pi*findgen(N)/(N-1)
phase_hong=fltarr(N)
for i=0,2 do $
  phase_hong=phase_hong + w(i)*(1.0-g(i)^2)/(4.0*!pi)/ $
    (1.0+g(i)^2-2.0*g(i)*cos(theta))^1.5

;Zook's phase function
a=0.195
K=0.4
phase_zook=(a/!pi)*exp(-K*(!pi-theta))

window,xs=501,ys=501
r2d=180.0/!pi
plot,theta*r2d,phase_hong,ylog=1,xstyle=1,charsize=1.3, $
  xtitle='scattering angle !7h!3    (degrees)',ystyle=1, $
  ytitle='phase function !7U(h)!3',yrange=[0.01,2.0]
oplot,theta*r2d,phase_zook,linestyle=2

