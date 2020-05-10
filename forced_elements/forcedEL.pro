;forcedEL.pro

;JC Liou provided these data which give the forced elements
;of small bodies orbiting between 0<r<5 AU. Angles are in degrees.

d=fltarr(5,969)
openr,1,'forcedEL.dat'
readf,1,d
close,1
r=reform(d(0,*))
e=reform(d(1,*))
w=reform(d(2,*))
i=reform(d(3,*))
O=reform(d(4,*))
window,xs=550,ys=550
plot,r,w,xrange=[0,3],yrange=[0,360],ystyle=1,charsize=1.3, $
  xtitle='r    (AU)',ytitle='longitude of periapse    (degrees)'
plot,r,O,xrange=[0,3],yrange=[0,360],ystyle=1,charsize=1.3, $
  xtitle='r    (AU)',ytitle='longitude of periapse    (degrees)'
plot_io,r,i,xrange=[0,3],yrange=[0.1,10],psym=10
