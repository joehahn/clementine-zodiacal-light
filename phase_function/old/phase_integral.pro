;phase_integral.pro
;plots the phase integral

d=fltarr(3,45,3)
openr,1,'phase_integral.dat'
readf,1,d
close,1
alpha=reform(d(0,*,*))
epsilon=reform(d(1,*,*))
II=reform(d(2,*,*))
r2d=180./!pi

window,xs=501,ys=501
plot,epsilon*r2d,II,psym=1,symsize=1,charsize=1.3, $
  xtitle='elongation angle !7e!3    (degrees)', $
  ytitle='phase integral',xrange=[-5,125],xstyle=1
for i=0,2 do begin $
  oplot,epsilon(*,i)*r2d,II(*,i),thick=2*i & $
  print,i,alpha(0,i),II(0,i) & $
endfor

plot,epsilon(*,1)*r2d,II(*,1),xrange=[0,30]
