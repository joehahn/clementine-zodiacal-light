;asteroid_stats.pro

pflag=0

;read Ted Bowell's list of asteroid orbit elements
read_astorb,number,name,H,diameter,a,e,i,O,w,M,data=data

;asteroid cumulative H distribution
p=where((a lt 6.0) and (number gt 1))
Hp=H(p)
Nbins=100
Haxis=30.0*findgen(Nbins)/(Nbins-1)
NH=fltarr(Nbins)
for ii=0,Nbins-1 do begin $
  j=where(Hp le Haxis(ii)) & $
  if (j(0) ne -1) then NH(ii)=n_elements(j) & $
endfor
window,xs=550,ys=550
plot_io,Haxis,NH,xrange=[0,20],yrange=[1,1.e6],charsize=1.3, $
  xtitle='absolute magnitude H',ytitle='N(mag<H)'

;plot (a,e,i) for objects brighter than H=15
p=where((a lt 6.0) and (H lt 15.0))
pp=where((a lt 6.0) and (H lt 15.0) and (i gt 20.0))
!p.multi=[0,1,2,0]
plot,a(p),e(p),psym=3,xrange=[0,6], $
  xtitle='a    (AU)',ytitle='e',yrange=[0,0.6]
oplot,a(pp),e(pp),psym=3,color=128
plot,a(p),i(p),psym=3,xrange=[0,6],xtitle='a    (AU)', $
  ytitle='i   (deg)',yrange=[0,50]
oplot,a(pp),i(pp),psym=3,color=128
oplot,[0,6],[20,20],color=128
!p.multi=0

;plot histogram of asteroid inclinations
histo,i(p),100,i_axis,Ni
thck=2*pflag+1
if (pflag eq 1) then setplot
plot,i_axis,Ni,psym=10,xrange=[0,40],yrange=[0,4500],ystyle=1, $
  charsize=1.3,xtitle='inclination i    (degrees)',thick=thck, $
  ytitle='asteroid inclination distribution N(i)', $
  xthick=thck,ythick=thck,charthick=thck
sigma=1.05*stdev(i(p))
r2d=180.0/!pi
coeff=0.0105*n_elements(p)/(sigma/r2d)^2
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'asteroid Sigma = ',sigma
pp=where(i(p) gt 20.0)
print,float(n_elements(pp))/float(n_elements(p))
;text for seminar figure
;xyouts,15,3900,'H<15 Asteroids',charsize=2.5,charthick=thck+2
;xyouts,23,3450,'!7r!3=6!9%!3',charsize=2.5,charthick=thck+2
if (pflag eq 1) then output_plot,'asteroid_i_dist.ps'
