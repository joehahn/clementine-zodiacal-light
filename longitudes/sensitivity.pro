;sensitivity.pro

;plot startracker sensitivity and determine effective
;bandpass

pflag=1

;wavelength & sensitivity for gain ID = 0,1,2
lambda = [ 4000, 5000, 6000, 7000, 8000, 9000, 10000 ]
s0     = [  0.0,  0.6,  1.6,  2.1,  2.3,  1.2,   0.0 ]
s1     = [  0.0,  1.6,  4.0,  5.0,  5.6,  3.4,   0.0 ]
s2     = [  0.0,  3.0,  8.0,  9.6, 11.6,  7.0,   0.0 ]

;normalize & average sensitivities
s0=s0/max(s0)
s1=s1/max(s1)
s2=s2/max(s2)
s=(s0+s1+s2)/3.0

if (pflag eq 0) then window,xs=501,ys=501
if (pflag eq 1) then setplot
thck=3
plot,lambda,s,psym=8,yrange=[0,1.1],charsize=1.3, $
  symsize=1,xtitle='wavelength    ('+string("305B)+')', $
  ystyle=1,ytitle='relative sensitivity', $
  xthick=thck,ythick=thck,charthick=thck
oplot,lambda,s,thick=3
;oplot,lambda,s0
;oplot,lambda,s1
;oplot,lambda,s2

;fit a spline to the data
lambda_spline=6000.0*findgen(1001)/1000.0+4000.0
s_spline=spline(lambda,s,lambda_spline)
;oplot,lambda_spline,s_spline,psym=1,symsize=0.3,color=128

;integrate sensitivity
I=int_tabulated(lambda_spline,s_spline)
idx=562
I_half=int_tabulated(lambda_spline(0:idx),s_spline(0:idx))
print,I,I_half/I

;mean wavelength
lambda_mean=lambda_spline(idx)
lambda_width=I
;plots,lambda_spline(idx)+[0,0],[0,1]

;effective filter
s_eff=s_spline*0.0
j=where(abs(lambda_spline-lambda_mean) le 0.5*lambda_width)
s_eff(j)=1.0
oplot,lambda_spline,s_eff,thick=5,color=150
I_eff=int_tabulated(lambda_spline(j),s_eff(j))
print,I,I_eff
plot,lambda,s,psym=8,yrange=[0,1.1],charsize=1.3, $
  symsize=1,xtitle='wavelength    ('+string("305B)+')',$
  ystyle=1,ytitle='relative sensitivity',noerase=1, $
  xthick=thck,ythick=thck,charthick=thck
oplot,lambda,s,thick=3
print,'effective bandpass (Angstroms) between ', $
  min(lambda_spline(j)),max(lambda_spline(j))
print,'mean wavelength (angstroms) = ',lambda_mean

if (pflag eq 1) then output_plot,'sensitivity.ps'
