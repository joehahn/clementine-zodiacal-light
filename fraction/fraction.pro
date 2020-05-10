;fraction.pro

pflag=1

;mass budget
date0 = [ 1955, 1967, 1976 ]	;date
f0    = [ 0.90, 0.90, 1.00 ]	;cometary fractio

;IR
date1 = [ 1994, 1995 ]
f1    = [ 0.60, 0.70 ]

;IDPs
date2 = [ 1993, 1993 ]
f2    = [ 0.20, 0.35 ]

;ecliptic ZL
date3 = [ 1993, 2001 ]
f3    = [ 0.35,  0.45 ]

;all-sky ZL
date4 = [ 1995, 2001 ]
f4    = [ 0.70, 0.88 ]

if (pflag eq 0) then window,xs=550,ys=550
if (pflag eq 1) then setplot
thck=2*pflag+1
plot,date0,f0,xrange=[1950, 2012],xstyle=1,ystyle=1, $
  yrange=[0,1.1],xtitle='Publication year', $
  ytitle='Cometary fraction',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck,thick=thck
oplot,date0,f0,psym=8,symsize=1.5
xyouts,1970,0.9,'mass-budget',charsize=1.3,charthick=thck

oplot,date1,f1,thick=thck
dx=0.1
dy=0.0015
oplot,date1,f1,psym=6,symsize=1.5
oplot,date1-dx,f1,psym=6,symsize=1.5
oplot,date1+dx,f1,psym=6,symsize=1.5
oplot,date1,f1-dy,psym=6,symsize=1.5
oplot,date1,f1+dy,psym=6,symsize=1.5
xyouts,1990.5,0.65,'IR',charsize=1.3,charthick=thck

oplot,date2,f2,thick=thck
oplot,date2-dx,f2,psym=5,symsize=1.5
oplot,date2+dx,f2,psym=5,symsize=1.5
oplot,date2,f2-dy,psym=5,symsize=1.5
oplot,date2,f2+dy,psym=5,symsize=1.5
xyouts,1987,0.25,'IDPs',charsize=1.3,charthick=thck

oplot,date3,f3,linestyle=2,thick=thck
plots,date3(1),f3(1),psym=7,symsize=1.5
plots,date3(1)-dx,f3(1),psym=7,symsize=1.5
plots,date3(1)+dx,f3(1),psym=7,symsize=1.5
plots,date3(1),f3(1)-dy,psym=7,symsize=1.5
plots,date3(1),f3(1)+dy,psym=7,symsize=1.5
xyouts,2003,0.90,'this',charsize=1.3,charthick=thck
xyouts,2003,0.865,'work',charsize=1.3,charthick=thck

oplot,date4,f4,linestyle=2,thick=thck
plots,date4(1)-dx,f4(1),psym=7,symsize=1.5
plots,date4(1)+dx,f4(1),psym=7,symsize=1.5
plots,date4(1),f4(1)-dy,psym=7,symsize=1.5
plots,date4(1),f4(1)+dy,psym=7,symsize=1.5
plots,date4(1),f4(1),psym=7,symsize=1.5
xyouts,2003,0.46,'this',charsize=1.3,charthick=thck
xyouts,2003,0.425,'work',charsize=1.3,charthick=thck

if (pflag eq 1) then output_plot,'fraction.ps'

