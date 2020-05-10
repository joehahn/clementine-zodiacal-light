;profiles.pro

nu_hong=1.45			;dust radial power law
nu_lamy=1.30			;dust radial power law
pflag=1				;plotting flag

;scattering angle
N=1001
scat_angle=!pi*findgen(N)/(N-1)
scat_angle(0)=0.5*scat_angle(1)

;get the phase law of Lamy & Perrin and Hong
phase_law_lamy=fltarr(N)
phase_law_hong=fltarr(N)
for j=0,N-1 do begin $
  phase_law_lamy(j)=phase_law(scat_angle(j),1.0     ) & $
  phase_law_hong(j)=phase_law(scat_angle(j),nu_hong ) & $
endfor
phase_law_lamy=phase_law_lamy+0.04/scat_angle^3
window,xs=550,ys=550
r2d=180.0/!pi
thck=2*pflag+1
if (pflag eq 1) then setplot
plot_io,scat_angle*r2d,phase_law_hong,yrange=[0.1,100],xstyle=1,$
  xtitle='scattering angle !9P!3    (degrees)', $
  ytitle='phase law !7w(!9P!3)',charsize=1.3, $
  thick=thck,xthick=thck,ythick=thck,charthick=thck
oplot,scat_angle*r2d,phase_law_lamy,thick=thck
xyouts,30,0.31,'Hong (1985)',charthick=thck,charsize=1
xyouts,17,10.0,'Lamy and Perrin (1986)',charthick=thck,charsize=1
if (pflag eq 1) then output_plot,'phase_function.ps'

stop

;compute ecliptic surface brightness profiles
N=101
elongation=0.5*!pi*findgen(N)/(N-1)
elongation(0)=0.5*elongation(1)
Z_lamy=fltarr(N)
Z_hong=fltarr(N)
lamy_intgrnd=phase_law_lamy*abs(sin(scat_angle))^nu_lamy
hong_intgrnd=phase_law_hong*abs(sin(scat_angle))^nu_hong
oplot,scat_angle*r2d,lamy_intgrnd,linestyle=1
oplot,scat_angle*r2d,hong_intgrnd,linestyle=2
for j=0,N-2 do begin $
  k=where(scat_angle ge elongation(j)) & $
  Z_lamy(j)=int_tabulated(scat_angle(k),lamy_intgrnd(k)) & $
  Z_hong(j)=int_tabulated(scat_angle(k),hong_intgrnd(k)) & $
endfor
Z_lamy=Z_lamy/sin(elongation)^(nu_lamy+1.0)
Z_hong=Z_hong/sin(elongation)^(nu_hong+1.0)
!p.multi=[0,1,2,0]
plot_oo,elongation*r2d,Z_lamy,xrange=[1,30],xstyle=1,ystyle=1, $
  yrange=[1,1.e5],ytitle='surface brightness Z(!7e!3)', $
  xtitle='elongation angle !7e!3    (degrees)',charsize=1.3
oplot,elongation*r2d,Z_hong,linestyle=2
j=where(Z_hong gt 0.0)
plot,elongation(j)*r2d,Z_lamy(j)/Z_hong(j),xrange=[0,30], $
  xtitle='elongation angle !7e!3    (degrees)',charsize=1.3, $
  ytitle='Z!lLamy!n/Z!lHong!n'
!p.multi=0




