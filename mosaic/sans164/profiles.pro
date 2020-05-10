;profiles.pro

;plot N,S,E,W profiles of mosaic

r2d=180.0/!pi
delta=10.0/r2d
cal=7.45e-14
pflag=1

restore,'mercator_sans164.dat'
sz=size(mosaic)
Nx=sz(1)
Ny=sz(2)
p_coords,coords,x_sun,y_sun,Nx,Ny
window,xs=Nx,ys=Ny

;get west profile
angle=abs(coords(*,*,1)-0.0/r2d)
rad=coords(*,*,0)
twopi=2.0*!pi
pix=where(((angle le delta) or (angle ge twopi-delta)) and $
  (good gt 0) and (rad ge 73))
tmplt_w=bytarr(Nx,Ny)+1
tmplt_w(pix)=0
;zap stars and planets
tmplt_w(534:540,350:356)=1
tmplt_w(607:613,354:360)=1
tmplt_w(634:640,354:360)=1
tmplt_w(641:647,350:356)=1
tmplt_w(670:676,405:410)=1
tmplt_w(755:761,364:370)=1
tmplt_ew=tmplt_w
tvscl,alog(mosaic*(1-tmplt_w)>1)
p_w=profile_idl(mosaic,x_sun,y_sun,r_w,err_w,tmplt_w)

;get north profile
angle=abs(coords(*,*,1)-90.0/r2d)
pix=where(((angle le delta) or (angle ge twopi-delta)) and $
  (good gt 0) and (rad ge 35))
tmplt_n=bytarr(Nx,Ny)+1
tmplt_n(pix)=0
tvscl,alog(mosaic*(1-tmplt_n)>1)
p_n=profile_idl(mosaic,x_sun,y_sun,r_n,err_n,tmplt_n)

;east profile
angle=abs(coords(*,*,1)-180.0/r2d)
pix=where(((angle le delta) or (angle ge twopi-delta)) and $
  (good gt 0) and (rad gt 61))
tmplt_e=bytarr(Nx,Ny)+1
tmplt_e(pix)=0
tmplt_ew(pix)=0
tvscl,alog(mosaic*(1-tmplt_e)>1)
p_e=profile_idl(mosaic,x_sun,y_sun,r_e,err_e,tmplt_e)

;south profile
angle=abs(coords(*,*,1)-270.0/r2d)
pix=where(((angle le delta) or (angle ge twopi-delta)) and $
  (good gt 0) and (rad ge 32))
tmplt_s=bytarr(Nx,Ny)+1
tmplt_s(pix)=0
;zap star
tmplt_s(389:395,156:162)=1
tvscl,alog(mosaic*(1-tmplt_s)>1)
p_s=profile_idl(mosaic,x_sun,y_sun,r_s,err_s,tmplt_s)
tvscl,alog(mosaic*(4-tmplt_w-tmplt_n-tmplt_e-tmplt_s)>1)

tvscl,alog(mosaic*(1-tmplt_ew)>1)

window,xs=550,ys=550
thck=2*pflag+1
clr=128
if (pflag eq 1) then setplot
plot_oo,r_e*ps_avg*r2d,p_e*cal,xrange=[2,35],thick=thck, $
  xstyle=1,yrange=[3.0e-13,3.0e-10], xthick=thck,ythick=thck, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3,ystyle=1, $
  ytitle='surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=10,xtickv=[2,3,4,5,6,7,8,9,10,20,30], $
  charthick=thck
xyouts,4.9,8.5e-11,'E',charthick=thck,charsize=1.3
oplot,r_w*ps_avg*r2d,p_w*cal,thick=thck,color=clr
xyouts,5.3,3.6e-11,'W',charthick=thck,color=clr,charsize=1.0
oplot,r_n*ps_avg*r2d,p_n*cal,thick=thck
xyouts,13.5,2.3e-12,'N',charthick=thck,charsize=1.3
oplot,r_s*ps_avg*r2d,p_s*cal,thick=thck,color=clr
xyouts,12.2,1.5e-12,'S',charthick=thck,color=clr,charsize=1.3

;overplot model E-W profile
asr=1.207e-8		;albedo*sigma_1*r_1, +/- 10%
Omega_sun=6.8e-5	;Sun solid angle, steradians
grid_factor=5
contrib   = [0.050, 0.450, 0.500]
nu_pop    = [2.000, 1.000, 1.45]
restore,'/net/lpi24/user2/hahn/clementine/zodiacal_light/fit/fit_nu_ast=1.0_nu_occ=2.0.dat'
sz=size(sin_epsilon)
Ngrid=sz(1)
Npop=n_elements(contrib)
model=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do begin $
  factor=asr*(Omega_sun/!pi)/cal/sin_epsilon^(nu_pop(pop)+1.0) &$
  model=model+contrib(pop)*factor*SB(*,*,pop) & $
endfor
p_coords,model_coords,Ngrid/2-1,Ngrid/2-1,Ngrid,Ngrid
angle=abs(model_coords(*,*,1)-180.0/r2d)
pix=where((angle le delta) or (angle ge twopi-delta))
tmplt_me=bytarr(Ngrid,Ngrid)+1
tmplt_me(pix)=0
;tvscl,enlarge(alog(model*(1-tmplt_me)>1),3)
p_me=profile_idl(model,Ngrid/2-1,Ngrid/2-1,r_me,err_me,tmplt_me)
r=r_me*grid_factor*ps_avg*r2d
;oplot,r,p_me*cal,linestyle=2
;fit power-law to e-w profile
e_min=2.5
e_max=30.0
p=where((r ge e_min) and (r le e_max))
e=r(p)
sin_e=sin(e/r2d)
ln_sin_e=alog(sin_e)
b=p_me(p)*cal
ln_b=alog(b)
coeff=svdfit(ln_sin_e,ln_b,2,yfit=ln_b_fit,sigma=coeff_err)
b_fit=exp(coeff(0))*sin_e^coeff(1)
;oplot,r(p),b,linestyle=2,thick=2
;oplot,r,b_fit
if (pflag eq 1) then output_plot,'profiles.ps'

;plot east-west profile
;generate a logarithmically increasing elongation-axis
Np=300
mn=alog(4.7/r2d/ps_avg)
mx=alog(30.0/r2d/ps_avg)
r_ew=exp((mx-mn)*findgen(Np)/(Np-1)+mn)
p_ew=profile_clem(mosaic,x_sun,y_sun,r_ew,1-tmplt_ew,err_ew)
plot_oo,r_ew*ps_avg*r2d,p_ew*cal,xrange=[4,35], $
  xstyle=1,yrange=[3.0e-13,1.0e-10], $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3,ystyle=1, $
  ytitle='ecliptic surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=8,xtickv=[4,5,6,7,8,9,10,20,30]

;plot errors
oplot,r_ew*ps_avg*r2d,(p_ew-err_ew)*cal
oplot,r_ew*ps_avg*r2d,(p_ew+err_ew)*cal

;fit power-law to e-w profile
e_min=5.0
e_max=27.0
p=where((r_ew*ps_avg*r2d ge e_min)and(r_ew*ps_avg*r2d le e_max))
e=r_ew(p)*ps_avg
sin_e=sin(e)
ln_sin_e=alog(sin_e)
b=p_ew(p)*cal
ln_b=alog(b)
coeff=svdfit(ln_sin_e,ln_b,2,yfit=ln_b_fit,sigma=coeff_err)
b_fit=exp(coeff(0))*sin_e^coeff(1)
oplot,e*r2d,b_fit,color=128,thick=5
oplot,r_ew*ps_avg*r2d,p_ew*cal,thick=2

;get errors
c0_err=0.06
plot_oo,r_ew*ps_avg*r2d,p_ew*cal,xrange=[4,35], $
  xstyle=1,yrange=[3.0e-13,1.0e-10],psym=3, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3,ystyle=1, $
  ytitle='ecliptic surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=8,xtickv=[4,5,6,7,8,9,10,20,30]
oplot,e*r2d,b_fit*(1.0+c0_err)
oplot,e*r2d,b_fit*(1.0-c0_err)
c1_err=0.04
plot_oo,r_ew*ps_avg*r2d,p_ew*cal,xrange=[4,35], $
  xstyle=1,yrange=[3.0e-13,1.0e-10],psym=3, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3,ystyle=1, $
  ytitle='ecliptic surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=8,xtickv=[4,5,6,7,8,9,10,20,30]
oplot,e*r2d,b_fit*sin_e^c1_err
oplot,e*r2d,b_fit/sin_e^c1_err
print,'b = ',exp(coeff(0)),' +/- ',exp(coeff(0))*c0_err, $
  ' / sin(e)^(1 + ',-coeff(1)-1.0,' +/- ',c1_err,' )'

;ratio of N-S profile
Np=150
r_ns_s=(12.2*findgen(Np)/(Np-1)+5.1)/ps_avg/r2d
p_n_s=spline(r_n,p_n,r_ns_s)
p_s_s=spline(r_s,p_s,r_ns_s)
;plot_oo,r_n,p_n,xstyle=1,ystyle=1
;oplot,r_ns_s,p_n_s,psym=1
;plot_oo,r_s,p_s,xstyle=1,ystyle=1
;oplot,r_ns_s,p_s_s,psym=1
;ratio of E-W profile
r_ew_s=(24.0*findgen(Np)/(Np-1)+6.0)/ps_avg/r2d
p_e_s=spline(r_e,p_e,r_ew_s)
p_w_s=spline(r_w,p_w,r_ew_s)
;plot_oo,r_e,p_e,xstyle=1,ystyle=1
;oplot,r_ew_s,p_e_s,psym=1
;plot_oo,r_w,p_w,xstyle=1,ystyle=1
;oplot,r_ew_s,p_w_s,psym=1
if (pflag eq 1) then setplot
plot,r_ew_s*ps_avg*r2d,p_e_s/p_w_s,yrange=[0.90,1.4], $
  ytitle='surface brightness ratio',charsize=1.3, $
  xtitle='elongation !7e!3    (degrees)',ystyle=1,xrange=[0,30],$
  xthick=thck,ythick=thck,charthick=thck,thick=thck
xyouts,10.0,1.035,'E/W',charthick=thck,charsize=1.5
oplot,r_ns_s*ps_avg*r2d,p_n_s/p_s_s,thick=thck,color=128
xyouts,4.0,1.11,'N/S',charthick=thck,charsize=1.5,color=128
if (pflag eq 1) then output_plot,'profile_ratios.ps'

