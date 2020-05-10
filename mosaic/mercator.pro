;mercator.pro

;generates mosaic image of the ZL in the mercator projection,
;which has ecliptic east to the left and north up.

;The planets have the (x,y) coordinates:
;Saturn (537,353) at 10.6 deg, Mars (610,357) at 16.2 deg,
;Mars(638,357)  at 18.3 deg, Saturn (644,353) at 18.7 deg,
;Mercury (758,367) at 27.4 deg,
;and an absent Venus has x~138 at -19.5 deg

;Set trim_image=1 to trim the edges of the mosaic
trim_image=1

;Fix datagap near Venus when set
fix_Venus=1

;data filenames
files=[ 'orbit66.dat' , 'orbit110.dat', 'orbit110_sunny.dat',  $
        'orbit164.dat', 'orbit193.dat', 'orbit206.dat'      ,  $
        'orbit253.dat' ]

;additional dc offset added to the data (counts/sec)
dc_correction=[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]

;range of uncertainty in dc_correction
dc_error     =[ 0.0, 2.0, 5.0, 3.0, 0.5, 3.0, 1.0 ]

;Read data
Nimages=n_elements(files)
for i=0,Nimages-1 do begin $
  print,'Reading ',files(i) & $
  restore,files(i) & $
  if (i eq 0) then begin $
    sz=size(counts) & $
    Nx=sz(1) & $
    Ny=sz(2) & $
    countss=fltarr(Nx,Ny,Nimages) & $
    exp_time=fltarr(Nx,Ny,Nimages) & $
    good_pixels=bytarr(Nx,Ny,Nimages) & $
    phii=fltarr(Nx,Ny,Nimages) & $
    thetaa=fltarr(Nx,Ny,Nimages) & $
    plate_scales=fltarr(Nimages) & $
    twist_angles=fltarr(Nimages) & $
  endif & $
  countss(0,0,i)=counts & $
  exp_time(0,0,i)=total_exp_time & $
  good_pixels(0,0,i)=good_pix & $
  phii(0,0,i)=phi & $
  thetaa(0,0,i)=theta & $
  plate_scales(i)=plate_scl & $
  twist_angles(i)=twist_angle & $
endfor
phi=phii
theta=thetaa
counts=countss

;use average plate-scale
ps_avg=avg(plate_scales)
ps_sd=stdev(plate_scales)
twopi=2.0*!pi
r2d=360./twopi
print,'average plate-scale = ',ps_avg,' radians/pixel'
print,'average plate-scale = ',ps_avg*r2d,' degrees/pixel'
print,'fractional 1-sigma error = ', $
  ps_sd/ps_avg/sqrt(n_elements(plate_scales))

;check images
flux=fltarr(Nx,Ny,Nimages)
px=where(good_pixels gt 0)
flux(px)=counts(px)/exp_time(px)

;resample image into ecliptic coordinate system,
;and smooth over a 2x2 box
x=-phi/ps_avg
y=theta/ps_avg
x_sun=-min(x)
x=x+x_sun
Nxs=round(max(x))+2
y_sun=-min(y)
y=y+y_sun
Nys=round(max(y))+2
total_counts=fltarr(Nxs,Nys)
total_exp_time=fltarr(Nxs,Nys)
good=intarr(Nxs,Nys)
for xi=0,Nx-1 do begin $
  for yi=0,Ny-1 do begin $
    for ti=0,Nimages-1 do begin $
      if (good_pixels(xi,yi,ti) eq 1) then begin $
        xx=x(xi,yi,ti) & $
        yy=y(xi,yi,ti) & $
        total_counts(xx,yy)=total_counts(xx:xx+1,yy:yy+1)+ $
          counts(xi,yi,ti)+dc_correction(ti)* $
          exp_time(xi,yi,ti) & $
        total_exp_time(xx,yy)=total_exp_time(xx:xx+1,yy:yy+1)+$
          exp_time(xi,yi,ti) & $
        good(xx,yy)=good(xx:xx+1,yy:yy+1)+1 & $
      endif & $
    endfor & $
  endfor & $
endfor
j=where(good gt 0)
mosaic=fltarr(Nxs,Nys)
mosaic(j)=total_counts(j)/total_exp_time(j)

;trim image
xmax=30.0 & $
ymax=30.0 & $
x=r2d*ps_avg*(indgen(Nxs)-x_sun)
y=r2d*ps_avg*(indgen(Nys)-y_sun)
if (trim_image eq 1) then begin $
  j=where(abs(x) le xmax) & $
  x=x(j) & $
  x0=min(j) & $
  x1=max(j) & $
  j=where(abs(y) le ymax) & $
  y=y(j) & $
  y0=min(j) & $
  y1=max(j) & $
  mosaic=mosaic(x0:x1,y0:y1) & $
  good=good(x0:x1,y0:y1) & $
  x_sun=x_sun-x0 & $
  y_sun=y_sun-y0 & $
  sz=size(mosaic) & $
  Nxs=sz(1) & $
  Nys=sz(2) & $
endif

;Sun
Rsun=4.64e-3			;angular radius of Sun, radians
xx=fltarr(Nxs,Nys)
xx1=findgen(Nxs)-x_sun
for i=0,Nys-1 do xx(0,i)=xx1
yy=fltarr(Nxs,Nys)
yy1=transpose(findgen(Nys)-y_sun)
for i=0,Nxs-1 do yy(i,0)=yy1
angle=acos(cos(xx*ps_avg)*cos(yy*ps_avg))	;angular distance from Sun
sun=where(angle le Rsun)

;cosmetic fix along Venus' row
if (fix_Venus eq 1) then begin $
  pix=where((good le 0) and (xx ge -x_sun) and $
    (xx le 200-x_sun) and (yy ge 350-y_sun) and $
    (yy le 400-y_sun)) & $
  for j=0,n_elements(pix)-1 do begin $
    x_pix=xx(pix(j))+x_sun & $
    y_pix=yy(pix(j))+y_sun & $
    y_pix2=-yy(pix(j))+y_sun & $
    mosaic(x_pix,y_pix)=mosaic(x_pix,y_pix2) & $
  endfor & $
  pix=where((good le 0) and (xx ge -x_sun) and $
    (xx le 200-x_sun) and (yy ge 350-y_sun) and $
    (yy le 400-y_sun)) & $
  for j=0,n_elements(pix)-1 do begin $
    x_pix=xx(pix(j))+x_sun & $
    y_pix=yy(pix(j))+y_sun & $
    y_pix2=-yy(pix(j))+y_sun & $
    mosaic(x_pix,y_pix)=mosaic(x_pix,y_pix2+3) & $
  endfor & $
  radial,rad,134,373,Nxs,Nys & $
  pix=where((good le 0) and (xx ge -x_sun) and $
    (xx le 200-x_sun) and (yy ge 350-y_sun) and $
    (yy le 400-y_sun) and (mosaic eq 0.0)) & $
  gd2=where((rad le 6.0) and (mosaic gt 0.0)) & $
  mosaic(pix)=avg(mosaic(gd2)) & $
endif

;total flux
gd=where(good ge 1)
Ngd=n_elements(gd)
total_flux=total(mosaic(gd))
print,'total flux (counts/sec) = ',total_flux & $
print,'Number of good pixels   = ',Ngd

;tv greyscale
m_min=5.0
m_max=max(mosaic)
lm0=alog(mosaic>m_min)
lm0_min=min(lm0)
lm1=lm0-lm0_min
lm1_max=max(lm1)
log_mosaic=254.0*lm1/lm1_max
log_mosaic(sun)=192.0
loadct,0
window,xs=Nxs,ys=Nys
tv,log_mosaic

;add color-coded intensity scale
log_mosaic2=fltarr(Nxs,Nys+20)
log_mosaic2(0,20)=log_mosaic
cal=7.5e-14
bx=254*findgen(Nxs)/(Nxs-1)
flux=cal*exp(bx*lm1_max/255.0+lm0_min)
flux_ticks=[5.0e-13, 1.6e-12, 5.0e-12, 1.6e-11, 5.0e-11, 1.6e-10]
flux_label=fix(flux_ticks/1.0e-13)
box=fltarr(Nxs,20)
for j=0,19 do box(0,j)=bx
log_mosaic2(0,0)=box
window,xs=Nxs,ys=Nys+20
loadct,39
tv,log_mosaic2
for j=0,n_elements(flux_ticks)-1 do begin $
  k=where(flux ge flux_ticks(j)) & $
  k=k(0) & $
  print,k,bx(k),flux(k) & $
  log_mosaic2(k-1:k+1,0:5)=0 & $
  ;check scale
  ;k=where(mosaic le flux(k)/cal) & $
  ;log_mosaic2(k+20*Nxs)=0 & $
  ;tv,log_mosaic2 & $
endfor
bx=bx*cal
tv,log_mosaic2

;contour map
loadct,0
log_mosaic(sun)=0.0
xrng=xmax*[-1,1]
yrng=xmax*[-1,1]
contour,alog(smooth(mosaic,4)>5),x,y,charsize=1.3, $
  xstyle=1,ystyle=1,xrange=xrng,yrange=yrng, $
  xtitle='ecliptic longitude    (degrees)',$
  ytitle='ecliptic latitude    (degrees)', $
  levels=findgen(25)*0.3
plots,[0,0],yrng,color=128
plots,xrng,[0,0],color=128
angl=twopi*findgen(101)/100.0
x_s=Rsun*r2d*cos(angl)
y_s=Rsun*r2d*sin(angl)
polyfill,x_s,y_s
one_arrow, 170.0, 800.0, 90.0, 'N'
one_arrow, 170.0, 800.0,180.0, 'E'

;tv colormap
loadct,39
tv,log_mosaic2
one_arrow, 80.0, 685.0, 90.0, 'N'
one_arrow, 80.0, 685.0,180.0, 'E'
set_plot,'ps'
device,xsize=17,ysize=18.0,xoffset=2,yoffset=7,/color, $
  /portrait,bits_per_pixel=8
tv,log_mosaic2
sz=30.0
one_arrow,1900.0,14900.0,90.0,'N',color=255,THICK=5.0, $
  CHARSIZE=1.5,ARROWSIZE=[sz*30.0,sz*9.0,35.0]
one_arrow,1900.0,14900.0,180.0,'E',color=255,THICK=5.0, $
  CHARSIZE=1.5,ARROWSIZE=[sz*30.0,sz*9.0,35.0]
xx=[-0.061, 0.12, 0.294, 0.477, 0.652, 0.833]
for j=0,n_elements(flux_ticks)-1 do begin $
  xyouts,xx(j),-22.0/Nys,string(flux_label(j)),normal=1, $
    color=0,charsize=1.3,charthick=3 & $
endfor
output_plot,'mosaic.ps'

;save as a 2xtiff file
loadct,40
log_mosaic3=log_mosaic2
radial,rad,x_sun,y_sun+20,Nxs,Nys+20
good3=fltarr(Nxs,Nys+20)
good3(0,20)=good
j=where((good3 eq 0) and (rad le 100.0))
log_mosaic3(j)=255.0
Nxs2=2*Nxs
Nys2=2*(Nys+20)
log_mosaic3=rebin(log_mosaic3,Nxs2,Nys2)
good3=rebin(good3<1,Nxs2,Nys2)
radial,rad,2*x_sun,2*(y_sun+20),Nxs2,Nys2
j=where((good3 eq 0) and (rad le 200.0))
log_mosaic3(j)=0.0
j=where(rad le 2.0*Rsun/ps_avg)
log_mosaic3(j)=192.0
log_mosaic3=transpose(rotate(log_mosaic3,1))
common colors,r1,g1,b1
write_tiff,'mosaic.tif',log_mosaic3,red=r1,green=g1,blue=b1
loadct,0

;save as a 1xtiff file
loadct,40
log_mosaic3=log_mosaic2
radial,rad,x_sun,y_sun+20,Nxs,Nys+20
good3=fltarr(Nxs,Nys+20)
good3(0,20)=good
log_mosaic3=transpose(rotate(log_mosaic3,1))
common colors,r1,g1,b1
write_tiff,'mosaic.tif',log_mosaic3,red=r1,green=g1,blue=b1
loadct,0

;save results
save,filename='mercator.dat',mosaic,good,ps_avg,angle, $
  x_sun,y_sun

stop

;ecliptic surface brightness profiles
loadct,0
prof=fltarr(Nxs)
for p=0,Nxs-1 do begin $
  im_sub=mosaic(p,y_sun-1:y_sun+1) & $
  gd_sub=good(p,y_sun-1:y_sun+1) & $
  j=where(gd_sub ne 0) & $
  if (j(0) ne -1) then prof(p)=median(im_sub(j)) & $
  if (j(0) eq -1) then prof(p)=-10.0 & $
endfor
prof_east=rotate(prof(0:x_sun),2)
prof_west=prof(x_sun:*)
;zap innermost value
j=where(prof_east gt -10.0)
prof_east(j(0))=prof_east(0)
j=where(prof_west gt -10.0)
prof_west(j(0))=prof_west(0)
;plot east profile
N=n_elements(prof_east)
elongation=ps_avg*findgen(N)
if (pflag eq 1) then setplot
plot,elongation*r2d,prof_east*cal,xrange=[3.5,35],xlog=1, $
  ylog=1,xstyle=1,yrange=[5.0e-13,3.0e-10],ystyle=1,thick=6, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3, $
  ytitle='ecliptic surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=8,xtickv=[4,5,6,7,8,9,10,20,30],xthick=5,ythick=5, $
  charthick=5
xyouts,20.0,1.4e-11,'Venus',charsize=1.3,charthick=4
;oplot,elongation*r2d,prof_west*cal,color=128
;resample profile on logarithmic interval
log_prof_east=alog(prof_east>1.0)
log_sin_e=alog(sin(elongation>elongation(1)))
xmin=-2.45
xmax=-1.3
Nx=101
x1=(xmax-xmin)*findgen(Nx)/(Nx-1)+xmin
p1=where((log_sin_e ge xmin) and (log_sin_e le xmax))
xmin=-1.0
xmax=-0.9
Nx=8
x2=(xmax-xmin)*findgen(Nx)/(Nx-1)+xmin
p2=where((log_sin_e ge xmin) and (log_sin_e le xmax))
x0=[x1,x2]
p0=[p1,p2]
y0=spline(log_sin_e(p0),log_prof_east(p0),x0)
coeff=poly_fit(x0,y0,1,yfit)
print,'east-west power law index = ',coeff(1) & $
print,'alpha = ',-coeff(1)-1.0 & $
print,'brightness factor = ',exp(coeff(0))*cal
fit=exp(coeff(0))*cal*sin(elongation>elongation(1))^coeff(1)
oplot,elongation*r2d,fit,color=128,thick=4
;oplot,elongation*r2d,fit*1.125*sin(elongation>elongation(1))^0.05,color=200
oplot,elongation*r2d,prof_east*cal
if (pflag eq 1) then output_plot,'profiles.ps'
;check fit
;plot,log_sin_e,log_prof_east,xrange=[-3,-0.5],yrange=[2,8], $
;  ystyle=1
;oplot,log_sin_e(p0),log_prof_east(p0),psym=1
;oplot,log_sin_e,coeff(0)+coeff(1)*log_sin_e,psym=3

;save results
save,filename='mercator.dat',mosaic,angle,cal,ps_avg,x,y,coeff, $
  x_sun,y_sun

stop

;Check dc_offsets

contour,alog(smooth(mosaic,4)>1), $
  levels=findgen(30)*0.25,xstyle=1,ystyle=1, $
  xrange=[000,900],yrange=[000,900]

;check image at Venus
a=mosaic(140:160,*)
for y=0,Nys-1 do a(0,y)=avg(a(*,y))
a=a(0,*)
plot,a,xrange=[300,450],yrange=[20,45]	;y=365

;check orbits 66 & 193
files=['orbit66.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0 ]
a=mosaic(*,170:190)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[450,650],yrange=[9,17]	;x=570
a=mosaic(410:430,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[50,200],yrange=[5,25]	;y=110

;check orbits 193 & 206
files=['orbit193.dat', 'orbit206.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(*,190:210)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,yrange=[10,20],xrange=[0,200]		;x=100
plot,a,xrange=[560,650]				;x=620
a=mosaic(*,340:370)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,yrange=[10,25],xrange=[0,200]		;x=100
a=mosaic(350:360,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[50,150],yrange=[5,20]		;y=100

;check orbits 164 & 193
files=['orbit164.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(190:210,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[200,300]			;y=260
a=mosaic(*,360:375)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[50,250],yrange=[15,75],ystyle=1	;x=150

;check orbits 193 & 253
files=['orbit193.dat', 'orbit253.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(80:100,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[100,250]				;y=170
plot,a,xrange=[300,500]				;y=390
a=mosaic(115:135,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[200,350],yrange=[20,50]		;y=260
a=mosaic(190:210,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[200,300]				;y=260

;check orbits 164 & 253
files=['orbit164.dat', 'orbit253.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(*,440:470)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[80,220]				;x=145
a=mosaic(*,270:300)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[50,220]				;x=170

;check orbits 66, 110, 193, & 253
files=['orbit66.dat','orbit110.dat','orbit193.dat' ,'orbit253.dat' ]
dc_correction=[ 0.0, 0.0, 0.0 , 0.0 ]
a=mosaic(*,490:510)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[500,750]			;x=650
a=mosaic(*,380:400)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[570,670],yrange=[15,60]	;x=640
a=mosaic(*,140:160)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[500,700],yrange=[0,20]	;x=620

;check orbits 66, 110, 110_sunny, 193, & 253
files=['orbit66.dat','orbit110.dat','orbit110_sunny.dat', $
  'orbit193.dat' ,'orbit253.dat' ]
dc_correction=[ 0.0, 0.0, 0.0, 0.0, 0.0 ]
a=mosaic(*,580:620)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[600,700],yrange=[0,20]	;x=640
plot,a,xrange=[720,825],yrange=[0,15]	;x=780
a=mosaic(670:700,*)
for i=0,Nys-1 do a(0,i)=median(a(*,i))
a=a(0,*)
plot,a,xrange=[200,400]			;y=270
