;mosaic.pro

;Set pflag=1 to output postscript files
pflag=0

;calibration constant coverts 1 ADU/sec/pixel to mean solar brightness
cal=1.007e-13

;Orbit 247 gives ugly results...don't use it.
files=['orbit66.dat',  'orbit110.dat', 'orbit110_sunny.dat',  $
       'orbit164.dat', 'orbit193.dat', 'orbit206.dat' ]
dc_correction=[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
;range of uncertainty in dc_correction
dc_error     =[ 0.0, 2.0, 4.0, 3.0, 0.0, 2.0 ]

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

;minus sign!
theta=-theta

;smooth orbit 110, if needed
;im=1
;gd=good_pixels(*,*,im)
;gd(182:*,*)=0
;good_pixels(0,0,im)=gd
;cnts=counts(*,*,im)*gd
;cnts=cnts(5:181,9:566)
;cnts=smooth(cnts,5)
;counts(5,9,im)=cnts

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
x=-sin(phi)*cos(theta)/ps_avg
y=-sin(theta)/ps_avg
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
        total_exp_time(xx,yy)=total_exp_time(xx:xx+1,yy:yy+1)+ $
          exp_time(xi,yi,ti) & $
        good(xx,yy)=good(xx:xx+1,yy:yy+1)+1 & $
      endif & $
    endfor & $
  endfor & $
endfor
j=where(good gt 0)
mosaic=fltarr(Nxs,Nys)
mosaic(j)=total_counts(j)/total_exp_time(j)

;cosmetic fix along Venus' row
save,filename='/tmp/mosaic.dat',mosaic
;for x=29,300 do begin $
;  y=where(good(x,365:375) eq 0) & $
;  y0=365+y(0) & $
;  y1=365+y(n_elements(y)-1) & $
;  if (y0 gt 364) then begin $
;    image(x,y0-1)=0.5*(image(x,y0-4:y0-2)+image(x,y0+2:y0+4)) & $
;    good(x,y0:y1)=1 & $
;  endif & $
;endfor

;trim image
x0=0
x1=770
y0=50
y1=730
mosaic=mosaic(x0:x1,y0:y1)
good=good(x0:x1,y0:y1)
x_sun=x_sun-x0
y_sun=y_sun-y0
sz=size(mosaic)
Nxs=sz(1)
Nys=sz(2)

;Sun
Rsun=4.64e-3			;angular radius of Sun, radians
radial,rad,x_sun,y_sun,Nxs,Nys
rad=rad*ps_avg
sun=where(rad le Rsun)

;tv greyscale
if (pflag eq 1) then setplot
if (pflag ne 1) then window,xs=Nxs,ys=Nys
log_mosaic=alog(mosaic>9)
log_mosaic=log_mosaic-min(log_mosaic)
log_mosaic=255.0*log_mosaic/max(log_mosaic)
log_mosaic(sun)=255.0
tv,log_mosaic
;one_arrow, 650.0, 650.0, 90.0, 'N'
;one_arrow, 650.0, 650.0,180.0, 'E'
;xyouts, -24.9, -0.25, '!9f!3',charsize=2
;xyouts,  23.2, -0.29, '!20T!3',charsize=2
;xyouts,  32.3, -0.28, '!9m!3',charsize=2
if (pflag eq 1) then output_plot,'grey_mosaic.ps'

;contour map
x=r2d*asin(ps_avg*(indgen(Nxs)-x_sun))
y=r2d*asin(ps_avg*(indgen(Nys)-y_sun))
if (pflag eq 1) then setplot
xrng=30.0*[-1,1] & yrng=30.0*[-1,1]
log_mosaic(sun)=0.0
contour,log_mosaic,x,y,nlevels=13,xstyle=1, $
  ystyle=1,charsize=1.3,xrange=xrng,yrange=yrng, $
  xtitle='ecliptic longitude    (degrees)', $
  ytitle='ecliptic latitude    (degrees)'
plots,[0,0],yrng,color=128
plots,xrng,[0,0],color=128
angle=twopi*findgen(101)/100.0
x_s=Rsun*r2d*cos(angle)
y_s=Rsun*r2d*sin(angle)
polyfill,x_s,y_s
one_arrow, 650.0, 650.0, 90.0, 'N'
one_arrow, 650.0, 650.0,180.0, 'E'
if (pflag eq 1) then output_plot,'contour.ps'

;tv colormap
if (pflag eq 1) then set_colorplot
loadct,40
log_mosaic(sun)=192.0
tv,log_mosaic
;one_arrow, 650.0, 650.0, 90.0, 'N'
;one_arrow, 650.0, 650.0,180.0, 'E'
;xyouts, -25.9, -0.25, '!9f!3 ',charsize=2
;xyouts,  22.2, -0.29, '!20T!3 ',charsize=2
;xyouts,  31.3, -0.28, '!9m!3 ',charsize=2
if (pflag eq 1) then output_plot,'color_mosaic.ps'
loadct,0

;ecliptic surface brightness profiles
if (pflag eq 1) then setplot
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
elongation=asin(ps_avg*findgen(N))
plot,elongation*r2d,prof_east*cal,xrange=[3.5,35],xlog=1, $
  ylog=1,xstyle=1,yrange=[5.0e-13,3.0e-10],ystyle=1,thick=3, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3, $
  ytitle='ecliptic surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=8,xtickv=[4,5,6,7,8,9,10,20,30]
;oplot,elongation*r2d,prof_west*cal,color=128
;resample profile on logarithmic interval
log_prof_east=alog(prof_east>1.0)
log_sin_e=alog(sin(elongation>elongation(1)))
xmin=-2.55
xmax=-1.3
Nx=101
x1=(xmax-xmin)*findgen(Nx)/(Nx-1)+xmin
p1=where((log_sin_e ge xmin) and (log_sin_e le xmax))
xmin=-1.0
xmax=-0.9
Nx=8
x2=(xmax-xmin)*findgen(Nx)/(Nx-1)+xmin
p2=where((log_sin_e ge xmin) and (log_sin_e le xmax))
x=[x1,x2]
p=[p1,p2]
y=spline(log_sin_e(p),log_prof_east(p),x)
coeff=poly_fit(x,y,1,yfit)
print,'east-west power law index = ',coeff(1) & $
print,'alpha = ',-coeff(1)-1.0 & $
print,'brightness factor = ',exp(coeff(0))*cal
fit=exp(coeff(0))*cal*sin(elongation>elongation(1))^coeff(1)
oplot,elongation*r2d,fit,color=128
;oplot,elongation*r2d, $
;  fit*1.1*sin(elongation>elongation(1))^0.05,color=200
if (pflag eq 1) then output_plot,'profiles.ps'
;check fit
;plot,log_sin_e,log_prof_east,xrange=[-3,-0.5],yrange=[2,8], $
;  ystyle=1
;oplot,log_sin_e(p),log_prof_east(p),psym=1
;oplot,log_sin_e,coeff(0)+coeff(1)*log_sin_e,psym=3

;save results
save,filename='mosaic.dat',mosaic,cal,coeff

stop

;Check dc_offsets

;check orbits 66 & 193
files=['orbit66.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(*,140:160)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[450,600],yrange=[9,15]	;x=530
a=mosaic(410:430,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[0,200],yrange=[5,20]	;y=100

;check orbits 193 & 206
files=['orbit193.dat', 'orbit206.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(*,70:100)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,yrange=[6,15],xrange=[0,300]		;x=70
plot,a,yrange=[10,15],xrange=[500,650]		;x=570
a=mosaic(*,50:60)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[0,300],yrange=[6,14]		;x=150
a=mosaic(*,340:370)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,yrange=[10,25],xrange=[20,120]		;x=65
a=mosaic(200:230,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[0,100],yrange=[5,15]		;y=30

;check orbits 164 & 193
files=['orbit164.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(140:160,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[350,450]			;y=380
a=mosaic(*,350:370)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[0,250]			;x=125
a=mosaic(190:210,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[200,300]			;y=250

;check orbits 66, 110, & 193
files=['orbit66.dat', 'orbit110.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0, 0.0 ]
a=mosaic(*,490:510)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[500,750]			;x=600
a=mosaic(*,380:400)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[550,650]			;x=570
a=mosaic(*,140:160)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[450,650]			;x=520

;check orbits 110 & 110_sunny
files=['orbit110.dat', 'orbit110_sunny.dat' ]
dc_correction=[ 0.0, 0.0 ]
a=mosaic(*,520:540)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[200,500],yrange=[0,20]	;x=300,430
a=mosaic(350:370,*)
for i=0,Nys-1 do a(0,i)=median(a(*,i))
a=a(0,*)
plot,a,xrange=[500,700]			;y=570
plot,a,xrange=[0,300]			;y~200
