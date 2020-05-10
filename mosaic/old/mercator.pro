;mercator.pro

;generates mosaic image of the ZL in the mercator projection,
;which has ecliptic east to the left and nort up.

;Set fix_image=1 to fill data-gaps
fix_image=0

;Set trim_image=1 to trim the edges of the mosaic
trim_image=0

;Set pflag=1 to output figures as postscript files
pflag=0

;calibration constant coverts 1 ADU/sec/pixel to mean solar brightness
cal=8.6e-14

;data filenames 
files=['orbit66.dat',  'orbit110.dat', 'orbit110_sunny.dat',  $
       'orbit164.dat', 'orbit193.dat', 'orbit206.dat',        $
       'orbit253.dat' ]

;additional dc offset added to the data (counts/sec)
dc_correction=[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]

;range of uncertainty in dc_correction
dc_error     =[ 0.0, 2.0, 5.0, 3.0, 0.5, 2.0, 1.0 ]

files=['orbit66.dat', 'orbit164.dat', 'orbit253.dat']
dc_correction=[ 0.0, 0.0, 0.0 ]

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
x=-phi/ps_avg
y=-theta/ps_avg
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

;fix image
if (fix_image eq 1) then begin $
  bad=where((angle lt 40.0/r2d) and (good eq 0)) & $
  Nbad=n_elements(bad) & $
  for i=0l,Nbad-1 do begin $
    pix=bad(i) & $
    x0=xx(pix) & $
    y0=yy(pix) & $
    x1=-x0 & $
    y1=y0 & $
    x2=-x0 & $
    y2=-y0 & $
    x3=x0 & $
    y3=-y0 & $
    x4=x_sun+[x0,x1,x2,x3] & $
    y4=y_sun+[y0,y1,y2,y3] & $
    gd=good(x4,y4) & $
    ;don't use rightmost images for repair
    j=where((gd gt 0) and (x4 le 580)) & $
    if (j(0) ne -1) then begin $
      mosaic(x0+x_sun,y0+y_sun)= $
        avg(mosaic(x4(j),y4(j))) & $
    endif & $
  endfor & $
  ;fix data-gap near Venus
  bad=where((good eq 0) and (xx le -330) and $
    (abs(yy) lt 5)) & $
  gd=where((good eq 1) and (xx le -330) and $
    (abs(yy) lt 5)) & $
  mosaic(bad)=avg(mosaic(gd)) & $
endif

;tv greyscale
if (pflag eq 1) then setplot
if (pflag ne 1) then window,xs=Nxs,ys=Nys
log_mosaic=alog(mosaic>5)
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

stop

;contour map
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
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=17,ysize=17.0,xoffset=3,yoffset=8,/color, $
    /portrait,bits_per_pixel=8 & $
endif
loadct,40
log_mosaic(sun)=192.0
tv,log_mosaic
loadct,0
;one_arrow, 750.0, 570.0, 90.0, 'N'
;one_arrow, 750.0, 570.0,180.0, 'E'
if (pflag eq 1) then output_plot,'color_mosaic.ps'

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

;check orbits 66 & 193
files=['orbit66.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0 ]
a=mosaic(*,170:190)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[450,650],yrange=[9,15]	;x=570
a=mosaic(410:430,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[50,200],yrange=[5,25]	;y=110

;check orbits 193 & 206
files=['orbit193.dat', 'orbit206.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(*,70:100)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,yrange=[6,15],xrange=[0,300]		;x=100
plot,a,yrange=[10,15],xrange=[500,700]		;x=630
a=mosaic(*,340:370)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,yrange=[10,25],xrange=[0,200]		;x=100
a=mosaic(200:230,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[0,150],yrange=[4,15]		;y=40

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

;check orbits 66, 110, & 193
files=['orbit66.dat', 'orbit110.dat', 'orbit193.dat' ]
dc_correction=[ 0.0, 0.0, 0.0 ]
a=mosaic(*,490:510)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[500,750]			;x=640
a=mosaic(*,380:400)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[550,650],yrange=[10,70]	;x=600
a=mosaic(*,140:160)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[500,700],yrange=[0,20]	;x=620

;check orbits 110 & 110_sunny
files=['orbit110.dat', 'orbit110_sunny.dat' ]
dc_correction=[ 0.0, 0.0 ]
a=mosaic(*,520:540)
for i=0,Nxs-1 do a(i,0)=median(a(i,*))
a=a(*,0)
plot,a,xrange=[200,500],yrange=[0,20]	;x=320,450
a=mosaic(350:370,*)
for i=0,Nys-1 do a(0,i)=median(a(*,i))
a=a(0,*)
plot,a,xrange=[500,700]			;y=590
plot,a,xrange=[0,300]			;y=210

;check orbits 193 & 253
files=['orbit193.dat', 'orbit253.dat' ]
dc_correction=[ 0.0, 0.0]
a=mosaic(70:100,*)
for i=0,Nys-1 do a(0,i)=avg(a(*,i))
a=a(0,*)
plot,a,xrange=[100,250]				;y=170
plot,a,xrange=[300,500]				;y=390
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
plot,a,xrange=[80,220]				;x=150
a=mosaic(*,270:300)
for i=0,Nxs-1 do a(i,0)=avg(a(i,*))
a=a(*,0)
plot,a,xrange=[60,240]				;x=190
