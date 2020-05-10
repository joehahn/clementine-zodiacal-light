;mercator_sans164.pro

;generates mosaic image of the ZL in the mercator projection,
;which has ecliptic east to the left and north up.
;Note: orbit 164 is not used!

;Set fix_image=1 to fill data-gaps
fix_image=0

;Set trim_image=1 to trim the edges of the mosaic
trim_image=1

;data filenames
files=[ '../orbit66.dat' , '../orbit110.dat', $
        '../orbit110_sunny.dat', '../orbit193.dat', $
        '../orbit206.dat', '../orbit253.dat' ]

;additional dc offset added to the data (counts/sec)
dc_correction=[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]

;range of uncertainty in dc_correction
dc_error     =[ 0.0, 2.0, 5.0, 0.5, 3.0, 1.0 ]

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

;total flux
gd=where(good ge 1)
Ngd=n_elements(gd)
total_flux=total(mosaic(gd))
print,'total flux (counts/sec) = ',total_flux & $
print,'Number of good pixels   = ',Ngd

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
log_mosaic=alog(mosaic>5)
log_mosaic=log_mosaic-min(log_mosaic)
log_mosaic=255.0*log_mosaic/max(log_mosaic)
log_mosaic(sun)=255.0
loadct,0
window,xs=Nxs,ys=Nys
tv,log_mosaic

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
loadct,39
log_mosaic(sun)=192.0
tv,log_mosaic
one_arrow, 80.0, 685.0, 90.0, 'N'
one_arrow, 80.0, 685.0,180.0, 'E'
tv,log_mosaic
sz=30.0
one_arrow,1900.0,14500.0,90.0,'N',color=255,THICK=5.0, $
  CHARSIZE=1.5,ARROWSIZE=[sz*30.0,sz*9.0,35.0]
one_arrow,1900.0,14500.0,180.0,'E',color=255,THICK=5.0, $
  CHARSIZE=1.5,ARROWSIZE=[sz*30.0,sz*9.0,35.0]
loadct,0

;save results
save,filename='mercator_sans164.dat',mosaic,good,ps_avg,angle, $
  x_sun,y_sun

