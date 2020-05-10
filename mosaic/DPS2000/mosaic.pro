;mosaic.pro

;Set pflag=1 to output postscript files
pflag=0

;calibration constant coverts 1 ADU/sec to mean solar brightness
cal=1.07e-13

;Filenames for each orbit
;Orbit 247 is ugly...don't use it.
files=['orbit66.dat',  'orbit110.dat', 'orbit164.dat', $
       'orbit193.dat', 'orbit206.dat' ]

;Read data
Nimages=n_elements(files)
for i=0,Nimages-1 do begin $
  print,'Reading ',files(i) & $
  restore,files(i) & $
  if (i eq 0) then begin $
    sz=size(avg_image) & $
    Nx=sz(1) & $
    Ny=sz(2) & $
    images=fltarr(Nx,Ny,Nimages) & $
    good_pixels=bytarr(Nx,Ny,Nimages) & $
    phii=fltarr(Nx,Ny,Nimages) & $
    thetaa=fltarr(Nx,Ny,Nimages) & $
    plate_scales=fltarr(Nimages) & $
    twist_angles=fltarr(Nimages) & $
  endif & $
  images(0,0,i)=avg_image & $
  good_pixels(0,0,i)=good_pix & $
  phii(0,0,i)=phi & $
  thetaa(0,0,i)=theta & $
  plate_scales(i)=plate_scl & $
  twist_angles(i)=twist_angle & $
endfor
phi=phii
theta=thetaa

;zap left, top & bottem rows of orbit 110
good_pixels(Nx-50:*,*,1)=0
good_pixels(*,0:120,1)=0
good_pixels(*,Ny-110:*,1)=0

;resample image into ecliptic coordinate system,
;and smooth over a 3x3 box
q=atan(tan(theta)/cos(phi))
eta=tan(theta)/cos(phi)
chi=cos(q)*tan(phi)/cos(q)
ps_avg=avg(plate_scales)
x=-round(chi/ps_avg)
xc=-min(x)+1
x=x+xc
Nxs=max(x)+2
y=round(eta/ps_avg)
yc=-min(y)+1
y=y+yc
Nys=max(y)+2
image=fltarr(Nxs,Nys)
good=intarr(Nxs,Nys)
for xi=0,Nx-1 do begin $
  for yi=0,Ny-1 do begin $
    for ti=0,Nimages-1 do begin $
      if (good_pixels(xi,yi,ti) eq 1) then begin $
        xx=x(xi,yi,ti) & $
        yy=y(xi,yi,ti) & $
        image(xx-1,yy-1)=image(xx-1:xx+1,yy-1:yy+1) $
          +images(xi,yi,ti) & $
        good(xx-1,yy-1)=good(xx-1:xx+1,yy-1:yy+1)+1 & $
      endif & $
    endfor & $
  endfor & $
endfor
j=where(good gt 0)
image(j)=image(j)/good(j)
good=good<1

;cosmetic fix along Venus' row
image_old=image
good_old=good
for x=29,300 do begin $
  y=where(good(x,365:375) eq 0) & $
  y0=365+y(0) & $
  y1=365+y(n_elements(y)-1) & $
  if (y0 gt 364) then begin $
    image(x,y0-1)=0.5*(image(x,y0-4:y0-2)+image(x,y0+2:y0+4)) & $
    good(x,y0:y1)=1 & $
  endif & $
endfor
Rsun=4.64e-3			;angular radius of Sun, radians
radial,rad,xc,yc,Nxs,Nys
rad=rad*ps_avg
sun=where(rad le Rsun)

if (pflag ne 1) then window,xs=Nxs,ys=Nys
print,'average plate-scale = ',ps_avg*180./!pi,' degrees/pixel'

;tv greyscale
if (pflag eq 1) then setplot
log_image=alog(image>1)
log_image=255.0*log_image/max(log_image)
log_image(sun)=255.0
tv,log_image
;one_arrow, 650.0, 650.0, 90.0, 'N'
;one_arrow, 650.0, 650.0,180.0, 'E'
;xyouts, -24.9, -0.25, '!9f!3',charsize=2
;xyouts,  23.2, -0.29, '!20T!3',charsize=2
;xyouts,  32.3, -0.28, '!9m!3',charsize=2
if (pflag eq 1) then output_plot,'grey_mosaic.ps'

;tv colormap
if (pflag eq 1) then set_colorplot
loadct,39
log_image(sun)=192.0
tv,log_image
;one_arrow, 650.0, 650.0, 90.0, 'N'
;one_arrow, 650.0, 650.0,180.0, 'E'
;xyouts, -25.9, -0.25, '!9f!3 ',charsize=2
;xyouts,  22.2, -0.29, '!20T!3 ',charsize=2
;xyouts,  31.3, -0.28, '!9m!3 ',charsize=2
if (pflag eq 1) then output_plot,'color_mosaic.ps'
loadct,0

;contour map
twopi=2.0*!pi
r2d=360.0/twopi
x=ps_avg*r2d*(indgen(Nxs)-xc)
y=ps_avg*r2d*(indgen(Nys)-yc)
if (pflag eq 1) then setplot
xrng=[-35,35]
yrng=[-35,35]
log_image(sun)=0.0
contour,log_image,x,y,nlevels=7,xstyle=1, $
  ystyle=1,charsize=1.3,xrange=xrng,yrange=yrng, $
  xtitle='ecliptic longitude    (degrees)', $
  ytitle='ecliptic latitude    (degrees)'
plots,[0,0],yrng,color=128
plots,xrng,[0,0],color=128
angle=twopi*findgen(101)/100.0
x_sun=Rsun*r2d*cos(angle)
y_sun=Rsun*r2d*sin(angle)
;polyfill,x_sun,y_sun
;one_arrow, 650.0, 650.0, 90.0, 'N'
;one_arrow, 650.0, 650.0,180.0, 'E'
if (pflag eq 1) then output_plot,'contour.ps'

;ecliptic surface brightness profiles
if (pflag eq 0) then begin $
  window,xs=600,ys=600 & $
  clr=255 & $
endif
if (pflag eq 1) then begin $
  set_colorplot & $
  clr=0 & $
endif
loadct,39
prof=fltarr(Nxs)
for p=0,Nxs-1 do begin $
  im_sub=image(p,yc-1:yc+1) & $
  gd_sub=good(p,yc-1:yc+1) & $
  j=where(gd_sub ne 0) & $
  if (j(0) ne -1) then prof(p)=median(im_sub(j)) & $
  if (j(0) eq -1) then prof(p)=-10.0 & $
endfor
prof_east=rotate(prof(0:xc),2)
prof_west=prof(xc:*)
x=ps_avg*r2d*findgen(n_elements(prof_east))
N=n_elements(x)
p=where(prof_east gt -9.0)
plot_oo,x(p(2):*),prof_east(p(2):*)*cal,xrange=[1,40], $
  xstyle=1,yrange=[2.0e-13,5.0e-10],ystyle=1,color=clr,thick=5, $
  xtitle='ecliptic latitude & longitude    (degrees)', $
  charsize=1.3,ytitle= $
  'surface brightness    (mean solar brightness B!l!9n!3!n)'
xyouts,4.7,1300*cal,'East',charsize=1,color=clr
p=where((prof_west gt -1.0) and (x le 26.0))
oplot,x(p(2:*)),prof_west(p(2:*))*cal,color=53,thick=5
xyouts,4.7,450*cal,'West',charsize=1,color=53
xyouts,20.0,150.0*cal,'Venus',color=clr
p=where((x ge 4.2) and (x le 17.0))
coeff=poly_fit(alog(x(p)),alog(prof_east(p)),1,yfit)
;oplot,x(p),exp(yfit)*cal,psym=1
print,'east profile power-law index = ',coeff(1)

;latitude surface brightness profiles
;if (pflag eq 1) then setplot
prof=fltarr(Nys)
for p=0,Nys-1 do begin $
  im_sub=image(xc-1:xc+1,p) & $
  gd_sub=good(xc-1:xc+1,p) & $
  j=where(gd_sub ne 0) & $
  if (j(0) ne -1) then prof(p)=median(im_sub(j)) & $
  if (j(0) eq -1) then prof(p)=-10.0 & $
endfor
prof_south=rotate(prof(0:yc),2)
prof_north=prof(yc:*)
y=ps_avg*r2d*findgen(n_elements(prof_north))
N=n_elements(x)
p=where(prof_south gt -9.0)
oplot,y(p(2):*),prof_south(p(2):*)*cal,thick=5,color=147
xyouts,1.4,2470*cal,'South',charsize=1,color=147
p=where(prof_north gt -1.0)
oplot,y(p(2):*),prof_north(p(2):*)*cal,color=254,thick=5
xyouts,1.4,1900*cal,'North',charsize=1,color=254
p=where((x ge 3.0) and (x le 20.0))
coeff=poly_fit(alog(y(p)),alog(prof_south(p)),1,yfit)
;oplot,y(p),exp(yfit)*cal,psym=1,color=196
print,'south profile power-law index = ',coeff(1)
p=where(prof_east gt -1.0)
if (pflag eq 1) then output_plot,'profiles.ps'
loadct,0

;profile ratios
p=where((prof_east gt 0.0) and (prof_south gt 0.0))
plot,x(p),prof_east(p)/prof_south(p),yrange=[2,5.5], $
  xtitle='angular distance    (degrees)',charsize=1.3, $
  ytitle='East/South adn West/North surface brightness ratios', $
  xrange=[0,30]
p=where((prof_west gt 0.0) and (prof_north gt 0.0))
oplot,x(p),prof_west(p)/prof_north(p)

;Distance of the ZL midplane from the ecliptic
if (pflag eq 1) then setplot
image=image_old
good=good_old
dx=1
dy=25
Ndy=2*dy+1
x=ps_avg*r2d*(indgen(Nxs)-xc)
y=fltarr(Nxs)
for p=0,Nxs-1 do begin $
  sub_image=image((p-dx)>0:(p+dx)<(Nxs-1),yc-dy:yc+dy) & $
  sub_good = good((p-dx)>0:(p+dx)<(Nxs-1),yc-dy:yc+dy) & $
  sub_med=fltarr(Ndy) & $
  for j=0,Ndy-1 do begin $
    row_image=sub_image(*,j) & $
    row_good=sub_good(*,j) & $
    k=where(row_good gt 0) & $
    if (k(0) ne -1) then sub_med(j)=median(row_image(k)) & $
    if (k(0) eq -1) then sub_med(j)=-2*Ndy & $
  endfor & $
  l=where(sub_med eq max(sub_med))-dy & $
  y(p)=l(0) & $
endfor
j=where(good(*,yc) eq 0)
y(j)=-2*Ndy
y=y*ps_avg
plot,x,y*r2d,psym=1,yrange=2.5*[-1,1],xrange=[-35,35], $
  ystyle=1,xstyle=1,xtitle='ecliptic longitude    (degrees)', $
  ytitle='ecliptic latitude    (degrees)',charsize=1.3, $
  symsize=0.7
if (pflag eq 1) then output_plot,'midplane.ps'
