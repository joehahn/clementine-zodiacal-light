;LPSC2001_mosaic.pro
;generates figures for 2001 LPSC abstract

;calibration constant coverts 1 ADU/sec to mean solar brightness
cal=1.007e-13

;Filenames for each orbit. Orbit 247 gives ugly results, not-used
files=['orbit66.dat',  'orbit110.dat', 'orbit110_sunny.dat',  $
       'orbit164.dat', 'orbit193.dat', 'orbit206.dat' ]

;dark-current corrections for flatfield
dc_correction=[ 0.0, 8.6, 15.0, 7.0, 1.3, 6.95 ]

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

;resample image into ecliptic coordinate system,
;and smooth over a 2x2 box
q=atan(tan(theta)/cos(phi))
eta=tan(theta)/cos(phi)
chi=cos(q)*tan(phi)/cos(q)
ps_avg=avg(plate_scales)
x=-chi/ps_avg
xc=-min(x)
x=x+xc
Nxs=round(max(x))+2
y=eta/ps_avg
yc=-min(y)
y=y+yc
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

;trim image
mosaic=mosaic(0:850,0:800)
good=good(0:850,0:800)
sz=size(mosaic)
Nxs=sz(1)
Nys=sz(2)

window,xs=Nxs,ys=Nys
print,'average plate-scale = ',ps_avg,' radians/pixel'
print,'average plate-scale = ',ps_avg*180./!pi,' degrees/pixel'
twopi=2.0*!pi
r2d=360./twopi
Rsun=4.64e-3			;angular radius of Sun, radians
radial,rad,xc,yc,Nxs,Nys
rad=rad*ps_avg
sun=where(rad le Rsun)

;tv greyscale
log_mosaic=alog(mosaic>1)
log_mosaic=255.0*log_mosaic/max(log_mosaic)
log_mosaic(sun)=255.0
tv,log_mosaic

stop

;tv colormap
set_plot,'ps'
device,xsize=16,ysize=14.1,xoffset=3,yoffset=8, $
  /color,/portrait,bits_per_pixel=8
loadct,39
log_mosaic(sun)=192.0
contour,log_mosaic,xstyle=5,ystyle=5,nlevels=5,nodata=1
tv,log_mosaic,!x.window(0),!y.window(0), $
  xsize=(!x.window(1)-!x.window(0)), $
  ysize=(!y.window(1)-!y.window(0)),/norm
log_mosaic(sun)=0.0
contour,mosaic,xstyle=5,ystyle=5, $
  levels=[37.0,111.,333.,1000.],noerase=1
contour,smooth(mosaic,5),xstyle=5,ystyle=5, $
  levels=[12.3],noerase=1
contour,smooth(mosaic,9),xstyle=5,ystyle=5, $
  levels=[6.15],noerase=1
x=r2d*atan(ps_avg*(indgen(Nxs)-xc))
dx=abs(10*fix(x/10)-x)
for i=0,Nxs-1 do begin $
  if (dx(i) lt 0.07) then begin $
    plots,[i,i],[0,-15] & $
    xyouts,i-140,-50,string(round(x(i))),/data & $
  endif & $
endfor
xyouts,255,-90,'ecliptic longitude    (degrees)'
y=r2d*atan(ps_avg*(indgen(Nys)-yc))
dy=abs(10*fix(y/10)-y)
for i=0,Nys-1 do begin $
  if (dy(i) lt 0.07) then begin $
    plots,[-15,0],[i,i] & $
    xyouts,-195,i-10,string(round(y(i))) & $
  endif & $
endfor
xyouts,-100,195,'ecliptic latitude    (degrees)',orient=90
output_plot,'color_mosaic_LPSC2001.ps'

;ecliptic surface brightness profiles
window,xs=600,ys=600
loadct,0
prof=fltarr(Nxs)
for p=0,Nxs-1 do begin $
  im_sub=mosaic(p,yc-1:yc+1) & $
  gd_sub=good(p,yc-1:yc+1) & $
  j=where(gd_sub ne 0) & $
  if (j(0) ne -1) then prof(p)=median(im_sub(j)) & $
  if (j(0) eq -1) then prof(p)=-10.0 & $
endfor
prof_east=rotate(prof(0:xc),2)
prof_west=prof(xc:*)
;zap innermost value
j=where(prof_east gt -10.0)
prof_east(j(0))=prof_east(0)
j=where(prof_west gt -10.0)
prof_west(j(0))=prof_west(0)
;generate average east-west profile
Nee=n_elements(prof_east)
Nw=n_elements(prof_west)
prof_ew=fltarr(max([Nee,Nw]))-10.0
for i=0,n_elements(prof_east)-1 do begin $
  if ((i lt Nee) and (i lt Nw)) then begin $
    if ((prof_east(i) gt -9.0) and (prof_west(i) gt -9.0)) then begin $
      prof_ew(i)=0.5*(prof_east(i)+prof_west(i)) & $
    endif else begin $
      if (prof_east(i) gt -9.0) then begin $
        prof_ew(i)=prof_east(i) & $
      endif else prof_ew(i)=prof_west(i) & $
    endelse & $
  endif & $
  if (i ge Nee) then prof_ew(i)=prof_west(i) & $
  if (i ge Nw) then prof_ew(i)=prof_east(i) & $
endfor
;replace Venus with east-profiles
x=r2d*atan(ps_avg*findgen(n_elements(prof_east)))
j=where((x gt 13) and (x lt 21.0))
prof_ew(j)=prof_west(j)
;replace Mercury with west-profile
j=where((x gt 27.0) and (x lt 27.7))
prof_ew(j)=prof_east(j)
plot_io,x,prof_east>1,xrange=[10,30]
oplot,x,prof_west
oplot,x,prof_ew,psym=1
;plot profiles
N=n_elements(x)
setplot
plot,x,prof_ew*cal,xrange=[3.5,35],xlog=1,ylog=1, $
  xstyle=1,yrange=[5.0e-13,3.0e-10],ystyle=1,thick=3, $
  xtitle='elongation !7e!3    (degrees)',charsize=1.3, $
  ytitle='ecliptic surface brightness Z(!7e!3)   (B!l!9n!3!n)', $
  xticks=8,xtickv=[4,5,6,7,8,9,10,20,30]
;oplot,x,prof_east*cal,color=128
;oplot,x,prof_west*cal,color=128
p=where((x ge 4.0) and (x le 13.0) and (prof_ew gt 0.0))
coeff=poly_fit(alog(sin(x(p)/r2d)),alog(prof_ew(p)),1,yfit)
;oplot,x(p),exp(yfit)*cal,color=128,psym=1
output_plot,'profiles_LPSC2001.ps'

print,'east-west power law index = ',coeff(1)
print,'alpha = ',-coeff(1)-1.0
print,exp(coeff(0))*cal
