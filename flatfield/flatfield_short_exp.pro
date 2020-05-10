;flatfield_short_exp.pro

;Computes a flatfield for a given value of dark_current.

;Coordinates for the optical center.
x_c=191.0      		
y_c=286.0

;subtract additional dark-current offset
dark_current=35.0

;read Bonnie's flatfields
flat0=readfits('/u2/hahn/clementine/flatfields/BFB32060.fits',h0)
flat1=readfits('/u2/hahn/clementine/flatfields/BFB32061.fits',h1)
flat2=readfits('/u2/hahn/clementine/flatfields/BFB32062.fits',h2)
;flat3=readfits('/u2/hahn/clementine/flatfields/BFB62060.fits',h3)
;flat4=readfits('/u2/hahn/clementine/flatfields/BFB62061.fits',h4)
;flat5=readfits('/u2/hahn/clementine/flatfields/BFB62062.fits',h5)
sz=size(flat0)
Nx=sz(1)
Ny=sz(2)
window,xs=2*Nx,ys=Ny

;form the median-average of flat0 & flat1 & flat2,
;and its standard-deviation. These are 30 millisecond exposures.
f012=[[[flat0]],[[flat1]],[[flat2]]]
f012_avg=average_image(f012)
sigma012=stdev_image(f012)
;subtract dark-current & destreak
flat=f012_avg-dark_current
flat=destreak(flat,0.03)

;form the median-average of flat3 & flat4 & flat5,
;and its standard-deviation. These are 60 millisecond exposures.
;f345=[[[flat3]],[[flat4]],[[flat5]]]
;f345_avg=average_image(f345)
;sigma345=stdev_image(f345)
;subtract dark-current & destreak
;flat=f345_avg-dark_current
;flat=destreak(flat,0.06)

;Rotate, transpose, and shift such
;that the blemish is aligned with the clementine data
flat=shift(rotate(transpose(flat),3),-1,2)

;angles
plate_scale=0.00131926	;plate scale (rad/pixel)
x=fltarr(Nx,Ny)
xx=findgen(Nx)-x_c
for i=0,Ny-1 do x(0,i)=xx
y=fltarr(Nx,Ny)
yy=transpose(findgen(Ny)-y_c)
for i=0,Nx-1 do y(i,0)=yy
theta=asin(plate_scale*y)
phi=asin(plate_scale*x/cos(theta))
angle=acos(cos(phi)*cos(theta))
;pixel solid angle
psa=1.0/cos(angle)

;zap outer rows & columns
good_pixels=bytarr(Nx,Ny)+1
good_pixels(*,0:1)=0
good_pixels(0,*)=0
good_pixels(Nx-2:*,*)=0
good=where(good_pixels ne 0)
bad=where(good_pixels eq 0)

;normalize so that inner 6degrees has a median value of 1.0
flatfield=fltarr(Nx,Ny)
r2d=180./!pi
pix=where((angle*r2d lt 2.5) and (good_pixels ne 0))
flatfield(good)=flat(good)/median(flat(pix))
flatfield(bad)=1.0
raw_flatfield=flatfield

;adjust outer part of flatfield (which was vignetted by filter)
;due to vignetting by lens housing and inefficiencies in
;fiber-optic at large angles.
correction=1.0-1.95*angle^2
pix=where(angle*r2d gt 20.0)
err_pix=where((angle*r2d lt 2.5) and (good_pixels eq 1))
err=stdev(flatfield(err_pix))
av=median(flatfield(err_pix))
print,'fractional error = ',err/av
flatfield(pix)=correction(pix)+err*randomn(seed,n_elements(pix))
flatfield(bad)=1.0
plot,angle(good)*r2d,flatfield(good),psym=3
oplot,angle(good)*r2d,correction(good),psym=3,color=128
tvscl,flatfield>0.5<1.05,Nx,0

;save & display flatfield
writefits, 'flatfield_short_exp.fits', flatfield
writefits, 'good_pix_short_exp.fits', good_pixels
restore,'orbit193.dat'
good=where(good_pix ne 0)
flat_avg_image=fltarr(Nx,Ny)
flat_avg_image(good)=avg_image(good)/flatfield(good)
tvscl,enlarge(avg_image(70:170,70:170),5)>14<32
tvscl,enlarge(shift(flatfield(70:170,70:170),0,0),5)
tvscl,enlarge(flat_avg_image(70:170,70:170),5)>14<32

;profile flatfield
r2d=180./!pi
dr=0.05/r2d
N=max(angle)/dr
prof=fltarr(N)
sd=fltarr(N)
r=findgen(N)*dr
for i=0,N-2 do begin $
  pix=where((angle ge r(i)) and (angle lt r(i+1)) and $
    (good_pixels eq 1)) & $
  if (pix(0) eq -1) then r(i)=-1.0 & $
  if (pix(0) ne -1) then $
    prof(i)=median(flatfield(pix)) & $
    if (n_elements(pix) gt 1) then  $
      sd(i)=stdev(flatfield(pix)) & $
endfor
!p.multi=[0,1,2,0]
plot,r*r2d,prof,xrange=[0,30],xstyle=1,charsize=1.3, $
  xtitle='angular distance    (degrees)', $
  ytitle='median flatfield value',yrange=[0,1.2]
oplot,angle*r2d,psa,psym=3,color=128
j=where(prof gt 0.0)
plot,r(j)*r2d,sd(j)/prof(j),xrange=[0,30],xstyle=1,charsize=1.3,$
  xtitle='angular distance    (degrees)', $
  ytitle='!7r!3/flatfield',yrange=[0,0.1]
!p.multi=0

;display raw and final flatfield
pflag=1
AB=fltarr(2*Nx,Ny)
AB(0,0)=raw_flatfield
AB(Nx,0)=flatfield
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=24,ysize=18,xoffset=2.0,yoffset=26, $
    /landscape,bits_per_pixel=8 & $
endif
AB=AB>0.5<1.05
AB=AB-min(AB)
AB=fix(AB*255.0/max(AB))
tv,AB
xyouts,-2.6,-0.015,'A',charsize=2,charthick=4,color=255
xyouts,14.6,-0.015,'B',charsize=2,charthick=4,color=255
if (pflag eq 1) then output_plot,'flatfield.ps'
write_tiff,'flatfield.tif',transpose(rotate(AB,1))

stop

;contour flatfield
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,smooth(flatfield,5),nlevels=7
output_plot

stop

;optimize by looping over numerous values of dark-current
N=61
flatfields=fltarr(Nx,Ny,N)
images=fltarr(Nx,Ny,N)
dc=fltarr(N)
pix=where((angle*r2d lt 2.5) and (good_pixels ne 0))
good=where(good_pixels eq 1)
bad =where(good_pixels ne 1)
window,xs=1.5*Nx,ys=Ny
for index=0,N-1,1 do begin $
  dc(index)=dark_current-(index-N/2) & $
  flat=f012_avg-dc(index) & $
  flat=destreak(flat,0.03) & $
  flat=shift(rotate(transpose(flat),3),-1,2) & $
  ff=fltarr(Nx,Ny) & $
  ff(good)=flat(good)/median(flat(pix)) & $
  ff(bad)=1.0 & $
  flatfields(0,0,index)=ff & $
  flatflux=avg_image/ff & $
  images(0,0,index)=flatflux & $
  print,index,dc(index) & $
endfor

;display results
z=''
for index=0,N-1,1 do begin $
  tvscl,enlarge(images(70:170,70:170,index),5)>14<32 & $
  print,index,dc(index) & $
endfor
