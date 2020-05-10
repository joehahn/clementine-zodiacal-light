;orbit66_alpha_aqr.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Set input parameters.				;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The leading part of the input fits-filenames.
lead_name='/u2/hahn/clementine/maybe66/conversions/lba'

;The range of image-numbers (e.g., the middle part
;of the fits filenames). Note---images 31-84 are skipped
;due to the camera's motion, but those are still useful images!
start=16
stop =28

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ -1 ]

;The trailing part of the fits filenames.
trail_name='v.066.fits'

;flatfield & good-pixel filenames
flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield_short_exp.fits'
good_file='/users/hahn/clementine/zodiacal_light/flatfield/good_pix_short_exp.fits'

;The (x,y) coordinates of 3+ stars that are used to co-register
;these images, and their (alpha,delta)=(RA,DEC) in radians,
;which is used to orient the images.
x2	=  124			;Sadalmelik = alpha Aqr
y2	=  89

;Sun's RA and DEC in radians
alpha_sun= 6.043145
delta_sun=-0.102731

;The half-width of the subarrays contain the alignment stars.
wx=10		;in x-direction
wy=25		;in y-direction

;If dark_file='', then use these data to determine the dark-current.
;But set dark_file = an appropriate filename in order to restore &
;use the dark-current info stored in another dataset.
dark_file=''

;If these data are to be used to measure the dark-current, then
;select the range of pixels (in the dark portion of the moon)
;where the dark-current will be measured. This is skipped if
;dark_file=''.
x_dark0=330
x_dark1=370
y_dark0=25
y_dark1=65

;orbit number
orbit=66

;On first entry, set lunar_limb=0 in order to run the
;proceedure that determines the lunar limb. Otherwise,
;set lunar_limb=1 to skip this step.
lunar_limb=1

;Zap the following pixels because charge bled away from 
;a saturated object. Just set zap = [-1] to skip this step.
zap_mask=intarr(384,576)+1
zap=[ -1 ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Main program.					;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;This loop reads each image and destreaks it.
;First, generate an array of image numbers, skipping
;numbers as appropriate.
image_num=indgen(stop-start+1)+start
for j=0,n_elements(image_num)-1 do begin $
  k=where(skip eq image_num(j)) & $
  if (k(0) ne -1) then image_num(j)=-1 & $
endfor
image_num=image_num(where(image_num ge 0))
Nimages=n_elements(image_num)
exp_time=fltarr(Nimages)		;exposure time (seconds)
headers=strarr(140,Nimages)
RA=fltarr(Nimages)			;right ascension of boresight (rad)
DEC=fltarr(Nimages)			;declination of boresight (rad)
RA_ret=fltarr(Nimages,4)		;right ascension of reticle points (rad)
DEC_ret=fltarr(Nimages,4)		;declination of reticle points (rad)
twist=fltarr(Nimages)			;twist angle (rad)
date=strarr(Nimages)			;date & time
temp=fltarr(Nimages)			;focal plane temperature in Kelvin
gain_id=strarr(Nimages)			;gain mode ID
offset_id=strarr(Nimages)		;offset mode ID
zero='0'

for t=0,Nimages-1 do begin $

  ;check for skips
  j=0 & $
  while (j(0) ne -1) do begin $
    j=where(skip eq image_num(t)) & $
    if (j(0) ne -1) then t=t+1 & $
  endwhile & $

  ;generate input filename
  middle=strtrim(string(image_num(t)),2) & $
  for j=1,4-strlen(middle) do middle=zero+middle & $
  file=lead_name+middle+trail_name & $

  ;read image and get size & header
  im=readfits(file,header,/silent) & $
  headers(0,t)=header & $
  sz=size(im) & $
  Nx=sz(1) & $
  Ny=sz(2) & $
  if (t eq 0) then begin $
    ;some setup on first entry only
    window,xs=2*Nx,ys=Ny & $
    images=intarr(Nx,Ny,Nimages) & $
  endif & $
  images(0,0,t)=im & $

  ;get exposure time (in seconds), reticle RA & DEC and twist angle,
  ;(in radians), and time/date (string) from the header.
  j=where(strmid(header,0,25) eq 'COMMENT EXPOSURE_DURATION') & $
  exp_time(t)=float(strmid(header(j(0)),27,10))/1000.0 & $
  print,file,'     exposure = ',exp_time(t),' sec' & $
  j=where(strmid(header,0,23) eq 'COMMENT RIGHT_ASCENSION') & $
  RA(t)=float(strmid(header(j(0)),26,9))*!pi/180.0 & $
  j=where(strmid(header,0,19) eq 'COMMENT DECLINATION') & $
  DEC(t)=float(strmid(header(j(0)),26,9))*!pi/180.0 & $
  j=where(strmid(header,0,24) eq 'COMMENT RETICLE_POINT_RA') & $
  RA_ret(t,0)=float(strmid(header(j(0)),29,6))*!pi/180.0 & $
  RA_ret(t,1)=float(strmid(header(j(0)),36,7))*!pi/180.0 & $
  RA_ret(t,2)=float(strmid(header(j(0)),44,7))*!pi/180.0 & $
  RA_ret(t,3)=float(strmid(header(j(0)),52,7))*!pi/180.0 & $
  j=where(strmid(header,0,33) eq 'COMMENT RETICLE_POINT_DECLINATION') & $
  DEC_ret(t,0)=float(strmid(header(j(0)),37,7))*!pi/180.0 & $
  DEC_ret(t,1)=float(strmid(header(j(0)),45,7))*!pi/180.0 & $
  DEC_ret(t,2)=float(strmid(header(j(0)),53,7))*!pi/180.0 & $
  DEC_ret(t,3)=float(strmid(header(j(0)),61,7))*!pi/180.0 & $
  j=where(strmid(header,0,19) eq 'COMMENT TWIST_ANGLE') & $
  twist(t)=float(strmid(header(j(0)),26,8))*!pi/180.0 & $
  j=where(strmid(header,0,18) eq 'COMMENT START_TIME') & $
  date(t)=strmid(header(j(0)),21,50) & $
  j=where(strmid(header,0,31) eq 'COMMENT FOCAL_PLANE_TEMPERATURE') & $
  temp(t)=float(strmid(header(j(0)),33,10)) & $
  j=where(strmid(header,0,20) eq 'COMMENT GAIN_MODE_ID') & $
  gain_id(t)=strmid(header(j(0)),28,5) & $
  j=where(strmid(header,0,22) eq 'COMMENT OFFSET_MODE_ID') & $
  offset_id(t)=strmid(header(j(0)),28,5) & $

  ;display input image
  ;im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<30 & $

endfor

;Find the saturated pixels, and then create a mask array
;having the same dimensions as the images array. Good rows &
;columns devoid of saturated pixels have mask=1, whereas bad
;rows and columns containing saturated pixels have mask=0
mask=bytarr(Nx,Ny,Nimages)+1
for t=0,Nimages-1 do begin $
  for y=0,Ny-1 do begin $
    j=where(images(*,y,t) ge 255) & $
    if (j(0) ne -1) then mask(*,y,t)=0 & $
  endfor & $
  for x=0,Nx-1 do begin $
    j=where(images(x,*,t) ge 255) & $
    if (j(0) ne -1) then mask(x,*,t)=0 & $
  endfor & $
endfor
;j=where((mask eq 1) and (images gt 160)) & plot,images(j),psym=3

;destreak the images & display original/destreaked images
images=float(images)
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  imd=destreak(im,exp_time(t)) & $
  images(0,0,t)=imd & $
  im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  imd=imd-median(imd(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<100 & $
  ;tvscl,imd*max(exp_time)*mask(*,*,t)/exp_time(t)>0<100,Nx,0 & $
endfor

;Subtract an average dark-current that is measured
;on a dark portion of the sky or moon. Be sure that there are
;not any cosmic-ray hits contaminating this measurement!
;plot,images(x_dark0:x_dark1,y_dark0:y_dark1,*),psym=3
dark=fltarr(Nimages)
for t=0,Nimages-1 do begin $
  dark(t)=avg(images(x_dark0:x_dark1,y_dark0:y_dark1,t)) & $
  im=images(*,*,t)-dark(t) & $
  images(0,0,t)=im & $
endfor

;flatfield the data
flatfield=readfits(flat_file)
good_flat_pix=readfits(good_file)
template=intarr(Nx,Ny)
template(x_dark0:x_dark1,y_dark0:y_dark1)=1
for t=0,Nimages-1 do begin $
  im=images(*,*,t)/flatfield & $
  images(0,0,t)=im & $
  im=im*mask(*,*,t) & $
  tvscl,(im+3.0*template)*max(exp_time)/exp_time(t)>0<20 & $
endfor

;Check the alignment of each image using the 3 stars.
;Determine the offsets (dx0,dy0), etc., for each frame
;using star0 to a precision (but not accuracy) of 1/MAG=0.2
;pixels. Choose a long-exposure image as the reference image.
register_image,images,6,dx2,dy2,MAG=5,X0=x2-wx,X1=x2+wx,Y0=y2-wy,Y1=y2+wy
;additional adjustments
dx2(9)=dx2(9)-1
dy2(3)=dy2(3)+1
dy2(5)=dy2(5)+1
dy2(7)=dy2(7)+1
dy2(9)=dy2(9)+1
dy2(11)=dy2(11)+1
;plot the 3 different offsets for dx & dy, average them,
;and then round off to do full-pixel shifts.
!p.multi=[0,1,2,0,0]
plot,dx2,psym=8,yrange=6*[-1,1],xrange=[-1,Nimages+1],xstyle=1, $
  ystyle=1,xtitle='image index',ytitle='!7D!3x    (pixels)'
dx=round(dx2)
oplot,dx,thick=3,psym=10
plot,dy2,psym=8,yrange=25*[-1,1],xrange=[-1,Nimages+1],xstyle=1, $
  ystyle=1,xtitle='image index',ytitle='!7D!3y    (pixels)'
dy=round(dy2)
oplot,dy,thick=3,psym=10
!p.multi=0
;Shift images and zap the outer edges
for t=0,Nimages-1 do begin $
  images(0,0,t)=shift(images(*,*,t),dx(t),dy(t)) & $
  mask(0,0,t)=shift(mask(*,*,t),dx(t),dy(t)) & $
endfor
ddx=max(abs(dx))+1
ddy=max(abs(dy))+1
zap_mask(0:ddx,*)=0
zap_mask(Nx-1-ddx:*,*)=0
zap_mask(*,0:ddy)=0
zap_mask(*,Ny-1-ddy:*)=0

;Coadd all good pixels (eg., where limb & mask both = 1)
limb_mask=mask
good_pix=bytarr(Nx,Ny)
counts=fltarr(Nx,Ny)		;summed counts in good pixels
total_exp_time=fltarr(Nx,Ny)	;summed exp. time in good pixels
for x=0,Nx-1 do begin $
  for y=0,Ny-1 do begin $
    j=where(limb_mask(x,y,*) eq 1) & $
    if (j(0) ne -1) then begin $
      counts(x,y)=total(images(x,y,j)) & $
      total_exp_time(x,y)=total(exp_time(j)) & $
      good_pix(x,y)=1 & $
    endif & $
  endfor & $
endfor

;zap any bad pixels
good_pix=good_pix*good_flat_pix*zap_mask
counts=counts*good_flat_pix*zap_mask
good=where(good_pix eq 1)
avg_image=fltarr(Nx,Ny)
avg_image(good)=counts(good)/total_exp_time(good)

;make pictures
im=alog(avg_image>1)
im=im-min(im)
im=byte(255*im/max(im))
window,xs=1.5*Nx,ys=Ny
tvscl,im

;get calibration star flux for alpha Aquarius
radius=5.0		;radius of aperture
l=12			;half-width of subarray centered at star
sub=avg_image(x2-l:x2+l,y2-l:y2+l)
radial,rad,l,l,2*l+1,2*l+1
plot,rad,sub,psym=8
prof=profile_idl(sub,l,l,r,err,sub*0)
pix=where((rad ge 6.0) and (rad le 10.0))
background=median(sub(pix)) 	;background in aperture
sub=sub-background
prof=profile_idl(sub,l,l,r,err,sub*0)
plot,r,prof,psym=10
oploterr,r,prof,err,3
oplot,r,r*0,color=128
star_flux=r*0			;integrated flux with distance
for i=0,n_elements(r)-1 do begin $
  pix=where(rad le r(i)) & $
  if (pix(0) ne -1) then star_flux(i)=total(sub(pix)) & $
endfor
plot,r,star_flux,xrange=[0,15],psym=10
j=where(r ge radius)
star_flux=star_flux(j(0))
plots,r(j(0)),star_flux,psym=8
print,'star-flux = ',star_flux,' counts/sec'
;adjust result for flatfield effects
x_c=191.0       ;Coordinates for the optical center.
y_c=286.0
t_c=bytarr(Nx,Ny)
t_c(x_c-radius:x_c+radius,y_c-radius:y_c+radius)=1
pix_c=where(t_c eq 1)
f_c=median(flatfield(pix_c))
t_star=bytarr(Nx,Ny)
t_star(x2-radius+30:x2+radius+30,y2-radius:y2+radius)=1
pix_star=where(t_star eq 1)
f_star=median(flatfield(pix_star))
tvscl,(avg_image+10.0*(t_c+t_star))>0<250
tvscl,(flatfield+0.3*(t_c+t_star)),Nx,0
plot,flatfield(pix_c),psym=1,yrange=[0.6,1.1]
oplot,flatfield(pix_star),psym=6
correction=f_star/f_c
print,'correction factor = ',correction
print,'fractional variation in flatfield = ', $
  stdev(flatfield(pix_star))/avg(flatfield(pix_star))

;analyze individual frames to get error estimate
flux=fltarr(Nimages)
for t=0,Nimages-1 do begin $
  sub=images(x2-l:x2+l,y2-l:y2+l,t)/exp_time(t) & $
  pix=where((rad ge 6.0) and (rad le 10.0)) & $
  background=avg(sub(pix)) & $
  sub=sub-background & $
  pix=where(rad le radius) & $
  flux(t)=total(sub(pix))*correction & $
endfor
error=stdev(flux)/sqrt(n_elements(flux)-1.0)
print,'corrrected star-flux = ',star_flux*correction,' counts/sec', $
  '+/- ',error

stop

;check the subarrays containing alignment stars.
erase
window,xs=1.5*Nx,ys=Ny
w=5
t=-1
t=t+1
sub=images(x2-w:x2+w,y2-w:y2+w,t)
contour,alog(sub>1),nlevels=9,xstyle=1,ystyle=1,title=string(t)
plots,[0,2*w+1],[w,w],color=128
plots,[w,w],[0,2*w+1],color=128

sub=avg_image(x2-w:x2+w,y2-w:y2+w)
contour,alog(sub>1),nlevels=9,xstyle=1,ystyle=1
plots,[0,2*w+1],[w,w],color=128
plots,[w,w],[0,2*w+1],color=128

for t=0,Nimages-1 do begin $
  sub=images(x2-wx:x2+wx,y2-wy:y2+wy,t) & $
  tvscl,enlarge(alog(sub>1),5) & $
  wait,0.1 & $
endfor
