;orbit66_beta_hyi.pro

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
start=1
stop =3

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
x0	=  311.0		;beta hydra
y0	=  462.0

;Sun's RA and DEC in radians
alpha_sun= 6.043145
delta_sun=-0.102731

;The half-width of the subarrays contain the alignment stars.
w=10

;If these data are to be used to measure the dark-current, then
;select the range of pixels (in the dark portion of the moon)
;where the dark-current will be measured. This is skipped if
;dark_file=''.
x_dark0=340
x_dark1=370
y_dark0=536
y_dark1=566

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

stop

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
;j=where((mask eq 1) and (images gt 70)) & plot,images(j),psym=3

;destreak the images & display original/destreaked images
images=float(images)
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  imd=destreak(im,exp_time(t)) & $
  images(0,0,t)=imd & $
  im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  imd=imd-median(imd(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<10 & $
  ;tvscl,imd*max(exp_time)*mask(*,*,t)/exp_time(t)>0<10,Nx,0 & $
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

;make the `point-source' flatfield
plate_scale=0.00131940	;plate scale (rad/pixel)
x_opt=191.0		;Coordinates for the optical center.
y_opt=286.0
xx=findgen(Nx)-x_opt
x=fltarr(Nx,Ny)
for i=0,Ny-1 do x(0,i)=xx
yy=transpose(findgen(Ny)-y_opt)
y=fltarr(Nx,Ny)
for i=0,Nx-1 do y(i,0)=yy
theta=asin(plate_scale*y)
phi=asin(plate_scale*x/cos(theta))
angle=acos(cos(phi)*cos(theta))
flatfield=readfits(flat_file)
flatfield=flatfield*cos(theta)
good_flat_pix=readfits(good_file)

;flatfield the data
template=intarr(Nx,Ny)
template(x_dark0:x_dark1,y_dark0:y_dark1)=1
for t=0,Nimages-1 do begin $
  im=images(*,*,t)/flatfield & $
  images(0,0,t)=im & $
  tvscl,(im+3.0*template)*max(exp_time)/exp_time(t)>0<10 & $
endfor

;align images
dx=[ 0,  11,  28 ]
dy=[ 0, -18, -52 ]
for t=0,Nimages-1 do begin $
  images(0,0,t)=shift(images(*,*,t),dx(t),dy(t)) & $
  tvscl,images(*,*,t)>0<8 & $
endfor
print,twist(0)*180./!pi
print,ra(0)*180./!pi
print,dec(0)*180./!pi

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
avg_image=fltarr(Nx,Ny)
good=where(good_pix gt 0)
avg_image(good)=counts(good)/total_exp_time(good)
tvscl,avg_image>0<50

;get calibration star flux for beta hydri
radius=5.0
l=12
sub=avg_image(x0-l:x0+l,y0-l:y0+l)
radial,rad,l,l,2*l+1,2*l+1
plot,rad,sub,psym=8
oplot,rad,rad*0,color=128
tvscl,enlarge(sub,5)>0<50
prof=profile_idl(sub,l,l,r,err,sub*0)
pix=where((rad ge 6.0) and (rad le 10.0))
background=avg(sub(pix)) 
sub=sub-background
prof=profile_idl(sub,l,l,r,err,sub*0)
plot,r,prof,psym=10,xrange=[0,15]
oploterr,r,prof,err,3
oplot,r,r*0,color=128
star_flux=r*0
for i=0,n_elements(r)-1 do begin $
  pix=where(rad le r(i)) & $
  if (pix(0) ne -1) then star_flux(i)=total(sub(pix)) & $
endfor
plot,r,star_flux,xrange=[0,15],psym=10
j=where(r ge radius)
star_flux=star_flux(j(0))
plots,r(j(0)),star_flux,psym=8
print,'star-flux = ',star_flux,' counts/sec'

;analyze individual frames to get error estimate
flux=fltarr(Nimages)
for t=0,Nimages-1 do begin $
  sub=images(x0-l:x0+l,y0-l:y0+l,t)/exp_time(t) & $
  pix=where((rad ge 6.0) and (rad le 10.0)) & $
  background=avg(sub(pix)) & $
  sub=sub-background & $
  pix=where(rad le radius) & $
  flux(t)=total(sub(pix)) & $
endfor
error=stdev(flux)/sqrt(n_elements(flux)-1.0)
print,'corrrected star-flux = ',star_flux,' counts/sec', $
  '+/- ',error

;check psf FWAHM
sub=images(211:221,50:60,0)
contour,sub
plot,sub(*,5)/max(sub),psym=10,yrange=[0,1.1],charsize=1.3, $
  xtitle='pixels',ytitle='normalized PSF flux'
oplot,sub(*,5)/max(sub),psym=8
oplot,sub(5,*)/max(sub),psym=10,linestyle=1
oplot,indgen(11),fltarr(11)+0.5,color=128

stop

;check the subarrays containing alignment stars.
erase
window,xs=1.5*Nx,ys=Ny
w=5
contour,alog(images(*,*,0)>1),xrange=x0+w*[-1,1], $
  yrange=y0+w*[-1,1],nlevels=9,xstyle=1,ystyle=1
plots,x0+w*[-1,1],y0+[0,0],color=128
plots,x0+[0,0],y0+w*[-1,1],color=128
contour,alog(images(*,*,1)>1),xrange=x0+w*[-1,1], $
  yrange=y0+w*[-1,1],nlevels=9,xstyle=1,ystyle=1
plots,x0+w*[-1,1],y0+[0,0],color=128
plots,x0+[0,0],y0+w*[-1,1],color=128
contour,alog(images(*,*,2)>1),xrange=x0+w*[-1,1], $
  yrange=y0+w*[-1,1],nlevels=9,xstyle=1,ystyle=1
plots,x0+w*[-1,1],y0+[0,0],color=128
plots,x0+[0,0],y0+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x0+w*[-1,1], $
  yrange=y0+w*[-1,1],nlevels=9,xstyle=1,ystyle=1
plots,x0+w*[-1,1],y0+[0,0],color=128
plots,x0+[0,0],y0+w*[-1,1],color=128
