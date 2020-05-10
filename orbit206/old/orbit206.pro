;orbit206.pro

;These data have problems: contours look funny.
;Double-check the dark-current & residual_flux.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Set input parameters.				;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The leading part of the input fits-filenames.
lead_name='/u2/hahn/clementine/new206/conversions/lba'

;The range of image-numbers (e.g., the middle part
;of the fits filenames). Note: images 39-40 are defective.
start= 1
stop = 38

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ -1 ]

;The trailing part of the fits filenames.
trail_name='v.206.fits'

;flatfield & good-pixel filenames
flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield_short_exp.fits'
good_file='/users/hahn/clementine/zodiacal_light/flatfield/good_pix_short_exp.fits'

;The (x,y) coordinates of 3+ stars that are used to co-register
;these images, and their (alpha,delta)=(RA,DEC) in radians,
;which is used to orient the images.
x0=94.3				;Theta Cetus
y0=200.55
alpha0	=  0.366621
delta0	= -0.142825
x1 =  86.0			;Mira = Omicron Cetus
y1 =  394.5
alpha1	=  0.608014
delta1	= -0.051970
x2=90.1				;Eta Cetus
y2=144.5
alpha2	=  0.299280
delta2	= -0.177714
x3	=  177.85		;alpha Pisces
y3	=  372.2
alpha3	= 0.532530
delta3	= 0.048237
x4	=  276.9		;omicron Pisces
y4	=  354.5
alpha4	= 0.459867
delta4	= 0.159833
x5	=  369.0		;eta Pisces
y5	=  344.5
alpha5	= 0.399171
delta5	= 0.267835
x6	=  131.9		;gamma Cetus
y6	=  497.5
alpha6	= 0.712531
delta6	= 0.056476
x7	=  214.5		;73 Cetus
y7	=  478.5
alpha7	= 0.646463
delta7	= 0.147656
x8	=  33.4			;Baten
y8	=  272.0
alpha8	= 0.486338
delta8	=-0.180380
x9	=  168.1		;iota cetus (not used)
y9	=  7.8
alpha9	= 0.084770
delta9	=-0.154006
x10	=  119.1		;Menkar (not used)
y10	=  556.4
alpha10	=  0.795344
delta10	=  0.071379

;Sun's RA and DEC in radians
alpha_sun=0.224238
delta_sun=0.096098

;The half-width of the subarrays contain the alignment stars.
w=7

;If dark_file='', then use these data to determine the dark-current.
;But set dark_file = an appropriate filename in order to restore &
;use the dark-current info stored in another dataset.
dark_file='/users/hahn/clementine/zodiacal_light/orbit193/dark.dat'

;If these data are to be used to measure the dark-current, then
;select the range of pixels (in the dark portion of the moon)
;where the dark-current will be measured. This is skipped if
;dark_file=''.
x_dark0=10
x_dark1=30
y_dark0=546
y_dark1=566

;The above aperture might not be completely dark (due to a faint
;zodiacal light or earthshine contribution), so the dark-current
;inferred in this aperture needs to be reduced by the following
;flux in counts/second. This residual can be determined using
;mosaic.pro
residual_flux=11.0

;orbit number
orbit=206

;On first entry, set lunar_limb=0 in order to run the
;proceedure that determines the lunar limb. Otherwise,
;set lunar_limb=1 to skip this step.
lunar_limb=1

;Zap the following pixels because charge bled away from 
;a saturated object; the value of the zapped pixels is
;replaced by the median value of the replacement_pix(els).
;Just set zap = [-1] to skip this step.
zap_mask=intarr(384,576)+1
zap_mask(382:*,*)=0		;zap right edge
;zap_mask(*,418:424)=0		;zap rows
zap_mask(287:291,*)=0		;  & columns near Venus
;radial,rad,289,421,384,576	;DONT zap Venus, at least
;j=where(rad le 30.)		;not for the DPS poster...
;zap_mask(j)=0
;replacement_pix=where((rad le 31.) and (rad ge 30.0))
zap=where(zap_mask eq 0)

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
  ;tvscl,im*max(exp_time)/exp_time(t)>0<100 & $

endfor

;Find the saturated pixels, and then create a mask array
;having the same dimensions as the images array. Good columns
;devoid of saturated pixels have mask=1, whereas bad columns
;containing saturated pixels have mask=0
mask=bytarr(Nx,Ny,Nimages)+1
for t=0,Nimages-1 do begin $
  for y=0,Ny-1 do begin $
    j=where(images(*,y,t) ge 255) & $
    ;if (j(0) ne -1) then mask(*,y,t)=0 & $
  endfor & $
  for x=0,Nx-1 do begin $
    j=where(images(x,*,t) ge 255) & $
    if (j(0) ne -1) then mask(x,*,t)=0 & $
  endfor & $
endfor
;j=where((mask eq 1) and (images gt 160)) & plot,images(j),psym=3

;subtract dark current & destreak images
coeff=[ 35.62, 13.33]		;anticipated dark-current
dark=coeff(0)+coeff(1)*exp_time
images_dd=fltarr(Nx,Ny,Nimages)
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  im=im-dark(t) & $
  imd=destreak(im,exp_time(t)) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<100 & $
  ;tvscl,imd*mask(*,*,t)*max(exp_time)/exp_time(t)>0<100,Nx,0 &$
  images_dd(0,0,t)=imd & $
endfor
plot,exp_time,exp_time,yrange=70*[-1,1],nodata=1,ystyle=1, $
  xtitle='exposure time    (sec)',charsize=1.3, $
  ytitle='residual dark current    (counts/sec)'
residual=fltarr(Nimages)
flatfield=readfits(flat_file)
sub_flat=flatfield(x_dark0:x_dark1,y_dark0:y_dark1)
for t=0,Nimages-1 do begin $
  sub=images_dd(x_dark0:x_dark1,y_dark0:y_dark1,t)/ $
    sub_flat & $
  oplot,exp_time(t)+sub*0,sub/exp_time(t),psym=3 & $
  j=where(exp_time eq exp_time(t)) & $
  residual(t)=avg( $
    images_dd(x_dark0:x_dark1,y_dark0:y_dark1,j)) & $
  plots,exp_time(t),residual(t)/exp_time(t),psym=6 & $
endfor
coeff_adj=poly_fit(exp_time,residual,1)
s=sort(exp_time)
oplot,exp_time(s),coeff_adj(0)/exp_time(s)+coeff_adj(1),thick=3
oplot,exp_time,exp_time*0+residual_flux,color=128
print,transpose(coeff)
print,transpose(coeff+coeff_adj-[0.0,residual_flux])
images=images_dd
flux=fltarr(Nx,Ny,Nimages)
for t=0,Nimages-1 do flux(0,0,t)=images(*,*,t)/exp_time(t)
plot,flux(x_dark0:x_dark1,y_dark0:y_dark1,*),psym=3, $
  xstyle=1
oplot,flux(x_dark0:x_dark1,y_dark0:y_dark1,*)*0.0+ $
  residual_flux,color=128,thick=3

;flatfield the data
flatfield=readfits(flat_file)*0.0+1.0
good_flat_pix=readfits(good_file)
template=intarr(Nx,Ny)
template(x_dark0:x_dark1,y_dark0:y_dark1)=1
for t=0,Nimages-1 do begin $
  im=images(*,*,t)/flatfield & $
  images(0,0,t)=im & $
  ;tvscl,(im+3.0*template)*max(exp_time)/exp_time(t)>0<20 & $
endfor

;check for any time-varying residual flux
time=fltarr(Nimages)
flux=fltarr(Nimages)
err=fltarr(Nimages)
for t=0,Nimages-1 do begin $
  time(t)=total(exp_time(0:t)) & $
  sub=images(x_dark0:x_dark1,y_dark0:y_dark1,t)/exp_time(t) & $
  flux(t)=avg(sub) & $
  err(t)=stdev(sub)/sqrt(n_elements(sub)) & $
endfor
plot,time,flux,psym=6,xtitle='time    (sec)',charsize=1.3, $
  ytitle='residual flux    (counts/sec)'
oploterr,time,flux,err,3
oplot,time,time*0+residual_flux,color=128
;remove time-varying residual flux
for t=0,Nimages-1 do begin $
  sub=images(x_dark0:x_dark1,y_dark0:y_dark1,t)/exp_time(t) & $
  a=avg(sub) & $
  images(0,0,t)=images(*,*,t)-(a-residual_flux)*exp_time(t) & $
endfor

;One first entry (e.g., if lunar_limb=0), click the mouse
;along four well-spaced points to determine the location
;of the lunar limb. Otherwise, restore the limb-data
;acquired during a previous iteration.
if (lunar_limb eq 0) then begin $
  xc=fltarr(Nimages) & $
  yc=fltarr(Nimages) & $
  R_moon=fltarr(Nimages) & $
  window,xs=Nx,ys=Ny & $
  for t=0,Nimages-1 do begin $
    im=images(*,*,t) & $
    im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
    im=im*0.7/exp_time(t)>0<saturate/2 & $
    limb_click,im,x,y,r & $
    xc(t)=x & $
    yc(t)=y & $
    R_moon(t)=r & $
    radial,rad,x,y,Nx,Ny & $
    pix=where(rad le r) & $
    lmb=bytarr(Nx,Ny)+1 & $
    lmb(pix)=0 & $
    tvscl,im*lmb & $
    print,t,exp_time(t) & $
  endfor & $
  save,filename='limb.dat',xc,yc,R_moon & $
endif

;Mask the moon. Also, increase the radius of the moon by
;a few pixels to compensate for any errors and to mask any
;lunar horizon glow.
restore,filename='limb.dat'
R_moon=R_moon+3
limb=bytarr(Nx,Ny,Nimages)
for t=0,Nimages-1 do begin $
  radial,rad,xc(t),yc(t),Nx,Ny & $
  pix=where(rad le R_moon(t)) & $
  lmb=bytarr(Nx,Ny)+1 & $
  lmb(pix)=0 & $
  limb(0,0,t)=lmb & $
  im=images(*,*,t) & $
  im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  ;tvscl,im*limb(*,*,t)*max(exp_time)/exp_time(t)>0<100 & $
endfor

;Check the alignment of each image using the 3 stars.
;Determine the offsets (dx0,dy0), etc., for each frame
;using star0 to a precision (but not accuracy) of 1/MAG=0.2
;pixels. Choose a long-exposure image as the reference image.
register_image,images,15,dx0,dy0,MAG=5,X0=x0-w,X1=x0+w,Y0=y0-w,Y1=y0+w
;using star1
register_image,images,15,dx1,dy1,MAG=5,X0=x1-w,X1=x1+w,Y0=y1-w,Y1=y1+w
;using star2
register_image,images,15,dx2,dy2,MAG=5,X0=x2-w,X1=x2+w,Y0=y2-w,Y1=y2+w
;plot the 3 different offsets for dx & dy, average them,
;and then round off to do full-pixel shifts.
!p.multi=[0,1,2,0,0]
plot,dx0,psym=8,yrange=3*[-1,1],xrange=[-1,Nimages+1],xstyle=1, $
  ystyle=1,xtitle='image index',ytitle='!7D!3x    (pixels)'
oplot,dx1,psym=6
oplot,dx2,psym=1,symsize=1.5	;don't use star2---too faint!
dx=round((dx0+dx1+dx2)/3.0)
oplot,dx,thick=3,psym=10
plot,dy0,psym=8,yrange=3*[-1,1],xrange=[-1,Nimages+1],xstyle=1, $
  ystyle=1,xtitle='image index',ytitle='!7D!3y    (pixels)'
oplot,dy1,psym=6
oplot,dy2,psym=1,symsize=1.5
dy=round((dy0+dy1+dy2)/3.0)
oplot,dy,thick=3,psym=10
!p.multi=0

;Coadd all good pixels (eg., where limb & mask both = 1)
limb_mask=limb*mask
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

;get plate-scale
x_opt=191.0	;Coordinates for the optical center.
y_opt=286.0
x_star=[ x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 ]	
y_star=[ y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10 ]
RA_star = [ alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha10 ]
DEC_star= [ delta0, delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8, delta9, delta10 ]
.run plate_scale
plate_scale, x_opt, y_opt, x_star, y_star, RA_star, DEC_star, $
  plate_scl, errors
print,'plate-scale = ',plate_scl,' radians/pixel'
print,min(errors)

;get RA & DEC
twist_angle=avg(twist)
RA_opt=avg(RA)
DEC_opt=avg(DEC)
.run RA_DEC
RA_DEC, plate_scl, x_opt, y_opt, x_star, y_star, RA_star, $
  DEC_star, twist_angle, RA_opt, DEC_opt, Nx, Ny, RA_array, $
  DEC_array

;get ecliptic coordinates, where phi is the longitude
;relative to the Sun, and theta is the ecliptic latitude
obl=0.4090928			;Earth's obliquity in radians	
x= cos(obl)*sin(RA_array)*cos(DEC_array)+sin(obl)*sin(DEC_array)
y=-sin(obl)*sin(RA_array)*cos(DEC_array)+cos(obl)*sin(DEC_array)
z= cos(RA_array)*cos(DEC_array)
phi=atan(x,z)
theta=asin(y)
phi_sun=acos(cos(alpha_sun)*cos(delta_sun))
if (alpha_sun gt !pi) then phi_sun=-phi_sun
phi=phi-phi_sun

;zap any bad pixels
if (zap(0) eq -1) then zap_mask=intarr(Nx,Ny)+1
good_pix=good_pix*good_flat_pix*zap_mask
counts=counts*good_flat_pix*zap_mask
good=where(good_pix eq 1)
avg_image=fltarr(Nx,Ny)
avg_image(good)=counts(good)/total_exp_time(good)

;make pictures
im=alog(avg_image>3<3000)
im=im-min(im)
im=byte(255*im/max(im))
d=sqrt(phi^2+theta^2)
sun_pix=where(d eq min(d))
im(sun_pix)=255
window,xs=Nx,ys=Ny
tvscl,im
contour,alog(smooth(avg_image,4)>3<3000),xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels',levels=findgen(30)*0.25+1
print,'Observation date = ',date(Nimages/2)

;save results
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,alog(smooth(avg_image,4)>3<3000),xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels',levels=findgen(9)*0.25+1
file='orbit'+strtrim(string(orbit)+'.ps',2)
output_plot,file
file='orbit'+strtrim(string(orbit)+'.dat',2)
save,filename=file,plate_scl,twist_angle,temp,gain_id, $
  offset_id,phi,theta,counts,total_exp_time,good_pix,avg_image

stop

;check the subarrays containing alignment stars.
erase
window,xs=1.5*Nx,ys=Ny
contour,alog(avg_image>1),xrange=x0+w*[-1,1], $
  yrange=y0+w*[-1,1],nlevels=9,xstyle=1,ystyle=1
plots,x0+w*[-1,1],y0+[0,0],color=128
plots,x0+[0,0],y0+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x1+w*[-1,1], $
  yrange=y1+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x1+w*[-1,1],y1+[0,0],color=128
plots,x1+[0,0],y1+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x2+w*[-1,1], $
  yrange=y2+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x2+w*[-1,1],y2+[0,0],color=128
plots,x2+[0,0],y2+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x3+w*[-1,1], $
  yrange=y3+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x3+w*[-1,1],y3+[0,0],color=128
plots,x3+[0,0],y3+w*[-1,1],color=128
contour,alog(avg_image>1<120),xrange=x4+w*[-1,1], $
  yrange=y4+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x4+w*[-1,1],y4+[0,0],color=128
plots,x4+[0,0],y4+w*[-1,1],color=128
contour,alog(avg_image>1<120),xrange=x5+w*[-1,1], $
  yrange=y5+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x5+w*[-1,1],y5+[0,0],color=128
plots,x5+[0,0],y5+w*[-1,1],color=128
contour,alog(avg_image>1<120),xrange=x6+w*[-1,1], $
  yrange=y6+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x6+w*[-1,1],y6+[0,0],color=128
plots,x6+[0,0],y6+w*[-1,1],color=128
contour,alog(avg_image>1<120),xrange=x7+w*[-1,1], $
  yrange=y7+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x7+w*[-1,1],y7+[0,0],color=128
plots,x7+[0,0],y7+w*[-1,1],color=128
contour,alog(avg_image>1<120),xrange=x8+w*[-1,1], $
  yrange=y8+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x8+w*[-1,1],y8+[0,0],color=128
plots,x8+[0,0],y8+w*[-1,1],color=128
contour,alog(avg_image>1<120),xrange=x9+w*[-1,1], $
  yrange=y9+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x9+w*[-1,1],y9+[0,0],color=128
plots,x9+[0,0],y9+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x10+w*[-1,1], $
  yrange=y10+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x10+w*[-1,1],y10+[0,0],color=128
plots,x10+[0,0],y10+w*[-1,1],color=128
for t=0,Nimages-1 do begin $
  tvscl,images(*,*,t)*max(exp_time)/exp_time(t)>0<50 & $
  sub0=images(x0-w>0:x0+w<Nx-1,y0-w>0:y0+w<Ny-1,t) & $
  tvscl,sub0*max(exp_time)/exp_time(t)>0<50,Nx,y0-w & $
  sub1=images(x1-w>0:x1+w<Nx-1,y1-w>0:y1+w<Ny-1,t) & $
  tvscl,sub1*max(exp_time)/exp_time(t)>0<50,Nx,y1-w & $
  sub2=images(x2-w>0:x2+w<Nx-1,y2-w>0:y2+w<Ny-1,t) & $
  tvscl,sub2*max(exp_time)/exp_time(t)>0<50,Nx+2*w,y2-w & $
endfor

;fluxes
flux=fltarr(Nx,Ny,Nimages)
for t=0,Nimages-1 do begin $
  f=fltarr(Nx,Ny) & $
  im=images(*,*,t) & $
  pix=where(limb_mask(*,*,t) ne 0) & $
  f(pix)=im(pix)/exp_time(t) & $
  flux(0,0,t)=f & $
endfor

;check rows & columns near Venus
y=425 & plot,avg_image(*,427)<40,thick=3 & oplot,avg_image(*,y),color=200
x=288 & plot,avg_image(292,*),xrange=[300,500],thick=3 & oplot,avg_image(x,*)


flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield_short_exp.fits'
flatfield=readfits(flat_file)
sz=size(flatfield)
Nx=sz(1)
Ny=sz(2)
contour,smooth(flatfield,5),nlevels=10
fr=transpose(rotate(flatfield,3))
radial,rad,30,300,Nx,Ny
j=where(rad le 150.0)
tmp=bytarr(Nx,Ny)
tmp(j)=1
contour,tmp,nlevels=1,noerase=1
flat2=flatfield
flat2(j)=fr(j)
contour,smooth(flat2,5),nlevels=10,noerase=1
im1=alog(smooth(avg_image,3)>3<3000)
im2=alog(smooth(avg_image*flatfield,3)>3<3000)
contour,im1,nlevels=20,color=128
contour,im2,nlevels=20,noerase=1

