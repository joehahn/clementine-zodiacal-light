;orbit110_sunny.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Set input parameters.				;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The leading part of the input fits-filenames.
lead_name='/u2/hahn/clementine/new110/conversions/lba'

;NEW!
;The range of image-numbers (e.g., the middle part
;of the fits filenames). Use only those frames containing
;Mercury that aren't moving too much.
start= 61
stop = 68

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ -1 ]

;The trailing part of the fits filenames.
trail_name='v.110.fits'

;NEW!
;flatfield & good-pixel filenames
flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield.fits'
good_file='/users/hahn/clementine/zodiacal_light/flatfield/good_pix.fits'

;NEW!
;The (x,y) coordinates of 3+ stars that are used to co-register
;these images, and their (alpha,delta)=(RA,DEC) in radians,
;which is used to orient the images.
x0=19.8				;Mercury *
y0=478.6
alpha0	=  5.747633
delta0	= -0.224184
x1=141.9			;Sadalmelik *
y1=368.9
alpha1	=  5.784821
delta1	= -0.005580
x2=124.2   			;Saturn *
y2=526.2
alpha2	=  5.896409
delta2	= -0.189457
x3=11.5	   		        ;Sadalsud *
y3=353.5
alpha3	=  5.635486
delta3	= -0.097234
x4=168.65			;Eniph *
y4=218.9
alpha4	=  5.690573
delta4	=  0.172351
x5=111.5	  	        ;Aqr-Theta4 *
y5=467.25
alpha5	=  5.833036
delta5	= -0.135845
x6=81.3			        ;Del Gamma 12 *
y6=42
alpha6	=  5.439537
delta6  =  0.281420
x7=175.5		        ;Sadachbia *
y7=411.6
alpha7  =  5.854082
delta7  = -0.024212
x8=47                           ; Rotanev
y8=39.25
alpha8  =  5.399831
delta8  =  0.254731

;Sun's RA and DEC in radians
alpha_sun= 6.191297
delta_sun=-0.039772

;The half-width of the subarrays contain the alignment stars.
w=7

;If dark_file='', then use these data to determine the dark-current.
;But set dark_file = an appropriate filename in order to restore &
;use the dark-current info stored in another dataset.
dark_file=''

;If these data are to be used to measure the dark-current, then
;select the range of pixels (in the dark portion of the moon)
;where the dark-current will be measured. This is skipped if
;dark_file is not =''.
x_dark0=10
x_dark1=30
y_dark0=10
y_dark1=30

;orbit number
orbit=110

;On first entry, set lunar_limb=0 in order to run the
;proceedure that determines the lunar limb. Otherwise,
;set lunar_limb=1 to skip this step.
lunar_limb=1

;NEW!
;Zap the following pixels because charge bled away from 
;a saturated object; the value of the zapped pixels is
;replaced by the median value of the replacement_pix(els).
;Just set zap = [-1] to skip this step.
zap_mask=intarr(384,576)+1
zap_mask(171:*,*)=0		;zap right half of image
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

;NEW!
;Find the saturated pixels, and then create a mask array
;having the same dimensions as the images array. Good columns
;devoid of saturated pixels have mask=1, whereas bad columns
;containing saturated pixels have mask=0
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
  ;tvscl,im*max(exp_time)/exp_time(t)>0<50 & $
  ;tvscl,imd*max(exp_time)/exp_time(t)>0<50,Nx,0 & $
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
    im=im*0.7/exp_time(t)>0<100 & $
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
  ;tvscl,im*limb(*,*,t)*max(exp_time)/exp_time(t)>0<10 & $
endfor

;NEW!
;Subtract an average dark-current that is measured
;on a dark portion of the sky or moon. Be sure that there are
;not any cosmic-ray hits contaminating this measurement!
;plot,images(x_dark0:x_dark1,y_dark0:y_dark1,*),psym=3
dark=fltarr(Nimages)
for t=0,Nimages-1 do $
  dark(t)=median(images(x_dark0:x_dark1,y_dark0:y_dark1,t))
dark_file_save='dark.dat'
dark_save=dark
exp_time_save=exp_time
save,dark_save,exp_time_save,filename=dark_file_save
plot,exp_time,dark,psym=1,xtitle='Exposure Time    (sec)', $
  ytitle='Dark Counts',yrange=[0.5*min(dark),1.5*max(dark)], $
  charsize=1.3,title='Orbit '+strtrim(string(orbit),2)
;compare to formula from Lawrence Livermore
dark_exp=(1000.0*exp_time+54.8)*0.0189* $
  exp(0.097*(median(temp)-273.15))+150.0-8.1*15.0
oplot,exp_time,dark_exp,color=128,psym=8
if (dark_file ne '') then restore,filename=dark_file 
for t=0,Nimages-1 do begin $
  j=where(exp_time_save eq exp_time(t)) & $
  dark(t)=avg(dark_save(j)) & $
  images(0,0,t)=images(*,*,t)-dark(t) & $
  plots,exp_time(t),dark(t),psym=6 & $
endfor 

;NEW!
;flatfield the data
flatfield=readfits(flat_file)
good_flat_pix=readfits(good_file)
template=intarr(Nx,Ny)
template(x_dark0:x_dark1,y_dark0:y_dark1)=1
for t=0,Nimages-1 do begin $
  im=images(*,*,t)/flatfield & $
  images(0,0,t)=im & $
  ;tvscl,(im+3.0*template)*max(exp_time)/exp_time(t)>0<3 & $
endfor

;Check the alignment of each image using the 3 stars.
;Determine the offsets (dx0,dy0), etc., for each frame
;using star0 to a precision (but not accuracy) of 1/MAG=0.2
;pixels. Choose a long-exposure image as the reference image.
register_image,images,1,dx0,dy0,MAG=5,X0=x0-w,X1=x0+w,Y0=y0-w,Y1=y0+w
;using star1
register_image,images,1,dx1,dy1,MAG=5,X0=x1-w,X1=x1+w,Y0=y1-w,Y1=y1+w
;using star2
register_image,images,1,dx2,dy2,MAG=5,X0=x2-w,X1=x2+w,Y0=y2-w,Y1=y2+w
;plot the 3 different offsets for dx & dy, average them,
;and then round off to do full-pixel shifts.
!p.multi=[0,1,2,0,0]
plot,dx0,psym=6,yrange=[-1,5],xrange=[-1,Nimages+1],xstyle=1, $
  ystyle=1,xtitle='image index',ytitle='!7D!3x    (pixels)'
oplot,dx1,psym=4
oplot,dx2,psym=1,symsize=1.5
dx=round((dx0+dx1+dx2)/3.0)
oplot,dx,thick=3,psym=10
plot,dy0,psym=6,yrange=3.5*[-1,1],xrange=[-1,Nimages+1],xstyle=1, $
  ystyle=1,xtitle='image index',ytitle='!7D!3y    (pixels)'
oplot,dy1,psym=6
oplot,dy2,psym=1,symsize=1.5
dy=round((dy0+dy1+dy2)/3.0)
oplot,dy,thick=3,psym=10
!p.multi=0
;Shift images and zap the outer edges
for t=0,Nimages-1 do begin $
  images(0,0,t)=shift(images(*,*,t),dx(t),dy(t)) & $
  limb(0,0,t)=shift(limb(*,*,t),dx(t),dy(t)) & $
  mask(0,0,t)=shift(mask(*,*,t),dx(t),dy(t)) & $
endfor
ddx=max(abs(dx))+1
ddy=max(abs(dy))+1
zap_mask(0:ddx,*)=0
zap_mask(Nx-1-ddx:*,*)=0
zap_mask(*,0:ddy)=0
zap_mask(*,Ny-1-ddy:*)=0

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
x_c=191.0	;Coordinates for the optical center.
y_c=286.0
x_star=[ x0, x1, x2, x3, x4, x5, x6, x7, x8]
y_star=[ y0, y1, y2, y3, y4, y5, y6, y7, y8] 
RA_star = [ alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8 ]
DEC_star= [ delta0, delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8 ]
.run plate_scale
plate_scale, x_c, y_c, x_star, y_star, RA_star, DEC_star, $
  plate_scl, errors
print,'plate-scale = ',plate_scl,' radians/pixel'
print,min(errors)

;get the twist-angle & (RA_c,DEC_c) coordinates of center.
twist_angle=avg(twist)
RA_c=avg(RA)-2.0*!pi
DEC_c=avg(DEC)
.run RA_DEC
RA_DEC, plate_scl, x_c, y_c, x_star, y_star, RA_star, $
  DEC_star, twist_angle, RA_c, DEC_c, Nx, Ny, RA_array, $
  DEC_array,ta=ta, rac=rac, decc=decc
r2d=180./!pi
print,'twist-angle = ',twist_angle*r2d,' degrees'
!p.multi=[0,1,3,0]
plot,ta*r2d
plot,rac*r2d
plot,decc*r2d
!p.multi=0
err0=sqrt((RA_array(x0,y0)-alpha0)^2+(DEC_array(x0,y0)-delta0)^2)/plate_scl
err1=sqrt((RA_array(x1,y1)-alpha1)^2+(DEC_array(x1,y1)-delta1)^2)/plate_scl
err2=sqrt((RA_array(x2,y2)-alpha2)^2+(DEC_array(x2,y2)-delta2)^2)/plate_scl
err3=sqrt((RA_array(x3,y3)-alpha3)^2+(DEC_array(x3,y3)-delta3)^2)/plate_scl
err4=sqrt((RA_array(x4,y4)-alpha4)^2+(DEC_array(x4,y4)-delta4)^2)/plate_scl
err5=sqrt((RA_array(x5,y5)-alpha5)^2+(DEC_array(x5,y5)-delta5)^2)/plate_scl
err6=sqrt((RA_array(x6,y6)-alpha6)^2+(DEC_array(x6,y6)-delta6)^2)/plate_scl
err7=sqrt((RA_array(x7,y7)-alpha7)^2+(DEC_array(x7,y7)-delta7)^2)/plate_scl
err8=sqrt((RA_array(x8,y8)-alpha8)^2+(DEC_array(x8,y8)-delta8)^2)/plate_scl
print,'Discrepancies between the observed & expected positions'
print,'of the 3 registration stars are LESS than
print,err0,err1,err2,err3,err4,err5,err6,err7,err8,' pixels.'

;get ecliptic coordinates, where phi is the longitude
;relative to the Sun, and theta is the ecliptic latitude
obl=0.4090928			;Earth's obliquity in radians	
x=cos(DEC_array)*cos(RA_array)
y=cos(obl)*cos(DEC_array)*sin(RA_array)+sin(obl)*sin(DEC_array)
z=-sin(obl)*cos(DEC_array)*sin(RA_array)+cos(obl)*sin(DEC_array)
phi=atan(y,x)
theta=atan(z,sqrt(x^2+y^2))
phi_sun=acos(cos(alpha_sun)*cos(delta_sun))
if (alpha_sun gt !pi) then phi_sun=-phi_sun
phi=phi-phi_sun

;zap any bad pixels
if (zap(0) eq -1) then zap_mask=intarr(Nx,Ny)+1
good_pix=good_pix*zap_mask
counts=counts*zap_mask
good=where(good_pix eq 1)
avg_image=fltarr(Nx,Ny)
avg_image(good)=counts(good)/total_exp_time(good)

;make pictures
im=alog(avg_image>1<1000)
im=im-min(im)
im=byte(255*im/max(im))
d=sqrt(phi^2+theta^2)
sun_pix=where(d eq min(d))
im(sun_pix)=255
window,xs=Nx,ys=Ny
tvscl,im
contour,smooth(im,3),nlevels=7,xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels'
print,'Observation date = ',date(Nimages/2)

;NEW!
;save results
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,smooth(im,3),nlevels=7,xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels'
file='orbit'+strtrim(string(orbit)+'_sunny.ps',2)
output_plot,file
file='orbit'+strtrim(string(orbit)+'_sunny.dat',2)
save,filename=file,plate_scl,twist_angle,temp,gain_id, $
  offset_id,phi,theta,counts,total_exp_time,good_pix

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
  yrange=y3+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x3+w*[-1,1],y3+[0,0],color=128
plots,x3+[0,0],y3+w*[-1,1],color=128

contour,alog(avg_image>1<120),xrange=x4+w*[-1,1], $
  yrange=y4+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x4+w*[-1,1],y4+[0,0],color=128
plots,x4+[0,0],y4+w*[-1,1],color=128

contour,avg_image>3<80,xrange=x5+w*[-1,1], $
  yrange=y5+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x5+w*[-1,1],y5+[0,0],color=128
plots,x5+[0,0],y5+w*[-1,1],color=128

contour,alog(avg_image>1<120),xrange=x6+w*[-1,1], $
  yrange=y6+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x6+w*[-1,1],y6+[0,0],color=128
plots,x6+[0,0],y6+w*[-1,1],color=128

contour,alog(avg_image>3<100),xrange=x7+w*[-1,1], $
  yrange=y7+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x7+w*[-1,1],y7+[0,0],color=128
plots,x7+[0,0],y7+w*[-1,1],color=128

contour,alog(avg_image>1<120),xrange=x8+w*[-1,1], $
  yrange=y8+w*[-1,1],nlevels=15,xstyle=1,ystyle=1
plots,x8+w*[-1,1],y8+[0,0],color=128
plots,x8+[0,0],y8+w*[-1,1],color=128

for t=0,Nimages-1 do begin $
  tvscl,images(*,*,t)*max(exp_time)/exp_time(t)>0<5 & $
  sub0=images(x0-w>0:x0+w<Nx-1,y0-w>0:y0+w<Ny-1,t) & $
  tvscl,sub0*max(exp_time)/exp_time(t)>0<5,Nx,y0-w & $
  sub1=images(x1-w>0:x1+w<Nx-1,y1-w>0:y1+w<Ny-1,t) & $
  tvscl,sub1*max(exp_time)/exp_time(t)>0<5,Nx,y1-w & $
  sub2=images(x2-w>0:x2+w<Nx-1,y2-w>0:y2+w<Ny-1,t) & $
  tvscl,sub2*max(exp_time)/exp_time(t)>0<5,Nx+2*w,y2-w & $
  wait,0.5 & $
endfor




