;orbit193.pro
;
;A non-flatfielded version of this data is used to check
;the flatfield. See ../flatfield/flatfield.pro
;
;Rows containing saturated pixels are NOT discarded.
;
;if fix_venus=1, then venus is replaced in each image with a
;local average proir to destreaking, then discarded. Set
;fix_venus=0 to generate the pretty picture w/Venus.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Set input parameters.				;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The leading part of the input fits-filenames.
lead_name='/u2/hahn/clementine/new193/conversions/lba'

;The range of image-numbers (e.g., the middle part
;of the fits filenames).
start=5866
stop =5905

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ -1 ]

;The trailing part of the fits filenames.
trail_name='z.193.fits'

;flatfield & good-pixel filenames
flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield_short_exp.fits'
good_file='/users/hahn/clementine/zodiacal_light/flatfield/good_pix_short_exp.fits'

;The (x,y) coordinates of 3+ stars that are used to co-register
;these images, and their (alpha,delta)=(RA,DEC) in radians,
;which is used to orient the images.
x0	=  61.0			;Theta Cetus
y0	=  309.0
alpha0	=  0.366621
delta0	= -0.142825
x1	=  55.5			;Eta Cetus
y1	=  252.2
alpha1	=  0.299280
delta1	= -0.177714
x2	=  58.1			;Mira = Omni Cetus
y2	=  500.1
alpha2	=  0.608014
delta2	= -0.051970
x3	=  339.3		;eta Pisces
y3	=  444.2
alpha3	= 0.399171
delta3	= 0.267835
x4	=  149.0		;Kaitain = alpha Pisces
y4	=  477.7
alpha4	= 0.532530
delta4	= 0.048237
x5	=  132.0		;iota Cetus
y5	=  110.5
alpha5	= 0.084770
delta5	=-0.154006
x6	=  187.2		;30 Pisces
y6	=  71.5
alpha6	= 0.008552
delta6	=-0.104965
x7	=  377.7		;beta Ares
y7	=  532.8
alpha7	= 0.500212
delta7	= 0.363169

;Sun's RA and DEC in radians
alpha_sun=0.184722
delta_sun=0.079425

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
x_dark0=330
x_dark1=370
y_dark0=10
y_dark1=50

;The above aperture might not be completely dark (due to a faint
;zodiacal light or earthshine contribution), so the dark-current
;inferred in this aperture needs to be reduced by the following
;flux in counts/second. This residual can be determined using
;mosaic.pro
residual_flux=1.3

;orbit number
orbit=193

;On first entry, set lunar_limb=0 in order to run the
;proceedure that determines the lunar limb. Otherwise,
;set lunar_limb=1 to skip this step.
lunar_limb=1

;Zap the following pixels because charge bled away from 
;a saturated object. Just set zap = [-1] to skip this step.
zap_mask=intarr(384,576)+1
zap_mask(258:264,*)=0		;zap columns near Venus
zap=where(zap_mask eq 0)

;set fix_venus=1 to replace pixels polluted by Venus with a
;local average. Set=0 to skip this when making the pretty picture
fix_venus=0

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
instrument=strarr(Nimages)		;instrument name
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
  j=where(strmid(header,0,23) eq 'COMMENT INSTRUMENT_NAME') & $
  instrument(t)=strmid(header(j(0)),28,23) & $

  ;display input image
  ;im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<100 & $

endfor

;store raw images
images_raw=images

;Find the saturated pixels, and then create a mask array
;having the same dimensions as the images array. Good columns
;devoid of saturated pixels have mask=1, whereas bad columns
;containing saturated pixels have mask=0
mask=bytarr(Nx,Ny,Nimages)+1
for t=0,Nimages-1 do begin $
  for y=0,Ny-1 do begin $
    j=where(images(*,y,t) ge 255) & $
    ;comment out this line when making the pretty picture.
    ;if (j(0) ne -1) then mask(*,y,t)=0 & $
  endfor & $
  for x=0,Nx-1 do begin $
    j=where(images(x,*,t) ge 255) & $
    if (j(0) ne -1) then mask(x,*,t)=0 & $
  endfor & $
endfor
;j=where((mask eq 1) and (images gt 160)) & plot,images(j),psym=3

;replace Venus-polluted pixels with a local average prior to
;destreaking the images.
if (fix_venus eq 1) then begin $
  xV=261 & $				;Venus x-coordinate
  yV=483 & $				;Venus y-coordinate
  rV=18	 & $				;Venus radius
  sub1=images(xV-rV:xV+rV,yV-3*rV:yV-rV,*) & $
  sub2=images(xV-rV:xV+rV,yV+rV:yV+3*rV,*) & $
  sub_avg=0.5*(sub1+sub2) & $
  radial,rad,rV,rV,2*rV+1,2*rV+1 & $
  pix=where(rad le rV) & $
  for t=0,Nimages-1 do begin $
    ;im=images(*,*,t) & $
    ;im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
    ;tvscl,im*max(exp_time)/exp_time(t)>0<60 & $
    sub_image=images(xV-rV:xV+rV,yV-rV:yV+rV,t) & $
    sub_sub_avg=sub_avg(*,*,t) & $
    sub_image(pix)=sub_sub_avg(pix) & $
    images(xV-rV,yV-rV,t)=sub_image & $
    ;im=images(*,*,t) & $
    ;im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
    ;tvscl,im*max(exp_time)/exp_time(t)>0<60,Nx,0 & $
  endfor & $
  radial,rad,xV,yV,Nx,Ny & $
  pix=where(rad le rV) & $
  zap_mask(pix)=0 & $
  zap=where(zap_mask eq 0) & $
endif

;iteratively determine the dark current for each image
plot,images(x_dark0:x_dark1,y_dark0:y_dark1,*),psym=3
dark=fltarr(Nimages)
sub_images=images(x_dark0:x_dark1,*,*)
flux=fltarr(Nimages)
flatfield=readfits(flat_file)
flat_sub=flatfield(x_dark0:x_dark1,*)
for iteration=0,15 do begin $
  for t=0,Nimages-1 do begin $
    im=sub_images(*,*,t)-dark(t) & $
    im=destreak(im,exp_time(t))/flat_sub & $
    flux(t)=avg(im(*,y_dark0:y_dark1)) & $
    dark(t)=dark(t)+flux(t)-residual_flux*exp_time(t) & $
    flux(t)=flux(t)/exp_time(t) & $
  endfor & $
  plot,exp_time,flux,psym=8,title=string(iteration), $
    xtitle='exposure time    (sec)',charsize=1.3, $
    ytitle='aperture flux    (counts/sec)' & $
  oplot,exp_time,exp_time*0+residual_flux,color=128 & $
endfor

;dark-subtract & destreak the images
images_dd=fltarr(Nx,Ny,Nimages)
for t=0,Nimages-1 do begin $
  im=images(*,*,t)-dark(t) & $
  im=destreak(im,exp_time(t)) & $
  images_dd(0,0,t)=im & $
endfor
images=images_dd

;flatfield the data
flatfield=readfits(flat_file)
;flatfield=fltarr(Nx,Ny)+1.0	;DON'T flatfield for flatfield_short_exp.pro
good_flat_pix=readfits(good_file)
template=intarr(Nx,Ny)
template(x_dark0:x_dark1,y_dark0:y_dark1)=1
flux=fltarr(Nimages)
for t=0,Nimages-1 do begin $
  im=images(*,*,t)/flatfield & $
  images(0,0,t)=im & $
  ;tvscl,(im+3.0*template)*max(exp_time)/exp_time(t)>0<30 & $
  flux(t)=avg(images(x_dark0:x_dark1,y_dark0:y_dark1,t)) $
    /exp_time(t) & $
endfor
;check flux in aperture
plot,exp_time,flux,psym=8,xtitle='exposure time    (sec)', $
    charsize=1.3,ytitle='aperture flux    (counts/sec)'
oplot,exp_time,exp_time*0+residual_flux,color=128

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
x_star=[ x0, x1, x2, x3, x4, x5, x6, x7 ]	
y_star=[ y0, y1, y2, y3, y4, y5, y6, y7 ]
RA_star = [ alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7 ]
DEC_star= [ delta0, delta1, delta2, delta3, delta4, delta5, delta6, delta7 ]
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
x=cos(obl)*sin(alpha_sun)*cos(delta_sun)+sin(obl)*sin(delta_sun)
z=cos(alpha_sun)*cos(delta_sun)
phi_sun=atan(x,z)
phi=phi-phi_sun

;print ecliptic longitude
if (phi_sun lt 0.0) then phi_sun=phi_sun+2.0*!pi
print,'ecliptic longitude = ',phi_sun*180./!pi,' degrees'
print,'Observation date = ',date(Nimages/2)

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
contour,im,nlevels=10,xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels'

;save results
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,alog(smooth(avg_image,4)>1),nlevels=30,xstyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels',ystyle=1
file='orbit'+strtrim(string(orbit)+'.ps',2)
output_plot,file
file='orbit'+strtrim(string(orbit)+'.dat',2)
save,filename=file,plate_scl,twist_angle,temp,gain_id, $
  offset_id,phi,theta,counts,total_exp_time,good_pix,avg_image

;make pretty picture
pretty=alog(avg_image>5<1000)
pretty=pretty-min(pretty)
pretty=255.0*pretty/max(pretty)
img=images(*,*,35)
radial,rad,xc(35),yc(35),Nx,Ny
j=where(rad le R_moon(35)-5)
pretty(j)=img(j)*2.0>0<255
sub=pretty(200:299,*)
x0=indgen(100)
x0=[x0(0:57),x0(65:*)]
x1=indgen(100)
for y=0,Ny-1 do begin $
  row=sub(*,(y-1)>0:(y+1)<(Ny-1)) & $
  for x=0,99 do row(x,0)=avg(row(x,*)) & $
  row=row(*,0) & $
  ;plot,x1,row,title=string(y) & $
  sub(0,y)=spline(x0,row(x0),x1) & $
  ;oplot,x1,sub(*,y),thick=3,color=128 & $
endfor
pretty(259,0)=sub(59:63,*)>0<255
left=pretty(254:258,*)
right=pretty(264:268,*)
a=0.5*(left+right)
pretty(259,325)=a(*,325:*)
smth=smooth(pretty,2)
pretty(259,325)=smth(259:263,325:*)
radial,rad_sub,61,483,100,Ny
radial,rad_pretty,261,483,Nx,Ny
j_sub   =where(rad_sub    le 10.0)
j_pretty=where(rad_pretty le 10.0)
pretty(j_pretty)=sub(j_sub)<255.0
tv,pretty

;Make the pretty picture only when Venus isn't fixed.
if (fix_Venus eq 0) then begin $
  save,filename='pretty_picture.dat',pretty & $
  set_plot,'ps' & $
  device,xsize=16,ysize=21,xoffset=3,yoffset=4.5, $
    /portrait,bits_per_pixel=8 & $
  tv,pretty & $
  output_plot,'pretty_picture.ps' & $
endif

stop

;show raw & reduced images for Icarus paper
pflag=0
t=26
mn=5.0
mx=55.0
raw=(images_raw(*,*,t)-dark(t))/exp_time(t)
raw=rotate(raw,1)>mn<mx
raw=255.0*raw/mx
ABC=fltarr(Ny,3*Nx)
ABC(0,2*Nx)=raw
log_avg_image=alog(rotate(avg_image,1)>5<1000)
log_avg_image=log_avg_image-min(log_avg_image)
log_avg_image=255.0*log_avg_image/max(log_avg_image)
ABC(0,Nx)=log_avg_image
restore,'pretty_picture.dat'
ABC(0,0)=rotate(pretty,1)
if (pflag eq 0) then window,xs=Ny,ys=3*Nx
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=11.5,ysize=23.0,xoffset=5.0,yoffset=2.5, $
    /portrait,bits_per_pixel=8 & $
endif
tv,ABC
xyouts,-75,395,'A',charsize=2,charthick=4,color=255
xyouts,-75,178,'B',charsize=2,charthick=4,color=255
xyouts,-75,-38,'C',charsize=2,charthick=4,color=255
if (pflag eq 1) then output_plot,'reduction_ABC.ps'

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
contour,alog(avg_image>1),xrange=x4+w*[-1,1], $
  yrange=y4+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x4+w*[-1,1],y4+[0,0],color=128
plots,x4+[0,0],y4+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x5+w*[-1,1], $
  yrange=y5+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x5+w*[-1,1],y5+[0,0],color=128
plots,x5+[0,0],y5+w*[-1,1],color=128
contour,alog(avg_image>1),xrange=x6+w*[-1,1], $
  yrange=y6+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x6+w*[-1,1],y6+[0,0],color=128
plots,x6+[0,0],y6+w*[-1,1],color=128
;contour,alog(avg_image>1),xrange=x7+w*[-1,1], $
;  yrange=y7+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
;plots,x7+w*[-1,1],y7+[0,0],color=128
;plots,x7+[0,0],y7+w*[-1,1],color=128
for t=0,Nimages-1 do begin $
  tvscl,images(*,*,t)*max(exp_time)/exp_time(t)>0<50 & $
  sub0=images(x0-w>0:x0+w<Nx-1,y0-w>0:y0+w<Ny-1,t) & $
  tvscl,sub0*max(exp_time)/exp_time(t)>0<50,Nx,y0-w & $
  sub1=images(x1-w>0:x1+w<Nx-1,y1-w>0:y1+w<Ny-1,t) & $
  tvscl,sub1*max(exp_time)/exp_time(t)>0<50,Nx,y1-w & $
  sub2=images(x2-w>0:x2+w<Nx-1,y2-w>0:y2+w<Ny-1,t) & $
  tvscl,sub2*max(exp_time)/exp_time(t)>0<50,Nx+2*w,y2-w & $
endfor

;check rows & columns near Venus
y=487 & plot,avg_image(*,488)<40,thick=3 & oplot,avg_image(*,y),psym=1 & print,max(avg_image(*,y))
x=257 & plot,avg_image(257,*),xrange=[400,550],thick=3 & oplot,avg_image(x,*)
