;orbit253.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Set input parameters.				;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The leading part of the input fits-filenames.
lead_name='/net/lpiip2/u2/hahn/clementine/new253/conversions/lba'

;The range of image-numbers (e.g., the middle part
;of the fits filenames).
start=5176
stop =5226

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ 5197, 5198, 5199 ]

;The trailing part of the fits filenames.
trail_name='z.253.fits'

;flatfield & good-pixel filenames
flat_file='/net/lpiip2/users/hahn/clementine/zodiacal_light/flatfield/flatfield.fits'
good_file='/net/lpiip2/users/hahn/clementine/zodiacal_light/flatfield/good_pix.fits'

;The (x,y) coordinates of 3+ stars that are used to co-register
;these images, and their (alpha,delta)=(RA,DEC) in radians,
;which is used to orient the images.

x0	=  278			;Venus - 0       
y0	=  55
alpha0	=  0.741803
delta0	=  0.284153

x1	=  126  		;Hamal - 1      
y1	=  118
alpha1	=  0.554899
delta1	=  0.409496

x2	=  137			;Sheretan       
y2	=  177
alpha2	=  0.500182             ;0.500212
delta2	=  0.363169


x3      =  110.25               ;Ari Lambda 9
y3      =  142.25
alpha3  =  0.514559
delta3  =  0.411829

;x3	=  339.3		;Ari Gamma 2-5   
;y3	=  444.2
;alpha3	=  0.495370
;delta3	=  0.336738

x4	=  37.5  		;Tri - alpha 2  
y4	=  113.5
alpha4	=  0.493413
delta4	=  0.516248

x5      =  107.5                ;Ari 7
y5      =  147.5
alpha5  =  0.505491
delta5  =  0.411500


;x5	=  120.75		;Psc 107
;y5	=  205.5
;alpha5	=  0.178889
;delta5	=  0.3537521

x6	=  77   		;Ari 35 -> #6
y6	=  40
alpha6	=  0.285277
delta6	=  0.483582

;x7	=  105
;y7	=  254
;alpha7	=  0.151311
;delta7	=  0.335809


x7	=  94	        	;Ari 41 -> #7
y7	=  23		
alpha7	=  0.741692             ;0.296677		
delta7	=  0.475785             ;0.475786

;Sun's RA and DEC in radians
alpha_sun=0.383151  ;0.18
delta_sun=0.160666  ;0.079425

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
x_dark1=50
y_dark0=10
y_dark1=50

;orbit number
orbit=253

;On first entry, set lunar_limb=0 in order to run the
;proceedure that determines the lunar limb. Otherwise,
;set lunar_limb=1 to skip this step.
lunar_limb=1

;Zap the following pixels because charge bled away from 
;a saturated object. Just set zap = [-1] to skip this step.
zap_mask=intarr(384,576)+1
zap_mask(382:*,*)=0		;zap right edge
zap_mask(*,480:486)=0		;zap rows
;zap_mask(259:263,*)=0		;  & columns near Venus
;radial,rad,261,483,384,576	;DONT zap Venus, at least
;j=where(rad le 26.)		;not for the DPS poster...
;zap_mask(j)=0
;replacement_pix=where((rad le 26.) and (rad ge 25.0))replacement_pix=where
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
  ;tvscl,im*max(exp_time)/exp_time(t) & $

endfor

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
j=where((mask eq 1) and (images gt 160)) & plot,images(j),psym=3

;stop

;subtract dark current & destreak images
coeff=[ 21.2, 13.7 ]		;anticipated dark-current
dark=coeff(0)+coeff(1)*exp_time
images_dd=fltarr(Nx,Ny,Nimages)
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  im=im-dark(t) & $
  imd=destreak(im,exp_time(t)) & $
;  tvscl,im*max(exp_time)/exp_time(t)>0<200 & $
;  tvscl,imd*mask(*,*,t)*max(exp_time)/exp_time(t)>0<200,Nx,0 &$
  images_dd(0,0,t)=imd & $
endfor
plot,exp_time,exp_time,yrange=200*[-1,1],nodata=1,ystyle=1, $
  xtitle='exposure time    (sec)',charsize=1.3, $
  ytitle='residual dark current    (counts/sec)'
dark_adj=fltarr(Nimages)

sub=images_dd(x_dark0:x_dark1,y_dark0:y_dark1,*)
sub_mask=mask(x_dark0:x_dark1,y_dark0:y_dark1,*)
sub_exp_time=sub*0

for t=0,Nimages-1 do begin $
  sub_exp_time(*,*,t) = exp_time(t) & $
  pix=where((sub_mask eq 1) and (sub_exp_time eq exp_time(t)) $ 
    and (sub_exp_time lt 0.35)) & $
  if(pix(0) ne -1) then begin & $
   oplot, sub_exp_time(pix), sub(pix)/exp_time(t),psym=3 & $
   dark_adj(t)=median(sub(pix)) & $
  endif & $
  plots,exp_time(t), dark_adj(t)/exp_time(t),psym=6 & $
endfor
data=where(exp_time lt 0.35)
coeff_adj=poly_fit(exp_time(data),dark_adj(data),1)
s=sort(exp_time)
oplot,exp_time(s),coeff_adj(0)/exp_time(s)+coeff_adj(1),thick=3
oplot,exp_time,exp_time*0,color=128
print,transpose(coeff)
print,transpose(coeff+coeff_adj)

images=images_dd
plot,images(x_dark0:x_dark1,y_dark0:y_dark1,*),psym=3


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
    limb_click,im,x,y,r & $
    xc(t)=x & $
    yc(t)=y & $
    R_moon(t)=r & $
    radial,rad,x,y,Nx,Ny & $
    pix=where(rad le r) & $
    ;tvscl,im*limb(*,*,t)*max(exp_time)/exp_time(t)>0<100 & $
    lmb=bytarr(Nx,Ny)+1 & $
    lmb(pix)=0 & $
    ;tvscl,im*lmb & $
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


;flatfield the data
flatfield=readfits(flat_file)
good_flat_pix=readfits(good_file)
template=intarr(Nx,Ny)
template(x_dark0:x_dark1,y_dark0:y_dark1)=1
for t=0,Nimages-1 do begin $
  im=images(*,*,t)/flatfield & $
  images(0,0,t)=im & $
  ;tvscl,(im+3.0*template)*max(exp_time)/exp_time(t)>0<100 & $
endfor

;Check the alignment of each image using the 3 stars.
;Determine the offsets (dx0,dy0), etc., for each frame
;using star0 to a precision (but not accuracy) of 1/MAG=0.2
;pixels. Choose a long-exposure image as the reference image.
register_image,images,18,dx0,dy0,MAG=5,X0=x0-w,X1=x0+w,Y0=y0-w,Y1=y0+w
;using star1
register_image,images,18,dx1,dy1,MAG=5,X0=x1-w,X1=x1+w,Y0=y1-w,Y1=y1+w
;using star2
register_image,images,18,dx2,dy2,MAG=5,X0=x2-w,X1=x2+w,Y0=y2-w,Y1=y2+w
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

;check linearity of flux in box
flux=fltarr(Nimages)
template=intarr(Nx,Ny)
template(150:180,360:390)=1
pix=where(template eq 1)
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  flux(t)=avg(im(pix)) & $
  im=im*(1.0+0.3*template) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<100 & $
endfor
plot,exp_time,flux,psym=1,xtitle='exposure time    (sec)', $
  ytitle='average flux in a box',charsize=1.3
coef=poly_fit(exp_time,flux,1,yfit)
s=sort(exp_time)
oplot,exp_time(s),yfit(s)
t=[0.0,exp_time(s)]
y=coef(0)+coef(1)*t
oplot,t,y,linestyle=1

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
x_star=[ x0, x1, x2, x3, x4, x5 ] 	
y_star=[ y0, y1, y2, y3, y4, y5 ] 
RA_star = [ alpha0, alpha1, alpha2, alpha3, alpha4, alpha5 ]  
DEC_star= [ delta0, delta1, delta2, delta3, delta4, delta5 ] 
.run plate_scale
plate_scale, x_c, y_c, x_star, y_star, RA_star, DEC_star, $
  plate_scl, errors
print,'plate-scale = ',plate_scl,' radians/pixel'
print,min(errors)

;get the twist-angle & (RA_c,DEC_c) coordinates of center.
twist_angle=avg(twist)
RA_c=avg(RA)
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

print,'Discrepancies between the observed & expected positions'
print,'of the registration stars are LESS than
print,err0,err1,err2,err3,err4,err5,' pixels.'

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
good_pix=good_pix*good_flat_pix*zap_mask
counts=counts*good_flat_pix*zap_mask
good=where(good_pix eq 1)
avg_image=fltarr(Nx,Ny)
avg_image(good)=counts(good)/total_exp_time(good)

stop

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
print,'Observation date = ',date(Nimages/2)

;save results
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,im,nlevels=10,xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels'
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
x0=[x0(0:57),x0(64:*)]
x1=indgen(100)
for y=0,Ny-1 do begin $
  row=sub(*,(y-1)>0:(y+1)<(Ny-1)) & $
  for x=0,99 do row(x,0)=avg(row(x,*)) & $
  row=row(*,0) & $
  ;plot,x1,row,title=string(y) & $
  sub(0,y)=spline(x0,row(x0),x1) & $
  ;oplot,x1,sub(*,y),thick=3,color=128 & $
endfor
pretty(258,0)=sub(58:63,*)>0<255
tv,pretty
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=4.5, $
  /portrait,bits_per_pixel=8
tv,pretty
output_plot,'pretty_picture.ps'
set_plot,'ps'
device,xsize=16,ysize=23,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,pretty,nlevels=7,xstyle=1,ystyle=1, $
  charsize=1.3,xtitle='pixels',ytitle='pixels'
output_plot,'pretty_contour.ps'

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
