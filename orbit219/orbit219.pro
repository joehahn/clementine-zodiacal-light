;orbit219.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;							;
;	Set input parameters.				;
;							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The leading part of the input fits-filenames.
lead_name='/u2/hahn/clementine/new219/conversions/lba'

;The range of image-numbers (e.g., the middle part
;of the fits filenames). Note: images 6220-6221 are
;problematic, and are ignored.
start=6183
stop =6219

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ 6187, 6192, 6217 ]

;The trailing part of the fits filenames.
trail_name='z.219.fits'

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
x0=25				;star0
y0=336
x1=71				;star1
y1=440
x2=359				;star2
y2=434

;The half-width of the subarrays contain the alignment stars.
w=7

;If these data are to be used to measure the dark-current, then
;select the range of pixels (in the dark portion of the moon)
;where the dark-current will be measured. This is skipped if
;dark_file=''.
x_dark0=10
x_dark1=30
y_dark0=546
y_dark1=566

;orbit number
orbit=219

;On first entry, set lunar_limb=0 in order to run the
;proceedure that determines the lunar limb. Otherwise,
;set lunar_limb=1 to skip this step.
lunar_limb=1

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
  ;tvscl,im*max(exp_time)/exp_time(t)>0<200 & $
  ;tvscl,imd*max(exp_time)*mask(*,*,t)/exp_time(t)>0<200,Nx,0 & $
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
  tvscl,im*limb(*,*,t)*max(exp_time)/exp_time(t)>0<100 & $
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

;Subtract an average dark-current that is measured
;on a dark portion of the sky or moon. Be sure that there are
;not any cosmic-ray hits contaminating this measurement!
;plot,images(x_dark0:x_dark1,y_dark0:y_dark1,*),psym=3
dark=fltarr(Nimages)
for t=0,Nimages-1 do $
  dark(t)=median(images(x_dark0:x_dark1,y_dark0:y_dark1,t))
dark_file_save='dark.dat'
exp_time_save=exp_time
dark_save=dark
save,dark_save,exp_time_save,filename=dark_file_save
plot,exp_time,dark,psym=1,xtitle='Exposure Time    (sec)', $
  ytitle='Dark Counts',yrange=[0.5*min(dark),max(dark)], $
  charsize=1.3,title='Orbit '+strtrim(string(orbit),2)
;compare to formula from Lawrence Livermore
dark_exp=(1000.0*exp_time+54.8)*0.0189* $
  exp(0.097*(median(temp)-273.15))+150.0-8.1*15.0
oplot,exp_time,dark_exp,color=128
if (dark_file ne '') then restore,filename=dark_file 
dark_save=dark
for t=0,Nimages-1 do begin $
  j=where(exp_time_save eq exp_time(t)) & $
  dark(t)=avg(dark_save(j)) & $
  images(0,0,t)=images(*,*,t)-dark(t) & $
  plots,exp_time(t),dark(t),psym=6 & $
endfor 

;Coadd all good pixels (eg., where limb & mask both = 1)
avg_image=fltarr(Nx,Ny)
limb_mask=limb*mask
good_pix=bytarr(Nx,Ny)
for x=0,Nx-1 do begin $
  for y=0,Ny-1 do begin $
    j=where(limb_mask(x,y,*) eq 1) & $
    if (j(0) ne -1) then begin $
      avg_image(x,y)=total(images(x,y,j))/total(exp_time(j)) & $
      good_pix(x,y)=1 & $
    endif & $
  endfor & $
endfor

;check linearity of flux in box
flux=fltarr(Nimages)
template=intarr(Nx,Ny)
template(150:180,360:390)=1
pix=where(template eq 1)
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  flux(t)=avg(im(pix)) & $
  im=im*limb_mask(*,*,t)*(1.0+0.3*template) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<100 & $
endfor
plot,exp_time,flux,psym=1,xtitle='exposure time    (sec)', $
  ytitle='average flux in a box',charsize=1.3
coef=poly_fit(exp_time,flux,1,yfit)
s=sort(exp_time)
oplot,exp_time(s),yfit(s)

;make pictures
im=alog(avg_image>3<3000)
im=im-min(im)
im=byte(255*im/max(im))
window,xs=Nx,ys=Ny
tvscl,im
contour,im,nlevels=10,xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels'

;save results
file='orbit'+strtrim(string(orbit)+'.gif',2)
write_gif,file,im
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,im,nlevels=10,xstyle=1,ystyle=1, $
  title='Orbit '+strtrim(string(orbit),2),charsize=1.3, $
  xtitle='pixels',ytitle='pixels'
file='orbit'+strtrim(string(orbit)+'.ps',2)
output_plot,file
file='orbit'+strtrim(string(orbit)+'.dat',2)
save,filename=file,avg_image,good_pix,RA_array,DEC_array, $
  plate_scl,twist_angle,dark,exp_time,temp,gain_id,offset_id
file='orbit'+strtrim(string(orbit)+'.fits',2)
writefits,file,avg_image

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
contour,alog(avg_image>1),xrange=x7+w*[-1,1], $
  yrange=y7+w*[-1,1],nlevels=25,xstyle=1,ystyle=1
plots,x7+w*[-1,1],y7+[0,0],color=128
plots,x7+[0,0],y7+w*[-1,1],color=128
for t=0,Nimages-1 do begin $
  tvscl,images(*,*,t)*max(exp_time)/exp_time(t)>0<50 & $
  sub0=images(x0-w>0:x0+w<Nx-1,y0-w>0:y0+w<Ny-1,t) & $
  tvscl,sub0*max(exp_time)/exp_time(t)>0<50,Nx,y0-w & $
  sub1=images(x1-w>0:x1+w<Nx-1,y1-w>0:y1+w<Ny-1,t) & $
  tvscl,sub1*max(exp_time)/exp_time(t)>0<50,Nx,y1-w & $
  sub2=images(x2-w>0:x2+w<Nx-1,y2-w>0:y2+w<Ny-1,t) & $
  tvscl,sub2*max(exp_time)/exp_time(t)>0<50,Nx+2*w,y2-w & $
endfor

;unsharp-mask
usm=avg_image-smooth(avg_image,5)
j=where(usm ge 50.0)
m=bytarr(Nx,Ny)
m(j)=1
tvscl,m
set_plot,'ps'
device,xsize=16,ysize=21,xoffset=3,yoffset=3.5, $
  /portrait,bits_per_pixel=8
contour,m,nlevels=1,xstyle=1,ystyle=1
output_plot




