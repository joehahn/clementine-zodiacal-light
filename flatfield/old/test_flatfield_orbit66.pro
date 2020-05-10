;test_flatfield_orbit66.pro
;
; Check the quality of the flatfield by monitoring
; the flux of various stars as the camera pans the
; sky in orbit66.
;
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
stop =84

;but skip these image-numbers, since these files don't exist.
;Set skip = [ -1 ] to avoid skipping any image-numbers.
skip=[ -1 ]

;The trailing part of the fits filenames.
trail_name='v.066.fits'

;flatfield & good-pixel filenames
flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield.fits'
good_file='/users/hahn/clementine/zodiacal_light/flatfield/good_pix.fits'

;If these data are to be used to measure the dark-current, then
;select the range of pixels (in the dark portion of the moon)
;where the dark-current will be measured. This is skipped if
;dark_file=''.
x_dark0=10
x_dark1=374
y_dark0=556
y_dark1=566

;On first entry, set star_clk=0 in order to run the
;proceedure that determines the target star;s
;(xc,yc) coordinates. Otherwise, set star_clk=1
;to skip this step.
star_clk=1
star_files=['orbit66_lower_left.dat','orbit66_upper_left.dat', $
            'orbit66_upper_right.dat', 'orbit66_left.dat' ]

;The plate-scale in radians per pixel
plate_scale=0.00137403

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
  ;tvscl,im*max(exp_time)/exp_time(t)>0<20 & $

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

;subtract dark current & destreak images
coeff=[ 36.9, 15.2 ]		;anticipated dark-current
dark=coeff(0)+coeff(1)*exp_time
for t=0,Nimages-1 do begin $
  im=images(*,*,t) & $
  im=im-dark(t) & $
  imd=destreak(im,exp_time(t)) & $
  ;tvscl,im*max(exp_time)/exp_time(t)>0<10 & $
  ;tvscl,imd*mask(*,*,t)*max(exp_time)/exp_time(t)>0<10,Nx,0 &$
  images(0,0,t)=imd & $
endfor

;flatfield the data
flatfield=readfits(flat_file)
good_pix=readfits(good_file)
for t=0,Nimages-1 do begin $
  images(0,0,t)=images(*,*,t)/flatfield & $
  ;tvscl,images(*,*,t)*mask(*,*,t)*max(exp_time)/ $
  ;  exp_time(t)>0<20 & $
endfor

;One first entry (e.g., if star_clk=0), click the mouse on the
;target star to get is (xc,yc) coordinates from each image.
star=3				;choose star!
window,xs=1.5*Nx,ys=Ny
if (star_clk eq 0) then begin $
  xc=fltarr(Nimages) & $
  yc=fltarr(Nimages) & $
  for t=0,Nimages-1 do begin $
    im=images(*,*,t) & $
    im=im*max(exp_time)/exp_time(t)>0<10 & $
    star_click,im,x,y & $
    xc(t)=x & $
    yc(t)=y & $
    radial,rad,x,y,Nx,Ny & $
    pix=where(rad le 5) & $
    lmb=bytarr(Nx,Ny) & $
    lmb(pix)=1 & $
    tvscl,im*(1.0+0.5*lmb)<10 & $
    print,t,x,y & $
  endfor & $
  save,filename=star_files(star),xc,yc & $
endif

;register each frame
tm=[ 20, 3, 71, 12]		;all images are registered to these frames for each star
Nstars=n_elements(star_files)
flux=fltarr(Nimages,Nstars)
av=fltarr(Nstars)
sd=fltarr(Nstars)
ro=fltarr(Nimages,Nstars)
for star=0,Nstars-1 do begin $
  restore,star_files(star) & $
  sh_images=fltarr(Nx,Ny,Nimages) & $
  sh_mask=bytarr(Nx,Ny,Nimages) & $
  sh_flat=fltarr(Nx,Ny,Nimages) & $
  flatfix=flatfield & $	;replace bubbles with nearby values
  flatfix(75,70)=flatfield(150:225,70:165) & $
  for t=0,Nimages-1 do begin $
    if (xc(t) ge 0) then begin $
      dx=-(xc(t)-xc(tm(star))) & $
      dy=-(yc(t)-yc(tm(star))) & $
      sh_images(0,0,t)=shift(images(*,*,t),dx,dy) & $
      msk=shift(good_pix*mask(*,*,t),dx,dy) & $
      if (dx gt 0) then msk(0:dx,*)=0 & $
      if (dx lt 0) then msk(Nx-1+dx:Nx-1,*)=0 & $
      if (dy gt 0) then msk(*,0:dy)=0 & $
      if (dy lt 0) then msk(*,Ny-1+dy:Ny-1)=0 & $
      sh_mask(0,0,t)=msk & $
      sh_flat(0,0,t)=shift(flatfix,dx,dy) & $
      ;tvscl,sh_images(*,*,t)*sh_mask(*,*,t)*max(exp_time)/ $
      ;  exp_time(t)>0<10 & $
    endif & $
  endfor & $
  ;compute star's flux
  radius=10 & $
  xo=xc(tm(star)) & $
  yo=yc(tm(star)) & $
  radial,rad,radius,radius,2*radius+1,2*radius+1 & $
  star_pix=where(rad le 4.0) & $
  for t=0,Nimages-1 do begin $
    if (xc(t) ge 0) then begin $
      sub=sh_images(xo-radius:xo+radius,yo-radius:yo+radius,t) & $
      sub_mask=sh_mask(xo-radius:xo+radius,yo-radius:yo+radius,t) & $
      b_pix=where((sub_mask gt 0) and (rad gt 5.0) and (rad lt 10.0)) & $
      star_pix=where((sub_mask gt 0) and (rad le 4.0)) & $
      sub=sub-median(sub(b_pix)) & $
      prof=profile_idl(sub,radius,radius,r,err,sub_mask-1) & $
      plot,r,prof,psym=10,yrange=[-5,30],ystyle=1 & $
      oplot,r,r*0,color=128 & $
      flux(t,star)=total(sub(star_pix))/exp_time(t) & $
      tvscl,enlarge(sub*sub_mask>(-2)<10,3),250,250 & $
    endif & $
  endfor & $
  ;adjust result for flatfield effects
  x_opt=191.0 & $      ;Coordinates for the optical center.
  y_opt=286.0 & $
  flat_opt=flatfield(x_opt-radius:x_opt+radius,y_opt-radius:y_opt+radius) & $
  flat_opt=median(flatfield) & $
  correction=fltarr(Nimages) & $
  for t=0,Nimages-1 do begin $
    if (xc(t) ge 0) then begin $
      flat_c=sh_flat(xo-radius:xo+radius,yo-radius:yo+radius,t) & $
      mask_c=sh_mask(xo-radius:xo+radius,yo-radius:yo+radius,t) & $
      pix=where(mask_c gt 0) & $
      correction(t)=median(flat_c(pix))/flat_opt & $
    endif & $
  endfor & $
  flux(0,star)=flux(*,star)*correction & $
  phi=atan(plate_scale*(xc-x_opt)) & $
  theta=atan(plate_scale*(yc-y_opt)*cos(phi)) & $
  ro(0,star)=acos(cos(phi)*cos(theta)) & $
  j=where((xc ge 0.0) and (ro(*,star) le 0.210)) & $
  av(star)=avg(flux(j,star)) & $
  sd(star)=stdev(flux(j,star)) & $
  flux(0,star)=flux(*,star)/av(star) & $
  sd(star)=sd(star)/av(star) & $
endfor

;plot normalized flux
r2d=180.0/!pi
star=3
plot,ro(*,star)*r2d,flux(*,star),psym=7, $
  xrange=[0,max(ro*r2d)], $
  yrange=[0,1.5],charsize=1.3,ytitle='normalized flux', $
  ystyle=1,xtitle='distance from optical center    (degrees)'
plots,[0,30],[1,1]
sym=[4,6,8,7]
loadct,3
for star=0,Nstars-2 do begin $
  j=where(flux(*,star) gt 0.0) & $
  j=j(sort(ro(j,star))) & $
  oplot,ro(j,star)*r2d,flux(j,star),psym=sym(star), $
    color=200-60*star & $
  oplot,ro(j,star)*r2d,flux(j,star),color=200-60*star, $
    thick=3 & $
  oploterr,ro(j,star)*r2d,flux(j,star), $
    ro(j,star)*0+sd(star),3 & $
endfor
loadct,0
