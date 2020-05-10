;correct_flatfield.pro
;
; Correct the flatfield by monitoring
; the flux of various stars as the camera pans the
; sky in orbit66. This correction simply adjusts the
; star's fluxes so that they match their fluxes near the
; optical center.
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
flat_file='/users/hahn/clementine/zodiacal_light/flatfield/flatfield_short_exp.fits'
good_file='/users/hahn/clementine/zodiacal_light/flatfield/good_pix_short_exp.fits'

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
star_files=[ 'star0.dat', 'star1.dat', 'star2.dat', $
             'star3.dat', 'star4.dat' ]

;The plate-scale in radians per pixel
plate_scale=0.00131926
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
  im=im-median(im(x_dark0:x_dark1,y_dark0:y_dark1)) & $
  tvscl,im*max(exp_time)/exp_time(t)>0<20 & $

endfor

stop

;Find the saturated pixels, and then create a mask array
;having the same dimensions as the images array. Good rows &
;columns devoid of saturated pixels have mask=1, whereas bad
;rows and columns containing saturated pixels have mask=0
mask=bytarr(Nx,Ny,Nimages)+1
for t=0,Nimages-1 do begin $
  ;for y=0,Ny-1 do begin $			;DON'T mask rows.
  ;  j=where(images(*,y,t) ge 255) & $
  ;  if (j(0) ne -1) then mask(*,y,t)=0 & $
  ;endfor & $
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

;adjust flatfield by adjust_flat so that each pixel 
;subtends the same solid angle.
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
adjust_flat=cos(angle)

;flatfield these data using the adjusted flatfield
flatfield=readfits(flat_file)
flatfield=flatfield*adjust_flat
good_pix=readfits(good_file)
good=where(good_pix eq 1)
bad=where(good_pix ne 1)
flatfield(bad)=1.0
for t=0,Nimages-1 do $
  images(0,0,t)=images(*,*,t)/flatfield

;One first entry (e.g., if star_clk=0), click the mouse on the
;target star to get is (xc,yc) coordinates from each image.
star=4				;choose star!
window,xs=1.5*Nx,ys=Ny
if (star_clk eq 0) then begin $
  xc=fltarr(Nimages) & $
  yc=fltarr(Nimages) & $
  for t=0,Nimages-1 do begin $
    im=images(*,*,t)*mask(*,*,t) & $
    ;im=images(*,*,t) & $	;for star=5
    im=alog(im>1<100) & $
    star_click,im,x,y & $
    xc(t)=x & $
    yc(t)=y & $
    radial,rad,x,y,Nx,Ny & $
    pix=where(rad le 5) & $
    lmb=bytarr(Nx,Ny) & $
    lmb(pix)=1 & $
    tvscl,im*(1.0+0.5*lmb)<50 & $
    print,t,x,y & $
  endfor & $
  save,filename=star_files(star),xc,yc & $
endif

;register each frame
tm=[ 4, 9, 62, 22, 39]	;all images are registered
Nstars=n_elements(star_files)	;to these frames for each star
flux=fltarr(Nimages,Nstars)
error=fltarr(Nimages,Nstars)
ro=fltarr(Nimages,Nstars)
sh_images=fltarr(Nx,Ny,Nimages)
sh_mask=bytarr(Nx,Ny,Nimages)
for star=0,Nstars-1 do begin $
  restore,star_files(star) & $
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
      ;tvscl,sh_images(*,*,t)*sh_mask(*,*,t)*max(exp_time)/ $
      ;  exp_time(t)>0<30 & $
    endif & $
  endfor & $
  ;compute star's flux
  radius=10 & $
  xo=xc(tm(star)) & $
  yo=yc(tm(star)) & $
  radial,rad,radius,radius,2*radius+1,2*radius+1 & $
  star_pix=where(rad le 5.0) & $
  for t=0,Nimages-1 do begin $
    if (xc(t) ge 0) then begin $
      sub=sh_images(xo-radius:xo+radius,yo-radius:yo+radius,t) &$
      sub_mask=sh_mask(xo-radius:xo+radius,yo-radius: $
        yo+radius,t) & $
      b_pix=where((sub_mask gt 0) and (rad gt 5.0) and $
       (rad lt 8.0)) & $
      star_pix=where((sub_mask gt 0) and (rad le 5.0)) & $
      ;subtract ZL background
      sub=sub-median(sub(b_pix)) & $
      prof=profile_idl(sub,radius,radius,r,err,sub_mask-1) & $
      plot,r,prof,psym=10,yrange=[-5,30],ystyle=1 & $
      oplot,r,r*0,color=128 & $
      flux(t,star)=total(sub(star_pix))/exp_time(t) & $
      tvscl,enlarge(sub*sub_mask>(-2)<40,3),250,250 & $
    endif & $
  endfor & $
  ;Get star's angular distance from optical center
  theta=asin(plate_scale*(yc-y_opt)) & $
  phi=asin(plate_scale*(xc-x_opt)/cos(theta)) & $
  ro(0,star)=acos(cos(phi)*cos(theta)) & $
endfor

;re-sort fluxes with distance
for star=0,Nstars-1 do begin $
  j=sort(ro(*,star)) & $
  ro(0,star)=ro(j,star) & $
  flux(0,star)=flux(j,star) & $
endfor

;plot normalized flux
r2d=180.0/!pi
sym=[4,6,8,7,1,5]
j=where(flux gt 0.0)
thck=3
setplot
plot,ro(j)*r2d,flux(j),psym=3,xrange=[0,30], $
  yrange=[0.5,1.2],nodata=1,ystyle=1, $
  charsize=1.3,ytitle='Normalized flux', $
  xtitle='Angular distance !7U!3!l0!n    (degrees)',$
  xthick=thck,ythick=thck,charthick=thck
;generate correction factor
a=0.5*findgen(1001)/1000.0
plots,[0,30],[1,1],color=128
for star=0,Nstars-1 do begin $
  j=where((ro(*,star)*r2d gt 12.5) and $
    (ro(*,star)*r2d lt 17.5)) & $
  factor=avg(flux(j,star))*1.0 & $
  j=where(flux(*,star) gt 0.0) & $
  oplot,ro(j,star)*r2d,flux(j,star)/factor,thick=thck & $
  oplot,ro(j,star)*r2d,flux(j,star)/factor,psym=sym(star), $
    thick=thck & $
endfor

;redo using non-flatfielded images
for t=0,Nimages-1 do $
  images(0,0,t)=images(*,*,t)*flatfield
;register each frame
tm=[ 4, 9, 62, 22, 39]	;all images are registered
Nstars=n_elements(star_files)	;to these frames for each star
flux2=fltarr(Nimages,Nstars)
error=fltarr(Nimages,Nstars)
ro=fltarr(Nimages,Nstars)
sh_images=fltarr(Nx,Ny,Nimages)
sh_mask=bytarr(Nx,Ny,Nimages)
for star=0,Nstars-1 do begin $
  restore,star_files(star) & $
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
      ;tvscl,sh_images(*,*,t)*sh_mask(*,*,t)*max(exp_time)/ $
      ;  exp_time(t)>0<30 & $
    endif & $
  endfor & $
  ;compute star's flux
  radius=10 & $
  xo=xc(tm(star)) & $
  yo=yc(tm(star)) & $
  radial,rad,radius,radius,2*radius+1,2*radius+1 & $
  star_pix=where(rad le 5.0) & $
  for t=0,Nimages-1 do begin $
    if (xc(t) ge 0) then begin $
      sub=sh_images(xo-radius:xo+radius,yo-radius:yo+radius,t) &$
      sub_mask=sh_mask(xo-radius:xo+radius,yo-radius: $
        yo+radius,t) & $
      b_pix=where((sub_mask gt 0) and (rad gt 5.0) and $
       (rad lt 8.0)) & $
      star_pix=where((sub_mask gt 0) and (rad le 5.0)) & $
      ;subtract ZL background
      sub=sub-median(sub(b_pix)) & $
      ;prof=profile_idl(sub,radius,radius,r,err,sub_mask-1) & $
      ;plot,r,prof,psym=10,yrange=[-5,30],ystyle=1 & $
      ;oplot,r,r*0,color=128 & $
      flux2(t,star)=total(sub(star_pix))/exp_time(t) & $
      ;tvscl,enlarge(sub*sub_mask>(-2)<40,3),250,250 & $
    endif & $
  endfor & $
  ;Get star's angular distance from optical center
  theta=asin(plate_scale*(yc-y_opt)) & $
  phi=asin(plate_scale*(xc-x_opt)/cos(theta)) & $
  ro(0,star)=acos(cos(phi)*cos(theta)) & $
endfor
;re-sort fluxes with distance
for star=0,Nstars-1 do begin $
  j=sort(ro(*,star)) & $
  ro(0,star)=ro(j,star) & $
  flux2(0,star)=flux2(j,star) & $
endfor
for star=0,Nstars-1 do begin $
  j=where((ro(*,star)*r2d gt 8.0) and $
    (ro(*,star)*r2d lt 12.0)) & $
  factor=avg(flux2(j,star))*1.0 & $
  j=where(flux2(*,star) gt 0.0) & $
  oplot,ro(j,star)*r2d,flux2(j,star)/factor,color=160, $
    thick=thck & $
  oplot,ro(j,star)*r2d,flux2(j,star)/factor,psym=sym(star), $
    color=160,thick=thck & $
endfor
output_plot,'flat_test.ps'
