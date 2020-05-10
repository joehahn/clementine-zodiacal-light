;flatfield.pro

;Computes a flatfield for a given value of dark_current.

;Coordinates for the optical center.
x_c=191.0      		
y_c=286.0

;subtract additional dark-current offset
dark_current=100.0

;read Bonnie's flatfields
flat0=readfits('/u2/hahn/clementine/flatfields/BFB32060.fits',h0)
flat1=readfits('/u2/hahn/clementine/flatfields/BFB32061.fits',h1)
flat2=readfits('/u2/hahn/clementine/flatfields/BFB32062.fits',h2)
flat3=readfits('/u2/hahn/clementine/flatfields/BFB62060.fits',h3)
flat4=readfits('/u2/hahn/clementine/flatfields/BFB62061.fits',h4)
flat5=readfits('/u2/hahn/clementine/flatfields/BFB62062.fits',h5)
sz=size(flat0)
Nx=sz(1)
Ny=sz(2)
window,xs=2*Nx,ys=Ny

;form the median-average of flat0 & flat1 & flat2,
;and its standard-deviation. These images were probably
;acquired at similar (but shorter) exposures.
;f012=[[[flat0]],[[flat1]],[[flat2]]]
;flat012=average_image(f012)
;sigma012=stdev_image(f012)

;form the median-average of flat3 & flat4 & flat5,
;and its standard-deviation. These images were probably
;acquired at similar (but longer) exposures.
f345=[[[flat3]],[[flat4]],[[flat5]]]
f345_avg=average_image(f345)
sigma345=stdev_image(f345)

;subtract dark-current & destreak
flat345=f345_avg-dark_current
flat345=destreak(flat345,0.06)

;factor out the filter's transmittance
n=1.53			;filter's index of refraction
plate_scale=0.00137396	;plate scale (rad/pixel)
x=fltarr(Nx,Ny)
xx=findgen(Nx)-x_c
for i=0,Ny-1 do x(0,i)=xx
y=fltarr(Nx,Ny)
yy=transpose(findgen(Ny)-y_c)
for i=0,Nx-1 do y(i,0)=yy
phi=atan(plate_scale*x)
theta=atan(plate_scale*y*cos(phi))
angle=acos(cos(phi)*cos(theta))
exponent=1.0/sqrt(1.0-(sin(angle)/n)^2)-1.0
trans=0.01^exponent
flat345=flat345/trans

;extrapolate into the vignetted corners
flat=flat345
rng=[0,350]
radial,rad,x_c,y_c,Nx,Ny
tmplt=intarr(Nx,Ny)
prof=profile_idl(flat,x_c,y_c,r,err,tmplt)
plot,rad,flat,psym=3,xrange=rng
oplot,r,prof,color=0,thick=7
Nzones=16
dt=2.0*!pi/Nzones
for i=0,Nzones-1 do begin $
  prof=profile_cone(flat,x_c,y_c,(i+0.501)*dt,dt/2,r,err, $
    tmplt,coords) & $
  if (i eq 14) then prof=profile_cone(flat,x_c,y_c, $
    5.69414,1.5*dt,r,err,tmplt,coords) & $
  pix=where(tmplt gt 0) & $
  good=where((rad ge 230) and (rad lt 270) and (tmplt eq 1)) & $
  bad=where((rad gt 275) and (tmplt eq 1)) & $
  if ((n_elements(good) gt 1) and (n_elements(bad) gt 1)) then begin $
      tmplt(good)=2 & $
      plot,rad(pix),flat(pix),psym=1,xrange=rng & $
      oplot,r,prof,color=64,thick=5 & $
      gd=where((r ge 196) and (r lt 270)) & $
      bad=where((tmplt gt 0) and (rad gt 275.0)) & $
      oplot,r(gd),prof(gd),thick=4,color=196 & $
      coeff=poly_fit(r(gd),prof(gd),2) & $
      prof_fit=coeff(0)+coeff(1)*r+coeff(2)*r^2 & $
      oplot,r,prof_fit,thick=2,color=128 & $
      noise=stdev(flat(good)) & $
      if (i eq 14) then noise=0.5*stdev(flat(good)) & $
      for j=0,n_elements(bad)-1 do begin $
        d=abs(r-rad(bad(j))) & $
        mn=min(d,idx) & $
        flat(bad(j))=prof_fit(idx)+noise*randomn(seed) & $
      endfor & $
      oplot,rad(bad),flat(bad),psym=1,color=64 & $
  endif & $
  tvscl,flat*(1+0.1*tmplt)>40 & $
  tvscl,flat>40,Nx,0 & $
endfor

;Rotate, transpose, and shift such
;that the blemish is aligned with the clementine data
flat=shift(rotate(transpose(flat),3),-1,2)
;zap outer rows & columns
good_pixels=fltarr(Nx,Ny)+1
good_pixels(*,0:1)=0
good_pixels(0,*)=0
good_pixels(Nx-2:*,*)=0
good=where(good_pixels ne 0)
bad=where(good_pixels eq 0)

;display & save flatfield
flatfield=fltarr(Nx,Ny)
flatfield(good)=flat(good)/avg(flat(good))
flatfield(bad)=1.0
writefits, 'flatfield.fits', flatfield
writefits, 'good_pix.fits', good_pixels
restore,'orbit193.dat'
flux=fltarr(Nx,Ny)
good=where(good_pix ne 0)
flux(good)=counts(good)/total_exp_time(good)
flatflux=fltarr(Nx,Ny)
flatflux(good)=flux(good)/flatfield(good)
tvscl,enlarge(flux(70:170,70:170),5)>14<32
tvscl,enlarge(shift(flatfield(70:170,70:170),0,0),5)
tvscl,enlarge(flatflux(70:170,70:170),5)>14<32

;loop over numerous values of dark-current
N=61
flatfields=fltarr(Nx,Ny,N)
images=fltarr(Nx,Ny,N)
dc=fltarr(N)
for index=0,N-1 do begin $
  dc(index)=dark_current-(index-N/2) & $
  flat=f345_avg-dc(index) & $
  flat=destreak(flat,0.06) & $
  flat=flat/trans & $
  for i=0,Nzones-1 do begin $
    prof=profile_cone(flat,x_c,y_c,(i+0.501)*dt,dt/2,r,err, $
      tmplt,coords) & $
    if (i eq 14) then prof=profile_cone(flat,x_c,y_c, $
      5.69414,1.5*dt,r,err,tmplt,coords) & $
    pix=where(tmplt gt 0) & $
    gd=where((rad ge 230) and (rad lt 270) and (tmplt eq 1)) & $
    bd=where((rad gt 275) and (tmplt eq 1)) & $
    if ((n_elements(gd) gt 1) and (n_elements(bd) gt 1)) then begin $
        tmplt(gd)=2 & $
        plot,rad(pix),flat(pix),psym=1,xrange=rng & $
        oplot,r,prof,color=64,thick=5 & $
        g=where((r ge 196) and (r lt 270)) & $
        b=where((tmplt gt 0) and (rad gt 275.0)) & $
        oplot,r(g),prof(g),thick=4,color=196 & $
        coeff=poly_fit(r(g),prof(g),2) & $
        prof_fit=coeff(0)+coeff(1)*r+coeff(2)*r^2 & $
        oplot,r,prof_fit,thick=2,color=128 & $
        noise=stdev(flat(gd)) & $
        if (i eq 14) then noise=0.5*stdev(flat(gd)) & $
        for j=0,n_elements(b)-1 do begin $
          d=abs(r-rad(b(j))) & $
          mn=min(d,idx) & $
          flat(b(j))=prof_fit(idx)+noise*randomn(seed) & $
        endfor & $
        oplot,rad(b),flat(b),psym=1,color=64 & $
    endif & $
    tvscl,flat*(1+0.1*tmplt)>40 & $
    tvscl,flat>40,Nx,0 & $
  endfor & $
  flat=shift(rotate(transpose(flat),3),-1,2) & $
  ff=fltarr(Nx,Ny) & $
  ff(good)=flat(good)/avg(flat(good)) & $
  ff(bad)=1.0 & $
  flatfields(0,0,index)=ff & $
  flatflux=fltarr(Nx,Ny) & $
  flatflux(good)=flux(good)/ff(good) & $
  images(0,0,index)=flatflux & $
  print,index,dc(index) & $
endfor

;display results
z=''
for index=0,N-1 do begin $
  tvscl,enlarge(images(70:170,70:170,index),5)>14<32 & $
  print,index,dc(index) & $
  ;read,z & $
endfor

