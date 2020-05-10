pro register_image, array, reference, x_shift, y_shift, $
  MAG=mag, X0=x0, X1=x1, Y0=y0, Y1=y1

;Determine the (x,y) shifts nedded to register each frame in the
;datacube array with respect to frame z=reference.
;If magnification MAG > 1, then each frame is magnified MAG times
;using the bicubic INTERPOLATE(); partial-pixel shifts are then
;computed to a precision 1/MAG. The crosscor() function
;is used to determine the needed shifts. Use SHIFT_IMAGE()
;to actually perform those shifts.
;
;INPUTS:
;	array = three dimensional datacube
;	reference = z-index of reference frame. All other images
;		    are regestered relative to this frame.
;
;OPTIONAL INPUTS:
;	MAG = zoom this image by a magnification factor MAG in
;	      order to calculate shifts to a precision of 1/MAG pixels.
;	X0,X1,Y0,Y1 = use the sub-cube array(X0:X1, Y0:Y1, *).
;OUTPUTS:
;	x_shift(z) = array indicating the shift (in pixels) along the
;		     x-axis needed to register frame z with the
;		     reference frame.
;	y_shift    = ditto for the y-axis.
;

;extract subarray if necessary
if (keyword_set(x0) eq 0) then begin $
  sub_array=array & $
endif else begin $
  sub_array=array(x0:x1,y0:y1,*) & $
endelse
sz=size(sub_array)
Nx=sz(1)
Ny=sz(2)
Nz=sz(3)

;magnify image if necessary
if (keyword_set(mag) eq 0) then mag=1
if (mag gt 1) then begin $
  Nx=Nx*mag & $
  Ny=Ny*mag & $
  x=findgen(Nx)/mag & $
  y=findgen(Ny)/mag & $
  sub_array_mag=fltarr(Nx,Ny,Nz) & $
  for z=0,Nz-1 do begin $
    sub_array_mag(0,0,z)=interpolate(sub_array(*,*,z), x, y, $
      grid=1, cubic=-0.5) & $
  endfor & $
  sub_array=sub_array_mag & $
endif

;setup
x_shift=fltarr(Nz)
y_shift=fltarr(Nz)
ref_image=sub_array(*,*,reference)

;loop over all images in the datacube
for z=0,Nz-1 do $
  if (z ne reference) then begin $
    other_image=sub_array(*,*,z) & $
    ;compute cross-correlation between reference and image z
    cc=crosscor(ref_image, other_image) & $
    ;find max value in cross-correlation map.
    mx=max(cc,index) & $
    x_shift(z)=(index mod Nx) - Nx/2 & $
    y_shift(z)=(index/Nx) - Ny/2 & $
  endif

;rescale shifts if necessary
if (mag gt 1) then begin $
  x_shift=x_shift/float(mag) & $
  y_shift=y_shift/float(mag) & $
endif

end
