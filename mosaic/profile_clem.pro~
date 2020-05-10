function profile_clem,array,xcenter,ycenter,r,tmplt,err

;Inputs: array, xcenter, ycenter, r, tmplt
;	use pixels with tmplt>0
;Outputs: err

N=n_elements(r)
prof=fltarr(N)
err=fltarr(N)

sz=size(array)
Nx=sz(1)
Ny=sz(2)
radial,rad,xcenter,ycenter,Nx,Ny

pix=where(tmplt gt 0)
Npix=n_elements(pix)
dist=rad(pix)
flux=array(pix)
s=sort(dist)
dist=dist(s)
flux=flux(s)

;get first point
mn=where(dist ge r(0))
mn=mn(0)
mx=where(dist ge r(1))
mx=mx(0)
Np=mx-mn+1
prof(0)=avg(flux(mn:mx))
;err(0)=stdev(flux(mn:mx))/sqrt(Np-1)
mom=moment(flux(mn:mx), sdev=stdev)
err(0)=stdev/sqrt(Np-1)

;get intermediate points
mn=where(dist ge r(0))
mn=mn(0)
mx=where(dist ge r(2))
mx=mx(0)-1
for i=1,N-2 do begin $
  prof(i)=avg(flux(mn:mx)) & $
  Np=mx-mn+1 & $
  ;err(i)=stdev(flux(mn:mx))/sqrt(Np-1) & $
  mom=moment(flux(mn:mx), sdev=stdev) & $
  err(i)=stdev/sqrt(Np-1) & $ 
  ;print,i,r(i-1),dist(mn),dist(mx),r(i+1) & $
  repeat mn=mn+1 until (dist(mn) gt r(i)) & $
  repeat mx=(mx+1)<(Npix-1) $
    until ((dist(mx) ge r((i+2)<(N-1))) or (mx eq Npix-1)) & $
  if (mx lt Npix-1) then mx=mx-1 & $
endfor

;get last point
mn=where(dist ge r(N-2))
mn=mn(0)
mx=Npix-1
Np=mx-mn+1
prof(N-1)=avg(flux(mn:mx))
;err(N-1)=stdev(flux(mn:mx))/sqrt(Np-1)
mom=moment(flux(mn:mx), sdev=stdev)
err(N-1)=stdev/sqrt(Np-1)

return,prof
end