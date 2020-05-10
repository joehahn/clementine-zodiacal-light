pro p_coords,coords,xcenter,ycenter,xdim,ydim

coords=fltarr(xdim,ydim,2)
x=findgen(xdim)-xcenter
y=findgen(ydim)-ycenter
for yi=0,ydim-1 do begin
  coords(0,yi,0)=sqrt(x^2+y(yi)^2)
  coords(0,yi,1)=atan(y(yi),x)
endfor
c=coords(*,*,1)
j=where(c lt 0.0)
if (j(0) ne -1) then c(j)=c(j)+2.0*!pi
coords(0,0,1)=c

;old code
;pi=acos(-1.)
;for y=0,ydim-1 do coords(0,y,0)=sqrt(x^2+float(y-ycenter)^2)
;for y=0,ycenter-1 do coords(0,y,1)=(2.0*pi-acos(x/coords(*,y,0)))
;for y=ycenter+1,ydim-1 do coords(0,y,1)=acos(x/coords(*,y,0))
;coords(0:xcenter-1,ycenter,1)=pi
;coords(xcenter,ycenter,1)=-3.5*pi
;c=coords(*,*,1)-pi/2.
;i=where(c lt 0.0)
;c(i)=c(i)+2.0*pi
;coords(0,0,1)=c
return
end
