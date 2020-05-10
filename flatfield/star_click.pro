pro star_click, image, xc, yc

;The image is displayed, the user clicks a star, and then
;a logarithmic contour map is displayed. Left-click on the 
;star to put it in the crosshairs, and the right-click the
;star to quit and return its (xc,yc) coordinates.
;Right-clicking the image immediatly returns (xc,yc)=(-1,-1).

;display image
tvscl,image

;get image size
sz=size(image)
Nx=sz(1)
Ny=sz(2)

;get star's initial (x,y) coordinates
cursor,xc,yc,wait=1,device=1
xc=round(xc)
yc=round(yc)

;quit on right-click
if (!err eq 4) then begin $
  xc=-1.0 & $
  yc=-1.0 & $
  return & $
endif

;get x-y range for contour plot
rad=7
xmin=(xc-rad)>0
xmax=(xc+rad)<(Nx-1)
ymin=(yc-rad)>0
ymax=(yc+rad)<(Ny-1)

;contour around star and get new (x,y) coordinates
!err=0
wait,0.2
while (!err ne 4) do begin $
  sub=alog(abs(image(xmin:xmax,ymin:ymax))>1) & $
  x_axis=indgen(xmax-xmin+1)+xmin & $
  y_axis=indgen(ymax-ymin+1)+ymin & $
  contour,sub,x_axis,y_axis,nlevels=9,xrange=[xmin,xmax], $
    yrange=[ymin,ymax],xstyle=1,ystyle=1, $
    xtitle='X    (pixels)',ytitle='Y    (pixels)' & $
  plots,[xmin,xmax],0.5*(ymin+ymax)+[0,0],color=128 & $
  plots,0.5*(xmin+xmax)+[0,0],[ymin,ymax],color=128 & $
  cursor,xc,yc,wait=1,data=1 & $
  xc=round(xc) & $
  yc=round(yc) & $
  xmin=(xc-rad)>0 & $
  xmax=(xc+rad)<(Nx-1) & $
  ymin=(yc-rad)>0 & $
  ymax=(yc+rad)<(Ny-1) & $
  wait,0.2 & $
endwhile

return
end
