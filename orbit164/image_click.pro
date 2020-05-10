pro image_click, image, xc, yc

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

return
end
