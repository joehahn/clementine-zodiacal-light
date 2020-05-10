pro limb_click, image, xc, yc, R

;This proceedure displays the image, and then the user clicks on
;4 sequential points along the lunar limb. The proceedure then
;determines the Moon's center coordinates (xc,yc) and radius R
;in pixels, assuming a spherical Moon.

;rescale & display image
img=image-min(image)
img=img*255.0/max(img)
tv,img

;get image size
sz=size(image)
Nx=sz(1)
Ny=sz(2)

x=fltarr(4)
y=fltarr(4)
for i=0,3 do begin $
  ;get point
  cursor,x0,y0,down=1,device=1 & $
  x(i)=x0 & $
  y(i)=y0 & $
  ;display +
  img((x(i)-3)>0<(Nx-1):(x(i)+3)>0<(Nx-1),y(i)>0<(Ny-1))=128 & $
  img(x(i)>0<(Nx-1),(y(i)-3)>0<(Ny-1):(y(i)+3)>0<(Ny-1))=128 & $
  tv,img & $
endfor

;Create a chord connecting points (x(0),y(0)) to (x(2),y(2))
;on the lunar limb. This chord is the line y=m02*x+b02 where:
m02=(y(2)-y(0))/(x(2)-x(0))
b02=y(2)-m02*x(2)

;Create a chord connecting points (x(1),y(1)) to (x(3),y(3))
;on the lunar limb. This chord is the line y=m13*x+b13 where:
m13=(y(3)-y(1))/(x(3)-x(1))
b13=y(3)-m13*x(3)

;The lines perpendicular to these chords have slopes
;-1/m02 and -1/m13 and y-intercepts:
b02_perp=0.5*(y(2)+y(0)+(x(2)+x(0))/m02)
b13_perp=0.5*(y(3)+y(1)+(x(3)+x(1))/m13)

;These perpendiculars intersect at the center of the moon at
xc=(b13_perp-b02_perp)/(1.0/m13-1.0/m02)
yc=-xc/m13+b13_perp

;The radius of the moon is
R=sqrt((x(1)-xc)^2+(y(1)-yc)^2)

end
