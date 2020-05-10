;longitudes.pro

;makes lines-of-sight map

pflag=1

;latitute and longitude of each line-of-sight in degrees
r2d=180.0/!pi
orbit  = [     66,    110,  164,   193,   206,   253 ]
lambda = [ 345.06, 354.18, 5.53, 11.51, 13.96, 23.72 ]/r2d
beta   = [   0.00,   0.00, 0.00,  0.00,  0.00,  0.00 ]

;zodiacal light symmetry axis, degrees
omega_zl=87.0/r2d		;ascending node
i_zl=10.0/r2d			;inclination (actually, its 3deg)

;setup display
mx=1.5
zmax=mx
zmin=-zmax
xmax=mx
xmin=-xmax
ymax=mx
ymin=-ymax
!x.s=[-xmin,1]/(xmax-xmin)
!y.s=[-ymin,1]/(ymax-ymin)
!z.s=[-zmin,1]/(zmax-zmin)
if (pflag eq 0) then window,xs=501,ys=501
erase
t3d,/reset,oblique=[-1.0, 30.0],scale=0.85*[1,1,1], $
  translate=[0.5,0.3,0.0]

;draw axes
mx=1.5
if (pflag eq 1) then setplot
plots,[0,mx],[0,0],[0,0],/T3D,/DATA,thick=7,color=128
plots,[0,0],[0,mx],[0,0],/T3D,/DATA,thick=7,color=128
plots,[0,0],[0,0],[0,mx],/T3D,/DATA,thick=7,color=128

;draw galactic north pole
lambda_g=180.02/r2d
beta_g=45.0/r2d				;actually = 29.81 deg
length=1.3
z=length*cos(lambda_g)*cos(beta_g)
x=length*sin(lambda_g)*cos(beta_g)
y=length*sin(beta_g)
plots,[0,x],[0,y],[0,z],/T3D,/DATA,linestyle=0,thick=2
xyouts,0.58,1.13,'GNP',charsize=1.5,charthick=3

;draw earth's orbit
angle=2.0*!pi*findgen(101)/100.0
x=sin(angle)
y=0.0
z=cos(angle)
plots,x,y,z,/T3D,/DATA,thick=3

;draw lines-of-sight
N_los=n_elements(orbit)
x_los=sin(lambda)*cos(beta)
y_los=sin(beta)
z_los=cos(lambda)*cos(beta)
f=1.6
for i=0,N_los-1 do begin $
  plots,[0,x_los(i)],[0,y_los(i)],[0,z_los(i)], $
    /T3D,/DATA,thick=3 & $
  plots,-f*[0,x_los(i)],-f*[0,y_los(i)],-f*[0,z_los(i)], $
    /T3D,/DATA,thick=3 & $
endfor

;draw zl midplane
y_zl=i_zl*sin(angle-omega_zl)
plots,x,y_zl,z,/T3D,/DATA,linestyle=2,thick=3

;text
xyouts,-0.115,-0.15,'!9n!3',charsize=2,charthick=4	;sun
xyouts,-1.27,-0.85,'!20x!3',charsize=2,charthick=2	;eqnx
xyouts,-0.14,1.24,'ENP',charsize=1.5,charthick=3	;pole
xyouts,-1.08, -0.57,'66',charsize=0.8
xyouts,-0.80, -0.61,'164',charsize=0.8
xyouts,-0.45, -0.57,'253',charsize=0.8

if (pflag eq 1) then output_plot,'longitudes.ps'