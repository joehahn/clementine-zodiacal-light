;geometry.pro

;draw viewing geometry

pflag=1
r2d=180.0/!pi
Po2=!pi/2.0

r1=1.0
D=2.1
theta=20.0/r2d
phi=75.0/r2d

eps=acos(cos(theta)*cos(phi))
r=sqrt(r1^2+D^2-2.0*r1*D*cos(eps))
beta=asin(D*sin(theta)/r)

;setup display
mx=1.3*r1
zmax=mx
zmin=-zmax
xmax=mx
xmin=-xmax
ymax=mx
ymin=-ymax
!x.s=[-xmin,1]/(xmax-xmin)
!y.s=[-ymin,1]/(ymax-ymin)
!z.s=[-zmin,1]/(zmax-zmin)
if (pflag eq 0) then window,xs=600,ys=600
if (pflag eq 1) then setplot
erase
t3d,/reset,oblique=[-1.0, 45.0],scale=0.85*[1,1,1], $
  translate=[0.3,0.5,0.0]

;draw axes
plots,[0,1.75],[0,0],[0,0],/T3D,/DATA,thick=3,color=128
plots,[0,0],[0,0.75],[0,0],/T3D,/DATA,thick=3,color=128
plots,[0,0],[0,0],[0,1.3],/T3D,/DATA,thick=3,color=128

;draw r1, r, D
g=r1^2+(r*cos(beta))^2-(D*cos(theta))^2
g=acos(g/2.0/r1/r/cos(beta))
;r
rv=r*[ cos(beta)*sin(g), sin(beta), cos(beta)*cos(g) ]
rvl=1.2*rv
;r1
r1v=r1*[ 0.0, 0.0, 1.0 ]
plots,[0,rv(0)],[0,0],[0,rv(2)],/T3D,/DATA,color=128,thick=5
plots,[rv(0),rv(0)],[0,rv(1)],[rv(2),rv(2)],/T3D,/DATA,color=128,thick=5
plots,[r1v(0),rv(0)],[0,0],[r1v(2),rv(2)],/T3D,/DATA,color=128,thick=5
plots,[0,r1v(0)],[0,r1v(1)],[0,r1v(2)],/T3D,/DATA,thick=7
plots,r1v(0),r1v(1),r1v(2),/T3D,/DATA,psym=8,symsize=2
plots,[0,rvl(0)],[0,rvl(1)],[0,rvl(2)],/T3D,/DATA,color=128,thick=5
plots,[0,rv(0)],[0,rv(1)],[0,rv(2)],/T3D,/DATA,thick=7
plots,rv(0),rv(1),rv(2),/T3D,/DATA,psym=8,symsize=2
plots,[r1v(0),rv(0)],[r1v(1),rv(1)],[r1v(2),rv(2)], $
  /T3D,/DATA,thick=7

;draw phi arc
N=101
angle=phi*findgen(N)/(N-1)
l=0.2*r1
x_arc=r1v(0)+l*sin(angle)
y_arc=r1v(1)+0*sin(angle)
z_arc=r1v(2)-l*cos(angle)
plots,x_arc,y_arc,z_arc,/T3D,/DATA,thick=3

;draw theta arc
angle=theta*findgen(N)/(N-1)
l=0.4*r1
x=l*cos(angle)
y=l*sin(angle)
z=0*sin(angle)
rot_angle=-(Po2-phi)
x_arc=r1v(0)+x*cos(rot_angle)-z*sin(rot_angle)
y_arc=r1v(1)+y
z_arc=r1v(2)+x*sin(rot_angle)+z*cos(rot_angle)
plots,x_arc,y_arc,z_arc,/T3D,/DATA,thick=3

;draw epsilon arc
angle=eps*findgen(N)/(N-1)
l=0.4*r1
x=l*sin(angle)
y=0*sin(angle)
z=-l*cos(angle)
rot_angle=-theta
x_arc=r1v(0)+x*cos(rot_angle)+y*sin(rot_angle)
y_arc=r1v(1)-x*sin(rot_angle)+y*cos(rot_angle)
z_arc=r1v(2)+z
plots,x_arc,y_arc,z_arc,/T3D,/DATA,thick=3

;draw scattering-angle arc
a0=Po2-phi
scat_angle=g+phi
angle=scat_angle*findgen(N)/(N-1)+a0
l=0.2*r1
x=l*cos(angle)
y=0*sin(angle)
z=l*sin(angle)
rot_angle=-theta
x_arc=rv(0)+x*cos(rot_angle)+y*sin(rot_angle)
y_arc=rv(1)-x*sin(rot_angle)+y*cos(rot_angle)
z_arc=rv(2)+z
plots,x_arc,y_arc,z_arc,/T3D,/DATA,thick=3

;draw beta arc
angle=beta*findgen(N)/(N-1)
l=0.4*r1
x=l*cos(angle)
y=l*sin(angle)
z=0*sin(angle)
rot_angle=(Po2-g)
x_arc=x*cos(rot_angle)-z*sin(rot_angle)
y_arc=y
z_arc=x*sin(rot_angle)+z*cos(rot_angle)
plots,x_arc,y_arc,z_arc,/T3D,/DATA,thick=3

;text
xyouts,-1.30,-0.40,'observer',charsize=1.5,charthick=3
xyouts,-0.53,0.178,'Sun',charsize=1.5,charthick=3
xyouts,0.80,0.62,'volume',charsize=1.5,charthick=3
xyouts,0.80,0.525,'element',charsize=1.5,charthick=3
xyouts,-0.74,-0.34,'!7u!3',charsize=1.2,charthick=3
xyouts,-0.51,-0.37,'!7h!3',charsize=1.2,charthick=3
xyouts,-0.52,-0.20,'!7e!3',charsize=1.2,charthick=3
xyouts,-0.025, 0.11,'!7b!3',charsize=1.2,charthick=3
xyouts,0.88, 0.27,'!9P!3',charsize=1.2,charthick=3
xyouts,-0.36,0.09,'!9n!3',charsize=2,charthick=3
xyouts,-1.15,-0.71,'x',charsize=1.5,charthick=3
xyouts,1.20,0.12,'y',charsize=1.5,charthick=3
xyouts,-0.33,0.81,'z',charsize=1.5,charthick=3
xyouts,-0.07,-0.16,'!7D!3',charsize=1.5,charthick=3
xyouts,0.24,0.32,'r',charsize=1.5,charthick=3
xyouts,-0.57,-0.03,'r',charsize=1.5,charthick=3
xyouts,-0.545,-0.045,'1',charsize=0.7,charthick=3

if (pflag eq 1) then output_plot,'geometry.ps'