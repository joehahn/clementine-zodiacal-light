;3D_phase_integral.pro

pflag=0

;restore observed image
if (keyword_set(mosaic) eq 0) then begin $
  restore, filename= $
    '/users/hahn/clementine/zodiacal_light/mosaic/mercator.dat'&$
  if (pflag eq 0) then window,xs=601,ys=601 & $
endif

;ZL model parameters
alpha=-coeff(1)		;brightness radial power law
k1=1.15			;ellipticity (smaller=more circular)
k2=1.00			;pointy-parameter (larger=rounder)

;sampling grid
theta_max=30.0			;maximum angle in degrees
delta_theta=0.5			;grid spacing
r2d=180./!pi
Ngrid=2.0*theta_max/delta_theta+1
theta=2.0*theta_max*(findgen(Ngrid)/(Ngrid-1)-0.5)/r2d
phi=theta
phase_integral=fltarr(Ngrid,Ngrid)

Npoints=31
d2r=!pi/180.0
for theta_index=0,Ngrid-1 do begin $
  for phi_index=0,Ngrid-1 do begin $
    thet=theta(theta_index) & $
    ph=phi(phi_index) & $
    if ((thet ne 0.0) or (ph ne 0.0)) then begin $
      epsilon=acos(cos(ph)*cos(thet)) & $
      angle=(!pi-epsilon)*findgen(Npoints)/(Npoints-1)+epsilon & $
      ;Zook's phase law
      K=0.77 & $
      phase_fn=exp(-K*(!pi-angle)) & $
      ;Leinart's latitude law
      sin_beta=sin(angle-epsilon)*sin(thet)/sin(epsilon) & $
      h=exp(-k1*(abs(sin_beta)^k2)) & $
      ;new fit
      beta=asin(abs(sin_beta)<1.0) & $
      h_fit=h-0.315 & $
      h_fit=h_fit/(h_fit) & $
      ;integrand
      integrand=phase_fn*h*(abs(sin(angle))^alpha) & $
      ;integrand=phase_fn*h_fit*(abs(sin(angle))^alpha) & $
      phase_integral(phi_index,theta_index)= $
        int_tabulated(angle,integrand)/(sin(epsilon)^alpha) & $
    endif & $
  endfor & $
endfor
N2=Ngrid/2
phase_integral(N2,N2)=max(phase_integral)
factor=0.87*exp(coeff(0))/phase_integral(N2+5,N2)/ $
  (sin(phi(N2+5))^alpha)
model=factor*phase_integral

;fix inner gap in mosaic & high-noise areas
sz=size(mosaic)
Nxs=sz(1)
Nys=sz(2)
mosaic_fix=mosaic
radial,rad,x_sun,y_sun,Nxs,Nys
j=where((rad le 150) and (mosaic eq 0.0))
mosaic_fix(j)=max(mosaic)
mosaic_fix(630:*,560:*)=0
mosaic_fix(630:*,0:110)=0

;compare model image to mosaic
thck=3
if (pflag eq 1) then setplot
contour,model,phi*r2d,theta*r2d,xstyle=1,ystyle=1, $
  levels=[15.6,31.3,62.5,125.,250.,500.,1000.],color=165, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
contour,mosaic_fix,x,y,xstyle=1,ystyle=1,noerase=1, $
  xrange=theta_max*[-1,1],yrange=theta_max*[-1,1], $
  levels=[31.3,62.5,125.,250.,500.,1000.],thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
contour,smooth(mosaic_fix,5),x,y,xstyle=1,ystyle=1,levels=[15.6], $
  xrange=theta_max*[-1,1],yrange=theta_max*[-1,1],noerase=1, $
  xtitle='ecliptic longitude !7u!3    (degrees)',thick=thck, $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3
if (pflag eq 1) then output_plot,'fit.ps'

stop

;profile
epsilon=phi(Ngrid/2+1:*)
profile=model(Ngrid/2+1:*,Ngrid/2)
plot_oo,epsilon*r2d,profile,xstyle=1
coeff=poly_fit(alog(sin(epsilon)),alog(profile),1,yfit)
;oplot,epsilon*r2d,exp(yfit),psym=8
j=where(x lt 0.0)
oplot,-x(j),mosaic(j,330),psym=1
print,'profile power-law index = ',coeff(1)

stop
