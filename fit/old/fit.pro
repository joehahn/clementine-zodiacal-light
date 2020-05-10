;fit.pro

;dust-model parameters
nu=1.45			;dust radial power law
asr=9.5e-9		;albedo*sigma_1*r_1, +/- 10%
Omega_sun=6.8e-5	;Sun solid angle, steradians

;calibration converts 1 count/sec/pixel into solar brightness
cal=7.45e-14

;sigma's for three different dust populations
sigma_pop = [0.000, 0.122, 0.576]

;ecliptic contributions by the three different dust populations
contrib   = [0.19, 0.27, 0.54]

;Calculate new surface-brightness integrals when <> 0
calc_SB=0

;generate postscript output if pflag=1
pflag=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;restore observed image
if (keyword_set(mosaic) eq 0) then begin $
  restore, filename= $
    '/users/hahn/clementine/zodiacal_light/mosaic/mercator.dat'&$
endif
;fix inner gap in mosaic & high-noise areas
sz=size(mosaic)
Nx=sz(1)
Ny=sz(2)
radial,rad,x_sun,y_sun,Nx,Ny
j=where((rad le 150) and (mosaic eq 0.0))
mosaic(j)=max(mosaic)
mosaic(610:*,560:*)=0
mosaic(630:*,0:110)=0
Nx=800
Ny=800
mosaic_fix=fltarr(Nx,Ny)
good_fix=fltarr(Nx,Ny)
dx=Nx/2-round(x_sun)+1
dy=Ny/2-round(y_sun)
mosaic_fix(dx,dy)=mosaic
good_fix(dx,dy)=good
x_sun=Nx/2
y_sun=Ny/2
x=ps_avg*(findgen(Nx)-x_sun)
y=ps_avg*(findgen(Ny)-y_sun)

;sky-coordinates (phi,theta)=(longitude,latitude)
grid_factor=5
Ngrid=Nx/grid_factor
phi_grid=max(x)*2.0*(findgen(Ngrid)/(Ngrid-1)-0.5)
theta_grid=phi_grid
Npop=n_elements(sigma_pop)
common share, beta, sigma
.run latitude_integrand.pro
r2d=180.0/!pi
Po2=!pi/2.0

if (calc_SB ne 0) then begin $

  ;Calculate latitude distribution h(beta) for diff' populations.
  Nbeta=201 & $
  beta_tab=Po2*findgen(Nbeta)/(Nbeta-1) & $
  h_tab=fltarr(Nbeta,Npop) & $
  print,'Calculating the latitude distribution h(beta)...' & $
  for j=0,Npop-1 do begin $ & $
    sigma=sigma_pop(j) & $
    for i=0,Nbeta-1 do begin $
      beta=beta_tab(i) & $
      h_tab(i,j)=qromo('latitude_integrand',beta,Po2,double=0, $
        eps=1.0e-5) & $
    endfor & $
    h_tab(0,j)=h_tab(*,j)/h_tab(0,j) & $    ;normalize to h(0)=1
    print,float(j+1)/Npop & $
  endfor & $
  ;some fixes
  h_tab(Nbeta-1,0)=h_tab(Nbeta-2,*) & $
  j=where(sigma_pop eq 0.0) & $
  h_tab(*,j)=1.0 & $
  ;total h
  h_total=fltarr(Nbeta) & $
  for j=0,Npop-1 do h_total=h_total+contrib(j)*h_tab(*,j) & $
  window,xs=550,ys=850 & $
  !p.multi=[0,1,2,0] & $
  plot,beta_tab*r2d,h_total,thick=1,xstyle=1,yrange=[0,1], $
    xtitle='latitude !7b!3    (degrees)',ytitle='h(!7b!3)', $
    charsize=1.3 & $
  for j=0,Npop-1 do $
    oplot,beta_tab*r2d,contrib(j)*h_tab(*,j) & $

  ;Calculate phase-law(scattering angle)
  scat_angle_tab=!pi*findgen(Nbeta)/(Nbeta-1) & $
  phase_law_tab=fltarr(Nbeta) & $
  for j=0,Nbeta-1 do $
    phase_law_tab(j)=phase_law(scat_angle_tab(j), nu) & $
  phase_law_sin_tab=phase_law_tab*abs(sin(scat_angle_tab))^nu & $
  plot,scat_angle_tab*r2d,phase_law_tab,xstyle=1, $
    xtitle='scattering angle !7h!3    (degrees)', $
    ytitle='!7w!3(!7h!3)',charsize=1.3 & $
  oplot,scat_angle_tab*r2d,phase_law_sin_tab,linestyle=2 & $
  !p.multi=0 & $

  ;calculate the surface-brightness integral SB
  print,'Calculating the surface-brightness integral...' & $
  SB=fltarr(Ngrid,Ngrid,Npop) & $
  sin_epsilon=fltarr(Ngrid,Ngrid) & $
  for pop=0,Npop-1 do begin $
    for i=0,Ngrid-1 do begin $
      for j=0,Ngrid-1 do begin $
        ;longitude phi, latitude theta, and elongation epsilon
        phi=phi_grid(i) & $
        theta=theta_grid(j) & $
        epsilon=acos(cos(phi)*cos(theta)) & $
        if (epsilon lt 1.0e-6) then epsilon=0.1*delta_phi & $
        sin_epsilon(i,j)=sin(epsilon)^(nu+1) & $
        ;solar latitude beta versus scattering angle scat_angle
        k=where(scat_angle_tab ge epsilon) & $
        b=asin(sin(scat_angle_tab(k)-epsilon)* $
          abs(sin(theta))/sin(epsilon)) & $
        s=sort(b) & $
        ss=sort(s) & $
        h=spline(beta_tab,h_tab(*,pop),b(s)) & $
        hs=h(ss) & $
        integrand=phase_law_sin_tab(k)*hs & $
        SB(i,j,pop)=int_tabulated(scat_angle_tab(k),integrand) &$
      endfor & $
      print,float(i+1)*(pop+1)/Ngrid/Npop & $
    endfor & $
  endfor & $
  save,filename='fit.dat',SB,sin_epsilon & $
endif

;Calculate dust model in units of counts/sec/pixel
restore,filename='fit.dat'
if (pflag eq 1) then setplot 
if (pflag eq 0) then window,xs=Nx,ys=Ny
model=fltarr(Ngrid,Ngrid)
models=fltarr(Ngrid,Ngrid,Npop)
factor=asr*(Omega_sun/!pi)/cal/sin_epsilon
for pop=0,Npop-1 do begin $
  models(0,0,pop)=factor*SB(*,*,pop) & $
  model=model+contrib(pop)*models(*,*,pop) & $
endfor

;contour mosaic & model
thck=2*pflag+1
rng=30.0*[-1,1]
lvls=rotate(1000.0*(2.0^(-findgen(7))),2)
contour,model,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
contour,mosaic_fix,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
contour,smooth(mosaic_fix,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
sz=30.0
xarrow=14200.0
yarrow=12400.0
one_arrow,xarrow,yarrow,90.0,'N',color=1,thick=4, $
  ARROWSIZE=[sz*30.0,sz*9.0,35.0],CHARSIZE=1.3
one_arrow,xarrow,yarrow,180.0,'E',color=1,thick=4, $
  ARROWSIZE=[sz*30.0,sz*9.0,35.0],CHARSIZE=1.3
if (pflag eq 1) then output_plot,'fit.ps'
print,'Faintest contour = ',lvls(0)*cal,' solar brightness units'

;plot each population's surface brightness
if (pflag eq 0) then window,xs=350,ys=950
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=8.5,ysize=24.0,xoffset=7,yoffset=2, $
    /portrait,bits_per_pixel=8 & $
endif
!p.multi=[0,1,3,0]
lvls2=rotate(16000.0*(2.0^(-findgen(14))),2)
contour,factor*sb(*,*,1),phi_grid*r2d,theta_grid*r2d, $
  xstyle=1,ystyle=1,levels=lvls2,xrange=rng,yrange=rng, $
  thick=thck,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.6, $
  xthick=thck,ythick=thck,charthick=thck
xyouts,-4.0,23.0,'LOW',charsize=1.0,charthick=thck
contour,factor*sb(*,*,2),phi_grid*r2d,theta_grid*r2d, $
  xstyle=1,ystyle=1,levels=lvls2,xrange=rng,yrange=rng, $
  thick=thck,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.6, $
  xthick=thck,ythick=thck,charthick=thck
xyouts,-2.0,24.5,'HI',charsize=1.0,charthick=thck
contour,factor*sb(*,*,0),phi_grid*r2d,theta_grid*r2d, $
  xstyle=1,ystyle=1,levels=lvls2,xrange=rng,yrange=rng, $
  thick=thck,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.6, $
  xthick=thck,ythick=thck,charthick=thck
xyouts,-8.5,22.0,'ISOTROPIC',charsize=1.0,charthick=thck
!p.multi=0
if (pflag eq 1) then output_plot,'isophotes.ps'

stop

;generate pseudo-error array for mosaic-fix data
prng=3
err=fltarr(Nx,Ny)
for i=0,Nx-1 do begin $
  for j=0,Ny-1 do begin $
    sub_mosaic=mosaic_fix((i-prng)>0:(i+prng)<(Nx-1), $
      (j-prng)>0:(j+prng)<(Ny-1)) & $
    sub_good =good_fix((i-prng)>0:(i+prng)<(Nx-1), $
      (j-prng)>0:(j+prng)<(Ny-1)) & $
    k=where(sub_good ge 1) & $ 
    Nk=n_elements(k) & $
    if (Nk gt 1) then $
      err(i,j)=stdev(sub_mosaic(k))/sqrt(Nk-1) & $
  endfor & $
endfor
radial,rad,Nx/2,Ny/2,Nx,Ny
gd=where((good_fix ge 1) and (mosaic_fix gt lvls(1)) and $
  (err gt 0.0) and (rad lt 285))
t=bytarr(Nx,Ny)
t(gd)=1
tvscl,t

;Searching for best fitting model via chi-squared minimization
print,'Scanning parameter space for best-fitting model...'
Nparams=101
newsize=grid_factor*Ngrid
ast=findgen(Nparams)/(Nparams-1)
htc=findgen(Nparams)/(Nparams-1)
chi2=fltarr(Nparams,Nparams)
for i=0,Nparams-1 do begin $
  for j=0,Nparams-1 do begin $
    as=ast(i) & $
    ht=htc(j) & $
    if (as+ht le 1.0) then begin $
      oc=1.0-as-ht & $
      contrb=[oc,as,ht] & $
      model=fltarr(Ngrid,Ngrid) & $
      for pop=0,Npop-1 do $
        model=model+contrb(pop)*models(*,*,pop) & $
      model=rebin(model,newsize,newsize) & $
      c2=total(((model(gd)-mosaic_fix(gd))/err(gd))^2) & $
      c2=c2/(n_elements(gd)-2) & $
      chi2(i,j)=c2 & $
      print,as,ht,oc,as+ht+oc,c2 & $
    endif & $
  endfor & $
endfor
i=where(chi2 le 0.0)
if (i(0) ne -1) then chi2(i)=max(chi2)
j=where(chi2 eq min(chi2))
jj=j/Nparams
ii=j-jj*Nparams
ast_best=ast(ii)
htc_best=htc(jj)
ooc_best=1.0-ast_best-htc_best
contrib=[ ooc_best, ast_best, htc_best ]

;contour chi-squared for the available parameter space
window,xs=550,ys=550
contour,chi2<2.0*min(chi2),ast,htc,nlevels=10,charsize=1.3, $
  xtitle='asteroid fraction',ytitle='HTC fraction'
contour,chi2,ast,htc,nlevels=20,charsize=1.3,noerase=1, $
  xtitle='asteroid fraction',ytitle='HTC fraction'
plots,ast_best,htc_best,psym=1,thick=3,symsize=2
htc_range=htc
ast_range=fltarr(Nparams)
for j=0,Nparams-1 do begin $
  i=where(chi2(*,j) eq min(chi2(*,j))) & $
  ast_range(j)=htc(i) & $
  plots,ast_range(j),htc_range(j),psym=1 & $
endfor

;vary parameters
delta_jj=15
htc_min=htc_range(jj-delta_jj)
ast_min=ast_range(jj-delta_jj)
plots,ast_min,htc_min,psym=8,symsize=2
contrib_min=[ 1.0-ast_min-htc_min ,ast_min, htc_min ]
model_min=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model_min=model_min+contrib_min(pop)*models(*,*,pop)
model=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model=model+contrib(pop)*models(*,*,pop)
htc_max=htc_range(jj+delta_jj)
ast_max=ast_range(jj+delta_jj)
plots,ast_max,htc_max,psym=8,symsize=2
contrib_max=[ 1.0-ast_max-htc_max ,ast_max, htc_max ]
model_max=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model_max=model_max+contrib_max(pop)*models(*,*,pop)
window,xs=1200,ys=400
!p.multi=[0,3,1,0]
contour,model_min,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[3,3,1,0]
contour,mosaic_fix,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[3,3,1,0]
contour,smooth(mosaic_fix,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[2,3,1,0]
contour,model,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[2,3,1,0]
contour,mosaic_fix,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[2,3,1,0]
contour,smooth(mosaic_fix,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[1,3,1,0]
contour,model_max,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[1,3,1,0]
contour,mosaic_fix,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[1,3,1,0]
contour,smooth(mosaic_fix,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=0

;Print results
print,'Oort Cloud comet    contribution: ', $
  contrib(0),' +/- ', contrib_min(0)-contrib(0) & $
print,'Asteroid + JF Comet contribution: ', $
  contrib(1),' +/- ', contrib_min(1)-contrib(1) & $
print,'Halley Family comet contribution: ', $
  contrib(2),' +/- ', contrib_max(2)-contrib(2)

