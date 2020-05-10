;fit_nu_occ=2.5.pro

;dust-model parameters
nu=1.45			;dust radial power law for phase function
asr=9.5e-9		;albedo*sigma_1*r_1, +/- 10%
Omega_sun=6.8e-5	;Sun solid angle, steradians

;calibration converts 1 count/sec/pixel into solar brightness
cal=7.45e-14

;sigma's for three different dust populations
sigma_pop = [0.000, 0.122, 0.576]

;ecliptic contributions by the three different dust populations
contrib   = [0.025, 0.375, 0.625]
min_chi2=568.446

;each population's radial power-law
nu_pop    = [2.50, 1.45, 1.45]

;Calculate new surface-brightness integrals when <> 0
calc_SB=1

;generate postscript output if pflag=1
pflag=1

;store results in this file
io_file='fit_nu_occ=2.5'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;restore observed image
if (keyword_set(mosaic) eq 0) then begin $
  ;rebin mosaic 
  grid_factor=5 & $
  restore, filename= $
    '/net/lpiip2/users/hahn/clementine/zodiacal_light/mosaic/mercator.dat'&$
  ;fix inner gap in mosaic & high-noise areas
  sz=size(mosaic) & $
  Nx=sz(1) & $
  Ny=sz(2) & $
  radial,rad,x_sun,y_sun,Nx,Ny & $
  j=where((rad le 150) and (mosaic eq 0.0)) & $
  mosaic(j)=max(mosaic) & $
  mosaic(610:*,560:*)=0 & $
  mosaic(630:*,0:110)=0 & $
  Nxs=800 & $
  Nys=800 & $
  mosaic_shift=fltarr(Nxs,Nys) & $
  good_shift=fltarr(Nxs,Nys) & $
  dx=Nxs/2-round(x_sun)+1 & $
  dy=Nys/2-round(y_sun) & $
  mosaic_shift(dx,dy)=mosaic & $
  good_shift(dx,dy)=good & $
  x_sun=Nxs/2 & $
  y_sun=Nys/2 & $
  x=ps_avg*(findgen(Nxs)-x_sun) & $
  y=ps_avg*(findgen(Nys)-y_sun) & $
  ;rebin mosaic and generate error array
  Nx=Nxs/grid_factor & $
  Ny=Nys/grid_factor & $
  mosaic_fix=fltarr(Nx,Ny) & $
  good_fix  =bytarr(Nx,Ny) & $
  err=fltarr(Nx,Ny) & $
  for i=0,Nx-1 do begin $
    for j=0,Ny-1 do begin $
      sub_mosaic=mosaic_shift( $
        (i-1)*grid_factor>0:(i+1)*grid_factor<(Nxs-1),   $
        (j-1)*grid_factor>0:(j+1)*grid_factor<(Nys-1)) & $
      sub_good =good_shift( $
        (i-1)*grid_factor>0:(i+1)*grid_factor<(Nxs-1),   $
        (j-1)*grid_factor>0:(j+1)*grid_factor<(Nys-1)) & $
      k=where(sub_good ge 1) & $ 
      Nk=n_elements(k) & $
      if (Nk gt 2*grid_factor^2) then begin $
        mosaic_fix(i,j)=avg(sub_mosaic(k)) & $
        good_fix(i,j)=1 & $
        err(i,j)=stdev(sub_mosaic(k))/sqrt(Nk-1) & $
      endif & $
    endfor & $
  endfor & $
  radial,rad,Nx/2,Ny/2,Nx,Ny & $
  gd=where((good_fix ge 1) and (mosaic_fix gt 31.25) and $
    (err gt 0.0) and (rad lt 325/grid_factor)) & $
  t=bytarr(Nx,Ny) & $
  t(gd)=1 & $
  ;zap star and earthshine-polluted regions.
  t(43:44,106:107)=0 & $
  t(82:110,90:*)=0 & $
  gd=where(t gt 0) & $
endif
window,xs=Nxs,ys=Nys
tvscl,alog(mosaic_shift>1)
tvscl,enlarge(alog(mosaic_fix*t>1),grid_factor)
tvscl,enlarge(alog(err*t>0.01),grid_factor)

;sky-coordinates (phi,theta)=(longitude,latitude)
Ngrid=Nx
phi_grid=max(x)*2.0*(findgen(Ngrid)/(Ngrid-1)-0.5)
theta_grid=phi_grid
Npop=n_elements(sigma_pop)
common share, beta, sigma
.run latitude_integrand.pro
r2d=180.0/!pi
Po2=!pi/2.0

if (calc_SB ne 0) then begin $

  ;Calculate latitude distribution h(beta) for diff' populations.
  Nbeta=101 & $;201 & $
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
  plot,scat_angle_tab*r2d,phase_law_tab,xstyle=1, $
    xtitle='scattering angle !7h!3    (degrees)', $
    ytitle='!7w!3(!7h!3)',charsize=1.3 & $
  oplot,scat_angle_tab*r2d, $
    phase_law_tab*abs(sin(scat_angle_tab))^nu,linestyle=2 & $
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
        sin_epsilon(i,j)=abs(sin(epsilon)) & $
        ;solar latitude beta versus scattering angle scat_angle
        k=where(scat_angle_tab ge epsilon) & $
        b=asin(sin(scat_angle_tab(k)-epsilon)* $
          abs(sin(theta))/sin(epsilon)) & $
        s=sort(b) & $
        ss=sort(s) & $
        h=spline(beta_tab,h_tab(*,pop),b(s)) & $
        hs=h(ss) & $
        integrand=phase_law_tab(k)*hs* $
          abs(sin(scat_angle_tab(k)))^nu_pop(pop) & $
        SB(i,j,pop)=int_tabulated(scat_angle_tab(k),integrand) &$
      endfor & $
      print,float(i+1)*(pop+1)/Ngrid/Npop & $
    endfor & $
  endfor & $
  save,filename=io_file+'.dat',SB,sin_epsilon & $
endif

;Calculate dust model in units of counts/sec/pixel
restore,filename=io_file+'.dat'
model=fltarr(Ngrid,Ngrid)
models=fltarr(Ngrid,Ngrid,Npop)
for pop=0,Npop-1 do begin $
  factor=asr*(Omega_sun/!pi)/cal/sin_epsilon^(nu_pop(pop)+1.0) &$
  models(0,0,pop)=factor*SB(*,*,pop) & $
  model=model+contrib(pop)*models(*,*,pop) & $
endfor

;plot each population's surface brightness
if (pflag eq 0) then window,xs=350,ys=950
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=8.5,ysize=24.0,xoffset=7,yoffset=2, $
    /portrait,bits_per_pixel=8 & $
endif
!p.multi=[0,1,3,0]
lvls2=rotate(16000.0*(2.0^(-findgen(14))),2)
rng=30.0*[-1,1]
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

;Searching for best fitting model via chi-squared minimization
print,'Scanning parameter space for best-fitting model...'
Nparams=41;101
ast=findgen(Nparams)/(Nparams-1)
htc=findgen(Nparams)/(Nparams-1)
occ=findgen(Nparams)/(Nparams-1)
chi2=fltarr(Nparams,Nparams,Nparams)
for i=0,Nparams-1 do begin $
  for j=0,Nparams-1 do begin $
    for k=0,Nparams-1 do begin $
      as=ast(i) & $
      ht=htc(j) & $
      oc=occ(k) & $
      contrb=[oc,as,ht] & $
      model=fltarr(Ngrid,Ngrid) & $
      for pop=0,Npop-1 do $
        model=model+contrb(pop)*models(*,*,pop) & $
      c2=total(((model(gd)-mosaic_fix(gd))/err(gd))^2) & $
      c2=c2/(n_elements(gd)-2) & $
      chi2(i,j,k)=c2 & $
    endfor & $
    print,as,ht,oc,c2 & $
  endfor & $
endfor
l=where(chi2 eq min(chi2))
kk=l/Nparams^2
jj=(l-kk*Nparams^2)/Nparams
ii=l-kk*Nparams^2-jj*Nparams
ast_best=ast(ii)
htc_best=htc(jj)
occ_best=occ(kk)
contrib=[ occ_best, ast_best, htc_best ]
min_chi2=min(chi2)
print,'Oort Cloud comet    contribution: ',contrib(0) & $
print,'Asteroid + JF Comet contribution: ',contrib(1) & $
print,'Halley Family comet contribution: ',contrib(2) & $
print,'Minimum chi^2 = ',min_chi2
;contour chi-squared for the available parameter space
window,xs=550,ys=550
z_min=2.0*min_chi2
contour,chi2(*,*,kk)<z_min,ast,htc,nlevels=10, $
  charsize=1.3,xtitle='asteroid fraction',ytitle='HTC fraction'
contour,chi2(*,*,kk)>z_min,ast,htc,nlevels=20,charsize=1.3, $
  noerase=1,xtitle='asteroid fraction',ytitle='HTC fraction'
plots,ast_best,htc_best,psym=1,thick=3,symsize=2
ast_range=ast
htc_range=fltarr(Nparams)
chi2_range=fltarr(Nparams)
for i=0,Nparams-1 do begin $
  j=where(chi2(i,*,kk) eq min(chi2(i,*,kk))) & $
  htc_range(i)=ast(j) & $
  chi2_range(i)=chi2(i,j,kk) & $
endfor
oplot,ast_range,htc_range,psym=1

;get range of low & high-inclination parameters
chi2_max=min_chi2+110.0
j=where(chi2_range le chi2_max)
htc_min=htc_range(j(0))
ast_min=ast_range(j(0))
plots,ast_min,htc_min,psym=8,symsize=2
contrib_min=[ occ_best, ast_min, htc_min ]
model_min=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model_min=model_min+contrib_min(pop)*models(*,*,pop)
model=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model=model+contrib(pop)*models(*,*,pop)
htc_max=htc_range(j(n_elements(j)-1))
ast_max=ast_range(j(n_elements(j)-1))
plots,ast_max,htc_max,psym=8,symsize=2
contrib_max=[ occ_best, ast_max, htc_max ]
model_max=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model_max=model_max+contrib_max(pop)*models(*,*,pop)
@lo_hi_plots
print,'Asteroid + JF Comet contribution: ', $
  contrib(1),' +/- ', 0.5*(contrib_max(1)-contrib_min(1)) & $
print,'Halley Family comet contribution: ', $
  contrib(2),' +/- ', 0.5*(contrib_min(2)-contrib_max(2))

;get range of OCC parameter
for k=0,Nparams-1 do chi2_range(k)=min(chi2(*,*,k))
k=where(chi2_range le chi2_max)
chi2_mn=min(chi2(*,*,k(0)))
l=where(chi2(*,*,k(0)) eq chi2_mn)
j=l/Nparams
i=l-j*Nparams
ast_min=ast(i)
htc_min=htc(j)
occ_min=occ(k(0))
contrib_min=[ occ_min ,ast_min, htc_min ]
model_min=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model_min=model_min+contrib_min(pop)*models(*,*,pop)
k=k(n_elements(k)-1)
chi2_mx=min(chi2(*,*,k))
l=where(chi2(*,*,k) eq chi2_mx)
j=l/Nparams
i=l-j*Nparams
ast_max=ast(i)
htc_max=htc(j)
occ_max=occ(k)
contrib_max=[ occ_max ,ast_max, htc_max ]
model_max=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model_max=model_max+contrib_max(pop)*models(*,*,pop)
@iso_plots
print,'Oort Cloud comet    contribution: ', $
  contrib(0),' +/- ', 0.5*(contrib_max(0)-contrib_min(0)) & $
print,'Minimum chi^2 = ',min_chi2

;contour mosaic & model
model=fltarr(Ngrid,Ngrid)
for pop=0,Npop-1 do $
  model=model+contrib(pop)*models(*,*,pop)
if (pflag eq 0) then window,xs=550,ys=550
if (pflag eq 1) then setplot
contour,model,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=1.3, $
  xthick=thck,ythick=thck,charthick=thck
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
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
if (pflag eq 1) then output_plot,io_file+'.ps'
print,'Faintest contour = ',lvls(0)*cal,' solar brightness units'
