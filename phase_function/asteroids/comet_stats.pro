;comet_stats.pro

pflag=1
q_max=2.5

;read the MPC's list of comet orbit elements
read_comets,name,q,e,i,O,w

;plot i vs a and compare to Levison & Duncan (1994) Fig. 6.
N=n_elements(q)
a=fltarr(N)
j=where(e lt 1.0)
a(j)=q(j)/(1.0-e(j))
k=where(e ge 1.0)
a(k)=-1.0
window,xs=550,ys=550
plot_oi,a(j),i(j),psym=8,xrange=[1,100],yrange=[0,180],ystyle=1,$
  charsize=1.3,xtitle='a    (AU)',ytitle='inclination i    (deg)'
plots,[7.4,7.4],[0,180],linestyle=1

;Tisserand parameter---compare to Fig. 8.
aj=5.2028
r2d=180.0/!pi
T=aj*(1.0-e)/q+2.0*sqrt(q*(1.0+e)/aj)*cos(i/r2d)
plot_oi,a(j),T(j),psym=8,xrange=[1,1.e2],charsize=1.3, $
  xtitle='a    (AU)',ytitle='Tisserand T'

;ALL comets' perihelion distribution
histo,q,50,q_axis,Nq
plot,q_axis,Nq,psym=10,thick=3,xtitle='perihelion q    (AU)', $
  ytitle='N(q)',charsize=1.3

;select JFC's with T>2 and q<2.5 AU
p_jfc=where((T ge 2.0) and (q le q_max))
histo,[min(q),q(p_jfc),max(q)],30,q_axis,Nq
Nq(0)=Nq(0)-1
Nq(n_elements(Nq)-1)=Nq(n_elements(Nq)-1)-1
oplot,q_axis,Nq,psym=10,thick=3,color=128

;plot histogram of JFC inclinations
if (pflag ne 1) then window,xs=550,ys=900
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=17.5,ysize=24.0,xoffset=2,yoffset=2, $
    /portrait,bits_per_pixel=8 & $
endif
thck=1+2*pflag
!p.multi=[0,1,3,0]
Nbins=18
bs=90.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i(p_jfc),min=0.0,max=90.0,binsize=bs)
plot,i_axis,Ni,psym=10,xrange=[0,50],yrange=[0,50],ystyle=1, $
  charsize=2.0,xtitle='inclination i    (degrees)',thick=thck, $
  ytitle='N(i)',xthick=thck,ythick=thck,charthick=thck,xticks=5,$
  xstyle=1,xminor=5
xyouts,28,41,'Jupiter-Family Comets',charsize=1.3,charthick=thck
sigma=8.0
r2d=180.0/!pi
coeff=1.4*n_elements(p_jfc)/n_elements(i_axis)/(sigma/r2d)^2
i_axis=0.5*findgen(181.0)
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'JFC Sigma = ',sigma

;select HTC's with T<2 and q<2.5 AU and a<35 AU
p_htc=where((T lt 2.0) and (q le q_max) and (a le 35.0) $
  and (a gt 0.0))
i_htc=i(p_htc)
j=where(i_htc gt 90.0)
if (j(0) ne -1) then i_htc(j)=180.0-i_htc(j)
Nbins=20
bs=180.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i_htc,min=0.0,max=180.0,binsize=bs)
plot,i_axis,Ni,psym=10,xrange=[0,90],yrange=[0,7],ystyle=1, $
  charsize=2.0,xtitle='inclination i    (degrees)',xstyle=1, $
  ytitle='N(i)',xticks=9,xthick=thck,ythick=thck,charthick=thck,$
  thick=thck,xminor=5
xyouts,54,5.6,'Halley-Type Comets',charsize=1.3,charthick=thck
sigma=33.0
coeff=4.3*n_elements(p_htc)/n_elements(i_axis)/(sigma/r2d)^2
i_axis=findgen(181.0)
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'HTC Sigma = ',sigma

;select OCC's with T<2 and q<2.5 AU and a>35 AU
p_occ=where((T lt 2.0) and (q le q_max) and (a gt 35.0))
i_occ=i(p_occ)
j=where(i_occ gt 90.0)
if (j(0) ne -1) then i_occ(j)=180.0-i_occ(j)
Nbins=25
bs=90.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i_occ,min=0.0,max=90.0,binsize=bs)
plot,i_axis,Ni,psym=10,xrange=[0,90],yrange=[0,25],ystyle=1, $
  charsize=2.0,xtitle='inclination i    (degrees)',xstyle=1, $
  ytitle='N(i)',xticks=9,xthick=thck,ythick=thck,charthick=thck,$
  thick=thck,xminor=5
xyouts,3,21,'Oort Cloud Comets',charsize=1.3,charthick=thck
i_axis=0.5*findgen(181.0)
Ni_exp=1.3*(2.0/!pi)*n_elements(p_htc)*sin(i_axis/r2d)
oplot,i_axis,Ni_exp,thick=thck
!p.multi=0
if (pflag eq 1) then output_plot,'comet_i_dist.ps'

;plot q vs i for various populations
window,xs=550,ys=550
plot,q,i,psym=8,xrange=[0,6],yrange=[0,180],ystyle=1,nodata=1, $
  xtitle='perihelion q    (AU)',charsize=1.3, $
  ytitle='inclination i    (degrees)'
oplot,q(p_jfc),i(p_jfc),psym=6
oplot,q(p_htc),i(p_htc),psym=4
oplot,q(p_occ),i(p_occ),psym=3
j=where(q gt q_max)
oplot,q(j),i(j),psym=8

;histogram 1/a for OCCs
x=(1.0-e)/q
p_lpc=where((a gt 34.2) or (a lt 0.0))
Nbins=200
bs=0.03/(Nbins-1)
a_axis=bs*findgen(Nbins)
Ni=histogram(x(p_lpc),min=0.0,max=0.03,binsize=bs)
plot,a_axis,Ni,psym=10,xrange=[0,0.005],yrange=[0,40], $
  charsize=1.3,xtitle='1/a    (AU!u-1!n)',xstyle=1,ytitle='N(i)'

stop

;cumulative inclination distribution & compare
;to Fig. 1 of Levison, Dones, & Duncan (2001)
Nbins=181
i_axis=findgen(Nbins)
f_htc=fltarr(Nbins)
N_htc=n_elements(p_htc)
i_htc=i(p_htc)
for j=0,Nbins-1 do begin $
  k=where(i_htc le i_axis(j)) & $
  if (k(0) ne -1) then f_htc(j)=n_elements(k) & $
endfor
f_htc=f_htc/n_elements(p_htc)
plot,i_axis,f_htc,psym=10,xrange=[0,180],xstyle=1
sigma=33.0
N_htc_exp=sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
f_htc_exp=fltarr(Nbins)
for j=0,Nbins-1 do f_htc_exp(j)=total(N_htc_exp(0:j))
f_htc_exp=f_htc_exp/f_htc_exp(Nbins-1)
oplot,i_axis,f_htc_exp,color=128


