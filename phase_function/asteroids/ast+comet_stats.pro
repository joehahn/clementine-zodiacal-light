;ast+comet_stats.pro

pflag=1

;ASTEROIDS
;read Ted Bowell's list of asteroid orbit elements
read_astorb,number,name,Ha,diameter,aa,ea,ia,Oa,wa,Ma;,data=data
;select asteroids objects brighter than H=15
p=where((aa lt 6.0) and (Ha lt 15.0))
;plot histogram of asteroid inclinations
histo,ia(p),100,i_axis,Ni
thck=2*pflag+1
if (pflag ne 1) then window,xs=550,ys=950
if (pflag eq 1) then begin $
  set_plot,'ps' & $
  device,xsize=17.5,ysize=24.0,xoffset=2,yoffset=2, $
    /portrait,bits_per_pixel=8 & $
endif
!p.multi=[0,1,4,0]
plot,i_axis,Ni,psym=10,xrange=[0,90],yrange=[0,4500],ystyle=1, $
  charsize=1.6,xtitle='inclination i    (degrees)',thick=thck, $
  ytitle='N(i)',xstyle=1,xticks=9, $
  xthick=thck,ythick=thck,charthick=thck
xyouts,54,3650,'Asteroids',charsize=1.3,charthick=thck
sigma=1.05*stdev(ia(p))
r2d=180.0/!pi
coeff=0.0105*n_elements(p)/(sigma/r2d)^2
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'asteroid Sigma = ',sigma
pp=where(ia(p) gt 20.0)
print,float(n_elements(pp))/float(n_elements(p))

;JFCs
q_max=2.5
;read the MPC's list of comet orbit elements
read_comets,name,q,e,i,O,w
N=n_elements(q)
a=fltarr(N)
j=where(e lt 1.0)
a(j)=q(j)/(1.0-e(j))
k=where(e ge 1.0)
a(k)=-1.0
;Tisserand parameter
aj=5.2028
T=aj*(1.0-e)/q+2.0*sqrt(q*(1.0+e)/aj)*cos(i/r2d)
;select JFC's with T>2 and q<2.5 AU
p_jfc=where((T ge 2.0) and (q le q_max))
;plot histogram of JFC inclinations
Nbins=18
bs=90.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i(p_jfc),min=0.0,max=90.0,binsize=bs)
plot,i_axis,Ni,psym=10,xrange=[0,90],yrange=[0,50],ystyle=1, $
  charsize=1.6,xtitle='inclination i    (degrees)',thick=thck, $
  ytitle='N(i)',xthick=thck,ythick=thck,charthick=thck,xticks=9,$
  xstyle=1,xminor=5
xyouts,54,41,'Jupiter-Family Comets',charsize=1.3,charthick=thck
sigma=8.0
r2d=180.0/!pi
coeff=1.4*n_elements(p_jfc)/n_elements(i_axis)/(sigma/r2d)^2
i_axis=0.5*findgen(181.0)
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'JFC Sigma = ',sigma

;HTCs
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
  charsize=1.6,xtitle='inclination i    (degrees)',xstyle=1, $
  ytitle='N(i)',xticks=9,xthick=thck,ythick=thck,charthick=thck,$
  thick=thck,xminor=5
xyouts,54,5.6,'Halley-Type Comets',charsize=1.3,charthick=thck
sigma=33.0
coeff=4.3*n_elements(p_htc)/n_elements(i_axis)/(sigma/r2d)^2
i_axis=findgen(181.0)
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'HTC Sigma = ',sigma

;OCCs
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
  charsize=1.6,xtitle='inclination i    (degrees)',xstyle=1, $
  ytitle='N(i)',xticks=9,xthick=thck,ythick=thck,charthick=thck,$
  thick=thck,xminor=5
xyouts,5,20,'Oort Cloud Comets',charsize=1.3,charthick=thck
i_axis=0.5*findgen(181.0)
Ni_exp=1.3*(2.0/!pi)*n_elements(p_htc)*sin(i_axis/r2d)
oplot,i_axis,Ni_exp,thick=thck
!p.multi=0




if (pflag eq 1) then output_plot,'ast+comet_i_dist.ps'
