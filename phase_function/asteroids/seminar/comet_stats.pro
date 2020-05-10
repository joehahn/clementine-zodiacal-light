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

;Tisserand parameter---compare to Fig. 8.
aj=5.2028
r2d=180.0/!pi
T=aj*(1.0-e)/q+2.0*sqrt(q*(1.0+e)/aj)*cos(i/r2d)

;select JFC's with T>2 and q<2.5 AU
p_jfc=where((T ge 2.0) and (q le q_max))
histo,[min(q),q(p_jfc),max(q)],30,q_axis,Nq
Nq(0)=Nq(0)-1
Nq(n_elements(Nq)-1)=Nq(n_elements(Nq)-1)-1

;plot histogram of JFC inclinations
if (pflag ne 1) then window,xs=550,ys=550
if (pflag eq 1) then setplot
thck=1+2*pflag
Nbins=18
bs=90.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i(p_jfc),min=0.0,max=90.0,binsize=bs)
plot,i_axis,Ni,psym=10,xrange=[0,50],yrange=[0,50],ystyle=1, $
  charsize=1.3,xtitle='inclination i    (degrees)',thick=thck, $
  ytitle='N(i)',xthick=thck,ythick=thck,charthick=thck,xticks=5,$
  xstyle=1,xminor=5
xyouts,14,44,'Jupiter-Family Comets',charsize=2.0, $
  charthick=thck+2
xyouts,29,39.5,'!7r!3=8!9%!3',charsize=2.5,charthick=thck+2
sigma=8.0
r2d=180.0/!pi
coeff=1.4*n_elements(p_jfc)/n_elements(i_axis)/(sigma/r2d)^2
i_axis=0.5*findgen(181.0)
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'JFC Sigma = ',sigma
if (pflag eq 1) then output_plot,'JFCs.ps'

;select HTC's with T<2 and q<2.5 AU and a<35 AU
p_htc=where((T lt 2.0) and (q le q_max) and (a le 35.0) $
  and (a gt 0.0))
Nbins=15
bs=180.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i(p_htc),min=0.0,max=180.0,binsize=bs)
if (pflag eq 1) then setplot
plot,i_axis,Ni,psym=10,xrange=[0,90],yrange=[0,9],ystyle=1, $
  charsize=1.3,xtitle='inclination i    (degrees)',xstyle=1, $
  ytitle='N(i)',xticks=9,xthick=thck,ythick=thck,charthick=thck,$
  thick=thck,xminor=5
xyouts,19,8.0,'Halley-Type Comets',charsize=2.5,charthick=thck+2
xyouts,45,7.0,'!7r!3=33!9%!3',charsize=2.5,charthick=thck+2
sigma=33.0
coeff=3.3*n_elements(p_htc)/n_elements(i_axis)/(sigma/r2d)^2
i_axis=findgen(181.0)
Ni_exp=coeff*sin(i_axis/r2d)*exp(-0.5*(i_axis/sigma)^2)
oplot,i_axis,Ni_exp,thick=thck
print,'HTC Sigma = ',sigma
if (pflag eq 1) then output_plot,'HTCs.ps'

;select OCC's with T<2 and q<2.5 AU and a>35 AU
p_occ=where((T lt 2.0) and (q le q_max) and (a gt 35.0))
i_occ=i(p_occ)
j=where(i_occ gt 90.0)
if (j(0) ne -1) then i_occ(j)=180.0-i_occ(j)
Nbins=25
bs=90.0/(Nbins-1)
i_axis=bs*findgen(Nbins)
Ni=histogram(i_occ,min=0.0,max=90.0,binsize=bs)
if (pflag eq 1) then setplot
plot,i_axis,Ni,psym=10,xrange=[0,90],yrange=[0,25],ystyle=1, $
  charsize=1.3,xtitle='inclination i    (degrees)',xstyle=1, $
  ytitle='N(i)',xticks=9,xthick=thck,ythick=thck,charthick=thck,$
  thick=thck,xminor=5
xyouts,5,22.5,'Oort Cloud Comets',charsize=2.5,charthick=thck+2
xyouts,11.5,20.0,'g(i)!9c!3sin(i)',charsize=2.5,charthick=thck+2
i_axis=0.5*findgen(181.0)
Ni_exp=1.3*(2.0/!pi)*n_elements(p_htc)*sin(i_axis/r2d)
oplot,i_axis,Ni_exp,thick=thck
!p.multi=0
if (pflag eq 1) then output_plot,'OCCs.ps'
