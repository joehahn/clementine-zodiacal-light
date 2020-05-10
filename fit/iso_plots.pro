;iso_plots.pro

!p.multi=[3,3,2,0]
contour,model_min,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1,$
  levels=lvls,color=165,xrange=rng,yrange=rng,charsize=2.0, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)', $
  xthick=thck,ythick=thck,charthick=thck, $
  title='!7v!3!u2!n='+strtrim(string(chi2_mn),2)
!p.multi=[3,3,2,0]
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1, $
  ystyle=1,noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[3,3,2,0]
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1,$
  ystyle=1,noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[2,3,2,0]
contour,model,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck, $
  title='!7v!3!u2!n='+strtrim(string(min(chi2)))
!p.multi=[2,3,2,0]
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[2,3,2,0]
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1, $
  ystyle=1,noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[1,3,2,0]
contour,model_max,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1,$
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck, $
  title='!7v!3!u2!n='+strtrim(string(chi2_mx),2)
!p.multi=[1,3,2,0]
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[1,3,2,0]
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=0
