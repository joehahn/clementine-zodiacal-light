window,xs=1200,ys=800
!p.multi=[0,3,2,0]
lvls=rotate(1000.0*(2.0^(-findgen(7))),2)
contour,model_min,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1,$
  levels=lvls,color=165,xrange=rng,yrange=rng,charsize=2.0, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)', $
  xthick=thck,ythick=thck,charthick=thck, $
  title='!7v!3!u2!n='+strtrim(string(chi2(ast_min(0)* $
    (Nparams-1),htc_min(0)*(Nparams-1),kk(0))),2)
!p.multi=[0,3,2,0]
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1, $
  ystyle=1,noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[0,3,2,0]
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1,$
  ystyle=1,noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[5,3,2,0]
contour,model,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1, $
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck, $
  title='!7v!3!u2!n='+strtrim(string(min(chi2)))
!p.multi=[5,3,2,0]
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[5,3,2,0]
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1, $
  ystyle=1,noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[4,3,2,0]
contour,model_max,phi_grid*r2d,theta_grid*r2d,xstyle=1,ystyle=1,$
  levels=lvls,color=165,xrange=rng,yrange=rng, $
  thick=7,xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck, $
  title='!7v!3!u2!n='+strtrim(string(chi2(ast_max(0)* $
    (Nparams-1),htc_max(0)*(Nparams-1),kk(0))),2)
!p.multi=[4,3,2,0]
contour,mosaic_shift,x*r2d,y*r2d,xstyle=1,ystyle=1,noerase=1, $
  xrange=rng,yrange=rng, $
  levels=lvls(1:*),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck
!p.multi=[4,3,2,0]
contour,smooth(mosaic_shift,5),x*r2d,y*r2d,xstyle=1,ystyle=1, $
  noerase=1,xrange=rng,yrange=rng, $
  levels=lvls(0),thick=thck, $
  xtitle='ecliptic longitude !7u!3    (degrees)', $
  ytitle='ecliptic latitude !7h!3    (degrees)',charsize=2.0, $
  xthick=thck,ythick=thck,charthick=thck

