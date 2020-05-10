;stats.pro

;extract some quantities from orbit 164 results

;read and display image
restore, 'orbit164.dat'
sz = size(avg_image)
Nx = sz[1]
Ny = sz[2]
window, xs=Nx, ys=Ny, retain=2
im = alog(avg_image>10.0)
tvscl, im

;find Sun
r = sqrt( phi^2 + theta^2)
mn = min(r, idx)
xy_sun = array_indices(im, idx)
x_sun = xy_sun[0]
y_sun = xy_sun[1]
im[x_sun, y_sun] = max(im)
tvscl, im

;;find elongation range along left-right axis
;r2d= 180.0/!pi
;phi_deg = phi[*, y_sun]*r2d
;zl_lr = avg_image[*, y_sun]
;plot, theta_deg, zl_lr
;j = where(zl_lr gt 0)
;phi_min = min(phi_deg[j])
;phi_max = max(phi_deg[j])
;print, 'left-right  elongation range (degrees) = ', phi_min, phi_max

;find elongation range along up-down axis
r2d= 180.0/!pi
theta_deg = theta[x_sun,*]*r2d
zl_ud = avg_image[x_sun, *]
plot, theta_deg, zl_ud
j = where(zl_ud gt 0)
theta_min = min(theta_deg[j])
theta_max = max(theta_deg[j])
print, 'up-down elongation range (degrees) = ', theta_min, theta_max

