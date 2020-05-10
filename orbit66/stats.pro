;stats.pro

;extract some quantities from orbit 66 results

;read and display image
restore, 'orbit66.dat'
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

;find elongation range
r2d= 180.0/!pi
phi_deg = phi[*, y_sun]*r2d
zl = avg_image[*, y_sun]
plot, phi_deg, zl
j = where(zl gt 0)
phi_min = min(phi_deg[j])
phi_max = max(phi_deg[j])
print, 'elongation range (degrees) = ', phi_min, phi_max

