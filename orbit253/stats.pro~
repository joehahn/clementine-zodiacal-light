;stats.pro

;extract some quantities from orbit 206 results

;read and display image
restore, 'orbit206.dat'
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

;find elongation range along the up-down phi axis
r2d= 180.0/!pi
angle = phi
angle_deg = angle[x_sun, *]*r2d
zl = avg_image[x_sun, *]
plot, angle_deg, zl
j = where(zl gt 0)
angle_min = min(angle_deg[j])
angle_max = max(angle_deg[j])
print, 'elongation range (degrees) = ', angle_min, angle_max
