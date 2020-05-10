function destreak, img, exp_time

;Herb Zook's destreaking algorithm. Inputs are the input
;image array img, and the exposure time exp_time in seconds.
;This function returns the destreaked array.

;do the math in double precision, but return a
;single-precision result
img=double(img)
exp_time=double(exp_time)
sz=size(img)
Nx=sz(1)
Ny=sz(2)

;pixel transfer time in seconds
tau=94.4d-6

;destreak
e=tau/exp_time
den=(1.d0+double(Ny-1)*e)*(1.d0-e)
ep=e/den
f=(1.d0+double(Ny-2)*e)/den
A=(f+ep)*img
for x=0,Nx-1 do A(x,0)=A(x,*)-ep*total(img(x,*))

return,float(A)
end
