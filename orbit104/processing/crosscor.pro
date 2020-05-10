function crosscor,a,b

;returns the cross-correlation of image a and b.
;note:reference point is at (s(1)/2,s(2)/2)

result=abs(fft((fft(a,-1)*conj(fft(b,-1))),1))    
s=size(a)
return,shift(result,s(1)/2,s(2)/2)
end                       
