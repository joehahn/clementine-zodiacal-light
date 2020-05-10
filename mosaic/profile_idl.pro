function profile_idl,array,xcenter,ycenter,r,err,tmplt

;use only those pixels where(tmplt eq 0)

sz=size(array)
radial,radius,xcenter,ycenter,sz(1),sz(2)
r_max=round(max(radius))+1
n_max=fix(4.1*sz(1))
r=fltarr(r_max)
n=fltarr(r_max)
p=fltarr(r_max,n_max)
for x=0,sz(1)-1 do $
  for y=0,sz(2)-1 do $
    if (tmplt(x,y) eq 0) then begin
      ri=round(radius(x,y))
      r(ri)=ri
      p(ri,n(ri))=array(x,y)
      n(ri)=n(ri)+1
    endif
pr=fltarr(r_max)
err=fltarr(r_max)
for ri=0,r_max-1 do begin
  if(n(ri) gt 0.0) then pr(ri)=avg(p(ri,0:n(ri)-1))
  ;if(n(ri) gt 1.0) then err(ri)=stdev(p(ri,0:n(ri)-1))/sqrt(n(ri))
  if(n(ri) gt 1.0) then begin 
    mom=moment(p(ri,0:n(ri)-1), sdev=stdev)
    err(ri)=stdev/sqrt(n(ri))
  endif
  ;print,ri,n(ri),reform(p(ri,0:n(ri)-1))
endfor

if ((xcenter eq fix(xcenter)) and (ycenter eq fix(ycenter)) and $
     (tmplt(xcenter,ycenter) eq 0)) then begin $
  r(1)=sqrt(2.0)
  arr=[array(xcenter-1,ycenter-1), array(xcenter+1,ycenter-1), $
       array(xcenter-1,ycenter+1), array(xcenter+1,ycenter+1)]
  pr(1)=avg(arr)
  err(1)=stdev(arr)/sqrt(4.0)
  r(0)=1.0
  arr=[array(xcenter,ycenter-1), array(xcenter,ycenter+1), $
       array(xcenter-1,ycenter), array(xcenter+1,ycenter)]
  pr(0)=avg(arr)
  err(0)=stdev(arr)/sqrt(4.0)
  r=[0,r]
  pr=[array(xcenter,ycenter),pr]
  err=[0.0,err]
endif

j=where((r ne 0.0) or (pr ne 0.0) or (err ne 0.0))
r=r(j)
pr=pr(j)
err=err(j)

return,pr
end
