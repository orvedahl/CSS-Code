;+^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  linrange.pro - creates a double precision array of n numbers starting
;                 with min and ending with max, in a linear fashion.
;
;  usage:  result=linrange(n,max,min)
;          where result = double-precision array with n elements
;                     n = number of points in the array
;                   max = highest number
;                   min = lowest number
;
;  notes:  returns -1 if error detected
;
;  M.DeRosa - 12 Apr 1995 - created
;             30 Nov 1998 - added usage message
;
;-^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function linrange,n,min,max

;  usage message
if n_params() lt 3 then begin
  print,'  result=linrange(n,min,max)'
  return,-1
endif

np=long(n)
mn=double(min)
mx=double(max)

vec=dindgen(np)
vec=vec*(mx-mn)/(np-1)
vec=vec+mn

return,vec
end

