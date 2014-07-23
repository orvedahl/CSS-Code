;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  weights_legendre.pro - This function returns the Legendre integration
;                         weights given an array of collocation points x.
;
;  usage:  weights = weights_legendre(x)
;       where weights = Legendre integration weights
;             x = array of collocation points, i.e. the zeroes of the Legendre
;                 function (of the first kind) of order nx (where nx =
;                 number of elements in the array x)
;
;  M.DeRosa - 13 Oct 2000 - cannibalized from spherical_transform.pro
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function weights_legendre,x

;  preliminaries
nx=n_elements(x)
weights=dblarr(nx)
costheta=double(x)
sintheta=sqrt(1-costheta^2)

;  set first two Legendre functions (evaluated at the collocation points)
Pm2=1
Pm1=x

;  iterate through the rest of the functions
for l=2,nx-1 do begin
  lr=1/double(l)
  P=(2-lr)*Pm1*costheta-(1-lr)*Pm2  ;  recursion relation Arfken 12.17a
  Pm2=Pm1
  Pm1=P
endfor

;  calculate dP_nx/dx evaluated at the collocation points
;    NOTE: P in the expression below is actually P_(nx-1), and in the
;    recursion relation listed below P_nx evaluated at the collocation points
;    is zero, by definition 
p_deriv=(nx*P)/sintheta^2  ;  recursion relation Arfken 12.26

;  calculate and then renormalize the weights
weights=2/(sintheta*p_deriv)^2
weights = weights * (2 * !dpi)

return,weights
end


