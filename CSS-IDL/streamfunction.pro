

function streamfunction, Vr, Vtheta, r, costheta, order=order, double=double

;------------------------------------------------------------
; This routine takes as input a divergenceless axisymmetric 
; vector field in spherical coordinates and computes from 
; it a streamfunction (a.k.a. a flux flunction).  The grid
; is decribed by r and costheta and can be non-uniform.
;------------------------------------------------------------
; INPUTS:
;
; Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
;              Dimensions are (N_Theta,N_R)
; r,costheta = the radius and cos(colatitude) of the grid.
;              r is assumed to vary from rmax to rmin and 
;              costheta from  1 to -1 (i.e. 90 degrees
;              to -90 degrees in latitude).
;              Dimensions are r(N_R), costheta(N_Theta)
; order      = If greater than zero, integration begins at the
;              outer shell and the north pole and proceeds
;              inward and southward.  If less than zero,
;              integration begins at the inner shell and 
;              south pole and proceeds upward and northward.
;              If equal to zero, both are done and an average
;              is taken.
;------------------------------------------------------------
; OUTPUTS:
;
; psi = the streamfunction
;------------------------------------------------------------

if n_elements(order) ne 1 then order = 0

sz=size(Vr) & N_Theta=sz(1) & N_R=sz(2)

theta = acos(costheta)
sintheta = sqrt(1.d0 - costheta^2)

dpsi_dr = fltarr(N_Theta,N_R) & dpsi_dtheta = dpsi_dr

for i=0,N_Theta-1 do begin
    dpsi_dr(i,*)     = - r*sintheta(i) * Vtheta(i,*)
    dpsi_dtheta(i,*) = r*r*sintheta(i) * Vr(i,*)
endfor

dtheta = dblarr(N_Theta) & dr = dblarr(N_R)

nt = N_Theta-1 & nr = N_R-1

if order ge 0 then begin
    ; double precision accumulation
    psi=dblarr(N_Theta,N_R)

    dtheta(1:nt) = theta(1:nt) - theta(0:nt-1)
    dr(1:nr) = r(1:nr) - r(0:nr-1)

    dtheta(0)=0 & dr(0)=0
    
    for i=1,nr do psi(1:nt,i) = psi(1:nt,i-1) + dpsi_dr(1:nt,i)*dr(i) 
    for i=1,nt do psi(i,1:nr) = psi(i-1,1:nr) + dpsi_dtheta(i,1:nr)*dtheta(i) 

endif 

if order le 0 then begin           
    psi2=dblarr(N_Theta,N_R)

    dtheta(0:nt-1) = theta(0:nt-1) - theta(1:nt)
    dr(0:nr-1) = r(0:nr-1) - r(1:Nr)

    dtheta(nt)=0 & dr(nr)=0

    for i=nr-1,0,-1 do psi2(0:nt-1,i) = psi2(0:nt-1,i+1) + dpsi_dr(0:nt-1,i)*dr(i)
    for i=nt-1,0,-1 do psi2(i,0:nr-1) = psi2(i+1,0:nr-1) + dpsi_dtheta(i,0:nr-1)*dtheta(i)

    if order lt 0 then begin
       if keyword_set(double) then return,double(psi2) else return,float(psi2)
    endif else begin
       psi = 0.5d0*(psi+psi2)
       print,'Doing both orders on streamfunction integration'
    endelse
endif

if keyword_set(double) then return,double(psi) else return,float(psi)


end
 


