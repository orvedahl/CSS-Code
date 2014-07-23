; idem sphe...MDrosa.pro but easier to call
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  spherical_transform.pro - This routine performs a spherical harmonic
;                            transform on a 2-D array.
;
;  usage:  B = spherical_transform(A,cp,lmax=lmax,period=period)
;       where B(lmax,lmax) = transformed array ordered (l,m)
;             A(n_phi,n_theta) = array to be transformed ordered (phi,theta)
;             cp = cosine of theta collocation points for theta grid
;             lmax = maximum l in expansion (default is (2*n_theta-1)/3)
;             period = periodicity factor in phi
;
;  notes: - All calculations are done in double precision.
;
;  M.Miesch - 15 Oct 1997 - acquired from Mark
;  M.DeRosa - 25 Oct 1999 - modified slightly (basic algorithm is in general
;                           the same as Mark's original version)
;             27 Oct 1999 - added period keyword
;             12 Sep 2000 - added weights keyword
;             12 Sep 2000 - changed normalization to something more sensible
;             14 Sep 2000 - fixed bug in m,l loop (near end of routine)
;             13 Oct 2000 - uses weights_legendre to compute theta integration
;                           weights, weights keyword now obsolete
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@weights_legendre.pro
function spherical_transform_new2,A,cp,lmax=lmax,period=period

;  This routine performs a spherical harmonic transformation on a 2-D
;  (N_phi,N_theta) array.  Currently all the work is done in idl, but future
;  versions may want to call C or Fortran routines for efficiency reasons.
;  The input parameters include the array A and the colocation points for the
;  grid (cos(theta)), which is of length N_theta, along with an optional
;  specification of lmax.

;  preliminaries
    costheta=double(cp)
    sintheta=sqrt(1-costheta^2)
    sz=size(A)
    n_phi=sz(1)
    n_theta=sz(2)
    if not keyword_set(lmax) then lmax=(2*n_theta-1)/3  else lmax=long(lmax(0))
    if not keyword_set(period) then period=1 else period=fix(period(0))

;  first compute the integration weights
    weights=weights_legendre(costheta)

;  next do the Fourier transform: phi -> m
    Bm=dcomplexarr(n_phi,n_theta,/noz)
    for i=0,n_theta-1 do Bm(*,i)=fft(A(*,i),/double)
    
;  finally the Legendre transform: theta -> l
    B=dcomplexarr(lmax+1,lmax/period+1) ;  only half of this array will be filled

;  Define N_mm such that Y_mm = N_mm sin^m(theta) exp(i m phi) i.e. it's the
;  normalization for the sectoral harmonics.  It will be useful below in
;  computing the spherical harmonics recursively.
    N_mm=dblarr(lmax+1)
    N_mm(0)=1/sqrt(4d0*!dpi)
    for m=1,lmax do N_mm(m)=-N_mm(m-1)*sqrt(1+1/double(2*m))

; first do m=0
    P_lm2=N_mm(0)
    P_lm1=P_lm2*costheta*sqrt(3d)
    B(0,0)=total(Bm(0,*)*P_lm2*weights,/double) ;  set l=0 m=0 term
    
    B(1,0)=total(Bm(0,*)*P_lm1*weights,/double) ;  set l=1 m=0 term

    for l=2,lmax do begin
        lr=double(l)
        c1=sqrt(4d0-1/lr^2)
        c2=-(1-1/lr)*sqrt((2*lr+1)/(2*lr-3))
        P_l=c1*costheta*P_lm1+c2*P_lm2
        B(l,0)=total(Bm(0,*)*P_l*weights,/double) ;  set m=0 term for all other l's

        P_lm2=P_lm1
        P_lm1=P_l
    endfor

;  note factor of 2 below accounts for the way IDL distributes power
;  in its fft, since only the l modes from 0 to lmax are used below
    Bm=2*Bm

; now the rest of the m's
    old_Pmm=N_mm(0)
    for m=1,lmax do begin

        P_lm2=old_Pmm*sintheta*N_mm(m)/N_mm(m-1)
        P_lm1=P_lm2*costheta*sqrt(double(2*m+3))
        old_Pmm=P_lm2
        if (m mod period) eq 0 then begin
            B(m,m/period)=total(Bm(m/period,*)*P_lm2*weights,/double) ;  set l=m term
            
            if m lt lmax then B(m+1,m/period)=$
              total(Bm(m/period,*)*P_lm1*weights,/double) ;  set l=m+1 term
            mr=double(m)
            
            for l=m+2,lmax do begin
                lr=double(l)
                c1=sqrt((4*lr^2-1)/(lr^2-mr^2))
                c2=-sqrt(((2*lr+1)*((lr-1)^2-mr^2))/((2*lr-3)*(lr^2-mr^2)))
                P_l=c1*costheta*P_lm1+c2*P_lm2
                B(l,m/period)=total(Bm(m/period,*)*P_l*weights,/double)
                
                P_lm2=P_lm1
                P_lm1=P_l
            endfor
        endif
    endfor

return,B
end            
