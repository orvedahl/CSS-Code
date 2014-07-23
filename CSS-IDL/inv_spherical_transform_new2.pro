; idem inv_.._MDRosa but easier to call!
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  inv_spherical_transform.pro - This routine performs an inverse spherical
;                                harmonic transform on a 2-D array.
;
;  usage:  A = inv_spherical_transform(B,cp,period=period,lmax=lmax,
;                                      mrange=mrange)
;       where B(lmax,lmax) = complex array to be transformed ordered (l,m)
;             A(n_phi,n_theta) = transformed array ordered (phi,theta)
;             cp = cosine of theta collocation points for theta grid
;             period = periodicity factor in phi, assumes input array
;                      contains m values which are integral multiples
;                      of period
;             lmax = set to max l value we want to use
;             mrange = set to be range of m values we want to use
;
;  notes: - All calculations are done in double precision.
;         - Will return an array of size (nphi,ntheta)=(2*ntheta,ntheta),
;           where ntheta is n_elements(cp)
;         - Routine is increasingly less accurate for higher l.  To see why,
;           look at table of Legendre functions (m=0 example) in Arfken - 
;           they are alternating series with increasingly larger numbers being
;           added and subtracted from each other
;
;  M.DeRosa - 12 Sep 2000 - created
;             24 Oct 2001 - fixed nasty bug related to sign of m=0 components
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function inv_spherical_transform_new2,B,cp,lmax=lmax,mrange=mrange,period=period

;  preliminaries
costheta=double(cp)
sintheta=sqrt(1-costheta^2)
ntheta=n_elements(cp)
nphi=2*ntheta
if keyword_set(period) then period=fix(period(0)) else period=1
if keyword_set(lmax) then lmax=long(lmax(0)) else lmax=n_elements(B(*,0))-1
case n_elements(mrange) of
  0:  mrange=[0,lmax]
  1:  mrange=[0,(fix(mrange(0))<lmax)]
  2:  mrange=[fix(mrange(0)),(fix(mrange(1))<lmax)]
endcase
Bamp=abs(B)

;  calculate array of phases
wh=where(Bamp ne 0,nwh)
B2=B
if nwh gt 0 then B2(wh)=B(wh)/Bamp(wh)
phase=acos(double(B2))
add=1d0*(imaginary(B2) lt 0)
mult=2d0*(imaginary(B2) lt 0)-1
phase=add*2*!dpi - mult*phase

;  set up array A
A=dblarr(nphi/period,ntheta)

;  take care of modes where m=0
CP_0_0=1/sqrt(4*!dpi)
if (mrange(0) eq 0) then begin

  ;  start with m=l=0 mode
  A=A+replicate(Bamp(0,0)*cos(phase(0,0))*CP_0_0,nphi,ntheta)

  ;  now do l=1 m=0 mode
  CP_1_0=sqrt(3d0)*costheta*CP_0_0
  Y=replicate(cos(phase(1,0)),nphi) # CP_1_0
  A=A+(Bamp(1,0)*Y)

  ;  do other l modes for which m=0
  CP_lm1_0=CP_0_0
  CP_l_0=CP_1_0
  for l=2,lmax do begin
    ld=double(l)
    CP_lm2_0=CP_lm1_0
    CP_lm1_0=CP_l_0
    c1=sqrt(4*ld^2-1)/ld
    c2=sqrt((2*ld+1)/(2*ld-3))*((ld-1)/ld)
    CP_l_0=c1*costheta*CP_lm1_0-c2*CP_lm2_0
    Y=replicate(cos(phase(l,0)),nphi) # CP_l_0
    A=A+(Bamp(l,0)*Y)
  endfor
endif

clearline = fifteenb()
;  loop through m's for m>0, and then loop through l's for each m
CP_m_m=CP_0_0
for m=1,mrange(1)-period do begin

  md=double(m)

  ;  do l=m mode first
  CP_mm1_mm1=CP_m_m
  CP_m_m=-sqrt(1+1/(2*md))*sintheta*CP_mm1_mm1
  if (mrange(0) le m) and ((m mod period) eq 0) then begin
    angpart=cos(md*2*!dpi*dindgen(nphi/period)/nphi + phase(m,m/period))

    A=A+Bamp(m,m/period)*(angpart#CP_m_m)

    ;  now do l=m+1 mode
    if lmax ge m+1 then begin
      CP_mp1_m=sqrt(2*md+3)*costheta*CP_m_m
      angpart=cos(md*2*!dpi*dindgen(nphi/period)/nphi + phase(m+1,m/period))
      A=A+Bamp(m+1,m/period)*(angpart#CP_mp1_m)
    endif

    ;  now do other l's
    if lmax ge m+2 then begin
      CP_lm1_m=CP_m_m
      CP_l_m=CP_mp1_m
      for l=m+2,lmax do begin
        ld=double(l)
        CP_lm2_m=CP_lm1_m
        CP_lm1_m=CP_l_m
        c1=sqrt((4*ld^2-1)/(ld^2-md^2))
        c2=sqrt(((2*ld+1)*((ld-1)^2-md^2))/((2*ld-3)*(ld^2-md^2)))
        CP_l_m=c1*costheta*CP_lm1_m-c2*CP_lm2_m
        angpart=cos(md*2*!dpi*dindgen(nphi/period)/nphi + phase(l,m/period))
        A=A+Bamp(l,m/period)*(angpart#CP_l_m)
      endfor
    endif
 endif

  per = strtrim(round((1.0d0*m)/double(mrange[1]-period)*100.0d0),2)
  lenp = strtrim(strlen(strtrim(per,2)),2)
  form="($,a"+lenp+",' % Completed',a,a)"
  print, form=form, per, '         ', clearline
endfor

return,A

end

function fifteenb

   return, string("15b)

end


