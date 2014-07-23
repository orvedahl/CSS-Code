@UniformDerivative6.pro
;Solves for a hydrostatic rho given g and dsdr and a constant Cp and gamma.

;Solves P = rho^gamma * exp[gamma * S/Cp] and dP/dr = -rho g => 
;drho/dr = -rho/Cp * dS/dr - rho^(2-gamma) * g * exp[-gamma * S/Cp]/gamma

;Builds the sixth order accuracy derivative matrix for either a uniform (css) or nonuniform grid (ash)
function build_matrix,n

   A = dblarr(n,n)
   A(*,*) = 0d0

   left  = [-147d0,360d0,-450d0,400d0,-225d0,72d0,-10d0]
   left1 = [-10d0,-77d0,150d0,-100d0,50d0,-15d0,2d0]
   left2 = [2d0,-24d0,-35d0,80d0,-30d0,8d0,-1d0]
   interior = [-1d0,9d0,-45d0,0d0,45d0,-9d0,1d0]
   rightn3 = [1d0,-8d0,30d0,-80d0,35d0,24d0,-2d0]
   rightn2 = [-2d0,15d0,-50d0,100d0,-150d0,77d0,10d0]
   rightn1 = [10d0,-72d0,225d0,-400d0,450d0,-360d0,147d0]
   A(0L,0L:6L) = left(*)
   A(1L,0L:6L) = left1(*)
   A(2L,0L:6L) = left2(*)
   for i=3L,n-4L do begin
       A(i,i-3L:i+3L) = interior(*)
   endfor
   A(n-3L,n-7L:n-1L)=rightn3(*)
   A(n-2L,n-7L:n-1L)=rightn2(*)
   ;A(n-1L,n-7L:n-1L)=rightn1(*)
   A = A/60d0
   A(n-1,n-1) = 1d0             ;Ensures dirichlet boundary at rmin

   return, transpose(A)

end

;Builds the b, in A##x+b = 0
function build_rhs,rho0,dsdrn,Sn,gravn,gamma,Cp,dlngdr,dlnCpdr

   rhs = (rho0)^(2d0-gamma)*gravn*exp(-gamma*Sn/Cp)/gamma+rho0*dsdrn/Cp;(for entropy version)
   rhs = rhs+dlngdr*rho0*(alog(rho0)+Sn/Cp)
   rhs = rhs+dlnCpdr*rho0*Sn/Cp
   
   return, transpose(rhs)
end

;Builds part of the Jacobian for A##x+b=0, namely in total 
;J_ij = A_ij + db_i/drho_j, this routine computes db_i/drho_j = delta_ij * { (2-gamma)*rho_j^(1-gamma)*g_j*exp(-gamma*S_j/Cp)/gamma+(dS/dr)_j/Cp } 
function build_Jac,n,rho0,dsdrn,Sn,gravn,gamma,Cp,dlngdr,dlnCpdr

   J = dblarr(n,n)
   for i=0L,n-2L do begin
       J(i,i) = (2d0-gamma(i))*(rho0(i))^(1d0-gamma(i))*gravn(i)*exp(-gamma(i)*Sn(i)/Cp(i))/gamma(i)+dsdrn(i)/Cp(i)
       J(i,i) = J(i,i)+dlngdr(i)*(alog(rho0(i))+1d0+Sn(i)/Cp(i))
       J(i,i) = J(i,i)+dlnCpdr(i)*Sn(i)/Cp(i)
   endfor
   J(n-1L,n-1L)=0d0 ;Ensures dirichlet boundary at rmin
   return, J

end

function hydrostatic_solver, imin, imax, r, rho, grav, dsdr, Cp, gamma, P, n1, $
                             dsdrn=dsdrn, Sn=Sn, Pn=Pn, Tn=Tn, gravn=gravn, stop=stop, nonconst=nonconst
    rmax = r(imax)
    rmin = r(imin)

   ;Set up a regularly gridded model for CSS
    n = n1
    dr = (rmax-rmin)/double(n-1L)
    If (n_elements(r) ne n1) Then Begin
        rn = dr*dindgen(n)+rmin
    Endif Else Begin
        rn=r
        rhon=rho
        gravn=grav
    EndElse

    If (not keyword_set(nonconst)) Then Begin
        Cpn = Cp+dindgen(n)*0d0
        gamn = gamma+dindgen(n)*0d0
        dlngdr = 0d0*Cpn
        dlnCpdr = 0d0*gamn
    Endif Else Begin
        gamn = interpol(gamma,r,rn,/SPLINE)
        dlngdr = UniformDerivative6(dr,alog(gamn),1)
        Cpn = interpol(Cp,r,rn,/SPLINE)
        dlnCpdr = UniformDerivative6(dr,alog(Cpn),1)
    Endelse

    If (n_elements(r) ne n1) Then Begin
        gravn = interpol(grav,r,rn,/SPLINE)
        rhon = interpol(rho,r,rn,/SPLINE)
                                ;Build entropy       
        dsdrn = interpol(dsdr,r,rn,/SPLINE)
        if (dsdrn(1) lt dsdrn(2)) then begin
            dsdrn(1) = (dsdrn(2)+dsdrn(0))/2.0d0
        endif
                                ;Integrate dsdr to get entropy
        Sn = dblarr(n)
        Sn(*)=0.0d0
                                ;forward
        for i=1L,n-1L do begin
            Sn(i) = int_tabulated(rn(0:i),dsdrn(0:i),/DOUBLE)
        endfor
                                ;backward
                                ;for i=n-2L,0L do begin
                                ;    Sn(i) = Sn(i) + int_tabulated(rn(i:n-1L),dsdrn(i:n-1L),/DOUBLE)
                                ;endfor
                                ;Sn(1L:n-2L) = Sn(1L:n-2L)/2.0d0
        
                                ; set the arbitrary constant
        ir=0L
        rho_o = rhon(ir)
        P_o = P(imin)
        Sn = Sn - Sn(ir) + Cpn(ir)*alog(P_o^(1d0/gamn(ir))/rho_o)
        
                                ;stop
    EndIf


    ;Now compute the density
    ;First build the derivative matrix
      
    A = build_matrix(n)
    A = A/dr
    A(n-1,n-1) = 1d0

    ;Iterate until within tolerance
    niter=0L
    nmax=20L
    eps = 1d-12
    err = 1d10
    rho0 = rhon
    grav0 = gravn
    drho = dblarr(n)
    drho(*) = 0d0
    dg = dblarr(n)
    dg(*) = 0d0
    while ((niter lt nmax) and (err gt eps)) do begin
       ;drho(*) = 0d0
       dg(*) = 0d0
       ;build the rhs
       rhs = build_rhs(rho0,dsdrn,Sn,grav0,gamn,Cpn,dlngdr,dlnCpdr)
       ;build the jacobian
       J = build_Jac(n,rho0,dsdrn,Sn,grav0,gamn,Cpn,dlngdr,dlnCpdr)
       ;solve system
       rhs = -reform(A##(rho0)+rhs)
       rhs(n-1) = 0d0
       J = A+J
       J(n-1,n-1) =1d0
       ludc,J,index
       drho = lusol(J,index,rhs)
       drho(n-1) = 0d0
       rho0 = rho0+drho
       rho0(n-1) = rhon(n-1)
       ;drho = drho+drhon
       ;rho0 = rho0+smooth(drho,20)
       ;Pn = rho0^(gamma)*exp(gamma*Sn/Cp)
       err = total(abs(drho))
       niter = niter+1L
       print, niter, err
    endwhile
       
    ;Now compute pressure and temperature
    Pn = rho0^(gamn)*exp(gamn*Sn/Cpn)
    Tn = gamn*Pn/rho0/Cpn/(gamn-1d0)    

;    window, 10
;    wset, 10
;    !P.MULTI=0
;    plot, r(imin:imax), P(imin:imax)
;    oplot, rn, Pn

    err = (A##Pn+rho0*grav0)/(rho0*grav0)
    err(n-1) = 0d0
    merr = mean(err)

    print,'hydrostatic error:', merr

    if keyword_set(stop) then stop

    gravn=grav0
    rhon=rho0

    return, rho0

end
