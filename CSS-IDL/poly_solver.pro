;Solves for a hydrostatic rho given g and dsdr and a constant Cp and gamma.

;Solves P = rho^gamma * exp[gamma * S/Cp] and dP/dr = -rho g => 
;drho/dr = -rho/Cp * dS/dr - rho^(2-gamma) * g * exp[-gamma * S/Cp]/gamma

;Builds the sixth order accuracy derivative matrix for either a uniform (css) or nonuniform grid (ash)
function build_matrix,nr

   A = dblarr(nr,nr)
   A(*,*) = 0d0
   left0 = [-147d0,360d0,-450d0,400d0,-225d0,72d0,-10d0]
   left1 = [-10d0,-77d0,150d0,-100d0,50d0,-15d0,2d0]
   left2 = [2d0,-24d0,-35d0,80d0,-30d0,8d0,-1d0]
   interior = [-1d0,9d0,-45d0,0d0,45d0,-9d0,1d0]
   rightn3 = [1d0,-8d0,30d0,-80d0,35d0,24d0,-2d0]
   rightn2 = [-2d0,15d0,-50d0,100d0,-150d0,77d0,10d0]
   rightn1 = [10d0,-72d0,225d0,-400d0,450d0,-360d0,147d0]
   A(0L,0L:6L) = left0(*)
   A(1L,0L:6L) = left1(*)
   A(2L,0L:6L) = left2(*)
   for i=3L,nr-4L do begin
       A(i,i-3L:i+3L) = interior(*)
   endfor
   A(nr-3L,nr-7L:nr-1L)=rightn3(*)
   A(nr-2L,nr-7L:nr-1L)=rightn2(*)
   A(nr-1L,nr-7L:nr-1L)=rightn1(*)
   A = A/60d0
   ;A(0,0) = 0d0                 ;Ensures dirichlet boundary at rmin
   ;A(nr-1L,nr-1L) = 0d0         ;Ensures Dirichlet boundary at rmax

   return, transpose(A)

end

;Builds part of the Jacobian for A##x+b=0, namely in total 
;J_ij = A_ij + db_i/drho_j, this routine computes db_i/drho_j = delta_ij * { (2-gamma)*rho_j^(1-gamma)*g_j*exp(-gamma*S_j/Cp)/gamma+(dS/dr)_j/Cp } 
function build_Jac,theta,nr,pidx

   J = diag_matrix((pidx-1d0)*theta^(pidx-2d0))
   J(0,0)=0d0 ;Ensures dirichlet boundary at rmin
   J(nr-1L,nr-1L)=0d0 ;Ensures dirichlet boundary at rmax
   return, J

end

pro poly_solver, xi, pidx, rho0, rho1, dtdr, xin=xin, theta=theta

    nr = n_elements(xi)
    ;Iterate until within tolerance
    eps = 1d-12
    xi0 = xi(0)
    xi1 = xi(nr-1L)
    xin = xi
    xi_iter=0L
    xi_err = 1d10
    t1 = (rho1/rho0)^(1d0/pidx)
    nmax=4L
    while ((xi_err gt eps) and (xi_iter lt nmax)) do begin
        niter=0L
        err = 1d10
        A = build_matrix(nr)
        dxi = xin(1)-xin(0)
        A = A/dxi        
        theta = make_array(nr,/double,value=1d0)
        u = theta*xin
        while ((niter lt nmax) and (err gt eps)) do begin
            ;J = build_Jac(theta,nr,pidx)
            ;D = A#A + diag_matrix(2d0/xin)#A + diag_matrix(theta^(pidx-1d0))
            ;J = A#A + diag_matrix(2d0/xin)#A + J
            D = A#A + diag_matrix((u/xin)^(pidx-1d0))
            J = A#A + diag_matrix(pidx*(u/xin)^(pidx-1d0))
            J(0,0) = 1d0
            J(nr-7L:nr-1L,nr-1L) = A(nr-7L:nr-1L,nr-1L)
            rhs = reform(D##u)
            rhs(0) = 0d0
            rhs(nr-1L) = dtdr
            ludc,J,index
            du = lusol(J,index,rhs)
            u(1L:*) = u(1L:*)-du(1L:*)
            D = A#A + diag_matrix((u/xin)^(pidx-1d0))
            err = total(abs(D##u))
            print, 'NR err = ', err
            niter = niter+1
        endwhile
        dudx = reform(A##u)
        du1 = dudx(nr-1L)
        u1 = u(nr-1L)
        xi1 = (u1+du1*xi1)/(t1+du1)
        print, 'xi1 = ', xi1
        xip1 = xin(nr-1L)
        xin = (xi1-xi0)*dindgen(nr)/(nr-1d0)+xi0
        xi_err = abs((xi1-xip1)/xi1)
        xi_iter = xi_iter+1L
    endwhile

    stop

end
