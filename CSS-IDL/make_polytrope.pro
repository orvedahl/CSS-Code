@poly_solver.pro
@NonUniformDerivative6.pro

pro make_polytrope, filename, poly_idx, new_nr, smod=smod, dencon=dencon, reverse=reverse, dsmin=dsmin, rmin=rmin, rmax=rmax, simple=simple, Ra=Ra

    ;======================; 
    ; Initialize variables ;
    ;======================;
    bigG = 6.67259d-8
    Lsun = 3.846d33
    imin = 0

    if(not keyword_set(print_idx)) then begin 
       print_idx=18
    endif

    If (not keyword_set(gamma)) then begin
       gamma = 1.6666666666666666666666666666666666667d0
    EndIf

    If keyword_set(smod) Then Begin
       filename='smod4.dat'
    EndIf

    openr, file_unit, filename, /get_lun

    Nr = long(0)

    if keyword_set(smod) then begin
       n1=Nr
       n2=Nr
       n3=Nr
       readf, file_unit, Nr, n1, n2, n3, format='(4I10)'
       data = dblarr(15)
       readf, file_unit, data, FORMAT='(5e16.9 / 5e16.9 / 5e16.9)'
       mass = data(0)
       rstar = data(1)
       L0 = data(2)
       print, data
    endif else begin
       readf, file_unit, Nr, format='(I10)'
       readf, file_unit, mass, rstar, L0,  FORMAT='(3e16.9)'
    endelse

    print, 'Number of Points ', ' Mass (g) ', ' Maximum Model Radius (cm) ', ' Luminosity (erg/cm^2/s)'
    print, Nr, mass, rstar, L0

    data = dblarr(20)
    radius = dblarr(Nr)
    dsdr = dblarr(Nr)
    gravity = dblarr(Nr)
    grav_est = dblarr(Nr)
    kappa_rad = dblarr(Nr)
    ddrkappa_rad = dblarr(Nr)
    L_conv = dblarr(Nr)
    L_rad = dblarr(Nr)
    L = dblarr(Nr)
    T = dblarr(Nr)
    P = dblarr(Nr)
    rho = dblarr(Nr)
    temp = dblarr(Nr)
    ;====================================;
    ; Read In Quantities from .rnd2 file ;
    ;====================================;
    if keyword_set(smod) then begin
       for i=Nr-1L,0L,-1 do begin
       
           readf, file_unit, data, FORMAT='(5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'

           radius(i) = data(0)
           if (i eq 0) then begin
              gravity(i) = 0.d0
           endif else begin
              gravity(i) = bigG*exp(data(1))*mass/radius(i)^2
              grav_est(i) = bigG*mass/radius(i)^2
           endelse

           T(i) = data(2)
           P(i) = data(3)
           rho(i) = data(4)
           kappa_rad(i) = data(7)
           L(i) = data(6)
           L_rad(i) = data(17)
           L_conv(i) = L(i)-L_rad(i)
           if keyword_set(reverse) then begin
              dsdr(i) = -data(19)
           endif else begin
              dsdr(i) = data(19)
           endelse

           if (print_idx eq 1) then begin
              temp(i) = gravity(i)
           endif else begin
              temp(i) = data(print_idx)
           endelse
       endfor
    endif else begin
       for i=Nr-1L,0,-1 do begin
           ;if (i eq Nr-1L) then begin 
           ;   readf, file_unit, data, FORMAT='(//5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'
           ;endif else begin
              readf, file_unit, data, FORMAT='(5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'
           ;endelse 

           radius(i) = data(0)
 
           if (i eq 0) then begin
              gravity(i) = 0.d0
           endif else begin
              gravity(i) = bigG*exp(data(1))*mass/radius(i)^2
              grav_est(i) = bigG*mass/radius(i)^2
           endelse
    
            T(i) = data(2)
            P(i) = data(3)
            rho(i) = data(4)
            kappa_rad(i) = data(7)
            L(i) = data(6)
            L_rad(i) = data(17)
            if keyword_set(reverse) then begin
               dsdr(i) = -data(19)
            endif else begin
               dsdr(i) = data(19)
            endelse

            if (print_idx eq 1) then begin
               temp(i) = gravity(i)
            endif else begin
               temp(i) = data(print_idx)
            endelse
       endfor
    endelse
    close, file_unit
    free_lun, file_unit
    Lmax = max(L)
    if keyword_set(dsmin) then begin
       imin = min(where(dsdr(nr/2:*) lt 0d0))+nr/2
    endif else begin
       imin = max(where(radius/rstar lt rmin)) 
    endelse
    if keyword_set(dencon) then begin
       imax = min(where(rho(imin)/rho ge dencon))
       rmax = radius(imax)/rstar
    endif else begin
       imax = min(where(radius/rstar gt rmax))
    endelse

   if (not keyword_set(dsmin)) then begin
       if (not keyword_set(rmax)) then begin
          rmax = 0.99
       endif
       if (not keyword_set(rmin)) then begin
          rmin = 0.85
       endif
       rmin = rstar*rmin
       rmax = rstar*rmax
    endif else begin
       rmin = radius(imin)
       rmax = radius(imax)
    endelse

    print, 'Array min and max', imin, imax
    print, 'Convection Zone radii: min and max', radius(imin), radius(imax)
    print, 'Top of convection zone is % of rmax',radius(imax)/radius(Nr-1)*100.
    print, 'The thickness of the convective zone is: ', rmax-rmin, 'cm or the outer', (radius(imax)-rmin)/rmax *100, '% of the star.'
    print, 'Density Contrast: ', rho(imin)/rho(imax), rho(imin), rho(imax)
    print, 'ds/dr values at the end points for boundary conditions: ', dsdr(imin), dsdr(imax)

    a = 7.56577d-15
    c = 2.99792458d10

    ;=======================================================;
    ; Build a hydrostatic polytropic background model (CSS) ;
    ;=======================================================;
    n = new_nr
    Cp=3.5d8
    gam1 = gamma-1.0d0
    T0 = T(imin)
    rho0 = rho(imin)
    rho1 = rho(imax)  
    P0 = gam1*Cp*rho0*T0/gamma

    If (not simple) Then Begin
        theta1 = (rho1/rho0)^(1d0/poly_idx)
        K = gam1*Cp*T0/rho0^(1d0/poly_idx)/gamma
        rn = ((poly_idx+1d0)*P0/(4d0*!PI*bigG))^(0.5d0)/rho0
        xi = rad/rn
        dtdr = NonUniformDerivative6(radius,T)
        dtdr = rn*dtdr(imin:*)/T(imin)
        dtdr1 = dtdr(imax-imin)
        poly_solver, xi, poly_idx, rho0, rho1, dtdr1, xin=xin, theta=theta
    EndIf Else Begin
        mass = dblarr(Nr)
        for i=1L,Nr-1L do begin
            mass(i) = 4d0*!PI*int_tabulated(radius(0:i),radius(0:i)^2*rho(0:i),/DOUBLE)
        endfor
        menc = mass(imin)
        x = 0.73d0
        y = 0.25d0
        kappa_e = 0.2*(1+x)
        kappa_g = 6.7d23*(x+y)*(1+x)

        rho1 = rho(imax)
        T1   = T(imax)
        rmax = radius(imax)
        Pr = 1d0
        mu0 = 1d10
        If (keyword_set(Ra)) Then Begin
            Ra1 = Ra
        Endif Else Begin
            Ra1 = 6d7
        EndElse

        delta = 10d0
        alpha = delta*rmin/(rmax-rmin)
        GM = gam1*(poly_idx+1d0)*Cp*T1*rmax*alpha/gamma
        K = gam1*Cp*T1/rho1^(1d0/poly_idx)/gamma

        dr = (rmax-rmin)/(n-1d0)
        rad = dr*dindgen(n)+rmin
        theta = 1d0+alpha*(rmax/rad-1d0)
        poly_rho = rho1*theta^poly_idx
        poly_T = T1*theta
        poly_P = K*poly_rho^(1d0+1d0/poly_idx)

        k0 = delta*rho1*Cp*sqrt(rmin*rmax*Cp*T1*gam1*(poly_idx+1d0)*(1d0-(poly_idx+1d0)*gam1/gamma)/gamma/Pr/Ra1*(1d0+alpha*(sqrt(rmax/rmin)-1d0))^(2d0*poly_idx-1d0))
        muexp = poly_idx
        mu = Pr*k0*theta^(muexp)/Cp
        kapexp = 0d0
        kr = k0*theta^(kapexp)
        kap = gamma*kr/poly_rho/Cp
        
        poly_Re = 2d5*poly_rho*(rmax-rmin)/mu

        ;poly_Ra = alpha*bigG*menc*(rmax-rmin)^4/rad^4*poly_rho/mu/theta
        ;poly_Ra = A^2*(rmax-rmin)^4*poly_rho*Cp*poly_T*gam1*Pr*(1d0-(poly_idx+1d0)*gam1/gamma)/gamma/mu/kap/rad^2
        Ck = k0/(rho1*Cp*(rmax-rmin)*sqrt(Cp*T1*gam1/gamma))
        poly_Ra = delta^2*(poly_idx+1d0)*(1-(poly_idx+1d0)*gam1/gamma)*rmax^2*rmin^2/rad^2/Pr/Ck^2*theta^(2*poly_idx-1d0-muexp-kapexp)/(rmax-rmin)^2
        grav = GM/rad^2

        !P.Charsize=2
        window,1
        plot, rad, poly_Ra
        window,2
        plot,rad,poly_Re

        th2 = 100d0
        th1 = 80d0
        dth = (th2-th1)/256d0*!PI/180d0
        dph = dth*sin(th1*!PI/180d0)
        dx = min([dr,rmin*dth,rmin*dph])
        taudk = min(dx^2/kap)
        taudmu = min(dx^2*poly_rho/mu)
        tauc = dx/sqrt(gam1*Cp*poly_T(0))

        print, taudmu, taudk, tauc
        print, 'rmin =', rmin
        print, 'rmax =', rmax
        print, 'rho1 =', rho1
        print, 'T1   =', T1
        print, 'mu0  =', mu0
        print, 'Pr   =', Pr
        print, 'Ra1  =', Ra1
        print, 'muexp  =', muexp
        print, 'kapexp =', kapexp
    EndElse
    

    stop

end
