;==========================================================================;
; This program reads in a .rnd (or.rnd2) file and determines the parameters
; necessary for a successful start of a new CSS case from a stellar
; model.
; 
; Kyle Augustson August 2006
;
;-------------------------------------------------------------------------
;Inputs:
; filename - the .rnd2 filename (and directory information if necessary)
; /stop - execute the stop at the end of the program
; print_idx - the quantity index in the .rnd2 file that the spare
;             window in the main window (window 0) will display
; /nomin - if you do not wish to compute the minimum index of the
;          convection zone
; /psf - output to postscript (not fully implemented)
; max_mach - the maximum Mach number that determines the
;            effective top of the convection zone.
;==========================================================================;
@docolor.pro
@NonUniformDerivative6.pro
@UniformDerivative6.pro
@hydrostatic_solver.pro

;=================================================;
; Mach number estimate (see Marc Derosa's thesis) ;
;=================================================;

function Mach, T, r, rho, L, c_p, gamma
   epsilon = (L/r^2/rho/!PI/4.)^(2./3.)/T/c_p
   return, (epsilon/(gamma-1.))^(1./2.)
end

function envelope, r, delta, two=two, rm=rm

   ;Typical delta is ~200d0
   r1 = min(r)
   r2 = max(r)
   If (keyword_set(rm)) Then Begin
       x = (rm-r)/(r2-r1)
       func = 1d0/(1d0+exp(delta*x))
   Endif Else Begin
       x = (r2-r)/(r2-r1)
       func = 2d0/(1d0+exp(delta*x^2))
       If (keyword_set(two)) Then Begin
           x = (r1-r)/(r2-r1)
           func = func+2d0/(1d0+exp(delta*x^2))
       Endif
   EndElse

   return, func

end

pro make_stellar_model, filename, stop=stop, print_idx=print_idx, rev=rev,rmin=rmin, rmax=rmax, $
                        new_nr=new_nr, poly_idx=poly_idx, gamma=gamma, smod=smod, dencon=dencon, dsmin=dsmin, $
                        ks_delta=ks_delta, Ra_top=Ra_top, rain=rain, nsqual=nsqual, constmu=constmu, $
                        constnu=constnu, Cp=Cp, nonuniform=nonuniform, dtmin=dtmin, nonconstgcp=nonconstgcp, $
                        fstar=fstar, surfcool=surfcool, tycho=tycho, midpt=midpt, shellavg=shellavg, evolS=evolS, $
                        botio=botio, fig1=fig1, cesam=cesam, mpt0=mpt0, nofrad=nofrad, Pr_top=Pr_top, evorad=evorad

    ;======================; 
    ; Initialize variables ;
    ;======================;
    bigG = 6.67259d-8
    Lsun = 3.846d33
    imin = 0

    if(not keyword_set(print_idx)) then begin 
       print_idx=18
    endif

    If (not keyword_set(Cp)) Then Begin
        Cp = 3.5d8
    EndIf

    If (not keyword_set(gamma)) then begin
       gamma = 1.6666666666666666666666666666666666667d0
    EndIf

    If keyword_set(smod) Then Begin
       filename='smod4.dat'
    EndIf

    If (not keyword_set(rain)) Then Begin
        rain =0
    EndIf
    
    if (not keyword_set(cesam)) then begin
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
        If (keyword_set(nonconstgcp)) Then gamma = dblarr(Nr)
        gravity = dblarr(Nr)
        grav_est = dblarr(Nr)
        kappa = dblarr(Nr)
        kappa_rad = dblarr(Nr)
        ddrkappa_rad = dblarr(Nr)
        L_conv = dblarr(Nr)
        L_rad = dblarr(Nr)
        L = dblarr(Nr)
        T = dblarr(Nr)
        P = dblarr(Nr)
        S = dblarr(Nr)
        rho = dblarr(Nr)
        temp = dblarr(Nr)
        all_qs = dblarr(Nr,20)
    endif
    ;====================================;
    ; Read In Quantities from .rnd2 file ;
    ;====================================;
    if (keyword_set(smod)) then begin
       for i=Nr-1L,0L,-1 do begin
       
           readf, file_unit, data, FORMAT='(5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'

           all_qs(i,*) = data

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
           L(i) = data(6)
           kappa(i) = data(7)
           If (keyword_set(nonconstgcp)) Then Begin
               gamma(i) = data(9)
           EndIf
           S(i) = data(16)
           L_rad(i) = data(17)
           L_conv(i) = L(i)-L_rad(i)
           if keyword_set(rev) then begin
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
       close, file_unit
       free_lun, file_unit
   endif

   if (keyword_set(tycho)) Then Begin
    ;Read in the tycho radial output header
       openr, funit, filename, /get_lun
       nhead = 39
       data = string(128)
    ;first line is not needed
       readf, funit, data, format='(a128)'
    ;second line is not needed
       readf, funit, data, format='(a128)'
    ;third line is number of lagrangian grid points
       readf, funit, data, format='(a128)'
       substrings = strsplit(data,' ',/EXTRACT)
       print, substrings

       nr = long(fix(strtrim(substrings(1))))
       nrt = nr+2L
    
    ;Skip a bunch of stuff
       for i=2,7 do begin
           readf, funit, data, format='(a128)'
       end

    ;tenth line is number of isotopes
       readf, funit, data, format='(a128)'
       substrings = strsplit(data,' ',/EXTRACT)
       print, substrings

       niso = fix(strtrim(substrings(1)))

       for i=9,42 do begin
           readf, funit, data, format='(a128)'
       end

    ;prepare to read in various quantities
       baseqs = 18
       elemqs = niso+1
       quantities = dblarr(baseqs+elemqs,nrt)
       x = nrt mod 5L
       for i=0,baseqs-1 do begin
           ir = 0L
           while (ir lt nrt-x-1L) do begin
               data = dblarr(5)
               readf, funit, data, format='(5e17.8)'
               quantities(i,ir:ir+4L) = data(*)
               ir = ir+5L
           endwhile
           if (x gt 0L) then begin
               data = dblarr(x)
               readf, funit, data, format='('+strtrim(string(x,format='(i2)'),2)+'e17.8)'
               quantities(i,nrt-x:nrt-1L) = data(*)
           endif
           data = string(128)
           readf, funit, data, format='(a128)'
       endfor 

       close,funit
       free_lun, funit
       nr = nr-1L
       radius = reform(quantities(0,0:nr))
       T = reform(quantities(4,0:nr))
       V = reform(quantities(5,0:nr))
       P = reform(quantities(6,0:nr))
       Mz = reform(quantities(7,0:nr))
       cv = reform(quantities(8,0:nr))
       M = reform(quantities(10,0:nr))
       kappa = reform(quantities(11,0:nr)) ;Actual frequency averaged opacities (OPAL2?)
       ss = reform(quantities(12,0:nr))
       snu = reform(quantities(13,0:nr))
       grav = reform(quantities(14,0:nr))
       grav(0) = 0.0d0
       S = reform(quantities(15,0:nr))
       del = reform(quantities(16,0:nr))
       delad = reform(quantities(17,0:nr))

       L = reform(quantities(9,0:nr))
       T(0) = T(1)    
       P(0) = P(1)

       rho = 1.0d0/V
       rho(0)=rho(1)

       eps = ss-snu
       eps(0) = eps(1)

       gamma_t = 1d0/(1d0-delad)
       gam1_t = gamma_t-1.0d0
    ;Pressure is total pressure, gas plus radiation, so to calculate Cp, remove radiative pressure
       Cp_t = gamma_t*P/T/rho/gam1_t

    ;Build entropy gradient, rho, and gravity to be smooth
    
       Hp = 0.25d0*(P+shift(P,1))*(shift(radius,1)-shift(radius,-1))/(P-shift(P,1))
       dsdr = -(del-delad)/Hp

       a = 7.56577d-15
       c = 2.99792458d10
       kappa_rad = 4.0d0*a*c*T^3/rho^2/Cp_t/kappa/3.0d0
       L_rad = -4d0*!dpi*radius^2*kappa_rad*rho*Cp_t*nud6(radius,T)
   endif

   If (keyword_set(cesam)) Then Begin
 ;Read in header
       openr, funit, filename, /get_lun
       data = string(128)
 ;first line is not needed
       readf, funit, data, format='(a128)'
 ;second line is not needed
       readf, funit, data, format='(a128)'
 ;third line is not needed
       readf, funit, data, format='(a128)'
 ;fourth line is not needed
       readf, funit, data, format='(a128)'

 ;Fifth line contains number of isotopes and elements and 
 ; their string names
       readf, funit, data, format='(a128)'
       substrings = strsplit(data,' ',/EXTRACT)
       print, substrings
       nelem = fix(substrings(0))
       elem_names = substrings(1:nelem-1)
       
 ;Sixth line contains number of points, number of global properties, number
 ; of local properties, and number of elements
       readf, funit, data, format='(a128)'
       substrings = strsplit(data,' ',/EXTRACT)
       nr_mod = fix(substrings(0))
       nglob = fix(substrings(1))
       ntot = fix(substrings(2))
       
 ;Seventh line is unknown... skip it
       readf, funit, data, format='(a128)'
       
 ;Lines 8-11 contain global quantities
       glbval = dblarr(nglob)
       readf, funit, glbval, FORMAT='(5e19.12 / 5e19.12 / 2e19.12)'

       mass = glbval(0)
       rstar = glbval(1)
       lstar = glbval(2)
       z0 = glbval(3)
       x0 = glbval(4)
       alpha_mlt = glbval(5)
       z1 = glbval(6)
       x1 = glbval(7)
       age = glbval(9)
       omega0 = glbval(10)

       msun = 1.98d33
       rsun = 6.9599d10
       lsun = 3.86d33

       Print, 'Stellar model properties of ', filename
       Print, ' mass=', mass, ' rstar=', rstar, ' lstar=', lstar, ' omega0=', omega0

 ;Everything else is local
       locval = dblarr(nr_mod,ntot)
       For i=0,nr_mod-1 do begin
           temp = dblarr(ntot)
           readf, funit, temp, FORMAT='(5e19.12 / 5e19.12 / 5e19.12 / 5e19.12 / 5e16.12 / 5e19.12 / 5e19.12 / 5e19.12 / 5e16.12 / 5e19.12 / 5e19.12 / 2e19.12)'
           locval(i,*) = temp
       endfor
    
       close, funit
       free_lun, funit

       radius = locval(*,0)
       lnm = locval(*,1)
       T = locval(*,2)
       P = locval(*,3)
       rho = locval(*,4)
       gradient1 = locval(*,5)
       gradient2 = locval(*,6)
       L = locval(*,7)
       kappa = locval(*,8)
       eps = locval(*,9) ;Nuclear+gravitational energy
       gamma1 = locval(*,10) ;Gamma 1
       delad = locval(*,11)
       del = locval(*,12) ;Delta?
       Cp_mod = locval(*,13)
       mueinv = locval(*,14)
       mu = locval(*,15)
       Vaissala = locval(*,16)
       Omega = locval(*,17)
       dlnkappa_dlnT = locval(*,18)
       dlnkappa_dlnrho = locval(*,19)
       deps_dlnT = locval(*,20)
       deps_dlnrho = locval(*,21)
       Ptot_Pgas = locval(*,22)
       delrad = locval(*,23)
       dgamma1_dlnP = locval(*,24)
       dgamma1_dlnT = locval(*,25)
       dgamma1_dlnY = locval(*,26)
       dP_drho = locval(*,27)
       dP_dT = locval(*,28)
       dP_dX = locval(*,29)
       de_drho = locval(*,30)
       de_dT = locval(*,31)
       de_dX = locval(*,32)
       ienergy = locval(*,33) ;Internal energy
       d2P_drho2 = locval(*,34)
       d2P_drhodT = locval(*,35)
       d2P_dT2 = locval(*,36)
       d2e_drho2 = locval(*,37)
       d2e_drhodT = locval(*,38)
       d2e_dT2 = locval(*,39)
       dK_dX = locval(*,40)
       d2K_dT2 = locval(*,41)
       deps_dX = locval(*,42)
       dX_dR = locval(*,43)
       JmB = locval(*,44)
       Edd = locval(*,45)
       
 ;isotopes
       isotopes = locval(*,46:*)

       mass_r = dblarr(nr_mod)
       mass_r(0) = 0d0
       for ir=1,nr_mod-1 do begin
           mass_r(ir) = 4d0*!dpi*int_tabulated(radius(0:ir),radius(0:ir)^2*rho(0:ir),/DOUBLE)
       endfor
       gamma_ad = 1d0/(1d0-delad)

       G = 6.67d-8
       gravity = G*mass_r/radius^2
       gravity(0)=0d0

       dedr = nud6(radius,ienergy)
       drhodr = nud6(radius,rho)
 ;From second law of thermo:
       ;dsdr = (dedr-P*drhodr/rho^2)/T

       dlnPdr = nud6(radius,P)/P
       dsdr = Cp_mod*(gradient2-delad)*dlnPdr

       nr = nr_mod
       a = 7.56577d-15
       c = 2.99792458d10
       kappa_rad = 4d0*a*c*T^3/(3d0*rho^2*Cp_mod*kappa)
       dtdr_mod = nud6(radius,T)
       L_rad = -4d0*!dpi*radius^2*kappa_rad*rho*Cp_mod*dtdr_mod
       L_conv = L-L_rad
   endif else begin
       full_data = dblarr(Nr,20)
       delad = dblarr(Nr)
       for i=Nr-1L,0,-1 do begin
           If (keyword_set(rev)) Then Begin
               readf, file_unit, data, FORMAT='(5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'
           Endif Else Begin
               if (i eq Nr-1L) then begin 
                   readf, file_unit, data, FORMAT='(//5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'
               endif else begin
                   readf, file_unit, data, FORMAT='(5e16.9 / 5e16.9 / 5e16.9 / 5e16.9)'
               endelse 
           Endelse

           all_qs(i,*) = data

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
            delad(i) = data(10)
            L_rad(i) = data(17)
            if keyword_set(rev) then begin
               dsdr(i) = -data(19)
            endif else begin
               dsdr(i) = data(19)
            endelse

            if (print_idx eq 1) then begin
               temp(i) = gravity(i)
            endif else begin
               temp(i) = data(print_idx)
           endelse
           full_data(i,*) = data
        endfor
        close, file_unit
        free_lun, file_unit
        gamma_mod = 1d0/(1d0-delad)
        Cp_mod = Cp*gamma_mod/gamma
    endelse

    Lmax = max(L)

    if (keyword_set(rmax)) Then begin
        imax = min(where(radius gt rmax*rstar))
    endif

    if keyword_set(dsmin) then begin
        imin = min(where(dsdr(nr/2:*) lt 0d0))+nr/2
    endif

    if keyword_set(dencon) then begin
        imin = max(where(rho/rho(imax) ge dencon))
    endif else begin
        imin = max(where(radius lt rmin*rstar)) 
    endelse

    rmin = radius(imin)
    rmax = radius(imax)
    
    machnum = Mach(T(imin:*), radius(imin:*), rho(imin:*), Lmax, Cp, gamma)

    print, 'Array min and max', imin, imax
    print, 'Convection Zone radii: min and max', radius(imin), radius(imax)
    print, 'Top of convection zone is % of rmax',radius(imax)/radius(Nr-1)*100.
    print, 'The thickness of the convective zone is: ', rmax-rmin, 'cm or the outer', (radius(imax)-rmin)/rmax *100, '% of the star.'
    print, 'Density Contrast: ', rho(imin)/rho(imax), rho(imin), rho(imax)
    print, 'ds/dr values at the end points for boundary conditions: ', dsdr(imin), dsdr(imax)
    drhodr = nud6(radius,rho)
    drhodr = drhodr(imin:imax)

    ;=================================================;
    ; Compute the radiative and convective luminosity ;
    ;=================================================;

    max_kappa = max(kappa_rad)
    temp = dblarr(nr)
    a = 7.56577d-15
    c = 2.99792458d10
    
    ;========================================;
    ; Plot the quantities in the main window ;
    ;========================================;

    !P.MULTI = [0,4,2,0,0]
    WINDOW, XSIZE=1280,YSIZE=924
    plot, radius(imin:imax)/rstar, kappa_rad(imin:imax), TITLE='kappa_rad',XTITLE='r/r_max',YTITLE='kappa_rad cm^2/g', charsize=2, ys=1, /ylog
    plot, radius(imin:imax)/rstar, gravity(imin:imax),TITLE='Gravity',XTITLE='r/r_max',YTITLE='g(r) cm/s^2',charsize=2,ys=1
    plot, radius(imin:imax)/rstar, T(imin:imax), TITLE='Temperature',XTITLE='r/r_max',YTITLE='T(r) K',charsize=2,ys=1
    plot, radius(imin:imax)/rstar, P(imin:imax), TITLE='Pressure',XTITLE='r/r_max',YTITLE='P(r) g/cm/s^2',charsize=2,ys=1
    plot, radius(imin:imax)/rstar, rho(imin:imax), TITLE='Density',XTITLE='r/r_max',YTITLE='rho(r) g/cm^3',charsize=2,ys=1
    plot, radius(imin:imax)/rstar, L_conv(imin:imax)/Lmax, TITLE='Convective Luminosity',XTITLE='r/r_max',YTITLE='Lc(r)/L0',charsize=2,ys=1
    plot, radius(imin:imax)/rstar, L_rad(imin:imax)/Lmax, TITLE='Radiative Luminosity',XTITLE='r/r_max',YTITLE='Lr(r)/L0',charsize=2,ys=1
    plot, radius(imin:imax)/rstar, abs(dsdr(imin:imax)), /ylog, TITLE='Entropy Gradient',XTITLE='r/r_max',YTITLE='log(d/dr S(r)) erg/K/cm',charsize=2,ys=1       

    ;============================================;
    ; Build a hydrostatic background model (CSS) ;
    ;============================================;
    n = new_nr
    ;Build entropy gradient
 
    ;calculate rmax based on desired density contrast
    rho0 = rho(imin)

    If (keyword_set(nonconstgcp)) Then Begin
        If (not keyword_set(smod)) Then Begin
                                ;Build entropy       
                                ;Integrate dsdr to get entropy
            S = dblarr(nr)
            S(*)=0d0
                                ;forward
            for i=1L,nr-1L do begin
                S(i) = int_tabulated(radius(0:i),dsdr(0:i),/DOUBLE)
            endfor
                                ; set the arbitrary constant
            ir=1150L
            rho_o = rho(ir)
            P_o = P(ir)
            Cp = 3.5d8
            gamma = 1.66666667d0
            S = S - S(ir) + Cp*alog(P_o^(1d0/gamma)/rho_o)
            
            gamma = (alog(P)+rho*T*S/P)/(alog(rho)+rho*T*S/P)
        EndIf
        Cp = gamma*P/(gamma-1d0)/rho/T
        ;stop
    EndIf

    gam1 = gamma-1d0

    If (keyword_set(nonuniform)) Then Begin
        nrd = 512L
        rnd = (rmax-rmin)*dindgen(nrd)/double(nrd-1L)+rmin
        Snrd = interpol(S,radius,rnd,/SPLINE)
        rhonrd = hydrostatic_solver(imin,imax,radius,rho,gravity,dsdr,Cp,gamma,P,nrd,Pn=Pnrd, Tn=Tnrd, gravn=gravnrd, dsdrn=dsdrnrd, Sn=Snrd,nonconst=nonconstgcp)
        ;cs = sqrt(gam1*Cp*Tnrd)
        cs = (rmax-rmin)*(1d0-exp(-2d0*dindgen(n)/double(n-1L)))/(1d0-exp(-2d0))+rmin
        ;cs = reverse(cs)

        ;regularly grid in cs
        print, 'nonuniform x1=', cs(0), ' x2=', cs(n-1L)
        xn = (cs(n-1L)-cs(0))*dindgen(n)/double(n-1L)+cs(0)
        rn = cs
        ;Find indices in r
        ;rn = dblarr(n)
        ;xn = dblarr(n)
        ;rn(0) = rmin
        ;rn(n-1) = rmax
        ;xn(0) = cs(0)
        ;xn(n-1) = cs(n-1L)
        ;irp=0
        ;for ic=1,n-2 do begin
        ;    ir = max(where(cs le x(ic)))
        ;    rn(ic) = rnd(ir)
        ;    xn(ic) = cs(ir)
        ;    irp = ir
        ;endfor
        dr = rn-shift(rn,1)
        dr(0) = rn(1)-rn(0)
        ;dr = min(dr)
        stop
    Endif Else Begin
        dr = dblarr(nr)
        dr(*) = (rmax-rmin)/(double(n)-1d0)
        rn = dr(1)*dindgen(n)+rmin
    EndElse
    ;kappar = 4.0d0*a*c*T^3/rho^2/Cp/kappa/3.0d0
    kapparn = interpol(kappa_rad,radius,rn,/SPLINE)
    ;Adjust kappar so that it is 1 solar luminosity at lower boundary (all flux carried by radiative)
    dtdr = nud6(radius,T)

    If (keyword_set(nonuniform)) Then Begin
        rhon = interpol(rhonrd,rnd,rn,/SPLINE)
        gravn = interpol(gravnrd,rnd,rn,/SPLINE)
        Pn = interpol(Pnrd,rnd,rn,/SPLINE)
        Tn = interpol(Tnrd,rnd,rn,/SPLINE)
        dsdrn = interpol(dsdrnrd,rnd,rn,/SPLINE)
        Sn = interpol(Snrd,rnd,rn,/SPLINE)
    Endif Else Begin
        rhon = hydrostatic_solver(imin,imax,radius,rho,gravity,dsdr,Cp,gamma,P,n, Pn=Pn, Tn=Tn, gravn=gravn, dsdrn=dsdrn, Sn=Sn,nonconst=nonconstgcp)
        inds = where(dsdrn eq 0d0)
        If (inds[0] ne -1) Then dsdrn(inds) = 1d-10
    Endelse

    If (keyword_set(shellavg)) Then Begin
        ;Read in the shell avg file
        read_shell_avg, shellavg, radius = r_sh, data = values_sh, QUANTITIES = quantities_sh
        nrec = n_elements(values_sh(0,0,*))
        If (n_elements(r_sh) eq new_nr) Then Begin
            rhon = reform(values_sh(*,0,nrec-1))
            Tn = reform(values_sh(*,1,nrec-1))
            Sn = reform(values_sh(*,2,nrec-1))
            Pn = rhon^(gamma)*exp(gamma*Sn/Cp)
            dsdrn = ud6(dr,Sn,1)
        Endif Else Begin
            Print, 'shell average radial points differs from new_nr'
            stop
        EndElse
    EndIf

    If (keyword_set(nonuniform)) Then Begin
        dtdr = nud6(rn,Tn)
        drhodr = nud6(rn,rhon)
    Endif Else Begin
        ;compact_fd6,dri,arr,b1,b2,ibc,dtype
        dri = 1d0/double(new_nr-1L)/dr
        dtdr = compact_fd6(dri,Tn,0,0,1,0)
        ;dtdr = ud6(dr,Tn,1)
        drhodr = compact_fd6(dri,rhon,0,0,1,0)
        ;drhodr = ud6(dr,rhon,1)
    EndElse
    model_deriv_r1 = dtdr(0)
    model_deriv_r2 = dtdr(new_nr-1)

    print, 'Model Derivatives; Numeric: ', model_deriv_r1, model_deriv_r2
    If keyword_set(poly_idx) Then Begin
       print, 'Analytic: ', A*rmax/rmin^2, A/rmax
    Endif

    ;kapparn = kapparn*rhon*Cp ;!!!!!!!!!!!!!Note the zero!!!!!!!!!!!!!
    L_rad = -4d0*!dpi*radius^2*rho*Cp_mod*kappa_rad*nud6(radius,T)
    Lradn = interpol(L_rad,radius,rn,/SPLINE)

    If (not keyword_set(Pr_top)) Then Pr_top = 1d0
    If (not keyword_set(Ra_top)) Then Ra_top = 1d6
    k_top  = rhon(n-1)*Tn(n-1)*(rmax-rmin)^2*sqrt(-gravn(n-1)*dsdrn(n-1)/Pr_top/Ra_top/Cp)
    mu_top = Pr_top*k_top/Tn(n-1)
    rho_top = rhon(n-1)
    T_top = Tn(n-1)

    If (keyword_set(constmu)) Then Begin
        mur = 0d0*rhon+mu_top
        kappas = k_top*Pr_top*Tn/T_top
    Endif Else If (keyword_set(constnu)) Then Begin
        mur = mu_top*(rhon/rho_top)
        kappas = k_top*(rhon/rho_top)*(Tn/T_top)
    EndIf Else Begin
        mu_min=0.5d0
        muexp = 60d0
        func = ((1d0-mu_min)*exp(-2d0*muexp*((rmin-rn)/(rmax-rmin))^2)+(1d0-mu_min)*exp(-muexp*((rn-rmax)/(rmax-rmin))^2)+mu_min)
        mur = mu_top*(rhon/rho_top)*func
        kappas = k_top*(Tn/T_top)*(rhon/rho_top)*func
    EndElse

    kapparn = -Lradn/(4d0*!dpi*rn^2*dtdr)
    If (kapparn(new_nr-1) lt 0) Then kapparn(new_nr-1) = kapparn(new_nr-2)
    Ra = -rhon^2*Tn*dsdrn*gravn*(rmax-rmin)^4/Pr_top/mur/kappas/Cp
    Re = 1d5*(rmax-rmin)*rhon/mur

    env = envelope(rn,1d3) ;,/two)
    r2 = max(rn)
    T2 = min(T)
    rho2 = min(rhon)
    kr2 = kapparn(new_nr-1)
    ds2 = dsdrn(new_nr-1)
    dt2 = dtdr(new_nr-1)
    ksfact = 1d-2
    k0 = -(max(L)/(4d0*!dpi*r2^2*ds2)+kr2*dt2/ds2+kappas(new_nr-1)*ksfact)
    kappa0 = k0*env+kappas*ksfact*(rho2/rhon)^3 ;/abs(dsdrn/ds2+1d0)
    If (keyword_set(nonuniform)) Then Begin
        epsn = nud6(rn,rn^2*kappa0*dsdrn)/rn^2/rhon/Tn
        Fepsn = -rhon*Tn*kappa0*dsdrn
    Endif Else Begin
        epsn0 = rn^2*kappa0*dsdrn ;+kapparn*dtdr*rn^2
        epsn = compact_fd6(dri,epsn0,0,0,1,0)/rn^2/rhon/Tn
        Fepsn = -kappa0*dsdrn ;-kapparn*dtdr
    EndElse

    !P.MULTI=[0,1,1]
    !P.Charsize=1.75
    window, 5
    Lrad = -4d0*!dpi*rn^2*kapparn*dtdr/max(L)
    rnorm = rn/max(radius)
    plot, rnorm, Lrad, yrange=[-0.1,1.3],/xs
    Leps =  -4d0*!PI*rn^2*kappa0*dsdrn/max(L)
    oplot, rnorm, Leps
    Ltot = Lrad+Leps
    oplot, rnorm, Ltot, linestyle=3,thick=2
    oplot, rnorm, 0d0*rn+1d0, linestyle=2

    th2 = 90d0
    th1 = 50d0
    ph2 = 30d0
    ph1 = 10d0
    ntheta = 512L
    nphi = 256L
    dth = (th2-th1)/ntheta*!PI/180d0
    dph = (ph2-ph1)/nphi*!PI/180d0

    kap = kappas/rhon/Tn
    dal = rn*dth
    daz = rn*dph*sin(th1*!dpi/180d0)
    
    ;If (keyword_set(nonuniform)) Then Begin
        taudk  = dblarr(new_nr)
        taudmu = dblarr(new_nr)
        tauc   = dblarr(new_nr)
        taudkr = dblarr(new_nr)
        taudk0 = dblarr(new_nr)
        for ir=0,new_nr-1 do begin
            dm = min([dr(ir),dal(ir),daz(ir)])
            taudk(ir)  = dm^2/kap(ir)
            taudmu(ir) = dm^2*rhon(ir)/mur(ir)
            tauc(ir)   = dm/sqrt(gam1*Cp*Tn(ir))
            taudkr(ir) = dm^2*rhon(ir)*Cp/kapparn(ir)
            If (keyword_set(surfcool)) Then Begin 
                taudk0(ir) = dm^2*rhon(ir)*Cp/abs(kappa0(ir))/Tn(ir)^3
            Endif Else Begin
                taudk0(ir)=0
            Endelse
        endfor
        print, min(dr), min(dal), min(daz)
        print, min(taudmu), min(taudk), min(tauc), min(taudkr), min(taudk0)
    ;Endif Else Begin
    ;    dx = min([dr,dal,daz])
    ;    taudkr = min(dx^2*rhon*Cp/kapparn)
    ;    If (keyword_set(surfcool)) Then Begin 
    ;        tauk0 = min(dx^2*rhon*Cp/abs(kappa0)/Tn^3)
    ;    Endif Else Begin
    ;        tauk0=0
    ;    EndElse
    ;    taudk = min(dx^2/kap)
    ;    taudmu = min(dx^2*rhon/mur)
    ;    tauc = min(dx/sqrt(gam1*Cp*Tn))
    ;    print, dr, dal, daz
    ;    print, taudmu, taudk, tauc, taudkr, tauk0
    ;EndElse

    ;hydrostatic check
    If (keyword_set(nonuniform)) Then Begin
        check = gam1*Cp*Tn*nud6(rn,alog(rhon))+gam1*Tn*dsdrn+gravn
    Endif Else Begin
        check = gam1*Cp*Tn*ud6(dr,alog(rhon),1)+gam1*Tn*dsdrn+gravn
    EndElse
    ;print, 'check = [', check, ']'

    mass = 4d0*!dpi*int_tabulated(rn,rhon*rn^2,/DOUBLE)

    ;write them out to a file
    openw, funit, 'background.dat', /get_lun
    printf, funit, format='(i5)', n
    ;printf, funit, format='(E16.9)', mass
    printf, funit, format='(E16.9)', gravn
    printf, funit, format='(E16.9)', rhon
    printf, funit, format='(E16.9)', Sn
    printf, funit, format='(E16.9)', dsdrn
    printf, funit, format='(E16.9)', mur
    If (keyword_set(nonuniform)) Then Begin
        dmudr = nud6(rn,mur)
    Endif Else Begin
        dmudr = compact_fd6(dri,mur,0,0,1,0)
        ;dmudr = ud6(dr,mur,1)
    EndElse
    printf, funit, format='(E16.9)', dmudr
    printf, funit, format='(E16.9)', kapparn
    If (keyword_set(nonuniform)) Then Begin
        dkappardr = nud6(rn,kapparn)
    Endif Else Begin
        dkappardr = compact_fd6(dri,kapparn,0,0,1,0)
        ;dkappardr = ud6(dr,kapparn,1)
    EndElse
    printf, funit, format='(E16.9)', dkappardr
    printf, funit, format='(E16.9)', kappas
    If (keyword_set(nonuniform)) Then Begin
        dkappasdr = nud6(rn,kappas)
    Endif Else Begin
        dkappasdr = compact_fd6(dri,kappas,0,0,1,0)
        ;dkappasdr = ud6(dr,kappas,1)
    EndElse
    printf, funit, format='(E16.9)', dkappasdr
    inds = where(abs(kappa0) lt 1d-10)
    If (inds[0] ne -1) Then kappa0(inds) = 0d0
    printf, funit, format='(E16.9)', kappa0
    If (keyword_set(nonuniform)) Then Begin
        dkappa0dr = nud6(rn,kappa0)
    Endif Else Begin
        dkappa0dr = compact_fd6(dri,kappa0,0,0,1,0)
    Endelse    
    inds = where(abs(dkappa0dr) lt 1d-10)
    If (inds[0] ne -1) Then dkappa0dr(inds) = 0d0
    inds = where(abs(epsn) lt 1d-10)
    epsn(inds) = 0d0
    printf, funit, format='(E16.9)', dkappa0dr
    printf, funit, format='(E16.9)', epsn
    printf, funit, format='(E16.9)', Lradn

    If (keyword_set(nonuniform)) Then Begin
        printf, funit, format='(E16.9)', rn
        dxdr=nud6(rn,xn)
        d2xdr2 = nud6(rn,dxdr)
        printf, funit, format='(E16.9)', dxdr
        printf, funit, format='(E16.9)', d2xdr2
    Endif

    If (keyword_set(evorad)) Then Begin
        rosskap = -(16d0*!dpi*a*c*Tn^3*rn^2*dtdr)/3d0/Lradn/rhon
        drkdr = compact_fd6(dri,rosskap,rosskap(0),rosskap(new_nr-1),1,0)
        printf, funit, format='(E16.9)', rosskap
        printf, funit, format='(E16.9)', drkdr
    EndIf
    close, funit
    free_lun, funit 

    If (keyword_set(fig1)) Then Begin
       ;Open a file dialog box to select (a) file(s).
        root_directory = '/freyr3/augustso/CSSDE/'
        files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
       ;Determine the number of files.
        n_files = n_elements(files)

       ;Handle the event that the user cancels.
        if(n_files eq 1) then begin
            if(files eq '') then return
        endif

       ;We assume the files read in are from the same run.
        read_shell_avg, files(0), radius = r, data = values, QUANTITIES = quantities

        num_records = n_elements(values(0,0,*))

       ;Determine the number of radial points.
        N_R = n_elements(r)

       ;Determine the number of quantities.
        N_Q = n_elements(quantities)+2
        quantities = [quantities,8,9]

        data = dblarr(N_R,N_Q,num_records*n_files)
        data(*,0:N_Q-3,0:num_records-1) = values

        if (n_files gt 1) then begin
            for i=1,n_files-1 do begin
                read_shell_avg, files(i), data = values
                num_rec = n_elements(values(0,0,*))
                data(*,0:N_Q-3,num_records:num_records+num_rec-1) = values
                num_records = num_records + num_rec
            endfor
        endif

        pres = dblarr(N_R,num_records)
        pres = gam1*Cp*data(*,0,*)*data(*,1,*)/gamma
        data(*,N_Q-2,*) = pres
        dr = r(1)-r(0)
        for i=0,num_records-1 do begin
            data(*,N_Q-1,i) = ud6(dr,data(*,2,i),1)
        endfor

        rho_m = dblarr(N_R)
        T_m = rho_m
        dsdr_m = rho_m
        for ir=0,N_R-1 do begin
            rho_m(ir) = mean(data(ir,0,*))
            T_m(ir) = mean(data(ir,1,*))
            dsdr_m(ir) = mean(data(ir,N_Q-1,*))
        endfor

        xsize=16.67d0
        ysize=10d0
        sunsym = sunsymbol()
        Set_Plot, "PS"
        docolor
        !P.MULTI=[0,1,1]
        !P.CHARSIZE=1
        Device, filename='rain_fig1.eps', /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize
        x0=0.1
        x1=0.48
        y0=0.12
        y1=0.97
        !P.POSITION=[x0,y0,x1,y1]
        !P.NOERASE=1 
        plot, r/max(radius), alog10(rho_m), thick=2, xtitle=textoidl('r/R_')+sunsym, ytitle=textoidl('log_{10} \rho')
        ;oplot, r/max(radius), alog10(rhon), thick=1.5, color=70
        oplot, radius/max(radius), alog10(rho), thick=2, color=115
        oplot, [r(0)/max(radius),r(0)/max(radius)],[-10,10]
        oplot, [r(n_r-1)/max(radius),r(n_r-1)/max(radius)],[-10,10]
        names = ['Sim.','1D']
        legend,names,COLORS=[0,120],PSYM=make_array(2,/INTEGER,value=0),pspacing=1.2,pos=[0.2,0.45],/normal,charsize=0.7,charthick=1

        x0=0.59
        x1=0.97
        !P.POSITION=[x0,y0,x1,y1]
        !P.NOERASE=1 
        plot, r/max(radius), alog10(T_m), thick=2, xtitle=textoidl('r/R_')+sunsym, ytitle=textoidl('log_{10} T'), yrange=[4,6]
        ;oplot, r/max(radius), alog10(Tn), thick=1.5, color=70
        oplot, radius/max(radius), alog10(T), thick=2, color=115
        oplot, [r(0)/max(radius),r(0)/max(radius)],[-10,10]
        oplot, [r(n_r-1)/max(radius),r(n_r-1)/max(radius)],[-10,10]

        DEVICE, /close
        SET_PLOT, 'x'
        print, 'Image ready.'
        SPAWN, 'kghostview rain_fig1.eps &'

        xsize=16.67d0
        ysize=10d0
        sunsym = sunsymbol()
        Set_Plot, "PS"
        docolor
        !P.MULTI=[0,1,1]
        !P.CHARSIZE=1
        Device, filename='nonrain_fig2.eps', /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize
        x0 = 0.1
        x1 = 0.97
        y0 = 0.12
        y1 = 0.97
        !P.POSITION=[x0,y0,x1,y1]
        plot, r/max(radius), dsdr_m, thick=2, xtitle=textoidl('r/R_')+sunsym, ytitle=textoidl('\partial S/\partial r')
        ;oplot, r/max(radius), alog10(rhon), thick=1.5, color=70
        oplot, radius/max(radius), dsdr, thick=2, color=115
        oplot, [r(0)/max(radius),r(0)/max(radius)],[-10,10]
        oplot, [r(n_r-1)/max(radius),r(n_r-1)/max(radius)],[-10,10]
        names = ['Sim.','1D']
        legend,names,COLORS=[0,120],PSYM=make_array(2,/INTEGER,value=0),pspacing=1.2,pos=[0.5,0.65],/normal,charsize=0.7,charthick=1

        DEVICE, /close
        SET_PLOT, 'x'
        print, 'Image ready.'
        SPAWN, 'kghostview nonrain_fig2.eps &'
    EndIf
    
 ;Figure out what the temperature and density are at r1-dr (for inflow
 ;outflow conditions).
    rnew = [rn(0)-dr,rn(*)]
    rhonew = interpol(rho,radius,rnew,/SPLINE)
    Tnew = interpol(T,radius,rnew,/SPLINE)
    print, 'Rho_rm1=', rhonew(0), 'T_rm1=', Tnew(0)

    If (keyword_set(stop)) Then Begin
       stop
    EndIf

end



