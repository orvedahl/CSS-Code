@NonUniformDerivative6.pro
@UniformDerivative6.pro
@hydrostatic_solver.pro

pro build_background, filename, imin=imin, newnr1=newnr1, newnr2=newnr2, dual_cheby=dual_cheby, rmid=rmid, dencon=dencon, docss=docss, doash=doash, stop=stop, compg=compg

    ;Constants
    msun = 1.98892d33
    rsun = 6.955d10
    lsun = 3.839d33
    
    a = 7.56577d-15
    c = 2.99792458d10
    Cp = 3.5d8
    gamma = 5.0d0/3.0d0
    openr, funit, filename, /get_lun

    ;read in first line of junk
    data = string(128)
    readf, funit, data, format='(a128)'
    ;The second line contains nzones
    nzones = lonarr(1)
    readf, funit, nzones, format='(i4)'
    nzones = long(reform(nzones))
    nzones = nzones(0)
    ;read in variable names
    readf, funit, data, format='(a128)'
    ;break up the string by blank spaces    
    qnames = strsplit(data,' ',/EXTRACT)
    print, qnames
  
    ;read in blank line
    readf, funit, data, format='(a128)'
  
    ;read in the data
    idata = intarr(1)
    data = dblarr(13)
    values = dblarr(nzones,13)
    zones = lonarr(nzones)
    for i=0L,nzones-1L do begin
        readf, funit, idata, data, format='(i4,13e16.9)'
        zones(i) = long(idata(0))
        values(i,*) = data(*)
    endfor

    close, funit
    free_lun, funit

    mr = msun*(1.0d0-values(*,0))
    r = rsun*10.0d0^(values(*,1))
    P = 10.0d0^(values(*,2))
    T = 10.0d0^(values(*,3))
    L = lsun*values(*,4)
    rho = 10.0d0^(values(*,5))
    kappa = 10.0d0^(values(*,6))
    eps = 10.0d0^(values(*,7))
    del = values(*,8)
    delad = values(*,9)
    delrad = values(*,10)
    x = values(*,11)
    y = values(*,12)

    ;Test convergence of model, i.e. is it hydrostatic?
    ;First compute the gravity
    bigG = 6.673d-8
    grav = bigG*mr/r^2

    if keyword_set(plots) then begin
       !P.CHARSIZE=2.0
       window, 1,xsize=1280, ysize=640
       wset, 1
       !P.MULTI=[0,3,2]
       plot, r, mr, /xs, yrange=[0.0d0,1.1d0*msun], title=textoidl('M_r')
       plot, r, P, /xs, yrange=[0.0d0,1.1d0*P(0)], title='P'
       plot, r, T, /xs, yrange=[0.0d0,1.1d0*T(0)], title='T'
       plot, r, rho, /xs, yrange=[0.0d0,1.1d0*rho(0)], title=textoidl('\rho')
       plot, r, L, /xs, yrange=[0.0d0,1.1d0*lsun], title='L'
       plot, r, kappa, /xs, title=textoidl('\kappa')
   
       window, 2,xsize=1280-fix(1280/3.0d0), ysize=(1280-fix(1280/3.0d0))/2.0d0
       wset,2
       !P.CHARSIZE=1.0
       !P.MULTI=[0,2,2]
       plot, r, eps, /xs, yrange=[0.0d0,1.1d0*eps(0)], title=textoidl('\epsilon')
       plot, r, del, /xs, yrange=[0.0d0,1.1d0*max(del)], title='Del'
       plot, r, delad, /xs, yrange=[0.0d0,1.1d0*max(delad)], title=textoidl('Del_{ad}')
       plot, r, delrad, /xs, yrange=[0.0d0,1.1d0*max(delrad)], title=textoidl('Del_{rad}')

       window, 3
       wset,3
       !P.MULTI=0
       !P.CHARSIZE=1.5
       plot, r, grav
       ;dpdr = (P-shift(P,1)+a*(T^4-shift(T^4,1))/3.0d0)/(r-shift(r,1))
       err = rho*grav + NonUniformDerivative6(r,P+a*T^4/3.0d0)

       window, 4
       wset,4
       !P.MULTI=0
       !P.CHARSIZE=1.5
       plot, r, alog10(abs(err/(rho*grav)))
    endif

    ;Set up a regularly gridded model for CSS
    if keyword_set(docss) then begin
       n = newnr1
       ;Build entropy gradient
       gam1 = gamma-1.0d0
       dsdr = Cp*(del-delad)*NonUniformDerivative6(r,alog(P))
       imin = min(where(dsdr lt 0.0d0))
       rmin = r(imin)

       ;calculate rmax based on desired density contrast
       rho0 = rho(imin)
       imax = min(where(rho0/rho gt dencon))
       rmax = r(imax)
      
       dr = (rmax-rmin)/(double(n)-1d0)
       rn = dr*dindgen(n)+rmin
       kappar = 4.0d0*a*c*T^3/rho^2/Cp/kappa/3.0d0
       kapparn = interpol(kappar,r,rn,/SPLINE)

       rhon = hydrostatic_solver(imin,imax,r,rho,grav,dsdr,Cp,gamma,P,n,/css, Pn=Pn, Tn=Tn, gravn=gravn, /compg)

       ;write them out to a file
       openw, funit, 'background.dat', /get_lun
       printf, funit, format='(i4)', n
       printf, funit, format='(E16.9)', gravn
       printf, funit, format='(E16.9)', rhon
       printf, funit, format='(E16.9)', Tn
       Cpn = dblarr(n)
       Cpn(*) = Cp
       printf, funit, format='(E16.9)', Cpn
       dcpndr = dblarr(n)
       dcpndr(*) = 0d0
       printf, funit, format='(E16.9)', dcpndr
       printf, funit, format='(E16.9)', kapparn
       dlnkappadr = UniformDerivative6(dr,alog(kapparn))
       printf, funit, format='(E16.9)', dlnkappadr
       close, funit
       free_lun, funit
  
       ;Check kappar and kappas
       frad = -rhon*Cp*kapparn*UniformDerivative6(dr,Tn)
       plot, rn, 4d0*!PI*rn^2*frad/lsun
       oplot, rn, make_array(n,value=1d0,/DOUBLE)

    endif

    ;Set up a Chebyshev gridded model for ASH
    if keyword_set(doash) then begin
       if keyword_set(dual_cheby) then begin
          n = newnr1+newnr2
       endif else begin
          n = newnr1
       endelse
       ;Build entropy gradient
       gam1 = gamma-1.0d0
       dsdr = Cp*(del-delad)*NonUniformDerivative6(r,alog(P))
       imin = min(where(dsdr lt 0.0d0))
       rmin = r(imin)

       ;calculate rmax based on desired density contrast
       rho0 = rho(imin)
       imax = min(where(rho0/rho gt dencon))
       rmax = r(imax)

       dr = (rmax-rmin)/(double(n)-1d0)
       rn = dr*dindgen(n)+rmin
       kappar = 4.0d0*a*c*T^3/rho^2/Cp/kappa/3.0d0
       kapparn = interpol(kappar,r,rn,/SPLINE)
       epsn = interpol(eps,r,rn,/SPLINE)
       Ln = interpol(L,r,rn,/SPLINE)

       if keyword_set(dual_cheby) then begin
          rhon = hydrostatic_solver(imin,imax,r,rho,grav,dsdr,Cp,gamma,P,newnr1,/ash,n2=nrnew2,dual_cheby=dual_cheby,rmid=rmid, Pn=Pn, Tn=Tn, Sn=Sn, dsdrn=dsdrn)
       endif else begin
          rhon = hydrostatic_solver(imin,imax,r,rho,grav,dsdr,Cp,gamma,P,newnr1,/ash, Pn=Pn, Tn=Tn, Sn=Sn, dsdrn=dsdrn)
       endelse

       ;write them out to a file
       openw, funit, 'background.dat', /get_lun
       printf, funit, format='(i4)', n
       for i=0L,n-1L do begin
           printf, funit, format='(10E16.9)', rn(i), gravn(i), rhon(i), Tn(i), Pn(i), Sn(i), Ln(i), dsdrn(i), kapparn(i), epsn(i)
       endfor
       close, funit
       free_lun, funit

       ;From HKT pg. 307
       ;alpha = 50.8d0/(1d-6*T(0))^(1.0d0/3.0d0)-2.0d0/3.0d0
       ;Compute total energy for unnormalized energy generation terms
       ;Leps = int_tabulated(r(0:imin),4.0d0*!PI*r(0:imin)^2*rho(0:imin)*(1d-6*T(0:imin))^alpha,/DOUBLE)
       ;print, 'The value by which the value computed in ASH is normalized by is', Leps
       ;dtdr = NonUniformDerivative6(r,T)
       ;print, 'The first energy generation norm is ', (L(imin)+4.0d0*!PI*r(imin)^2*4.0d0*a*c*T(imin)^3/rho(imin)/kappa(imin)/3.0d0*dtdr(imin))/Leps
    endif

    if keyword_set(stop) then stop

end
