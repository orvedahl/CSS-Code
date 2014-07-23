@UniformDerivative6.pro
@NonUniformDerivative6.pro

pro build_flow_stats, iter, css_case, dir0=dir0, omega0=omega0, nonuniform=nonuniform, Cp=Cp, gamma=gamma, precomp=precomp, rmid=rmid

    if (not keyword_set(dir0)) then begin
       dir0 = '/freyr3/augustso/CSSDE/'+css_case+'/Checkpoints/'
    endif    

    k=1000000L
    If (iter lt k) Then Begin
        zero_char = '0'
        i=6L
        pretemp = ' '
        stemp = ' '
        While ((iter lt k) and (i gt 0L)) Do Begin
            short_temp = string(i,format='(i1)')
            format = '(i'+strtrim(short_temp,2)+')'
            fmt = strtrim(format,2)
            stemp = string(iter,format=fmt)
            pretemp = strtrim(pretemp+zero_char,2)
            i=i-1L
            k=k/10L
        Endwhile
    Endif Else Begin
        pretemp = ' '
        stemp = string(iter, format='(i7)')
    Endelse

    dir0 = strtrim(dir0,2)+strtrim(pretemp+stemp+'/',2)

    fname = strtrim(dir0+'header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format                                                                              
    temp = lonarr(5)
    readf, file_unit, temp,format='(5i9)'

    nr   = fix(temp(0))
    nth  = fix(temp(1))
    nph  = fix(temp(2))
    nq   = fix(temp(3))
    iter = long(temp(4))

    print, temp

    temp = dblarr(8)
    readf, file_unit, temp, format='(8e14.7)'
    dt   = temp(0)
    time = temp(1)
    r1   = temp(2)
    r2   = temp(3)
    th1  = temp(4)
    th2  = temp(5)
    ph1  = temp(6)
    ph2  = temp(7)

    if (not keyword_set(omega0)) then begin
        omega0 = 2.666d-6
    endif

;Read in background file
    openr, funit, dir0+'background.dat', /get_lun
    readf, funit, format='(i5)', n
    gravn = dblarr(n)
    rhon = dblarr(n)
    Sn = dblarr(n)
    dsdrn = dblarr(n)
    mur = dblarr(n)
    dmudr = dblarr(n)
    opacity = dblarr(n)
    dopdlnd = dblarr(n)
    dopdlnT = dblarr(n)
    kappas = dblarr(n)
    dkappasdr = dblarr(n)
    kappa0 = dblarr(n)
    dkappa0dr = dblarr(n)
    epsn = dblarr(n)
    Lradn = dblarr(n)
    readf, funit, format='(E16.9)', gravn
    readf, funit, format='(E16.9)', rhon
    readf, funit, format='(E16.9)', Sn
    readf, funit, format='(E16.9)', dsdrn
    readf, funit, format='(E16.9)', mur
    readf, funit, format='(E16.9)', dmudr
    readf, funit, format='(E16.9)', opacity
    readf, funit, format='(E16.9)', dopdlnd
    readf, funit, format='(E16.9)', dopdlnT
    readf, funit, format='(E16.9)', kappas
    readf, funit, format='(E16.9)', dkappasdr
    readf, funit, format='(E16.9)', kappa0
    readf, funit, format='(E16.9)', dkappa0dr
    readf, funit, format='(E16.9)', epsn
    readf, funit, format='(E16.9)', Lradn
    If (keyword_set(nonuniform)) Then Begin
        readf, funit, format='(E16.9)', rn
        dxdr=nud6(rn,xn)
        d2xdr2 = nud6(rn,dxdr)
        readf, funit, format='(E16.9)', dxdr
        readf, funit, format='(E16.9)', d2xdr2
    EndIf
    close, funit
    free_lun, funit 

    qnames = ['rho','u','v','w','s']

    fname = strtrim(dir0+'checkpoint.',2)

    data = dblarr(nr,nth,nph,nq)
    temp = dblarr(nth,nph)
 ;read in files                                                                                 
    for iq=0,nq-1 do begin
        tempname = strtrim(fname+qnames(iq),2)
        openr, file_unit, tempname, /get_lun, /swap_if_big_endian
        for ir=0L,nr-1L do begin
            temp(*,*) = 0d0
            readu, file_unit, temp
            data(ir,*,*,iq) = temp
        endfor
        close, file_unit
        free_lun, file_unit
    endfor

    ur = reform(data(*,*,*,1))
    ut = reform(data(*,*,*,3))
    up = reform(data(*,*,*,2))
    rho = reform(data(*,*,*,0))
    S   = reform(data(*,*,*,4))
    data = 0d0

    gam1 = gamma-1d0
    Cv = Cp/gamma

    enter = 0d0
    leave = 0d0
    tiny  = 1d-16
    rhoin = 0d0
    mfout = 0d0
    ir    = 174
    for it=0,nth-1 do begin
        for ip=0,nph-1 do begin
            leave = max([0d0,ur(ir,it,ip)])/max([tiny,ur(ir,it,ip)])
            enter = 1d0-leave
            rhoin = rhoin+enter*rho(ir,it,ip)
            mfout = mfout+leave*rho(ir,it,ip)*ur(ir,it,ip)
        endfor
    endfor

    uin = -mfout/rhoin

    s_in = 0d0
    Fke = 0.5d0*total(rho*ur*(ur^2+ut^2+up^2))
    Fenout = 0d0
    rhoin2 = 0d0
    lstar = 1.5640284d34
    Ftot = lstar*nth*nph/(4d0*!dpi*r2^2)
    T = rho^gam1*exp(S/Cv)/gam1/Cv
    Tn = rhon^gam1*exp(Sn/Cv)/gam1/Cv
    for it=0,nth-1 do begin
        for ip=0,nph-1 do begin
            leave = max([0d0,ur(ir,it,ip)])/max([tiny,ur(ir,it,ip)])
            enter = 1d0-leave
            rhoin2 = rhoin2 + enter*rho(ir,it,ip)^gamma
            Fenout = Fenout + leave*ur(ir,it,ip)*Cp*rho(ir,it,ip)*(T(ir,it,ip)-Tn(ir))
        endfor
    endfor

    temp1 = (Ftot-Fke-Fenout)*rhon(ir)^gam1/(uin*Cp*Tn(ir)*rhoin2)
    temp2 = rhon(ir)^gam1*rhoin/rhoin2
    s_in = Cv*alog(temp1+temp2)

stop

;Step two compute horizontal averages of Re, Ro, etc...
;Compute curl to find Rossby

    dr = (r2-r1)/double(nr-1)
    dth = (th2-th1)/double(nth-1)
    dph = (ph2-ph1)/double(nph-1)

    r = dr*dindgen(nr)+r1
    theta = dth*dindgen(nth)+th1
    phi = dph*dindgen(nph)+ph1

;Build up a list of statistics on the flow (snapshot)
    ir=nr-17
    hist = spherical_histogram(ur(ir,*,*), theta, locations=bins, /latitude, nbins=100)
    for ir=nr-16,nr-8 do begin
        hist = hist+spherical_histogram(ur(ir,*,*), theta, locations=bins, /latitude, nbins=100)
    endfor
    plot, bins, hist
    upos = (abs(ur)+ur)/2d0
    uneg = (ur-abs(ur))/2d0
    mpos = mean(upos)
    sdpos = stddev(upos,/double)
    mneg = mean(uneg)
    sdneg = stddev(uneg,/double)
    print, 'mpos=', mpos, ' sdpos=', sdpos
    print, 'mneg=', mneg, ' sdneg=', sdneg

    sint = sin(theta)
    cost = cos(theta)

    if (keyword_set(precomp)) Then Begin
        restore, dir0+'dsdr.sav'
        restore, dir0+'vort.sav'
        restore, dir0+'enst.sav'
        restore, dir0+'T.sav'
    Endif Else Begin
        dsdr = dblarr(nr,nth,nph)
        dutdr = dblarr(nr)
        dupdr = dblarr(nr)
    
        durdt = dblarr(nth)
        dupdt = dblarr(nth)

        durdp = dblarr(nph)
        dutdp = dblarr(nph)

        vort = dblarr(nr,nth,nph,3)
        T = dblarr(nr,nth,nph)

;Radial derivatives

        for ip=0,nph-1 do begin
            for it=0,nth-1 do begin
                dsdr(*,it,ip) = ud6(dr,S(*,it,ip),1)
                dutdr = ud6(dr,ut(*,it,ip),1)
                dupdr = ud6(dr,up(*,it,ip),1)
                vort(*,it,ip,1) = -dupdr-up(*,it,ip)/r
                vort(*,it,ip,2) = dutdr+ut(*,it,ip)/r
                T(*,it,ip) = rho(*,it,ip)^gam1*exp(S(*,it,ip)/Cv)/gam1/Cv
            endfor
        endfor
    
        print, 'radial derivatives done'
        
;Theta derivatives
        for ip=0,nph-1 do begin
            for ir=0,nr-1 do begin
                durdt = ud6(dth,ur(ir,*,ip),1)
                dupdt = ud6(dth,up(ir,*,ip),1)
                vort(ir,*,ip,0) = dupdt/r(ir)+cost*up(ir,*,ip)/sint/r(ir)
                vort(ir,*,ip,2) = vort(ir,*,ip,2)-durdt/r(ir)
            endfor
        endfor
        
        print, 'theta derivatives done'
        
;Phi derivatives
        for it=0,nth-1 do begin
            for ir=0,nr-1 do begin
                durdp = ud6(dph,ur(ir,it,*),1)
                dutdp = ud6(dph,ut(ir,it,*),1)
                vort(ir,it,*,0) = vort(ir,it,*,0)-dutdp/r(ir)/sint(it)
                vort(ir,it,*,1) = vort(ir,it,*,1)-durdp/r(ir)/sint(it)
            endfor
        endfor
        
        print, 'phi derivatives done'
        
;vorticity
        enst = sqrt(vort(*,*,*,0)^2+vort(*,*,*,1)^2+vort(*,*,*,2)^2)
        
 ;Save fields
        save, dsdr, filename=dir0+'dsdr.sav'
        save, vort, filename=dir0+'vort.sav'
        save, enst, filename=dir0+'enst.sav'
        save, T, filename=dir0+'T.sav'
    endelse

    print, 'vorticity done'

    enst_avg = dblarr(nr)
    Ro = dblarr(nr)
    cRo = dblarr(nr)
    Re = dblarr(nr)
    Ra = dblarr(nr)
    Re_fluc = dblarr(nr)
    Pr = dblarr(nr)
    Ek = dblarr(nr)
    ur_avg = dblarr(nr)
    ut_avg = dblarr(nr)
    up_avg = dblarr(nr)
    urms = dblarr(nr)
    urrms = dblarr(nr)
    urmin = dblarr(nr)
    urmax = dblarr(nr)
    utrms = dblarr(nr)
    uprms = dblarr(nr)
    ufluc_rms = dblarr(nr)
    ufluc = dblarr(nr,nth,nph,3)
    for ir=0,nr-1 do begin
        ur_avg = mean(ur(ir,*,*))
        ut_avg = mean(ut(ir,*,*))
        up_avg = mean(up(ir,*,*))
    endfor
    
    for ip=0,nph-1 do begin
        for it=0,nth-1 do begin
            ufluc(*,it,ip,0) = ur(*,it,ip)-ur_avg
            ufluc(*,it,ip,1) = ut(*,it,ip)-ut_avg
            ufluc(*,it,ip,2) = up(*,it,ip)-up_avg
        endfor
    endfor

    for ir=0,nr-1 do begin
        urrms(ir) = sqrt(total(ur(ir,*,*)^2)/nth/nph)
        urmin(ir) = min(ur(ir,*,*))
        urmax(ir) = max(ur(ir,*,*))
        utrms(ir) = sqrt(total(ut(ir,*,*)^2)/nth/nph)
        uprms(ir) = sqrt(total(up(ir,*,*)^2)/nth/nph)
        urms(ir) = sqrt(total(ur(ir,*,*)^2+ut(ir,*,*)^2+up(ir,*,*)^2)/nth/nph)
        ufluc_rms(ir) = sqrt(total(ufluc(ir,*,*,0)^2+ufluc(ir,*,*,1)^2+ufluc(ir,*,*,2)^2)/nth/nph)
    endfor

    nur = mur/rhon
    If (keyword_set(rmid)) Then Begin
        D = r2-rmid
    Endif Else Begin
        D = r2-r1
    EndElse
    for ir=0,nr-1 do begin
        enst_avg(ir) = mean(reform(enst(ir,*,*)))
        Ro(ir) = enst_avg(ir)/omega0/2d0
        Re(ir) = urms(ir)*D/nur(ir)
        Re_fluc(ir) = ufluc_rms(ir)*D/nur(ir)
        Ra(ir) = max([-mean(rho(ir,*,*)*T(ir,*,*)*dsdr(ir,*,*))*gravn(ir)*D^4/kappas(ir)/nur(ir)/Cp,1d0])
        Ek(ir) = nur(ir)/2d0/omega0/D^2
        Pr(ir) = nur(ir)*mean(rho(ir,*,*)*T(ir,*,*))/kappas(ir)
        cRo(ir) = sqrt(abs(Ra(ir))*Ek(ir)^2*Pr(ir))
    endfor

;Compute rms Ra, Re, Ro, etc...

    Ra_mean = mean(Ra)
    Re_mean = mean(Re)
    Re_fluc_mean = mean(Re_fluc)
    Ro_mean = mean(Ro)
    cRo_mean = mean(cRo)
    Ek_mean = mean(Ek)
    Pr_mean = mean(Pr)

    print, '<Ra>', Ra_mean, ' <Re>', Re_mean, " <Re'>", Re_fluc_mean, ' <Ro>', Ro_mean
    print, ' <conv Ro>', cRo_mean, '<Pr>', Pr_mean, '<Ek>', Ek_mean

;Write out results
    
    openw, file_unit, dir0+'stats.dat', /get_lun
    printf, funit, format='(E16.9)', Ra_mean
    printf, funit, format='(E16.9)', Re_mean
    printf, funit, format='(E16.9)', Re_fluc_mean
    printf, funit, format='(E16.9)', Ro_mean
    printf, funit, format='(E16.9)', cRo_mean
    printf, funit, format='(E16.9)', enst_avg
    printf, funit, format='(E16.9)', Ro
    printf, funit, format='(E16.9)', cRo
    printf, funit, format='(E16.9)', Re
    printf, funit, format='(E16.9)', Ra
    printf, funit, format='(E16.9)', Re_fluc
    printf, funit, format='(E16.9)', Pr
    printf, funit, format='(E16.9)', Ek
    printf, funit, format='(E16.9)', ur_avg
    printf, funit, format='(E16.9)', ut_avg
    printf, funit, format='(E16.9)', up_avg
    printf, funit, format='(E16.9)', urms
    printf, funit, format='(E16.9)', urrms
    printf, funit, format='(E16.9)', urmin
    printf, funit, format='(E16.9)', urmax
    printf, funit, format='(E16.9)', utrms
    printf, funit, format='(E16.9)', uprms
    printf, funit, format='(E16.9)', ufluc_rms
    close, file_unit
    free_lun, file_unit

    stop

;Done!

end
