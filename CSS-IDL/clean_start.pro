pro clean_start, iter, Cp, css_case=css_case

    ;Remove perturbations from rho, make each surface constant pressure by
    ;setting drho = -rho*ds/cp

    ;Open s and rho checkpoints
    dir = '/freyr3/augustso/CSSDE/'+css_case+'/Checkpoints/'+iter+'/'

    fname = strtrim(dir+'header',2)

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

    print, temp

    close, file_unit
    free_lun, file_unit

    frho = dir+'checkpoint.rho'
    fs = dir+'checkpoint.s'

    temp = dblarr(nth,nph)
    rho = dblarr(nr,nth,nph)
    openr, file_unit, frho, /get_lun, /swap_if_big_endian
    for ir=0L,nr-1L do begin
        temp(*,*) = 0d0
        readu, file_unit, temp
        rho(ir,*,*) = temp
    endfor
    close, file_unit
    free_lun, file_unit

    temp = dblarr(nth,nph)
    s = dblarr(nr,nth,nph)
    openr, file_unit, fs, /get_lun, /swap_if_big_endian
    for ir=0L,nr-1L do begin
        temp(*,*) = 0d0
        readu, file_unit, temp
        s(ir,*,*) = temp
    endfor
    close, file_unit
    free_lun, file_unit

    smean = dblarr(nr)
    rhomean = dblarr(nr)
    ds = dblarr(nr,nth,nph)
    for ir=0,nr-1L do begin
        smean(ir) = mean(s(ir,*,*))
        rhomean(ir) = mean(rho(ir,*,*))
        rho(ir,*,*) = rhomean(ir) ;zero out perturbations
        ds(ir,*,*) = s(ir,*,*)-smean(ir)
    endfor

    rho = rho*(1d0-ds/Cp) ;Add in constant pressure perturbations

    ;write out new rho
    openw, file_unit, frho+'_mod', /get_lun, /swap_if_big_endian

    plane = dblarr(nth,nph)
    for ir=0L,nr-1L do begin
        plane = reform(rho(ir,*,*))
        writeu, file_unit, plane
    endfor
    close, file_unit
    free_lun, file_unit
stop
end
