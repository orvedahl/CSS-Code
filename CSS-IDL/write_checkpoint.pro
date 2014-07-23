pro write_checkpoint, iter, css_case=css_case, data=data, rad=rad, theta=theta, phi=phi, time=time, tag=tag

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
        EndWhile
    EndIf Else Begin
        pretemp = ' '
        stemp = string(iter, format='(i7)')
    EndElse
    
    iterstr = strtrim(pretemp+stemp,2)
    write_dir = '/freyr3/augustso/CSSDE/'+css_case+'/Checkpoints/'+iterstr+tag+'/'

    fname = strtrim(write_dir+'header',2)
    openw, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = lonarr(5)
    nr = n_elements(rad)
    nth = n_elements(theta)
    nph = n_elements(phi)
    nq = n_elements(data(0,0,0,*))

    temp(0) = nr
    temp(1) = nth
    temp(2) = nph
    temp(3) = nq
    temp(4) = iter

    printf, file_unit, temp,format='(5i9)'

    temp = dblarr(8)
    dt = time(0)
    t  = time(1)
    r1 = min(rad)
    r2 = max(rad)
    th1 = min(theta)
    th2 = max(theta)
    ph1 = min(phi)
    ph2 = max(phi)

    temp(0) = dt
    temp(1) = t
    temp(2) = r1
    temp(3) = r2
    temp(4) = th1
    temp(5) = th2
    temp(6) = ph1
    temp(7) = ph2

    printf, file_unit, temp, format='(8e14.7)'

    close, file_unit
    free_lun, file_unit

    qnames = ['rho','w','v','u','s','br','bt','bp']
    files = write_dir+'checkpoint.'+qnames
    data = dblarr(nr,nth,nph,nq)
    for iv=0,nq-1 do begin
        temp = dblarr(nth,nph)
        openw, file_unit, files(iv), /get_lun, /swap_if_big_endian
        for ir=0L,nr-1L do begin
            temp(*,*) = 0d0
            writeu, file_unit, temp
            data(ir,*,*,iv) = temp
        endfor
        close, file_unit
        free_lun, file_unit
    endfor

stop

end
