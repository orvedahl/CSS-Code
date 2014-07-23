pro read_checkpoint, iter, css_case=css_case, data=data, rad=rad, theta=theta, phi=phi, time=time

  pretemp = ' '
  stemp = string(iter, format='(i7.7)')

    
    iterstr = strtrim(pretemp+stemp,2)
    read_dir = '/freyr1/augustso/CSSDE/'+css_case+'/Checkpoints/'+iterstr+'/'

    fname = strtrim(read_dir+'header',2)
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
    t    = temp(1)
    time = dblarr(2)
    time(0) = dt
    time(1) = t
    r1   = temp(2)
    r2   = temp(3) 
    th1  = temp(4)
    th2  = temp(5)
    ph1  = temp(6)
    ph2  = temp(7)

    print, temp

    close, file_unit
    free_lun, file_unit

    qnames = ['rho','w','v','u','s','bt','bp','br']
    files = read_dir+'checkpoint.'+qnames
    data = dblarr(nth,nph,nr,nq)
    for iv=0,nq-1 do begin
        openr, file_unit, files(iv), /get_lun, /swap_if_big_endian
        a = assoc(file_unit,dblarr(nth,nph,nr))
        data[*,*,*,iv] = a[0]
        close, file_unit
        free_lun, file_unit
    endfor

    rad = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
    phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1

stop

end
