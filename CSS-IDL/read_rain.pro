pro read_rain, iter, css_case

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

    ;read in old checkpoint
    read_dir = '/freyr1/augustso/CSSDE/'+css_case+'/Checkpoints/'+iterstr+'/'

    fname = strtrim(read_dir+'header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = lonarr(6)
    readf, file_unit, temp,format='(6i9)'

    nr   = fix(temp[0])
    nth  = fix(temp[1])
    nph  = fix(temp[2])
    nq   = fix(temp[3])
    iter = long(temp[4])
    nsqs = long(temp[5])

    print, temp

    temp = dblarr(8)
    readf, file_unit, temp, format='(8e14.7)'
    dt   = temp[0]
    time = temp[1]
    r1   = temp[2]
    r2   = temp[3]
    th1  = temp[4]
    th2  = temp[5]
    ph1  = temp[6]
    ph2  = temp[7]

    print, temp

    close, file_unit
    free_lun, file_unit

    fname = strtrim(read_dir+'checkpoint.rain_squalls',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian
    a = assoc(file_unit, dblarr(13))
    readbuf = dblarr(13,nsqs)
    for ir=0,nsqs-1 do begin
       readbuf[*,ir] = a[ir]
    endfor
    close, file_unit
    free_lun, file_unit

    duration = readbuf[0,*]
    sigmaT = readbuf[1,*]
    start_time = readbuf[2,*]
    rainsize = readbuf[3,*]
    dels = readbuf[4,*]
    thetaloc = readbuf[5,*]
    philoc = readbuf[6,*]
    scaling = readbuf[7,*]
    iloct = fix(readbuf[9,*])
    ilocp = fix(readbuf[10,*])
    vel = fix(readbuf[11,*])

    stop

end
