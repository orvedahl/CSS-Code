pro extract3d, iter, vdfname, css_case=css_case

    dir1 = '/freyr1/augustso/CSSDE/'+css_case+'/Checkpoints/'
    dir2 = dir1

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
   
    dir1 = strtrim(dir1,2)+strtrim(pretemp+stemp+'/',2)
    dir2 = strtrim(dir2,2)+strtrim(pretemp+stemp+'/',2)+'Spherical/'

    fname = strtrim(dir1+'header',2) 

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

    qnames = ['rho','u','v','w','s','br','bt','bp','rho_fluc','s_fluc']
    names = qnames+'.dat'

    fname = strtrim(dir1+'checkpoint.',2)

    data = dblarr(nth,nph,nr,nq+2)
    temp = dblarr(nth,nph,nr)
    ;read in files
    for iq=0,nq-1 do begin
        tempname = strtrim(fname+qnames[iq],2)
        openr, file_unit, tempname, /get_lun, /swap_if_big_endian
        readu, file_unit, temp
        data[*,*,*,iq] = temp
        close, file_unit
        free_lun, file_unit
    endfor

    stop

    rho_avg = dblarr(nr)
    s_avg = dblarr(nr)

    for ir=0,nr-1 do begin
        rho_avg[ir] = mean(data[*,*,ir,0])
        s_avg[ir] = mean(data[*,*,ir,4])
        data[*,*,ir,nq] = data[*,*,ir,0]-rho_avg[ir]
        data[*,*,ir,nq+1] = data[*,*,ir,4]-s_avg[ir]
    endfor

    for iq=0,nq+1 do begin        
        fout = strtrim(dir1+names(iq),2)
        ;Write out the spherical data
        openw, file_unit, fout, /get_lun, /swap_if_big_endian
        writeu, file_unit, transpose(data[*,*,*,iq],[1,0,2])
        close, file_unit
        free_lun, file_unit 
    endfor

    th1 = th1*180d0/!PI-90d0
    th2 = th2*180d0/!PI-90d0
    ph1 = ph1*180d0/!PI
    ph2 = ph2*180d0/!PI

    extstr = strtrim(string(ph1,format='(f8.2)'),2)+':'+strtrim(string(th1,format='(f8.2)'),2)+':'+strtrim(string(r1/r2,format='(f8.2)'),2)+':' $
             +strtrim(string(ph2,format='(f8.2)'),2)+':'+strtrim(string(th2,format='(f8.2)'),2)+':'+strtrim(string(1.0,format='(f8.2)'),2)

    sizestr = strtrim(string(nph,format='(i4)'),2)+'x'+strtrim(string(nth,format='(i4)'),2)+'x'+strtrim(string(nr,format='(i4)'),2)

    cmd = 'vdfcreate -coordsys "spherical" -dimensio ' + strtrim(sizestr,2) + ' -numts 1 -level 2 -order 0:1:2 -extents ' $
           + strtrim(extstr,2) + ' -varnames rho:u:v:w:s:br:bt:bp:rho_fluc:s_fluc ' + strtrim(dir2+vdfname,2)

    print, cmd
    spawn, cmd

    for iq=0,nq+1 do begin 
        cmd = 'raw2vdf -dbl -varname ' + strtrim(qnames[iq],2) + ' ' +strtrim(dir2+vdfname,2) + ' ' + strtrim(dir1+names[iq],2)
        print, cmd
        spawn, cmd
    endfor

end
