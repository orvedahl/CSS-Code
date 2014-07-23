pro read_checkpoint_file, iter, css_case, iq, dir0=dir0

    if (not keyword_set(dir0)) then begin
       dir0 = '/freyr3/augustso/CSSDE/'+css_case+'/Checkpoints/'
    endif    

    pretemp = ' '
    stemp = string(iter, format='(i7.7)')

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

    qnames = ['rho','w','v','u','s','bt','bp','br']

    fname = strtrim(dir0+'checkpoint.',2)

    data = dblarr(nr,nth,nph)
    temp = dblarr(nth,nph)
 ;read in files                                                                                 

    tempname = strtrim(fname+qnames(iq),2)
    openr, file_unit, tempname, /get_lun, /swap_if_big_endian
    for ir=0L,nr-1L do begin
        temp(*,*) = 0d0
        readu, file_unit, temp
        data(ir,*,*) = temp
    endfor
    close, file_unit
    free_lun, file_unit
stop
end
