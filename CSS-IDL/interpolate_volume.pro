pro interpolate_volume, istep, dir, new_nth, new_nph

    zero_char = '0'
    i=5L
    k=100000L
    pretemp = ' '
    stemp = ' '
    While ((istep lt k) and (i gt 0L)) Do Begin
       short_temp = string(i,format='(i1)')
       format = '(i'+strtrim(short_temp,2)+')'
       fmt = strtrim(format,2)
       stemp = string(istep,format=fmt)
       pretemp = strtrim(pretemp+zero_char,2)
       i=i-1L
       k=k/10L 
    Endwhile
    
    fname = strtrim(dir+'checkpoint_'+strtrim(pretemp+stemp,2)+'_header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = intarr(5)
    readf, file_unit, temp,format='(5i6)'

    nth  = fix(temp(0))
    nph  = fix(temp(1))
    nr   = fix(temp(2))
    nq   = fix(temp(3))
    iter = fix(temp(4))

    print, temp

    temp = fltarr(8)
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

    fname = strtrim(dir+'checkpoint_'+strtrim(pretemp+stemp,2),2)
    
    openr, file_unit, fname, /get_lun, /swap_if_big_endian
    temp = assoc(file_unit,fltarr(nth,nph,nr,nq))
    volumes = fltarr(nth,nph,nr,nq)
    volumes(*,*,*,*) = temp(*,*,*,*,0)
    close, file_unit
    free_lun, file_unit

    clearline = fifteenb()

    ;Loop over r and interpolate theta-phi planes
    new_vol = fltarr(new_nth,new_nph,nr,nq)
    for iq=0L,nq-1L do begin
        for ir=0L,nr-1L do begin
            new_vol(*,*,ir,iq) = congrid(volumes(*,*,ir,iq),new_nth,new_nph,CUBIC=-0.5)
            per = strtrim(round((1.0d0*ir)/(1.0d0*(nr-1.0d0))*100.0d0),2)
            lenp = strtrim(strlen(strtrim(per,2)),2)
            form="($,a"+lenp+",' % Completed',a,a)"
            print, form=form, per, '         ', clearline
        endfor
        print, 'iq = ', iq, 'of', nq-1
    endfor

    theta = (th2-th1)*findgen(new_nth)/(new_nph-1d0)+th1
    plot, theta, new_vol(*,new_nph/2,nr/2,1)

    ;Write out the spherical data
    openw, file_unit, fname+'_new', /get_lun, /swap_if_big_endian

    a = assoc(file_unit,fltarr(new_nth,new_nph,nr,nq))
    a(0) = new_vol
    close, file_unit
    free_lun, file_unit 
end

function fifteenb

   return, string("15b)

end
