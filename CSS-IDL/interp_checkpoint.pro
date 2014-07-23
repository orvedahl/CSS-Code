;Interpolate onto a new grid
;Only does this in the horizontal
;If starting a new case where nr needs to be changed, will need to
;extract fluctuations from background and apply to new checkpoint

@interp.pro

pro interp_checkpoint, iter, newbounds, css_case1=css_case1, keep_old=keep_old

    If (not keyword_set(css_case2)) Then css_case2='solar_mag4'

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
    read_dir = '/freyr1/augustso/CSSDE/'+css_case1+'/Checkpoints/'+iterstr+'/'
    write_dir = '/freyr1/augustso/CSSDE/'+css_case1+'/Checkpoints/'+iterstr+'_interp/'
    ;write_dir = '/freyr3/augustso/CSSDE/'+css_case2+'/Checkpoints/'+iterstr+'/'

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

    qnames = ['rho','w','v','u','s','br','bt','bp']
    files = read_dir+'checkpoint.'+qnames
    data = dblarr(nth,nph,nr,nq)
    temp = dblarr(nth,nph,nr)
    for iv=0,nq-1 do begin
        openr, file_unit, files(iv), /get_lun, /swap_if_big_endian
        readu, file_unit, temp
        data[*,*,*,iv] = temp
        close, file_unit
        free_lun, file_unit
    endfor
    stop

    nth2  = fix(newbounds(0))
    nph2  = fix(newbounds(1))
    If (keyword_set(keep_old)) Then Begin
        newth1  = th1
        newth2  = th2
        newph1  = ph1
        newph2  = ph2
    Endif Else Begin
        newth1  = newbounds(2)
        newth2  = newbounds(3)
        newph1  = newbounds(4)
        newph2  = newbounds(5)
    Endelse
    nr2 = fix(newbounds(6)) ;Cannot change r1, r2
    dph = (ph2-ph1)/double(nph-1L)
    oldphi = dph*dindgen(nph)+ph1
    dth = (th2-th1)/double(nth-1L)
    oldtheta = dth*dindgen(nth)+th1
    dr = (r2-r1)/double(nr-1L)
    oldrad = dr*dindgen(nr)+r1

    dph = (newph2-newph1)/double(nph2-1L)
    newphi = dph*dindgen(nph2)+newph1
    dth = (newth2-newth1)/double(nth2-1L)
    newtheta = dth*dindgen(nth2)+newth1
    dr2 = (r2-r1)/double(nr2-1)
    newrad = dr2*dindgen(nr2)+r1

    data2 = dblarr(nth2,nph2,nr2,nq)

    If (nr2 ne nr) Then Begin
       for iq=0,nq-1 do begin
          Print, 'on quantity:', qnames[iq]
          interp3d, data[*,*,*,iq], oldphi, oldtheta, oldrad, newphi, newtheta, newrad, data2[*,*,*,iq]
       endfor
    Endif Else Begin
       for iq=0,nq-1 do begin
          Print, 'on quantity:', qnames[iq]
          for ir=0,nr-1 do begin
             interp2d, data[*,*,ir,iq], oldphi, oldtheta, newphi, newtheta, data2[*,*,ir,iq]
          endfor
       endfor
    EndElse

    ;Write out to new checkpoint
    fname = strtrim(write_dir+'header',2)

    openw, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = lonarr(5)
    temp(0) = nr2
    temp(1) = nth2
    temp(2) = nph2
    temp(3) = nq
    temp(4) = iter
    printf, file_unit, temp,format='(5i9)'

    temp = dblarr(8)
    temp(0) = dt
    temp(1) = time
    temp(2) = r1
    temp(3) = r2
    temp(4) = newth1
    temp(5) = newth2
    temp(6) = newph1
    temp(7) = newph2
    printf, file_unit, temp, format='(8e14.7)'

    close, file_unit
    free_lun, file_unit

    ;Ensure BC
    ;Radial velocity
    data2[*,*,0,3] = 0d0
    data2[*,*,nr2-1,3] = 0d0

    ;Radial field zero on r, perfect conductor
    data2[*,*,0,7] = 0d0
    data2[*,*,nr2-1,7] = 0d0

    ;Theta velocity
    data2[0,*,*,1] = 0d0
    data2[nth2-1,*,*,1] = 0d0

    ;Radial,Phi Field on theta
    data2[0,*,*,6:7] = 0d0
    data2[nth2-1,*,*,6:7] = 0d0

    files = write_dir+'checkpoint.'+qnames
    for iv=0,nq-1L do begin
        openw,file_unit,files(iv),/get_lun,/swap_if_big_endian        
        temp = data2[*,*,*,iv]
        writeu, file_unit, temp
        close, file_unit
        free_lun, file_unit
    endfor
stop
end
