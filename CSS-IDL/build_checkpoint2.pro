pro build_checkpoint2, iter, nr, nth, nphi, nnr=nnr, nnt=nnt, nnp=nnp, read_dir=read_dir, write_dir=write_dir, reiter=reiter

    If (not keyword_set(read_dir)) Then read_dir='/freyr3/augustso/CSSDE/dir_n1t3_ld/Checkpoints/'

    If (not keyword_set(write_dir)) Then write_dir='/freyr3/augustso/CSSDE/dir_n1t3_ld/Checkpoints/Rebuilt/'

    If (keyword_set(reiter)) Then Begin
        k=1000000L
        If (reiter lt k) Then Begin
            zero_char = '0'
            i=6L
            pretemp = ' '
            stemp = ' '
            While ((iter lt k) and (i gt 0L)) Do Begin
                short_temp = string(i,format='(i1)')
                format = '(i'+strtrim(short_temp,2)+')'
                fmt = strtrim(format,2)
                stemp = string(reiter,format=fmt)
                pretemp = strtrim(pretemp+zero_char,2)
                i=i-1L
                k=k/10L
            EndWhile
        EndIf Else Begin
            pretemp = ' '
            stemp = string(reiter, format='(i7)')
        EndElse

        reiter_string = strtrim(pretemp+stemp,2)
        write_dir = strtrim(write_dir,2)+strtrim(reiter_string+'/',2)
    Endif Else Begin
        write_dir = strtrim(write_dir,2)+strtrim(iter_string+'/',2)
    Endelse

    k=100000L
    If (iter lt k) Then Begin
        zero_char = '0'
        i=5L
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
        stemp = string(iter, format='(i6)')
    EndElse

    iter_string = strtrim(pretemp+stemp,2)

    read_dir = strtrim(read_dir,2)+strtrim(iter_string+'/',2)

    exists = file_test(read_dir, /directory)
    If (exists eq 0) Then Begin 
        Print, 'Read directory does not exist.' 
        stop 
    Endif



    exists = file_test(write_dir, /directory)
    If (exists eq 0) Then Begin 
        Print, 'Write directory does not exist.' 
        Print, 'Making write directory...'
        file_mkdir, write_dir
        file_chmod, write_dir, '777'o
    EndIf

    ;Read in first header
    fname = strtrim(read_dir+'checkpoint_'+strtrim(iter_string,2)+'.',2)
    fhead = strtrim(read_dir+'checkpoint_'+strtrim(iter_string,2)+'_header.',2)
    filename = strtrim(fhead+'00000',2)
    openr, file_unit, filename, /get_lun, /swap_if_big_endian


    temp = intarr(11)
    readf, file_unit, temp,format='(11i6)'
;    temp = intarr(5)
;    readf, file_unit, temp,format='(5i6)'

    mynr = fix(temp(0))
    mynth  = fix(temp(1))
    mynphi  = fix(temp(2))
    nv   = fix(temp(3))
    myr1 = fix(temp(4))
    myr2 = fix(temp(5))
    myth1 = fix(temp(6))
    myth2 = fix(temp(7))
    myphi1 = fix(temp(8))
    myphi2 = fix(temp(9))

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
    
    bounds = [nr,nth,nphi,r1,r2,th1,th2,ph1,ph2]

    print, 'bounds = [', bounds, ']'

    full_data = dblarr(nr,nth,nphi,nv)

    nnodes = nnr*nnt*nnp

    For in=0L,nnodes-1L Do Begin

        ;Read in headers
        k=10000
        If (in lt k) Then Begin
            i=0L
            zero_char = '0'
            i=4L
            pretemp = ' '
            While ((in lt k) and (i gt 0L)) Do Begin
                short_temp = string(i,format='(i1)')
                format = '(i'+strtrim(short_temp,2)+')'
                fmt = strtrim(format,2)
                stemp = string(in,format=fmt)
                pretemp = strtrim(pretemp+zero_char,2)
                i=i-1L
                k=k/10L
            EndWhile
        Endif Else Begin
            pretemp = ' '
            stemp = string(in, format='(i5)')
        EndElse
       
                                ;Read in first header
        fhead = strtrim(read_dir+'checkpoint_'+strtrim(iter_string,2)+'_header.'+strtrim(pretemp+stemp,2),2)
        openr, file_unit, fhead, /get_lun, /swap_if_big_endian
 
        temp = intarr(11)
        readf, file_unit, temp,format='(11i6)'
        
        mynr = fix(temp(0))
        mynth  = fix(temp(1))
        mynphi  = fix(temp(2))
        nv   = fix(temp(3))
        myr1 = fix(temp(4))
        myr2 = fix(temp(5))
        myth1 = fix(temp(6))
        myth2 = fix(temp(7))
        myphi1 = fix(temp(8))
        myphi2 = fix(temp(9))

        close, file_unit
        free_lun, file_unit

        ;Read in data

        filename = strtrim(fname,2)+strtrim(pretemp+stemp,2)
        openr, file_unit, filename, /get_lun, /swap_if_big_endian

        read_buff = dblarr(mynth,mynphi)
        For iv=0,nv-1 Do Begin
            For ir=0,mynr-1 Do Begin
                read_buff(*,*) = 0d0
                readu,file_unit, read_buff
                full_data(myr1+ir-1L,myth1-1L:myth2-1L,myphi1-1L:myphi2-1L,iv) = read_buff
            Endfor
        Endfor
        close, file_unit
        free_lun, file_unit

    Endfor

    ;Write out header

    fhead = strtrim(strtrim(write_dir)+'header',2)

    openw, file_unit, fhead, /get_lun, /swap_if_big_endian

    temp = lonarr(5)
    temp(0) = nr
    temp(1) = nth
    temp(2) = nphi
    temp(3) = nv
    If (keyword_set(reiter)) Then Begin
        temp(4) = reiter
    Endif Else Begin
        temp(4) = iter
    EndElse

    printf, file_unit, temp, format='(5i9)'

    temp = dblarr(8)
    temp(0) = dt
    temp(1) = time
    temp(2) = r1
    temp(3) = r2
    temp(4) = th1
    temp(5) = th2
    temp(6) = ph1
    temp(7) = ph2

    printf, file_unit, temp, format='(8e14.7)'

    close, file_unit
    free_lun, file_unit

    ;Write out data
    names = ['rho','w','v','u','s','br','bt','bp']
    for iv=0L,nv-1L do begin
        fname = strtrim(strtrim(write_dir)+'checkpoint.'+strtrim(names(iv),2),2)
        
        openw, file_unit, fname, /get_lun, /swap_if_big_endian

        plane = dblarr(nth,nphi)
        for ir=0L,nr-1L do begin
            plane = reform(full_data(ir,*,*,iv))
            writeu, file_unit, plane
        endfor
        close, file_unit
        free_lun, file_unit
    endfor

end
