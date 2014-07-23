pro read_scalar, filename, values=values, iters=iters, ie=ie, mkuniq=mkuniq, mag=mag, imag=imag

    niters = file_lines(filename)

    If (keyword_set(imag)) Then Begin
       niters = niters-imag
    EndIf

    openr, funit, filename, /get_lun, /swap_if_big_endian

    print, 'Reading niters=', niters 

    iters = lindgen(niters)+1L
    ltmp = 0L
    If (keyword_set(mag)) Then Begin
       values = dblarr(9,niters)
       If (not keyword_set(imag)) Then Begin
          temp = dblarr(10)
          for it=0L,niters-1L do begin
             temp(*) = 0d0
             readf, funit, temp, format='(i10,9E15.7)'
             iters(it) = long(temp(0))
             values(*,it) = temp(1:9)
          endfor
       Endif Else Begin
          print, 'Skipping ', imag, ' lines.'
          temp = dblarr(8,imag)
          readf, funit, temp, format='(i10,7E15.7)'

          temp = dblarr(10)
          for it=0L,niters-1L do begin
             temp(*) = 0d0
             readf, funit, temp, format='(i10,9E15.7)'
             iters(it) = long(temp(0))
             values(*,it) = temp(1:9)
          endfor
       EndElse
    Endif Else Begin
       values = dblarr(7,niters)
       temp = dblarr(8)

       for it=0L,niters-1L do begin
          temp(*) = 0d0
          readf, funit, temp, format='(i10,7E15.7)'
          iters(it) = long(temp(0))
          values(*,it) = temp(1:7)
       endfor
    EndElse 
    close, funit
    free_lun, funit

    If (keyword_set(mkuniq)) Then Begin
        iteruniq = uniq(iters, sort(iters))
        values = values(*,iteruniq)
        iters = iters(iteruniq)
        niters = n_elements(iters)
        ;Write out fixed file:
        openw, funit, filename, /get_lun, /swap_if_big_endian
        If (keyword_set(mag)) Then Begin
           for it=0L,niters-1L do begin
              temp = values(*,it)
              printf, funit, iters(it), temp, format='(i10,9E15.7)'
           endfor
        EndIf Else Begin
           for it=0L,niters-1L do begin
              temp = values(*,it)
              printf, funit, iters(it), temp, format='(i10,7E15.7)'
           endfor
        EndElse
        close,funit
        free_lun,funit
    EndIf

end
