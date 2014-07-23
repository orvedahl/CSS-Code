@get_flux_names.pro
pro read_flux_bal, file, quantities=quantities, q_names=q_names, radius=radius, data=flux_bal_values, qoi=qoi, times=times

    openr, file_unit, file, /get_lun, /swap_if_big_endian

    temp = fltarr(2)
    readu, file_unit, temp

    nr = fix(temp(0))
    nq = fix(temp(1))

    averages = assoc(file_unit,fltarr(nr,nq+1))

    r1 = averages(2,0,0)
    r2 = averages(3,0,0)
    radius = (r2-r1)*dindgen(nr)/(nr-1.0d0)+r1

    quantities = indgen(nq)+1

    If keyword_set(qoi) Then Begin
       iq = where(quantities eq qoi)+1
       If (iq eq 0) Then Begin
          Print, 'Quantity ', iq,' not available'
          return
       EndIf
    EndIf

    q_names = get_flux_names(quantities)

    record = 0
    while not(eof(file_unit)) do begin
        temp = reform(averages(0:nr-1,0:nq,record))
        record = record+1
    endwhile
    
    if (record eq 0) then begin
        num_records = 1
    endif else begin
        num_records = record
    endelse

    if keyword_set(qoi) then begin
       flux_bal_values = dblarr(nr,num_records)
    endif else begin
       flux_bal_values = dblarr(nr,nq,num_records)
    endelse

    times = dblarr(num_records)

    for i=0,num_records-1 do begin
        if keyword_set(qoi) then begin
           flux_bal_values(*,i) = averages(*,iq+1,i)
           times(i) = averages(5,0,i)
        endif else begin
           flux_bal_values(*,*,i) = averages(*,1:nq,i)
           times(i) = averages(5,0,i)
        endelse
    endfor    
    close, file_unit
    free_lun, file_unit
end
