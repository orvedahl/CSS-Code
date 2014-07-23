pro read_shell_avg, file, quantities=quantities, q_names=q_names, radius=radius, data=shell_avg_values, oldfile=oldfile

    openr, file_unit, file, /get_lun, /swap_if_big_endian

    temp = dblarr(2)
    readu, file_unit, temp

    nr = fix(temp(0))
    nq = fix(temp(1))

    print, nr, nq
    If (keyword_set(oldfile)) Then Begin
        averages = assoc(file_unit,dblarr(nr,nq+1))
    Endif Else Begin
        averages = assoc(file_unit,dblarr(nr,nq+2))
    EndElse

    r1 = averages(2,0,0)
    r2 = averages(3,0,0)
    If (keyword_set(oldfile)) Then Begin
        radius = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    Endif Else Begin
        radius = reform(averages(*,1,0))
    EndElse
    quantities = fix(averages(4:nq+3,0,0))
    print, quantities

    If keyword_set(qoi) Then Begin
       iq = where(quantities eq qoi)+1
       If (iq lt 0) Then Begin
          Print, 'Quantity ', iq,' not available'
          return
       EndIf
    EndIf

    q_names = get_quantity_name(quantities)

    record = 0
    while not(eof(file_unit)) do begin
        If (keyword_set(oldfile)) Then Begin
            temp = reform(averages(0:nr-1,0:nq,record))
        Endif Else Begin
            temp = reform(averages(0:nr-1,0:nq+1,record))
        EndElse
        record = record+1
    endwhile
    
    if (record eq 0) then begin
        num_records = 1
    endif else begin
        num_records = record
    endelse

    if keyword_set(qoi) then begin
       shell_avg_values = dblarr(nr,num_records)
    endif else begin
       shell_avg_values = dblarr(nr,nq,num_records)
    endelse

    If (keyword_set(oldfile)) Then Begin
        for i=0,num_records-1 do begin
            if keyword_set(qoi) then begin
                shell_avg_values(*,i) = averages(*,iq+1,i)
            endif else begin
                shell_avg_values(*,*,i) = averages(*,1:nq,i)
            endelse
        endfor    
    Endif Else Begin
        for i=0,num_records-1 do begin
            if keyword_set(qoi) then begin
                shell_avg_values(*,i) = averages(*,iq+2,i)
            endif else begin
                shell_avg_values(*,*,i) = averages(*,2:nq+1,i)
            endelse
        endfor    
    EndElse
    close, file_unit
    free_lun, file_unit
end
