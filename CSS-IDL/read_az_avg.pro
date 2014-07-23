pro read_az_avg, file, qoi=qoi, quantities=quantities, q_names=q_names, radius=radius, theta=theta, data=az_avg_values, num_recs=num_recs
    
    ;It is assumed that the filename is a single string with no directory information
    ;if byte(1,0) then begin
        openr, file_unit, file, /swap_if_big_endian, /get_lun ;open the file in read-only mode for little endian format
    ;endif else begin
    ;     openr, file_unit, file, /get_lun ;open the file in read-only mode for big endian format
    ;endelse

    temp = dblarr(3)
    readu, file_unit, temp

    nth = fix(temp(0))
    nr = fix(temp(1))
    nq = fix(temp(2))

    print, temp
    close, file_unit
    free_lun, file_unit

    openr, file_unit, file, /get_lun, /swap_if_big_endian
    averages = assoc(file_unit,dblarr(nth,nr+2,nq))
    r1 = averages[3,0,0,0]
    r2 = averages[4,0,0,0]
    th1 = averages[5,0,0,0]
    th2 = averages[6,0,0,0]
 
    theta = (th2-th1)*dindgen(nth)/(nth-1.0d0)+th1
    quantities = fix(averages[7:nq+6,0,0,0])
    radius = averages[0:nr-1,1,0]

    If keyword_set(qoi) Then Begin
       iq = where(quantities eq qoi)+1
       If (iq lt 0) Then Begin
          Print, 'Quantity ', iq,' not available'
          return
       EndIf
    EndIf

    print, quantities
    q_names = get_quantity_name(quantities)
    print, q_names

    record = 0
    while not(eof(file_unit)) do begin
        temp = averages[0,1,0,record]
        record = record+1
    endwhile

    if (record eq 0) then begin
        num_records = 1
    endif else begin
        num_records = record
    endelse

    num_recs = num_records
    print, 'num_records = ', num_records

    if keyword_set(qoi) then begin
       az_avg_values = dblarr(nth,nr,num_records)
    endif else begin
       az_avg_values = dblarr(nth,nr,nq,num_records)
    endelse

    for i=0,num_records-1 do begin
        if keyword_set(qoi) then begin
           az_avg_values[*,*,i] = averages[*,2:nr+1,iq,i]
        endif else begin
           az_avg_values[*,*,*,i] = averages[*,2:nr+1,*,i]
        endelse
    endfor    
    close, file_unit
    free_lun, file_unit
end
