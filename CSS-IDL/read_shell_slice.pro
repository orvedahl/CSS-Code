pro read_shell_slice, file, qoi=qoi, quantities=quantities, q_names=q_names, radii=radii, theta=theta, phi=phi, data=shell_slice_values, nonuni=nonuni, recmax=recmax

    openr, file_unit, file, /get_lun, /swap_if_big_endian

    temp = dblarr(5)
    readu, file_unit, temp

    nradii = fix(temp(0))
    nr = fix(temp(1))
    nth = fix(temp(2))
    nphi = fix(temp(3))
    nq = fix(temp(4))

    print, nth, nphi, nradii, nq
    close, /all
    free_lun, file_unit
    
    openr, file_unit, file, /get_lun, /swap_if_big_endian
    averages = assoc(file_unit,dblarr(nth,nphi+1,nradii,nq))

    r1 = averages(5,0,0,0,0)
    r2 = averages(6,0,0,0,0)
    th1 = averages(7,0,0,0,0)
    th2 = averages(8,0,0,0,0)
    phi1 = averages(9,0,0,0,0)
    phi2 = averages(10,0,0,0,0)
 
    phi = (phi2-phi1)*dindgen(nphi)/(nphi-1.0d0)+phi1
    theta = (th2-th1)*dindgen(nth)/(nth-1.0d0)+th1
   
    quantities = fix(averages(11:nq+10,0,0,0,0))
    If (keyword_set(nonuni)) Then Begin
        Case (nonuni) of
            1: begin
                radius = (r2-r1)*(1d0-exp(-2d0*dindgen(nr)/double(nr-1L)))/(1d0-exp(-2d0))+r1
            end
            Default: begin
                radius = (r2-r1)*dindgen(nr)/(nr-1.0d0)+r1
            end
        EndCase
    Endif Else Begin
        radius = (r2-r1)*dindgen(nr)/(nr-1.0d0)+r1
    EndElse
    radii = radius(fix(averages(nq+11:nradii+nq+10,0,0,0,0))-1)

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
    If (keyword_set(recmax)) Then Begin
       If (recmax ge 1) Then record = recmax
    Endif Else Begin
       while not(eof(file_unit)) do begin
          temp = averages(0,1:nphi,0,0,record)
          record = record+1
       endwhile
    EndElse
    
    if (record eq 0) then begin
        num_records = 1
    endif else begin
        num_records = record
    endelse

    print, 'num_records = ', num_records

    if keyword_set(qoi) then begin
       shell_slice_values = dblarr(nth,nphi,nradii,num_records)
    endif else begin
       shell_slice_values = dblarr(nth,nphi,nradii,nq,num_records)
    endelse

    for i=0,num_records-1 do begin
        if keyword_set(qoi) then begin
           shell_slice_values(*,*,*,i) = averages(*,1:nphi,*,iq,i)
        endif else begin
           shell_slice_values(*,*,*,*,i) = averages(*,1:nphi,*,*,i)
        endelse
    endfor    
    close, file_unit
    free_lun, file_unit
end
