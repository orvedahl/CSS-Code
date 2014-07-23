@write_bulk_cartesian.pro
pro make_cartesian_chunk, iter, q, css_case, triangle_file=triangle_file, resize=resize, trick=trick, float=float

    dir='/freyr1/augustso/CSSDE/'+css_case+'/Checkpoints/'+iter+'/'
    fname = strtrim(dir+'header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = lonarr(5)
    readf, file_unit, temp,format='(5i9)'

    nr  = fix(temp(0))
    nth  = fix(temp(1))
    nph   = fix(temp(2))
    nq   = fix(temp(3))
    it = fix(temp(4))

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

    fname = strtrim(dir+'checkpoint.'+q,2)
 
    ;Set up grid (assumes regular r, theta, phi grid (so pre-interpolate!))
    radius = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
    phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1

    openr, file_unit, fname, /get_lun, /swap_if_big_endian
    
    if keyword_set(float) then begin
       temp = assoc(file_unit, fltarr(nth,nph,nr))
    endif else begin
       temp = assoc(file_unit, dblarr(nth,nph,nr))
    endelse

    bulk_data = reform(temp(0))

    close, file_unit
    free_lun, file_unit

    bounds = [nr,nth,nph,r1,r2,th1,th2,ph1,ph2]

    print, 'bounds ', bounds

    filename=q+'.cart'

    write_bulk_cartesian, bulk_data, bounds, filename, triangle_file=triangle_file, dir=dir, resize=resize

end
