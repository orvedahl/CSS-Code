pro load_data, dir, irng, istep
common blk, data, radius, phi, times

    nfs = long((irng[1]-irng[0])/istep)
    iters = lindgen(nfs)*istep+irng[0]

    ;Get file names
    root_dir = dir+'/Checkpoints/'

    pretemp = ' '
    stemp = string(irng[0], format='(i7.7)')
    iterstr = strtrim(pretemp+stemp,2)

    read_dir = root_dir+iterstr+'/'

    fname = strtrim(read_dir+'header',2)
    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = lonarr(4)
    readf, file_unit, temp,format='(5i9)'

    nr   = fix(temp(0))
    nph  = fix(temp(1))
    nth  = fix(temp(2))
    nq   = fix(temp(3))
    iter = long(temp(4))

    print, temp

    temp = dblarr(8)
    readf, file_unit, temp, format='(8e14.7)'
    dt   = temp(0)
    t    = temp(1)
    r1   = temp(2)
    r2   = temp(3) 
    th1  = temp(4)
    th2  = temp(5)
    ph1  = temp(6)
    ph2  = temp(7)

    print, temp

    close, file_unit
    free_lun, file_unit

    qnames = ['rho','w','v','u','s','bt','bp','br']
    files = read_dir+'checkpoint.'+qnames
    data = dblarr(nth,nph,nr,nq,nfs)
    for iv=0,nq-1 do begin
        openr, file_unit, files(iv), /get_lun, /swap_if_big_endian
        a = assoc(file_unit,dblarr(nth,nph,nr))
        data[*,*,*,iv,0] = a[0]
        close, file_unit
        free_lun, file_unit
    endfor

    radius = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1
    times = dblarr(nfs)
    times[0] = t
    for it=1L,nfs-1L do begin
       pretemp = ' '
       stemp = string(iters[it], format='(i7.7)')
       iterstr = strtrim(pretemp+stemp,2)
       
       read_dir = root_dir+iterstr+'/'
       
       fname = strtrim(read_dir+'header',2)
       openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format
       
       temp = lonarr(4)
       readf, file_unit, temp,format='(4i9)'
       
       nr   = fix(temp(0))
       nph  = fix(temp(1))
       nq   = fix(temp(2))
       iter = long(temp(3))
       
       temp = dblarr(7)
       readf, file_unit, temp, format='(7e14.7)'
       dt   = temp(0)
       t    = temp(1)
       r1   = temp(2)
       r2   = temp(3) 
       th1  = temp(4)
       ph1  = temp(5)
       ph2  = temp(6)
       
       close, file_unit
       free_lun, file_unit
       
       files = read_dir+'checkpoint.'+qnames
       for iv=0,nq-1 do begin
          openr, file_unit, files(iv), /get_lun, /swap_if_big_endian
          a = assoc(file_unit,dblarr(nph,nr))
          data[*,*,iv,it] = a[0]
          close, file_unit
          free_lun, file_unit
       endfor
       
       times[it] = t
    endfor

end

pro checkpoint_movie, compute_vort=compute_vort
common blk, data, radius, phi, times

   nfs = n_elements(times)
   nr = n_elements(radius)
   nphi = n_elements(phi)
   r1 = min(radius)
   r2 = max(radius)
   If (keyword_set(compute_vort)) Then Begin
      vort = dblarr(nphi,nr,nfs)
      for it=0,nfs-1 do begin
         vp = reform(data[*,*,2,it])
         vr = reform(data[*,*,3,it])
         dpi = 1d0/(max(phi)-min(phi))
         dri = 1d0/(r2-r1)
         for ir=0,nr-1 do begin
            tmp = reform(vr[*,ir])
            tmp2 = compact_fd6(dpi,tmp,0d0,0d0,1,dtype=3)/radius[ir]
            vort[*,ir,it] = tmp2
         endfor
         for ip=0,nphi-1 do begin
            tmp = reform(vp[ip,*])
            b1 = vp[ip,0]/r1
            b2 = vp[ip,nr-1]/r2
            tmp2 = compact_fd6(dri,tmp,b1,b2,4,dtype=0)
            vort[ip,*,it] = vort[ip,*,it]-tmp2-tmp/radius
         endfor
         print, 'on it=', it
      endfor
   EndIf
   stop

end
