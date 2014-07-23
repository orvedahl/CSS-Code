pro load_data
common vccom, data, radius, theta, phi
    dir1 = '/freyr3/augustso/CSSDE/solar_test/Checkpoints/2260000/'

    fname = strtrim(dir1+'header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = lonarr(5)
    readf, file_unit, temp,format='(5i9)'

    nr   = fix(temp[0])
    nth  = fix(temp[1])
    nph  = fix(temp[2])
    nq   = fix(temp[3])
    iter = long(temp[4])

    print, temp

    temp = dblarr(8)
    readf, file_unit, temp, format='(8e14.7)'
    dt   = temp[0]
    time = temp[1]
    r1   = temp[2]
    r2   = temp[3]
    th1  = temp[4]
    th2  = temp[5]
    ph1  = temp[6]
    ph2  = temp[7]

    print, temp

    close, file_unit
    free_lun, file_unit

    radius = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
    phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1

    qnames = ['rho','u','v','w','s']

    fname = strtrim(dir1+'checkpoint.',2)

    data = dblarr(nth,nph,nr,nq)
    temp = dblarr(nth,nph,nr)
    ;read in files
    for iq=0,nq-1 do begin
        tempname = strtrim(fname+qnames[iq],2)
        print, 'reading ', tempname
        openr, file_unit, tempname, /get_lun
        temp[*,*,*] = 0d0
        readu, file_unit, temp
        data[*,*,*,iq] = temp
        close, file_unit
        free_lun, file_unit
    endfor

    
end

pro velocity_characterization
common vccom, data, radius, theta, phi

    rsun = 6.96d10
    nr = n_elements(radius)
    nth = n_elements(theta)
    nph = n_elements(phi)
    docolor
    choose_color, /ash_color, /quiet
    vrslice = congrid(reform(data[*,*,nr/2,1]),1024,1024,cubic=-0.5)
    tvscl, transpose(vrslice)

    vr_rms = dblarr(nr)
    vt_rms = dblarr(nr)
    vp_rms = dblarr(nr)
    vh_rms = dblarr(nr)
    for ir=0,nr-1 do begin
        ntot = double(nth)*double(nph)
        tmp = reform(data[7:nth-8,*,ir,1])
        vr_rms[ir] = sqrt(total(tmp^2)/ntot)
        tmp = reform(data[7:nth-8,*,ir,2])^2
        vt_rms[ir] = sqrt(total(tmp)/ntot)
        tmp = reform(data[7:nth-8,*,ir,3])^2
        vp_rms[ir] = sqrt(total(tmp)/ntot)
        tmp = reform(data[7:nth-8,*,ir,2])^2+reform(data[7:nth-8,*,ir,3])^2
        vh_rms[ir] = sqrt(total(tmp)/ntot)
    endfor

    vrp_rms = dblarr(nr)
    vhp_rms = dblarr(nr)
    for ir=0,nr-1 do begin
        ntot = double(nth)*double(nph)
        tmp = reform(data[7:nth-8,*,ir,1])
        vrp_rms[ir] = sqrt(total((tmp-vr_rms[ir])^2)/ntot)
        tmp = reform(data[7:nth-8,*,ir,2]-vt_rms[ir])^2+reform(data[7:nth-8,*,ir,3]-vp_rms[ir])^2
        vhp_rms[ir] = sqrt(total(tmp-vh_rms[ir])/ntot)
    endfor

    xsize = 15d0
    ysize = 15d0

    nth = n_elements(theta)
   
    filename = 'css_rms_velocities.eps'
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, /TT_FONT, SET_FONT='Helvetica'

    !P.MULTI=0
    !P.FONT=1
    !P.NOERASE=1
    !P.CHARSIZE=2
    x0 = 0.175d0
    x1 = 0.95d0
    y1 = 0.95d0
    y0 = 0.16d0
    !P.POSITION=[x0,y0,x1,y1]
    plot, radius/rsun, vhp_rms/1d2, charsize=2, xtitle='r/R'+sunsymbol(), ytitle=textoidl('v_{rms} (m s^{-1})'), $
      thick=3, xthick=3, ythick=3
    oplot, radius/rsun, vrp_rms/1d2, linestyle=1, thick=3

    ;y1 = y0-0.07d0
    ;y0 = y1-0.25d0
    ;!P.POSITION=[x0,y0,x1,y1]
    ;plot, radius/rsun, vh_rms/1d2, charsize=2, xtitle='r/R'+sunsymbol(), ytitle=textoidl('v_{rms} (m s^{-1})'), $
    ;  thick=3, xthick=3, ythick=3

    ;y1 = y0-0.07d0
    ;y0 = y1-0.25d0
    ;!P.POSITION=[x0,y0,x1,y1]
    ;plot, radius/rsun, vh_rms/vr_rms, charsize=2, xtitle='r/R'+sunsymbol(), ytitle=textoidl('v_{h}^{rms}/v_{r}^{rms}'), yrange=[0d0,3d0]

    DEVICE, /close
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

stop

end
