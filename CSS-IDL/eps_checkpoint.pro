@display_css_shell_slice.pro
pro load_fields, iter, css_case, dir=dir
    common epscheck, data, phi, theta, r, nr, nth, nph, nq

    If (keyword_set(dir)) Then Begin
        dir1 = strtrim(dir,2)
    Endif Else Begin
        dir1 = '/freyr1/augustso/CSSDE/'+css_case+'/Checkpoints/'

        stemp = string(iter, format='(i7.7)')
        dir1 = strtrim(dir1,2)+strtrim(stemp+'/',2)
    EndElse

    fname = strtrim(dir1+'header',2)

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

    phi=(ph2-ph1)*dindgen(nph)/(nph-1d0) ; +ph1
    theta=(th2-th1)*dindgen(nth)/(nth-1d0)+th1
    r=(r2-r1)*dindgen(nr)/(nr-1d0)+r1

    fname = strtrim(dir1+'checkpoint.',2)
    qnames = ['rho','w','v','u','s','br','bt','bp']
    data = dblarr(nr,nth,nph,nq)
    temp = dblarr(nth,nph)
    ;read in files
    for iq=0,nq-1 do begin
        tempname = strtrim(fname+qnames(iq),2)
        openr, file_unit, tempname, /get_lun, /swap_if_big_endian
        for ir=0L,nr-1L do begin
            temp(*,*) = 0d0
            readu, file_unit, temp
            data(ir,*,*,iq) = temp
        endfor
        close, file_unit
        free_lun, file_unit
    endfor

    print, 'files loaded'
    print, 'size of data', size(data)

end

pro eps_checkpoint, filename, quan, rad, rstar=rstar, sclmn=sclmn, sclmx=sclmx, ortho=ortho, latchop=latchop
    common epscheck, data, phi, theta, r, nr, nth, nph, nq

    If (not keyword_set(rstar)) Then rstar = 6.96d10

    y0 = cos(theta(0))
    y1 = cos(theta(nth-1))
    x0 = sin(theta(nth/2))*cos(phi(0))
    x1 = sin(theta(nth/2))*cos(phi(nph-1))
    ;aspect = (y1-y0)/(x1-x0);y/x
    aspect = 1d0
    print, 'aspect =', aspect
    space=0.01d0    
    If (aspect gt 1d0) Then Begin
        aspect = 0.5d0*aspect
        ysize = 4d0*2.54d0
        xsize = ysize/aspect
    Endif Else Begin
        aspect = 1d0 ;aspect
        xsize = 4d0*2.54d0
        ysize = xsize*aspect
    EndElse
    csize = 1.0

    Set_Plot, "PS"
    !P.MULTI=0
    !P.FONT=1
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica', /TT_FONT
    choose_color, /ash_color, /quiet

    css_lon = phi*180d0/!dpi
    If (keyword_set(latchop)) Then Begin
        ;css_lat = reverse(90d0-theta[0:latchop]*!radeg)
        css_lat = 90d0-theta[0:latchop]*!radeg
  
        If (keyword_set(sclmn)) Then Begin
            mnv = sclmn*min(data(rad,0:latchop,*,quan))
        Endif Else Begin
            mnv = 0.5*min(data(rad,0:latchop,*,quan))
        Endelse

        If (keyword_set(sclmx)) Then Begin
            mxv = sclmx*max(data(rad,0:latchop,*,quan))
        Endif Else Begin
            mxv = 0.5*max(data(rad,0:latchop,*,quan))
        EndElse
        draw = transpose(reform(data(rad,0:latchop,*,quan)))
    Endif Else Begin
        ;css_lat = reverse(90d0-theta*!radeg)
        css_lat = 90d0-theta*!radeg
        If (keyword_set(sclmn)) Then Begin
            mnv = sclmn*min(data(rad,*,*,quan))
        Endif Else Begin
            mnv = 0.5*min(data(rad,*,*,quan))
        Endelse

        If (keyword_set(sclmx)) Then Begin
            mxv = sclmx*max(data(rad,*,*,quan))
        Endif Else Begin
            mxv = 0.5*max(data(rad,*,*,quan))
        EndElse
        draw = transpose(reform(data(rad,*,*,quan)))
    EndElse

    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]

 ;Figure out bounds from size of domain
    x0 = 0.05d0
    y0 = 0.05d0
    x1 = 0.95d0
    y1 = 0.95d0
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits
    If (not keyword_set(ortho)) Then Begin
        ortho = [0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))]
    EndIf
    display_css_shell_slice, css_lat, css_lon, draw, /publication, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox
    x0 = 0.01d0
    y0 = 0.01d0
    x1 = 0.99d0
    y1 = 0.99d0
    !P.POSITION=[x0,y0,x1,y1]
    radstr = strtrim(string(r(rad)/rstar,format='(f4.2)')+textoidl(' R_{*}'))
    print, 'radstr=',radstr
    xyouts, /NORMAL, 0.44, 0.965, radstr,charsize=csize,charthick=1,color=0
    xyouts, /NORMAL, x0-0.5d0*space, 0.5-space, 'Eq.',charsize=csize,charthick=1,color=0

    ;xyouts, /NORMAL, 0.25, 0.02, textoidl('-60\circ'),charsize=csize,charthick=1,color=0
    ;xyouts, /NORMAL, 0.25, 0.96, textoidl('+60\circ'),charsize=csize,charthick=1,color=0
    xyouts, /NORMAL, x0+4d0*space, 0.03, textoidl('-30\circ'),charsize=csize,charthick=1.0,color=0,wid=wid
    xyouts, /NORMAL, x1-wid-4d0*space, 0.03, textoidl('+30\circ'),charsize=csize,charthick=1.0,color=0
    ;xyouts, /NORMAL, x0+space, 0.96, '(a)',charsize=csize,charthick=1.0,color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'evince '+filename+' &'
    stop
end
