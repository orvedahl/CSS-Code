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

    openr, file_unit, fname, /get_lun

    temp = lonarr(5)
    readf, file_unit, temp,format='(6i9)'

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
    data = dblarr(nth,nph,nr,nq)
    ;read in files
    for iq=0,nq-1 do begin
        tempname = strtrim(fname+qnames(iq),2)
        openr, file_unit, tempname, /get_lun, /swap_if_big_endian
        a = assoc(file_unit,dblarr(nth,nph,nr))
        data[*,*,*,iq] = a[0]

        close, file_unit
        free_lun, file_unit
    endfor

    print, 'files loaded'
    print, 'size of data', size(data)

end

pro ASP_fourpanel, filename, quan, rad, sclmn=sclmn, sclmx=sclmx
    common epscheck, data, phi, theta, r, nr, nth, nph, nq

    rstar = 6.96d10

    y0 = cos(theta(0))
    y1 = cos(theta(nth-1))
    x0 = sin(theta(nth/2))*cos(phi(0))
    x1 = sin(theta(nth/2))*cos(phi(nph-1))

    space = 0.01d0    
    xsize = 7.5d0
    ysize = xsize/3.5d0

    csize = 0.7

    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Times-Roman', /CMYK, /INCHES

    !P.MULTI=0
    !P.FONT=0
    !P.NOERASE=1
    !P.THICK=2.5
    !P.CHARSIZE=csize
    docolor
    choose_color, /ash_color, /quiet

    css_lon = phi*180d0/!dpi
    css_lat = 90d0-theta*!radeg
    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]
    ortho = [0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))]
    highlightbox = limits

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[0]*min(data(*,*,rad[0],quan))
    Endif Else Begin
       mnv = 0.5*min(data(*,*,rad[0],quan))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[0]*max(data(*,*,rad[0],quan))
    Endif Else Begin
       mxv = 0.5*max(data(*,*,rad[0],quan))
    EndElse
    draw = transpose(reform(data(*,*,rad[0],quan)))

 ;Figure out bounds from size of domain
    dx = 0.225
    x0 = space
    x1 = x0+dx

    y0 = space
    y1 = 1d0-space
    !P.POSITION=[x0,y0,x1,y1]

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    xyouts, /NORMAL, x0+0.12d0*dx, space, textoidl('-60\circ'),color=0, charsize=0.5
    xyouts, /NORMAL, x0+0.12d0*dx, 1d0-4d0*space, textoidl('+60\circ'),color=0, charsize=0.5
    ;xyouts, /NORMAL, x0-space, 0.5-8d0*space, textoidl('0\circ'),color=0, charsize=0.5
    ;xyouts, /NORMAL, x1, 0.5-8d0*space, textoidl('90\circ'),color=0, charsize=0.5
    xyouts, /NORMAL, x0, 0.95, '(a)',color=0

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[1]*min(data(*,*,rad[1],quan))
    Endif Else Begin
       mnv = 0.5*min(data(*,*,rad[1],quan))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[1]*max(data(*,*,rad[1],quan))
    Endif Else Begin
       mxv = 0.5*max(data(*,*,rad[1],quan))
    EndElse
    draw = transpose(reform(data(*,*,rad[1],quan)))

    x0 = x1 + space
    x1 = x0 + dx
    !P.POSITION=[x0,y0,x1,y1]

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    xyouts, /NORMAL, x0, 0.95, '(b)',color=0

    x0 = x1 + space
    x1 = x0 + dx
    x0p = x0
    x1p = x1

    x0 = x1 + 0.5d0*space
    x1 = 1d0 - 0.5d0*space
    !P.POSITION=[x0,y0,x1,y1]

    read_jpeg, 'AZ_Omega_Labeled_Howe_2009_shrink.jpg', image, ctable, COLORS=256, DITHER=1, /TWO_PASS
    tvlct, ctable

    dx = x1-x0
    dy = dx*xsize/ysize
    yy0 = 0.5d0-dy/2d0
    tv, image, x0, yy0, xsize=dx, ysize=dy, /normal

    choose_color, /ash_color, /quiet

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[2]*min(data(*,*,rad[2],quan))
    Endif Else Begin
       mnv = 0.5*min(data(*,*,rad[2],quan))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[2]*max(data(*,*,rad[2],quan))
    Endif Else Begin
       mxv = 0.5*max(data(*,*,rad[2],quan))
    EndElse
    draw = transpose(reform(data(*,*,rad[2],quan)))

    x0 = x0p
    x1 = x1p
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    xyouts, /NORMAL, x0, 0.95, '(c)',color=0
    xyouts, /NORMAL, x1-1.5d0*space, 0.95, '(d)',color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'evince '+filename+' &'
    stop
end
