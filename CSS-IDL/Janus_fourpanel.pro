pro load_fields, numrecs=numrecs, filename=filename
    common epscheck, data1, data2, data3, phi1, theta1, r1, quan1, phi2, theta2, r2, quan2

 ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr1/augustso/CSSDE/'
    file = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select Shell Slice File to Read')
       
 ;Determine the number of files.
    n_files = n_elements(file)
       
    print, 'loading file:',file
 ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(file eq '') then return
    endif
       
 ;Load sld case
    read_shell_slice, file, phi=phi1, theta = theta1, radii = r1, data = data1, QUANTITIES = quan1
    print, 'file loaded'
    print, 'size of data', size(data)

 ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr1/augustso/CSSDE/'
    file = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select Shell Slice File to Read')
       
 ;Determine the number of files.
    n_files = n_elements(file)
       
    print, 'loading file:',file
 ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(file eq '') then return
    endif
       
 ;Load sld case
    read_shell_slice, file, phi=phi2, theta = theta2, radii = r2, data = data2, QUANTITIES = quan2
    print, 'file loaded'
    print, 'size of data', size(data2)
end

pro Janus_fourpanel, filename, rad, quan, sclmn=sclmn, sclmx=sclmx
    common epscheck, data1, data2, data3, phi1, theta1, r1, quan1, phi2, theta2, r2, quan2

    nr = n_elements(r1)
    nth = n_elements(theta1)
    nphi = n_elements(phi1)
    nq = n_elements(quan1)

    rstar = 6.96d10

    space = 0.01d0    
    xsize = 7.5d0
    ysize = xsize/3.75d0

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

    css_lon = phi1*180d0/!dpi
    css_lat = 90d0-theta1*!radeg
    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]
    ortho = [0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))]
    highlightbox = limits

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[0]*min(data1(*,*,rad[0],quan[0],0))
    Endif Else Begin
       mnv = 0.5*min(data1(*,*,rad[0],quan[0],0))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[0]*max(data1(*,*,rad[0],quan[0],0))
    Endif Else Begin
       mxv = 0.5*max(data1(*,*,rad[0],quan[0],0))
    EndElse
    draw = transpose(reform(data1(*,*,rad[0],quan[0],0)))

 ;Figure out bounds from size of domain
    dx = 0.9d0/4d0
    x0 = 2d0*space
    x1 = x0+dx

    y0 = space
    y1 = 1d0-space
    !P.POSITION=[x0,y0,x1,y1]

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    plots, [x0,x1], [y1,y1]-7d0*space, thick=2, /normal
    xyouts, /NORMAL, x0, 2d0*space, textoidl('-2.5\circ'),color=0, charsize=0.5
    xyouts, /NORMAL, x1-2.5d0*space, 2d0*space, textoidl('+2.5\circ'),color=0, charsize=0.5
    xyouts, /NORMAL, x0, 0.95, '(a)',color=0

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[1]*min(data1(*,*,rad[1],quan[1],0))
    Endif Else Begin
       mnv = 0.5*min(data1(*,*,rad[1],quan[1],0))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[1]*max(data1(*,*,rad[1],quan[1],0))
    Endif Else Begin
       mxv = 0.5*max(data1(*,*,rad[1],quan[1],0))
    EndElse
    draw = transpose(reform(data1(*,*,rad[1],quan[1],0)))

    x0 = x1 + 2d0*space
    x1 = x0 + dx
    !P.POSITION=[x0,y0,x1,y1]

    choose_color, /ryb, /quiet
    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    choose_color, /ash_color, /quiet
    plots, [x0,x1], [y1,y1]-7d0*space, thick=2, /normal
    xyouts, /NORMAL, x0, 0.95, '(b)',color=0

    css_lon = phi2*180d0/!dpi
    css_lat = 90d0-theta2*!radeg
    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]
    ortho = [0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))]
    highlightbox = limits

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[2]*min(data2(*,*,rad[2],quan[2],0))
    Endif Else Begin
       mnv = 0.5*min(data2(*,*,rad[2],quan[2],0))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[2]*max(data2(*,*,rad[2],quan[2],0))
    Endif Else Begin
       mxv = 0.5*max(data2(*,*,rad[2],quan[2],0))
    EndElse
    draw = transpose(reform(data2(*,*,rad[2],quan[2],0)))

    x0 = x1 + 2d0*space
    x1 = x0 + dx
    !P.POSITION=[x0,y0,x1,y1]

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    choose_color, /ash_color, /quiet
    plots, [x0,x1], [y1,y1]-7d0*space, thick=2, /normal
    xyouts, /NORMAL, x0, 2d0*space, textoidl('-2.5\circ'),color=0, charsize=0.5
    xyouts, /NORMAL, x1-2.5d0*space, 2d0*space, textoidl('+2.5\circ'),color=0, charsize=0.5
    xyouts, /NORMAL, x0, 0.95, '(c)',color=0

    If (keyword_set(sclmn)) Then Begin
       mnv = sclmn[3]*min(data2(*,*,rad[3],quan[3]))
    Endif Else Begin
       mnv = 0.5*min(data2(*,*,rad[3],quan[3]))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx[3]*max(data2(*,*,rad[3],quan[3]))
    Endif Else Begin
       mxv = 0.5*max(data2(*,*,rad[3],quan[3]))
    EndElse
    draw = transpose(reform(data2(*,*,rad[3],quan[3])))

    x0 = x1 + 2d0*space
    x1 = x0 + dx
    !P.POSITION=[x0,y0,x1,y1]
    choose_color, /ryb, /quiet
    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=ortho, latlon_window=limits, highlightbox=highlightbox

    choose_color, /ash_color, /quiet
    plots, [x0,x1], [y1,y1]-7d0*space, thick=2, /normal
    xyouts, /NORMAL, x0, 0.95, '(d)',color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'evince '+filename+' &'
    stop
end
