pro srt2010_fig1ab, filename, ir=ir, iq=iq

    root_directory = '/freyr3/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select CSS Shell Slice File to Read')
 
    n_files = n_elements(files)

    if(n_files eq 1) then begin
       if(files eq '') then return
    endif

    read_shell_slice, files(0), phi=phi, theta = theta, radii = r, data = temp, QUANTITIES = quantities
    numrec = n_elements(temp(0,0,0,0,*))
    nr = n_elements(r)
    nth = n_elements(theta)
    nph = n_elements(phi)
    nq = n_elements(quantities)
  
    quantity_names = strarr(nq)
    for i=0,nq-1 do begin
        quantity_names(i) = get_quantity_name(quantities(i))
    endfor

    If (not keyword_set(iq)) Then Begin
        iq=3
    EndIf

    If ((not keyword_set(ir)) and (ir ne 0)) Then Begin 
        ir=nr-1
    EndIf

    ur = dblarr(nth,nph)
    print, ir, iq
    ur(*,*) = temp(*,*,ir,iq,numrec-1)

    lon = phi*180d0/!dpi
    lat = reverse(90d0-theta*!radeg)
    
    limits=[min(lat),min(lon),max(lat),max(lon)]
    
    If (iq eq 3) Then Begin
        mnv = min(ur(*,*))
        mxv = max(ur(*,*))
    Endif Else begin
        ;ur = ur-mean(ur)
        mnv = min(ur(*,*))
        mxv = max(ur(*,*))
    EndElse
    ;inds = where(ur lt 0d0)
    ;mneg = mean(ur(inds))
    ;inds = where(ur lt mneg,complement=cinds)
    ;ur(inds) = 0d0
    ;ur(cinds) = 1d0
    draw = transpose(ur)

    xsize=10d0
    ysize=10d0
    choose_color, /ash_color, /quiet
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, /COLOR, BITS_PER_PIXEL=8, XSIZE=xsize, YSIZE=ysize

    x0 = 0.05
    y0 = 0.05
    x1 = 0.95
    y1 = 0.95
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits
    display_shell_slice, lat, lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, orthographic=[0.5d0*(min(lon)+max(lon)),0.5d0*(min(lat)+max(lat))], latlon_window=limits, highlightbox=highlightbox, /nolines
    ;plots, [0,x1],[y1-0.0375,y1-0.0375],/normal, thick=4, color=255
    ;plots, [x0+0.02,x1-0.02],[y1-0.037,y1-0.037],/normal, thick=3, color=0
    ;plots, [0,x1],[y1-0.037,y1-0.037],/normal, thick=3, color=255
    ;plots, [x0+0.005,x1-0.005],[y1-0.037,y1-0.037],/normal, thick=2, color=0

    ;plots, [0,x1],[y0+0.0355,y0+0.0355],/normal, thick=4, color=255
    ;plots, [x0+0.02,x1-0.02],[y0+0.036,y0+0.036],/normal, thick=3, color=0
    ;plots, [0,x1],[y0+0.0355,y0+0.0355],/normal, thick=4, color=255
    ;plots, [x0+0.005,x1-0.005],[y0+0.036,y0+0.036],/normal, thick=3, color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'
end
