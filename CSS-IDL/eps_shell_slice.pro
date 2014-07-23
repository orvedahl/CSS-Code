pro load_fields, numrecs=numrecs, filename=filename
    common epshell, data, phi, theta, r, num_records, N_R, N_Theta, N_Phi, N_Q

    If (keyword_set(filename)) Then Begin
       restore, filename
       N_R = n_elements(rn)
       N_Theta = n_elements(theta)
       N_Phi = n_elements(phi)
       N_Q = 4
       num_records = 1
       r = rn
       data = dblarr(N_Theta,N_Phi,N_R,4,1)
       data[*,*,*,0:2] = curl
       data[*,*,*,3] = ens
       curl = 0
       ens = 0
    Endif Else Begin
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
       read_shell_slice, file, phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities
       print, 'file loaded'
       print, 'size of data', size(data)
       num_records = n_elements(data(0,0,0,0,*))
       N_R = n_elements(r)
       N_Theta = n_elements(theta)
       N_Phi = n_elements(phi)
       N_Q = n_elements(quantities)
    EndElse
end

pro eps_shell_slice, filename, quan, rad, rec, rstar=rstar, sclmn=sclmn, sclmx=sclmx, aspect=aspect
    common epshell, data, phi, theta, r, num_records, N_R, N_Theta, N_Phi, N_Q

    If (not keyword_set(rstar)) Then rstar = 6.96d10

    y0 = cos(theta(0))
    y1 = cos(theta(N_Theta-1))
    x0 = sin(theta(N_Theta-1))*cos(phi(0))
    x1 = sin(theta(N_Theta-1))*cos(phi(N_Phi-1))
    If (not keyword_set(aspect)) Then Begin
        aspect = (y1-y0)/(x1-x0) ;y/x
    Endif
        
    print, 'aspect =', aspect
    space=0.01d0    
    ;If (aspect gt 1d0) Then Begin
    ;    aspect = 0.7*aspect
    ;    ysize = 4d0*2.54d0
    ;    xsize = ysize/aspect
    ;Endif Else Begin
    ;    aspect = 1.2d0*aspect
    ;    xsize = 4d0*2.54d0
    ;    ysize = xsize*aspect
    ;Endelse
    xsize = 5*2.54d0
    ysize = aspect*xsize
    csize = 1.5
    cs2 = 2

    !P.MULTI=0
    !P.FONT=1
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    choose_color, /ash_color, /quiet
    css_lon = phi*180d0/!dpi
    css_lat = reverse(90d0-theta*!radeg)

    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]
    
    If (keyword_set(sclmn)) Then Begin
        mnv = sclmn*min(data(*,*,rad,quan,rec))
    Endif Else Begin
        mnv = 0.5*min(data(*,*,rad,quan,rec))
    Endelse

    If (keyword_set(sclmx)) Then Begin
       mxv = sclmx*max(data(*,*,rad,quan,rec))
    Endif Else Begin
       mxv = 0.5*max(data(*,*,rad,quan,rec))
    Endelse
    print, mnv, mxv
    draw = transpose(reform(data(*,*,rad,quan,rec)))
    
 ;Figure out bounds from size of domain
    x0 = 0.06d0
    y0 = 0.05d0
    x1 = 0.99d0
    y1 = 0.95d0
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox
    x0 = 0.02d0
    y0 = 0.05d0
    x1 = 0.99d0
    y1 = 0.95d0

    !P.POSITION=[x0,y0,x1,y1]
    sunsym = "!D!9!Z(6E)!X"
    radstr = strtrim(string(r(rad)/rstar,format='(f4.2)'))+'!U !NR'
    ;print, 'radstr=',radstr
    xyouts, /NORMAL, 0.81, 0.95, radstr, charsize=cs2,charthick=1,color=0,width=wid
    xyouts, /NORMAL, 0.81+wid, 0.95, sunsym, charsize=csize,charthick=1,color=0,font=-1

    xyouts,/NORMAL,0.07,0.04,textoidl('0\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,0.92,0.04,textoidl('+60\circ'),charsize=csize,charthick=1.0,color=0

    xyouts,/NORMAL,0.01,0.07,textoidl('-20\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,0.01,0.905,textoidl('+20\circ'),charsize=csize,charthick=1.0,color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'evince '+filename+' &'
    stop
end
