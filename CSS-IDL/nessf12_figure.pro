pro load_fields
    common epshell, data, phi, theta, r, num_records, N_R, N_Theta, N_Phi, N_Q

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
end

pro nessf12_shell_slice, filename, quan, rad, rec
    common epshell, data, phi, theta, r, num_records, N_R, N_Theta, N_Phi, N_Q

    If (not keyword_set(rstar)) Then rstar = 6.96d10

    space=0.01d0    
    xsize = 7.5d0*2.54d0
    ysize = xsize*(25d0/20d0)/3d0
    csize = 1.0
    csize2 = 1.25*csize
    Set_Plot, "PS"
    !P.MULTI=0
    !P.FONT=1
    !P.CHARSIZE=1.0
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT

    nth = n_elements(theta)
    nph = n_elements(phi)
    ;ntnew = 2L*nth
    ;npnew = 2L*nph

    ;dp = max(phi)-min(phi)
    ;phi2 = 2d0*dp*dindgen(npnew)/double(npnew-1)+min(phi)

    ;dt = max(theta)-min(theta)
    ;theta2 = 2d0*dt*dindgen(ntnew)/double(ntnew-1)+min(theta)

    css_lon = phi*180d0/!dpi
    css_lat = reverse(90d0-theta*!radeg)

    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]
    
    slice = reform(data[*,*,rad,quan[0],rec])
    ;datanew = dblarr(ntnew,npnew)
    ;datanew[0:nth-1,0:nph-1]=slice
    ;datanew[0:nth-1,nph:npnew-1]=slice
    ;datanew[nth:ntnew-1,0:nph-1]=reverse(slice,1)
    ;datanew[nth:ntnew-1,nph:npnew-1]=reverse(slice,1)
    ;mnv = 0.2*min(slice)
    mxv = 7d4
    mnv = -mxv

    draw = transpose(reform(slice))

 ;Figure out bounds from size of domain
    rsz = xsize/3d0*0.85d0
    dx = rsz/xsize
    x0 = 0.13d0/2.5d0
    x1 = x0+dx
    y0 = 0.05d0
    y1 = 0.95d0
    !P.NOERASE=1
    !P.POSITION=[x0,y0,x1,y1]
    print, !P.POSITION
    highlightbox = limits
    choose_color, /ash_color, /quiet
    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox

    ;colourbar,[x0+0.06,0.02,x0+0.14,0.05],mx=mxv,mn=mnv,textmn='-1',textmx='1',text_char_size=csize,$
    ;           plttype='ps',textu=textoidl('kms^{-1}'),text_thick=1,/pmval ;mpstxt

    xyouts,/NORMAL,x0+0.10,0.025,textoidl('\pm')+'!U !N700!U !Nm!E !Ns!U-1!N',charsize=csize,charthick=1.0,color=0
    ;radstr = strtrim(string(r[rad]/rstar,format='(f4.2)'),2)+'!U !NR'
    radstr = textoidl('v_{r}')
    ;print, 'radstr=',radstr
    xyouts, /NORMAL, (x1+x0)/2d0, 0.95, radstr, charsize=csize2,charthick=1,color=255,width=wid
    xyouts, /NORMAL, (x1+x0)/2d0-wid/2, 0.95, radstr, charsize=csize2,charthick=1,color=0
    xyouts, /NORMAL, x0-0.02, 0.93, 'a', charsize=csize2,charthick=1,color=0
    ;xyouts, /normal, 0.445/3d0+wid, 0.945, "!9!Z(6E)!X", font=-1, charsize=0.75

    xyouts,/NORMAL,0.05d0,0.03,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0,width=wid
    xyouts,/NORMAL,x1-0.01/3d0-wid,0.03,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    xyouts,/NORMAL,0.015,0.15,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,0.02,0.49,textoidl('Eq.'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,0.0075,0.83,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    slice = reform(data[*,*,rad,quan[1],rec])
    ;datanew = dblarr(ntnew,npnew)
    ;datanew[0:nth-1,0:nph-1]=slice
    ;datanew[0:nth-1,nph:npnew-1]=slice
    ;datanew[nth:ntnew-1,0:nph-1]=reverse(slice,1)
    ;datanew[nth:ntnew-1,nph:npnew-1]=reverse(slice,1)
    mnv = 0.2d0*min(slice)
    mxv = -mnv ;0.5d0*max(slice)

    draw = transpose(reform(slice))

    x0 = x1+0.04
    x1 = x0+dx
    !P.POSITION=[x0,y0,x1,y1]
    choose_color, /mag_color, /quiet
    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox

    choose_color, /ash_color, /quiet
    colourbar2,[x0-0.03,0.38,x0-0.0175,0.62],mx=mxv,mn=mnv,plttype='ps'
    xyouts,/NORMAL,x0+0.12,0.025,textoidl('\pm')+'!U !N1!U !NkG',charsize=csize,charthick=1.0,color=0
    ;radstr = strtrim(string(r[rad]/rstar,format='(f4.2)'),2)+'!U !NR'
    radstr = textoidl('B_{r}')
    ;print, 'radstr=',radstr
    xyouts, /NORMAL, (x1+x0)/2d0, 0.95, radstr, charsize=csize2,charthick=1,color=255,width=wid
    xyouts, /NORMAL, (x1+x0)/2d0-wid/2, 0.95, radstr, charsize=csize2,charthick=1,color=0
    xyouts, /NORMAL, x0-0.02, 0.93, 'b', charsize=csize2,charthick=1,color=0
    ;xyouts, /normal, x0+(0.445-0.13)/3d0+wid, 0.945, "!9!Z(6E)!X", font=-1, charsize=0.75

    xyouts,/NORMAL,x0-0.01/3d0,0.03,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0,width=wid
    xyouts,/NORMAL,x1-0.01/3d0-wid,0.03,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    slice = reform(data[*,*,rad,quan[2],rec])
    ;datanew = dblarr(ntnew,npnew)
    ;datanew[0:nth-1,0:nph-1]=slice
    ;datanew[0:nth-1,nph:npnew-1]=slice
    ;datanew[nth:ntnew-1,0:nph-1]=reverse(slice,1)
    ;datanew[nth:ntnew-1,nph:npnew-1]=reverse(slice,1)
    mnv = 0.33d0*min(slice)
    mxv = -mnv ;0.5d0*max(slice)

    draw = transpose(reform(slice))

    x0 = x1+0.04
    x1 = x0+dx
    !P.POSITION=[x0,y0,x1,y1]
    choose_color, /mag_color, /quiet
    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox

    colourbar2,[x0-0.03,0.38,x0-0.0175,0.62],mx=mxv,mn=mnv,plttype='ps'
    xyouts,/NORMAL,x0+0.12,0.025,textoidl('\pm')+'!U !N1!U !NkG',charsize=csize,charthick=1.0,color=0

    choose_color, /ash_color, /quiet
    ;radstr = strtrim(string(r[rad]/rstar,format='(f4.2)'),2)+'!U !NR'
    radstr = textoidl('B_{\phi}')
    ;print, 'radstr=',radstr
    xyouts, /NORMAL, (x1+x0)/2d0, 0.95, radstr, charsize=csize2,charthick=1,color=255,width=wid
    xyouts, /NORMAL, (x1+x0)/2d0-wid/2, 0.95, radstr, charsize=csize2,charthick=1,color=0
    xyouts, /NORMAL, x0-0.02, 0.93, 'c', charsize=csize2,charthick=1,color=0
    ;xyouts, /normal, x0+(0.445-0.13)/3d0+wid, 0.945, "!9!Z(6E)!X", font=-1, charsize=0.75

    xyouts,/NORMAL,x0-0.01/3d0,0.03,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0,width=wid
    xyouts,/NORMAL,x1-0.01/3d0-wid,0.03,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'
    stop
end
