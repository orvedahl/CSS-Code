pro load_fields
common fields, css_ur, surf_ur, css_lon, css_lat, surf_lon, surf_lat
    rstar = 6.96d10
    root_directory = '/freyr3/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select CSS Shell Slice File to Read')
    
    n_files = n_elements(files)
    
    if(n_files eq 1) then begin
        if(files eq '') then stop
    endif

    read_shell_slice, files(0), phi=phi, theta = theta, radii = r, data = temp
    print, 'css_radii=', r/rstar
    numrec = n_elements(temp(0,0,0,0,*))
    nr = n_elements(r)
    nth = n_elements(theta)
    nph = n_elements(phi)
    nq = n_elements(quantities)

    css_ur = dblarr(nth,nph)
    css_ur(*,*) = reform(temp(*,*,nr-1,4,numrec-1))
    
    css_lon = phi*180d0/!dpi
    css_lat = reverse(90d0-theta*!radeg)

    ;Resort Surface convection data
    restore, '/scratch/trampeda/scratch/bob/96x96x20/Uy.728.4.save' ;'/freyr3/augustso/SteinData/Uy.318.6-462.0.save'

    surf_ur_tmp = -reform(uy(*,265,*)) ;40 is about photosphere pick out 4 degree patch
    surf_ur = dblarr(2500,1250)
    surf_ur(0:999,0:999) = surf_ur_tmp
    surf_ur(1000:1999,0:999) = surf_ur_tmp
    surf_ur(2000:2499,0:999) = surf_ur_tmp(0:499,*)
    surf_ur(*,1000:1249) = surf_ur(*,0:249)

    surf_ur = surf_ur*1d8/60d0
    
    ;Build fake coordinates based up knowledge that it is 4 degrees square, assume uniform grid
    nz = n_elements(surf_ur(*,0))
    nx = n_elements(surf_ur(0,*))
    dlon = 10d0/double(nx-1)
    dlat = 20d0/double(nz-1)
    surf_lon = dlon*dindgen(nx)-5d0
    surf_lat = dlat*dindgen(nz)

    indsn = where(surf_ur lt 0d0, complement=indsp)
    meansn = Mean(surf_ur(indsn))
    meansp = Mean(surf_ur(indsp))

    indcn = where(css_ur lt 0d0, complement=indcp)
    css_ur(indcn) = css_ur(indcn)*min(surf_ur(indsn))/min(css_ur(indcn))
    meancn = Mean(css_ur(indcn))
    meancp = Mean(css_ur(indcp))    
    ;css_ur(indcn) = css_ur(indcn)*meansn/meancn
    css_ur(indcp) = css_ur(indcp)*meansp/meancp
    ;css_ur = css_ur*max(abs(surf_ur))/max(abs(css_ur))

    ;stop
end

pro css_msc_compare, filename
common fields, css_ur, surf_ur, css_lon, css_lat, surf_lon, surf_lat

    space=0.01d0

    rstar = 6.96d10
    xsize = 6d0*2.54d0
    ysize = xsize

    !P.MULTI=0
    !P.FONT=1
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    choose_color, /ash_color, /quiet

    x0 = space
    y0 = 0.1
    x1 = x0+0.5d0-3d0*space
    y1 = 0.9
    
    limits=[min(surf_lat),min(surf_lon),max(surf_lat),max(surf_lon)]

    mxv = 0.9*max(surf_ur)
    mnv = 0.3d0*min(surf_ur)
    draw = transpose(surf_ur)

    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = [min(surf_lat),min(surf_lon),min(surf_lat)+8d0,min(surf_lon)+8d0]
    display_css_shell_slice, surf_lat, surf_lon, draw, /publication,  minv=mnv, maxv=mxv, res_scale=0.3, /nolines, $
      orthographic=[0.5d0*(min(surf_lon)+max(surf_lon)),0.5d0*(min(surf_lat)+max(surf_lat))], latlon_window=limits, highlightbox=highlightbox

    ;xyouts,/NORMAL,space,0.925,'a',charsize=1.5,charthick=1.0,color=0    
    xyouts,/NORMAL,x0+space,0.05,textoidl('-5\circ'),charsize=2,charthick=1.0,color=0
    xyouts,/NORMAL,x1-8d0*space,0.05,textoidl('+5\circ'),charsize=2,charthick=1.0,color=0

    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]    

    mxv = 0.75d0*max(css_ur)
    mnv = 0.15d0*min(css_ur)
    draw = transpose(css_ur)

 ;Figure out bounds from size of domain
    x0 = x1+3d0*space
    y0 = 0.1d0
    x1 = x0+0.5d0-3d0*space
    y1 = 0.9d0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, minv=mnv, maxv=mxv, res_scale=0.3, /nolines, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox

    ;xyouts,/NORMAL,x0,0.925,'b',charsize=1.5,charthick=1.0,color=0
    xyouts,/NORMAL,x0+space,0.05,textoidl('-5\circ'),charsize=2,charthick=1.0,color=0
    xyouts,/NORMAL,x1-8d0*space,0.05,textoidl('+5\circ'),charsize=2,charthick=1.0,color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'
    stop
end
