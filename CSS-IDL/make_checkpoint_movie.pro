@display_shell_slice.pro
@read_checkpoint_file.pro

pro make_checkpoint_movie, iter, css_case, iq, xsize=xsize, ysize=ysize, dir=dir, ct=ct

    if not keyword_set(xsize) then begin
        xsize = 512
    endif

    if not keyword_set(ysize) then begin
        ysize=xsize
    endif

    if not keyword_set(ct) then begin
        ct = 'ct_az_avg.bin'
    endif

    rsun = 6.96d10
 ;We assume the files read in are from the same run.
    read_checkpoint_file, iter, css_case, iq, phi=phi, theta = theta, rad= r, data = data
    
 ;Determine the number of radial points.
    N_R = n_elements(r)
    steps = (N_R+95L)
 ;Determine the number of meridional points.
    N_Theta = n_elements(theta)
    
    N_Phi = n_elements(phi)
       
    lon = phi*180d0/!dpi
    lat = reverse(90d0-theta*!radeg)
    
    print, 'rendering...'
    ctr=bytarr(256) & ctg=ctr & ctb=ctr
    close,3 & openr,3,ct
    readu,3,ctr,ctg,ctb
    tvlct,ctr,ctg,ctb
    close,3
    
    set_plot, 'Z', /Copy
    device, set_resolution=[xsize,ysize],Z_Buffer=0
    
    print, 'Percent complete:'
    k=0L
    time = 0d0
    For ir=N_R-2,0,-1 Do Begin
 ;We assume the files read in are from the same run.
       limits=[min(lat),min(lon),max(lat),max(lon)]
       erase
       !P.MULTI=0
       !P.NOERASE=1
       fill = bytarr(xsize,ysize)
       fill[*,*] = 255b
       tv,fill,/device
       
       draw = Transpose(Reform(data[*,*,ir]))
       nlevs = 15
       levs=get_shell_slice_levels(draw,lat,nlevs,mn=mnv,mx=mxv,pct=0.02, /uneven_range)
       
       print, mnv, mxv
       
       !P.POSITION=[0.02,0.02,0.98,0.98]
       r1 = r[ir]
       display_shell_slice, lat, lon, draw, /publication, minv=mnv, maxv=mxv, /nolines, $
                            res_scale=0.3, orthographic=[0.5d0*(min(lon)+max(lon)),0.5d0*(min(lat)+max(lat))], latlon_window=limits, /white
       
       xyouts, 10, 10, strtrim(string(r1/rsun,FORMAT='(G10.2)'),2)+'R'+sunsymbol(), /device, charthick = 2, charsize=1.5, color=0
       
       image3D = TVRD()
       If ((ir eq 15) or (ir eq 63) or (ir eq 111)) Then Begin
          for it=0,31 do begin 
             write_png, strtrim(string(k,form='(i7.7)'),2)+'.png', image3D,ctr,ctg,ctb
             k = k+1L
          endfor
          print, ' '
          print, strtrim(string((1.0d2*k)/(1.0d0*steps-1.0d0),FORMAT='(G10.4)'),2)+'%'
          print, ' '
       Endif Else Begin
          write_png, strtrim(string(k,form='(i7.7)'),2)+'.png', image3D,ctr,ctg,ctb
          print, ' '
          print, strtrim(string((1.0d2*k)/(1.0d0*steps-1.0d0),FORMAT='(G10.4)'),2)+'%'
          print, ' '
          k = k+1L
       EndElse          
    EndFor
 end
