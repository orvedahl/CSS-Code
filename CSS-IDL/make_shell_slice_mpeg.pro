@display_shell_slice.pro

pro make_shell_slice_mpeg, filename, iq, ir, xsize=xsize, ysize=ysize, frames=frames, dir=dir, ct=ct, dt=dt, nstep=nstep, exist=exist, twodepth=twodepth, jpeg=jpeg, png=png, minmax=minmax, rmmean=rmmean

    If (not keyword_set(dir)) Then Begin
        dir = '/freyr3/augustso/CSSDE/'
    Endif

    if not keyword_set(xsize) then begin
        xsize = 700
    endif
    if not keyword_set(ysize) then begin
        ysize=xsize
    endif
    if not keyword_set(frames) then begin
        frames = 10
    endif

    if not keyword_set(ct) then begin
        ct = 'ct_az_avg.bin'
    endif

    if not keyword_set(dt) then begin
        dt = 6.0684d0
    endif

    if not keyword_set(nstep) then begin
       nstep=100
    endif

    rsun = 6.96d10

    If (not keyword_set(exist)) Then Begin
    files = DIALOG_PICKFILE(PATH=dir,TITLE='Select Shell Slice Files to Read',/MULTIPLE_FILES)
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(file eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_slice, files(0), phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities

    ;Determine the number of radial points.
    N_R = n_elements(r)

    ;Determine the number of meridional points.
    N_Theta = n_elements(theta)

    N_Phi = n_elements(phi)

    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

    num_records = n_elements(data[0,0,0,0,*])

    lon = phi*180d0/!dpi
    lat = reverse(90d0-theta*!radeg)

    nlevs = 15
    If (keyword_set(twodepth)) Then Begin
        mnv = dblarr(2)
        mxv = dblarr(2)
        levs=get_shell_slice_levels(transpose(reform(data[*,*,ir[0],iq,0])),lat,nlevs,mn=mnt,mx=mxt,pct=0.02, /uneven_range)
        mnv[0]=mnt
        mxv[0]=mxt
        levs=get_shell_slice_levels(transpose(reform(data[*,*,ir[1],iq,0])),lat,nlevs,mn=mnt,mx=mxt,pct=0.02, /uneven_range)
        mnv[1]=mnt
        mxv[1]=mxt
    Endif Else Begin
       If (keyword_set(rmmean)) Then Begin
          levs=get_shell_slice_levels(transpose(reform(data[*,*,ir,iq,0]-mean(data[*,*,ir,iq,0]))),lat,nlevs,mn=mnv,mx=mxv,pct=0.02, /uneven_range)
       Endif Else Begin
          levs=get_shell_slice_levels(transpose(reform(data[*,*,ir,iq,0])),lat,nlevs,mn=mnv,mx=mxv,pct=0.02, /uneven_range)
       EndElse
    EndElse

    print, mnv, mxv
    If (keyword_set(minmax)) Then Begin
       mnv = minmax[0]
       mxv = minmax[1]
    EndIf

    n_files_recs = (1L*num_records)*(1L*n_files)
    int_steps = lindgen(n_files_recs)
    step_text = ash_12_renum(int_steps, ash_12_numbering_starts=1, /quiet)

    print, 'rendering...'
    ctr=bytarr(256) & ctg=ctr & ctb=ctr
    close,3 & openr,3,ct
    readu,3,ctr,ctg,ctb
    tvlct,ctr,ctg,ctb
    close,3

    set_plot, 'Z', /Copy
    device, set_resolution=[xsize,ysize],Z_Buffer=0

    print, 'reading Shell files from ', dir
    print, 'Percent complete:'
    k=0L
    time = 0d0
    For ifile=0,n_files-1 Do Begin
                                ;We assume the files read in are from the same run.
        read_shell_slice, files(ifile), data = data
        num_records = n_elements(data[0,0,0,0,*])
        For irec=0,num_records-1 Do Begin
            limits=[min(lat),min(lon),max(lat),max(lon)]
            time = time+dt*nstep/3600d0
            
;            mnv = min(data(*,*,ir,iq,irec))
;            mxv = max(data(*,*,ir,iq,irec))

            erase
            !P.MULTI=0
            !P.NOERASE=1
            fill = bytarr(xsize,ysize)
            fill[*,*] = 255b
            tv,fill,/device

            If (keyword_set(twodepth)) Then Begin
                !P.POSITION=[0.04,0.04,0.475,0.96]
                r1 = r[ir[0]]
                draw = Transpose(Reform(data[*,*,ir[0],iq,irec]))
                display_shell_slice, lat, lon, draw, /publication, minv=mnv[0], maxv=mxv[0], /nolines, $
                  res_scale=0.3, orthographic=[0.5d0*(min(lon)+max(lon)),0.5d0*(min(lat)+max(lat))], latlon_window=limits, /white
            Endif Else Begin
                !P.POSITION=[0.05,0.05,0.95,0.95]
                r1 = r[ir]
                If (keyword_set(rmmean)) Then Begin
                   draw = Transpose(Reform(data[*,*,ir,iq,irec]-mean(data[*,*,ir,iq,irec])))
                Endif Else Begin
                   draw = Transpose(Reform(data[*,*,ir,iq,irec]))
                EndElse
                display_shell_slice, lat, lon, draw, /publication, minv=mnv, maxv=mxv, /nolines, $
                  res_scale=0.3, orthographic=[0.5d0*(min(lon)+max(lon)),0.5d0*(min(lat)+max(lat))], latlon_window=limits, /white
            EndElse


            xyouts, 5, 5, strtrim(string(time,FORMAT='(G10.3)'),2)+' hours', /dev, charthick = 2, charsize=1.5, color=0
            xyouts, 0.515/2d0-0.03, 0.975, strtrim(string(r1/rsun,FORMAT='(G10.2)'),2)+'R'+sunsymbol(), /normal, charthick = 2, charsize=1.5, color=0

            If (keyword_set(twodepth)) Then Begin
                !P.POSITION=[0.525,0.04,0.96,0.96]
                draw = Transpose(Reform(data[*,*,ir[1],iq,irec]))
                display_shell_slice, lat, lon, draw, /publication, minv=mnv[1], maxv=mxv[1], /nolines, $
                  res_scale=0.3, orthographic=[0.5d0*(min(lon)+max(lon)),0.5d0*(min(lat)+max(lat))], latlon_window=limits, /white
                xyouts, 5, 5, strtrim(string(time,FORMAT='(G10.3)'),2)+' hours', /dev, charthick = 2, charsize=1.5, color=0
                r2 = r[ir[1]]
                xyouts, (0.525+0.96)/2d0-0.03, 0.975, strtrim(string(r2/rsun,FORMAT='(G10.2)'),2)+'R'+sunsymbol(), /normal, charthick = 2, charsize=1.5, color=0
            EndIf

            image3D = TVRD()
            If (keyword_set(jpeg)) Then Begin
                image = bytarr(3,xsize,ysize)
                image[0,*,*] = ctr[image3D]
                image[1,*,*] = ctg[image3D]
                image[2,*,*] = ctb[image3D]
                write_jpeg, step_text(k)+'.jpg', image,True=1,Quality=100
            Endif

            If (keyword_set(png)) Then Begin
                write_png, step_text(k)+'.png', image3D,ctr,ctg,ctb
            EndIf

            print, ' '
            print, strtrim(string((1.0d2*k)/(1.0d0*n_files*num_records-1.0d0),FORMAT='(G10.4)'),2)+'%'
            print, ' '
            k = k+1L

        EndFor
    EndFor
EndIf

bitrate=fix(double(xsize)*double(ysize)*double(frames)/1250.0d0)

if (keyword_set(jpeg)) Then Begin
    spawn, "mencoder mf://'*.jpg' -mf" + ' w=' + strtrim(string(xsize),2)+':h=' + strtrim(string(ysize),2) + ':fps='+strtrim(string(frames),2)+':type=jpg -ovc lavc -lavcopts vcodec=ffv1:vbitrate='+strtrim(string(bitrate),2)+':autoaspect=0 -nosound -o ' + filename

    print, "mencoder mf://'*.jpg' -mf" + ' w=' + strtrim(string(xsize),2)+':h=' + strtrim(string(ysize),2) + ':fps='+strtrim(string(frames),2)+':type=jpg -ovc lavc -lavcopts vcodec=ffv1:vbitrate='+strtrim(string(bitrate),2)+':autoaspect=0 -nosound -o ' + filename
endif

if (keyword_set(png)) Then Begin
    spawn, "mencoder mf://'*.png' -mf" + ' w=' + strtrim(string(xsize),2)+':h=' + strtrim(string(ysize),2) + ':fps='+strtrim(string(frames),2)+':type=png -ovc lavc -lavcopts vcodec=ffv1:vbitrate='+strtrim(string(bitrate),2)+':autoaspect=0 -nosound -o ' + filename

    print, "mencoder mf://'*.png' -mf" + ' w=' + strtrim(string(xsize),2)+':h=' + strtrim(string(ysize),2) + ':fps='+strtrim(string(frames),2)+':type=png -ovc lavc -lavcopts vcodec=ffv1:vbitrate='+strtrim(string(bitrate),2)+':autoaspect=0 -nosound -o ' + filename
endif

end
