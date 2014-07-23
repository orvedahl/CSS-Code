pro bsd2011_omega, filename, omega0

 ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr1/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius = r, data = values, QUANTITIES = quantities ;, /oldfile

    num_records = n_elements(values(0,0,*))
    ;Determine the number of radial points.
    N_R = n_elements(r)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

    data = dblarr(N_R,N_Q,num_records*n_files)
    data(*,0:N_Q-1,0:num_records-1) = values

    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_avg, files(i), data = values ;, /oldfile
            num_rec = n_elements(values(0,0,*))
            data(*,0:N_Q-1,num_records:num_records+num_rec-1) = values
            num_records = num_records + num_rec
        endfor
    endif

    ;Get vphiavg
    vphi = dblarr(N_R)
    for ir=0,N_R-1 do begin
        vphi(ir) = mean(data(ir,6,*))
    endfor
    
    obar = (omega0+vphi/r)*1d9/(2d0*!dpi)+27d0

;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr1/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius = r2, data = values, QUANTITIES = quantities ;, /oldfile

    num_records = n_elements(values(0,0,*))
    ;Determine the number of radial points.
    N_R = n_elements(r2)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

    data = dblarr(N_R,N_Q,num_records*n_files)
    data(*,0:N_Q-1,0:num_records-1) = values

    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_avg, files(i), data = values ;, /oldfile
            num_rec = n_elements(values(0,0,*))
            data(*,0:N_Q-1,num_records:num_records+num_rec-1) = values
            num_records = num_records + num_rec
        endfor
    endif

    ;Get vphiavg
    vphi = dblarr(N_R)
    for ir=0,N_R-1 do begin
        vphi(ir) = mean(data(ir,6,*))
    endfor
    
    obar2 = (omega0+0.5d0*vphi/r2)*1d9/(2d0*!dpi)+25d0

    omnew = 2.666d-6
    conang = 1d9*omnew*(min(r)/r)^2/2d0/!dpi+55d0

    rsun = 6.96d10
    rmid = 0.94d0*rsun
    func = 1d0/(1d0+exp(-50d0*(r2-rmid)/(max(r2)-min(r2))))
    solarang = 470d0-7d3*(abs(1d0-r2/rmid)^2)*func-3d0*(abs(r2-rmid)/(max(r2)-min(r2)))^1.5*(1d0-func)

    xsize=4.5d0*2.54d0
    ysize=xsize
    !P.FONT=0
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Times-Roman'
    choose_color, /rainbow, /quiet
    x0 = 0.2
    y0 = 0.14
    x1 = 0.9
    y1 = 0.9
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    !P.CHARSIZE=1.0
    !P.THICK=4
    !X.THICK=4
    !Y.THICK=4
    rmin = min(r2)
    rmax = max(r2)
 ;Create a plot of omega at slices
 ;degree angles (from equator = 0)
    max_omega=max([max(obar2),max(conang)])
    min_omega=min([min(obar2),min(conang)])
    
    linestyles = [0,0,0,2];,0,0]
    colors = [0,40,200,200];,190,0]
    thickn = [3,3,3,3] ;,3,3,5]
    rsun = 6.96d10
    
    plot, r2/rsun, solarang, xtitle=textoidl('r/R'), ytitle=textoidl('\Omega (nHz)'), $
     linestyle=linestyles[0], xs=1, charthick=1, color=0,  yrange=[400d0,480d0], /ys ;[min_omega,max_omega]
    oplot, r2/rsun, obar2, linestyle=linestyles[1], color=colors[1]
    oplot, r/rsun, obar, linestyle=linestyles[2], color=colors[2]
    oplot, r/rsun, conang, linestyle=linestyles[3], color=colors[3]
    xyouts, /normal, 0.5725, 0.045, "!9!Z(6E)!X", font=-1, charsize=0.75
    names = ['Sun','Coll.','No Coll.','C. Ang.']
    legend,names,COLORS=[colors], LIN=[linestyles],/bottom,/left,charthick=1
    
    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'evince '+filename+' &'

end
