pro stability_regions, filename, omega0, Cp, stop=stop

 ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr3/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius = r, data = values, QUANTITIES = quantities

    num_records = n_elements(values(0,0,*))
    ;Determine the number of radial points.
    N_R = n_elements(r)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

    data = dblarr(N_R,N_Q,num_records*n_files)
    data(*,0:N_Q-1,0:num_records-1) = values

    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_avg, files(i), data = values
            num_rec = n_elements(values(0,0,*))
            data(*,0:N_Q-1,num_records:num_records+num_rec-1) = values
            num_records = num_records + num_rec
        endfor
    endif

    ;Get vphiavg
    vphi = dblarr(N_R)
    dsdr = dblarr(N_R)
    for ir=0,N_R-1 do begin
        vphi(ir) = mean(data(ir,6,*))
        dsdr(ir) = mean(data(ir,8,*))
    endfor
    bigG = 6.667d-8
    mass = 1.98d33
    grav = bigG*mass/r^2

    n2 = grav*dsdr/Cp

    obar = (omega0+vphi/r)*1d9/(2d0*!dpi)+40d0
    conang = obar(N_R/2)*(r(N_R/2)/r)^2+6d0
    ;obar = (obar-470d0)*50d0/(max(obar)-min(obar))+475d0
    onorm = 1d9*omega0
    rsun = 6.96d10
    rmid = 0.94d0*rsun
    func = 1d0/(1d0+exp(-50d0*(r-rmid)/(max(r)-min(r))))
    solarang = 470d0-7d3*(abs(1d0-r/rmid)^2)*func-3d0*(abs(r-rmid)/(max(r)-min(r)))^1.5*(1d0-func)

    Rwca = r*deriv(r,conang^2)
    Rwsa = r*deriv(r,solarang^2)

    n2 = n2/omega0^2
    Rwca = Rwca/onorm^2
    Rwsa = Rwsa/onorm^2

    xsize=10d0
    ysize=10d0
    !P.FONT=1
    Set_Plot, "PS"
    Device, filename='stabreg.eps', /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    choose_color, /rainbow, /quiet
    x0 = 0.2
    y0 = 0.14
    x1 = 0.9
    y1 = 0.9
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    r1 = min(r)
    r2 = max(r)

    ymx = max([max(Rwsa),max(Rwca)])
    ymn = min([min(Rwsa),min(Rwca)])

    xmx = max(n2)
    xmn = min(n2)
    ;xmx = 0d0
    ;xmn = -100d0
    
    linestyles = [2,0]
    colors = [0,0]
    thickn = [3,3] 
    sunsym = sunsymbol()
    rsun = 6.96d10
    
    plot, n2, Rwsa, xtitle=textoidl('N^2/Omega^2'), ytitle=textoidl('R_w/\Omega^2'), title='Stability Regions', $
      yrange=[ymn,ymx], xrange=[xmn,xmx], linestyle=linestyles[0], /xs, /ys, thick=2, xthick=3, ythick=3, charthick=3, color=0
    oplot, n2, Rwca, linestyle=linestyles[1], thick=thickn[1], color=colors[1]
    oplot, n2, -n2, thick=thickn[1],color=30
    ;names = ['Sun','C. Ang.']
    ;legend,names,COLORS=[colors], LIN=[linestyles],/bottom,/left,charthick=2,thick=thickn
    
    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview stabreg.eps &'

stop

end
