pro load_fields
    common sohofig, urnwtt, ursldt, urnwtb, ursldb, phi, theta, r

    ;Select files to load
    ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr3/augustso/CSSDE/solar_wvh_bio/Shell_Slices'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
        if(files eq '') then return
    endif

    read_shell_slice, files(0), phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities
    print, 'file loaded'
    print, 'size of data', size(data)
    num_records = n_elements(data(0,0,0,0,*))
    N_R = n_elements(r)
    N_Theta = n_elements(theta)
    N_Phi = n_elements(phi)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)
    ;pick out radial velocities
    urnwtt = dblarr(N_Phi,N_Theta,long(num_records)*long(n_files))
    urnwtt(*,*,0:num_records-1) = reform(data(*,*,N_R-1,4,*))
    urnwtb = dblarr(N_Phi,N_Theta,long(num_records)*long(n_files))
    urnwtb(*,*,0:num_records-1) = reform(data(*,*,N_R/2-1,4,*))
    nrec = num_records
    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_slice, files(i),data = data
            urnwtt(*,*,nrec:nrec+num_records-1) = reform(data(*,*,N_R-1,4,*))
            urnwtb(*,*,nrec:nrec+num_records-1) = reform(data(*,*,N_R/2-1,4,*))
            nrec = nrec+num_records
        endfor
    endif

    ;Select files to load
    ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr3/augustso/CSSDE/hyper_test4/Shell_Slices'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
        if(files eq '') then return
    endif

    read_shell_slice, files(0), phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities
    print, 'file loaded'
    print, 'size of data', size(data)
    num_records = n_elements(data(0,0,0,0,*))
    N_R = n_elements(r)
    N_Theta = n_elements(theta)
    N_Phi = n_elements(phi)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)
    ;pick out radial velocities
    ursldt = dblarr(N_Phi,N_Theta,long(num_records)*long(n_files))
    ursldt(*,*,0:num_records-1) = reform(data(*,*,N_R-1,4,*))
    ursldb = dblarr(N_Phi,N_Theta,long(num_records)*long(n_files))
    ursldb(*,*,0:num_records-1) = reform(data(*,*,N_R/2-1,4,*))
    nrec = num_records
    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_slice, files(i),data = data
            ursldt(*,*,nrec:nrec+num_records-1) = reform(data(*,*,N_R-1,4,*))
            ursldb(*,*,nrec:nrec+num_records-1) = reform(data(*,*,N_R/2-1,4,*))
            nrec = nrec+num_records
        endfor
    endif
end

pro soho2010_fig1a
    common sohofig, urnwtt, ursldt, urnwtb, ursldb, phi, theta, r

    nnwt = n_elements(urnwtt(0,0,*))
    nsld = n_elements(ursldt(0,0,*))
    nth = n_elements(theta)

    urnwtm = 0d0*urnwtt
    for in=0L,nnwt-1L do begin
        urnwtm(*,*,in) = urnwtt(*,*,in)-mean(urnwtt(*,*,in))
    endfor

    ursldm = 0d0*ursldt
    for in=0L,nsld-1L do begin
        ursldm(*,*,in) = ursldt(*,*,in)-mean(ursldt(*,*,in))
    endfor

    ph1 = phi(0)
    ph2 = max(phi)
    nph = n_elements(phi)
    dph = (ph2-ph1)/double(nph-1)
    rx = (ph2-ph1)/(dindgen(nph/2)+1d0)/1d8

    ;Compute spectra
    urnwtfft = dcomplexarr(nph,nth)
    ursldfft = dcomplexarr(nph,nth)
    urnwtfft(*,*) = dcomplex(0d0,0d0)
    ursldfft(*,*) = dcomplex(0d0,0d0)
    for in=0L,nnwt-1L do begin
        tfft = fft(reform(urnwtm(*,*,in)),/double)
        urnwtfft = urnwtfft+tfft
    endfor
    urnwtfft = urnwtfft/double(nnwt)

    tfa = total(urnwtfft,2)
    tfa = tfa(0:nph/2L-1L)+reverse(tfa(nph/2L:*))
    urnwta = abs(tfa)

    for in=0L,nsld-1L do begin
        tfft = fft(reform(ursldm(*,*,in)),/double)
        ursldfft = ursldfft+tfft
    endfor
    ursldfft = ursldfft/double(nsld)

    tfa = total(ursldfft,2)
    tfa = tfa(0:nph/2L-1L)+reverse(tfa(nph/2L:*))
    urslda = abs(tfa)

    norm = max([max(urnwta),max(urslda)])
    urnwta = urnwta/max(urnwta);norm
    urslda = urslda/max(urslda);norm

    rx1 = rx*r(3)
    nr = n_elements(r)
    rx2 = rx*r(nr-1)

    xsize=16d0
    ysize=16d0
    !P.MULTI=0
    !P.FONT=1
    !P.CHARSIZE=2
    Set_Plot, "PS"
    filename = 'soho2010_power.eps'
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    !P.POSITION=[0.15,0.15,0.95,0.95]
    choose_color, /ash_color, /quiet

    plot, rx2, 0d0*urslda, /xlog, /xs, xtitle='k (Mm)', ytitle='Power', thick=1, xthick=3, ythick=3
    oplot, rx2, smooth(urnwta,5), color=125, thick=3
    oplot, rx2, smooth(urslda,5), thick=3

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

    stop

end

pro soho2010_fig1b
    common sohofig, urnwtt, ursldt, urnwtb, ursldb, phi, theta, r

    nnwt = n_elements(urnwtt(0,0,*))
    nsld = n_elements(ursldt(0,0,*))
    nth = n_elements(theta)

    urnwttm = 0d0*urnwtt
    for in=0L,nnwt-1L do begin
        urnwttm(*,*,in) = urnwtt(*,*,in)-mean(urnwtt(*,*,in))
    endfor

    ursldtm = 0d0*ursldt
    for in=0L,nsld-1L do begin
        ursldtm(*,*,in) = ursldt(*,*,in)-mean(ursldt(*,*,in))
    endfor

    urnwtbm = 0d0*urnwtb
    for in=0L,nnwt-1L do begin
        urnwtbm(*,*,in) = urnwtb(*,*,in)-mean(urnwtb(*,*,in))
    endfor

    ursldbm = 0d0*ursldb
    for in=0L,nsld-1L do begin
        ursldbm(*,*,in) = ursldb(*,*,in)-mean(ursldb(*,*,in))
    endfor

    minn = min([min(urnwttm),min(ursldtm)])
    maxn = max([max(urnwttm),max(ursldtm)])
    nph = n_elements(phi)
    nth = n_elements(theta)
    ntot = double(nth)*double(nph)
    urnwtth = histogram(urnwttm(*,*,0),nbins=200,loc=ust,min=minn,max=maxn)/ntot
    for in=1,nnwt-1L do begin
        urnwtth = urnwtth + histogram(urnwttm(*,*,in),nbins=200,loc=ust,min=minn,max=maxn)/ntot
    endfor
    urnwtth = urnwtth/double(nnwt)+1d-7

    ursldth = histogram(ursldtm(*,*,0),nbins=200,loc=ust,min=minn,max=maxn)/ntot
    for in=1,nsld-1L do begin
        ursldth = ursldth + histogram(ursldtm(*,*,in),nbins=200,loc=ust,min=minn,max=maxn)/ntot
    endfor
    ursldth = ursldth/double(nsld)+1d-7
    ust = ust/1d5

    minn = min([min(urnwtbm),min(ursldbm)])
    maxn = max([max(urnwtbm),max(ursldbm)])
    nph = n_elements(phi)
    nth = n_elements(theta)
    ntot = double(nth)*double(nph)
    urnwtbh = histogram(urnwtbm(*,*,0),nbins=200,loc=usb,min=minn,max=maxn)/ntot
    for in=1,nnwt-1L do begin
        urnwtbh = urnwtbh + histogram(urnwtbm(*,*,in),nbins=200,loc=usb,min=minn,max=maxn)/ntot
    endfor
    urnwtbh = urnwtbh/double(nnwt)+1d-7

    ursldbh = histogram(ursldbm(*,*,0),nbins=200,loc=usb,min=minn,max=maxn)/ntot
    for in=1,nsld-1L do begin
        ursldbh = ursldbh + histogram(ursldbm(*,*,in),nbins=200,loc=usb,min=minn,max=maxn)/ntot
    endfor
    ursldbh = ursldbh/double(nsld)+1d-7
    usb = usb/1d5

    xsize=16d0
    ysize=16d0
    !P.MULTI=0
    !P.FONT=1
    !P.CHARSIZE=3
    Set_Plot, "PS"
    filename = 'soho2010_pdf.eps'
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    !P.POSITION=[0.2,0.21,0.95,0.96]
    choose_color, /ash_color, /quiet

    plot, ust, 0d0*ursldth, /xs, xtitle=textoidl('u_r (km s^{-1})'), thick=1, xthick=5, ythick=5, /ylog, $
      yrange=[1d-5,1d-1], xrange=[-2.5d0,1d0], xticks=3, xtickv=[-2,-1,0,1]
    oplot, ust, urnwtth, color=125, thick=10
    oplot, ust, ursldth, thick=10
    oplot, usb, urnwtbh, color=125, thick=10, linestyle=2
    oplot, usb, ursldbh, thick=10, linestyle=2

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

    stop

end
