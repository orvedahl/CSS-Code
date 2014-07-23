pro flux_bal_timeavg, css_case=css_case, Lstar=Lstar, no_fac=no_fac, hyper=hyper, psf=psf, rstar=rstar, $
                      yrange=yrange, dtavg=dtavg, botio=botio, topio=topio, stop=stop

    If (not keyword_set(Lstar)) Then Lstar = 3.86d33 ;solar
    If (not keyword_set(rstar)) Then rstar = 6.96d10 ;solar
    If (not keyword_set(dtavg)) Then dtavg = 1d0

    ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr1/augustso/CSSDE/'
    If (keyword_set(css_case)) Then Begin
        root_directory=root_directory+strtrim(css_case)
    EndIf
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius = r, data = values, QUANTITIES = quantities

    num_records = 1L*n_elements(values(0,0,*))

    ;Determine the number of radial points.
    N_R = n_elements(r)

    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

    data = dblarr(N_R,N_Q,num_records*n_files)
    data(*,0:N_Q-1,0:num_records-1) = values

    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_avg, files(i), data = values
            num_rec = 1L*n_elements(values(0,0,*))
            data(*,0:N_Q-1,num_records:num_records+num_rec-1) = values
            num_records = num_records + num_rec
        endfor
    endif

    quantity_names = strarr(N_Q)
    for i=0,N_Q-1 do begin
        quantity_names(i) = get_quantity_name(quantities(i))
    endfor

    If (not keyword_set(hyper)) Then Begin
        fksavg = dblarr(N_R)
        fksavg(*) = 0d0
        fenavg = fksavg*0d0
        fkeavg = fksavg*0d0

        iks = where(quantities eq 15)
        ien = where(quantities eq 16)
        ike = where(quantities eq 17)

        for ir=0L,N_R-1L do begin
            fksavg(ir) = mean(data(ir,iks,*))
            fenavg(ir) = mean(data(ir,ien,*))
            fkeavg(ir) = mean(data(ir,ike,*))
        endfor

        Lks = 4d0*!dpi*r^2*fksavg
        Len = 4d0*!dpi*r^2*fenavg
        Lke = 4d0*!dpi*r^2*fkeavg
        fvdavg = fksavg*0d0
        facavg = fksavg*0d0
        ivd = where(quantities eq 19)
        iac = where(quantities eq 18)
        for ir=0L,N_R-1L do begin
            fvdavg(ir) = mean(data(ir,ivd,*))
            facavg(ir) = mean(data(ir,iac,*))
        endfor
        If (keyword_set(no_fac)) Then facavg = 0d0*facavg
        Lvd = 4d0*!dpi*r^2*fvdavg
        Lac = 4d0*!dpi*r^2*facavg
        Ltot = Len+Lke+Lvd

        minv = min([min(Lks),min(Len),min(Lke),min(Lvd),min(Lac)])
        maxv = max([max(Lks),max(Len),max(Lke),max(Lvd),max(Lac)])
        ;r = r/rstar*100d0
        docolor
        loadct, 13
        colors = [80,110,140,170,200,230]
        plot, r, 0d0*r+1d0, yrange=[minv/Lstar,maxv/Lstar]
        oplot, r, Ltot/Lstar, color=colors[0]
        oplot, r, Lks/Lstar, color=colors[1]
        oplot, r, Len/Lstar, color=colors[2]
        oplot, r, Lke/Lstar, color=colors[3]
        oplot, r, Lvd/Lstar, color=colors[4]
        oplot, r, Lac/Lstar, color=colors[5]

        If (keyword_set(psf)) Then Begin
            xsize=16d0
            ysize=16d0
            !P.MULTI=0
            !P.FONT=1
            !P.CHARSIZE=3
            colors = [30,80,150,225,240]
            Set_Plot, "PS"
            filename = 'flux_bal_timeavg_visc.eps'
            Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, $
              SET_FONT='Helvetica Italic', /TT_FONT
            !P.POSITION=[0.2,0.21,0.95,0.96]
            plot, r, 0d0*r, yrange=[minv/Lstar,maxv/Lstar], /ys, ytitle='L/L!Dsun', xtitle='r/R!Dsun!N (%)', $
              thick=5, xthick=5, ythick=5, color=0, xrange=[double(round(r(0))),100d0], xticks=4, xtickv=100d0*[0.91d0,0.93d0,0.95d0,0.97d0,0.99d0], $
              yticks=5, ytickv=[-2d0,-1d0,0d0,1d0,2d0,3d0], linestyle=2
            oplot, r, 0d0*r+1d0, color=0, linestyle=2, thick=10
            oplot, r, 1d0+0d0*Ltot/Lstar, color=colors[0], thick=10
            oplot, r, Lks/Lstar, color=colors[1], thick=10
            oplot, r, Len/Lstar-(Ltot/Lstar-1d0), color=colors[2], thick=10
            oplot, r, Lke/Lstar, color=colors[3], thick=10
            oplot, r, Lac/Lstar, color=colors[4], thick=10
            DEVICE, /close
            SET_PLOT, 'x'
            print, 'Image ready.'
            SPAWN, 'kghostview '+filename+' &'
        EndIf

    Endif Else Begin
        fksavg = dblarr(N_R)
        fksavg(*) = 0d0
        fenavg = fksavg*0d0
        fkeavg = fksavg*0d0
        facavg = fksavg*0d0

        iks = where(quantities eq 16)
        ien = where(quantities eq 17)
        ike = where(quantities eq 18)
        iac = where(quantities eq 20)

        for ir=0L,N_R-1L do begin
            fksavg(ir) = mean(data(ir,iks,*))
            fenavg(ir) = mean(data(ir,ien,*))
            fkeavg(ir) = mean(data(ir,ike,*))
            facavg(ir) = mean(data(ir,iac,*))
        endfor

        If (keyword_set(botio)) Then Begin
            fksavg(0L:N_R/2-1L) = 0d0
        EndIf

        If (keyword_set(topio)) Then Begin
            fksavg(N_R/2:*) = 0d0
        EndIf

        Lks = 4d0*!dpi*r^2*fksavg
        Lke = 4d0*!dpi*r^2*fkeavg
        Len = 4d0*!dpi*r^2*fenavg
        Lac = 4d0*!dpi*r^2*facavg
        Ltot = Lks+Len+Lke
        minv = min([min(Lks),min(Len),min(Lke),min(Lac)])
        maxv = max([max(Lks),max(Len),max(Lke),max(Lac)])
        r = r /rstar*100d0
        docolor
        loadct, 13
        colors = [255,80,150,140,170]
        If (keyword_set(yrange)) Then Begin
            plot, r, 0d0*r+1d0, yrange=yrange
        Endif Else Begin
            plot, r, 0d0*r+1d0, yrange=[minv/Lstar,maxv/Lstar]
        EndElse
        oplot, r, Ltot/Lstar, color=colors[0]
        oplot, r, Lks/Lstar, color=colors[1]
        oplot, r, Len/Lstar, color=colors[2]
        oplot, r, Lke/Lstar, color=colors[3]
        oplot, r, Lac/Lstar, color=colors[4]

        If (keyword_set(psf)) Then Begin
            xsize=16d0
            ysize=16d0
            !P.MULTI=0
            !P.FONT=1
            !P.CHARSIZE=3
            colors = [30,80,150,225,240]
            Set_Plot, "PS"
            filename = 'flux_bal_timeavg_hyper.eps'
            Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, $
              SET_FONT='Helvetica Italic', /TT_FONT
            !P.POSITION=[0.2,0.21,0.95,0.96]
            plot, r, 0d0*r, ytitle='L/L!Dsun', xtitle='r/R!Dsun!N (%)', $
              thick=5, xthick=5, ythick=5, color=0, xrange=[double(round(r(0))),100d0], xticks=4, xtickv=100d0*[0.91d0,0.93d0,0.95d0,0.97d0,0.99d0], $
              yticks=5, ytickv=[-2d0,-1d0,0d0,1d0,2d0,3d0], linestyle=2
            oplot, r, 0d0*r+1d0, color=0, linestyle=2, thick=10
            oplot, r, Ltot/Lstar, color=colors[0], thick=10
            oplot, r, Lks/Lstar, color=colors[1], thick=10
            oplot, r, Len/Lstar, color=colors[2], thick=10
            
            oplot, r, Lke/Lstar, color=colors[3], thick=10
            oplot, r, Lac/Lstar, color=colors[4], thick=10
            DEVICE, /close
            SET_PLOT, 'x'
            print, 'Image ready.'
            SPAWN, 'kghostview '+filename+' &'
        EndIf

    EndElse

    print, 'averaged over ', double(num_records)*100d0*dtavg/3600d0,' hours'
If (keyword_set(stop)) Then stop
end
