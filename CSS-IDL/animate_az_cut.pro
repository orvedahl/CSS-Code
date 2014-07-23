@image_polar.pro

pro animate_az_cut, iter, mpgname, qs=qs, dir=dir, xsize=xsize, ysize=ysize, frames=frames, ct=ct

    if (not keyword_set(qs)) then begin
        qs = [0]
    endif

    if (not keyword_set(dir)) then begin
       dir = '/freyr3/augustso/CSSDE/Vapor/big_case_constnu/Checkpoints/'
    endif

    if not keyword_set(frames) then begin
        frames = 30
    endif

    if not keyword_set(ct) then begin
        ct = 'ct_az_avg.bin'
    endif

    k=100000L
    If (iter lt k) Then Begin
        zero_char = '0'
        i=5L
        pretemp = ' '
        stemp = ' '
        While ((iter lt k) and (i gt 0L)) Do Begin
            short_temp = string(i,format='(i1)')
            format = '(i'+strtrim(short_temp,2)+')'
            fmt = strtrim(format,2)
            stemp = string(iter,format=fmt)
            pretemp = strtrim(pretemp+zero_char,2)
            i=i-1L
            k=k/10L
        EndWhile
    EndIf Else Begin
        pretemp = ' '
        stemp = string(iter, format='(i6)')
    EndElse
    
    dir = strtrim(dir,2)+strtrim(pretemp+stemp+'/',2)

    fname = strtrim(dir+'checkpoint_'+strtrim(pretemp+stemp,2)+'_header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = intarr(5)
    readf, file_unit, temp,format='(5i6)'

    nr   = fix(temp(0))
    nth  = fix(temp(1))
    nph  = fix(temp(2))
    nq   = fix(temp(3))
    iter = fix(temp(4))

    print, temp

    temp = dblarr(8)
    readf, file_unit, temp, format='(8e14.7)'
    dt   = temp(0)
    time = temp(1)
    r1   = temp(2)
    r2   = temp(3)
    th1  = temp(4)
    th2  = temp(5)
    ph1  = temp(6)
    ph2  = temp(7)

    print, temp

    close, file_unit
    free_lun, file_unit

    fname = strtrim(dir+'checkpoint_'+strtrim(pretemp+stemp,2),2)
    
    openr, file_unit, fname, /get_lun, /swap_if_big_endian
    volumes = assoc(file_unit,dblarr(nr,nth,nph,nq))
    data = dblarr(nr,nth,nph,nq+2)
    data(*,*,*,0:nq-1) = reform(volumes(*,*,*,*,0))
    close, file_unit
    free_lun, file_unit

    rho_avg = dblarr(nr)
    s_avg = dblarr(nr)

    for ir=0,nr-1 do begin
        rho_avg(ir) = mean(data(ir,*,*,0))
        s_avg(ir) = mean(data(ir,*,*,4))
        data(ir,*,*,nq) = data(ir,*,*,0)-rho_avg(ir)
        data(ir,*,*,nq+1) = data(ir,*,*,4)-s_avg(ir)
    endfor

    nqs = n_elements(qs)

    if not keyword_set(xsize) then begin
        xsize = 2*nr*nqs
    endif

    if not keyword_set(ysize) then begin
        ysize=fix(1.04d0*nth*r2*cos(th1)/(!dpi/2d0-th1)/r1)
        dy = 2*r2*cos(th1)
        dx = r2-r1*sin(th1)
        xsize=fix(2.2*ysize*dx/dy*nqs) ;fix(1.02d0*nqs*(r2-r1*sin(th1))*nr/(r2-r1))
    endif

    print, 'xsize=',xsize, ' ysize=',ysize

    r = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
    costheta = cos(theta)

    mins = dblarr(nqs)
    maxs = dblarr(nqs)
    for iq=0,nqs-1 do begin
        mins(iq) = min(data(*,*,*,qs(iq)))
        maxs(iq) = max(data(*,*,*,qs(iq)))
    endfor

    int_steps = lindgen(nph)+1
    step_text = ash_12_renum(int_steps, ash_12_numbering_starts=0, /quiet)

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
    For ip=0,nph-1 Do Begin
        dx = 1d0/nqs
        x0=0d0
        y0=0.02
        x1=x0+dx
        y1=0.98
        erase
        !P.MULTI=0
        !P.NOERASE=1
        For iq=0,nqs-1 Do Begin
            !P.POSITION=[x0,y0,x1,y1]
            nlevs=10            ;  number of contour levels
            image_polar, data(*,*,ip,qs(iq)), r, costheta, mini=mnv, maxi=mxv, spacing=0.5d0*(r2-r1)/(nr-1d0), /noborder, /nocontour;, /window
            x0=x1
            x1=x0+dx
        EndFor
        image3D = TVRD()
        image = bytarr(3,xsize,ysize)
        image(0,*,*) = ctr[image3D]
        image(1,*,*) = ctg[image3D]
        image(2,*,*) = ctb[image3D]
        write_jpeg, step_text(k)+'.jpg', image,True=1,Quality=100
        print, ' '
        print, strtrim(string((1d2*k)/(1d0*nph-1d0),FORMAT='(G10.4)'),2)+'%'
        print, ' '
        k = k+1L
    EndFor


    bitrate=fix((1d0*xsize*ysize*frames)/1250d0)

    spawn, "mencoder mf://'*.jpg' -mf" + ' w=' + strtrim(string(xsize),2)+':h=' + strtrim(string(ysize),2) + ':fps='+strtrim(string(frames),2)+':type=jpg -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate='+strtrim(string(bitrate),2)+':autoaspect=0 -nosound -o ' + mpgname

    print, "mencoder mf://'*.jpg' -mf" + ' w=' + strtrim(string(xsize),2)+':h=' + strtrim(string(ysize),2) + ':fps='+strtrim(string(frames),2)+':type=jpg -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate='+strtrim(string(bitrate),2)+':autoaspect=0 -nosound -o ' + mpgname

end
