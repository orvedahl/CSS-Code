@read_scalar.pro
pro plot_scalar, css_case, filename=filename, is=is, ie=ie, log=log, psf=psf, stop=stop, mkuniq=mkuniq, mag=mag, imag=imag

    If (not keyword_set(filename)) Then filename='/home8/ryor5023/Programs/CSS_Code/Runs/'+css_case+'/scalar.dat'

    read_scalar, filename, values=data, iters=iters, mkuniq=mkuniq, mag=mag, imag=imag

    If (not keyword_set(is)) Then is=2L
    If (not keyword_set(ie)) Then ie=n_elements(iters)-1L

    times = data(0,is:ie)
    dt = data(1,is:ie)
    ke = data(2,is:ie)
    cke = data(3,is:ie)
    mcke = data(4,is:ie)
    drke = data(5,is:ie)
    mach = data(6,is:ie)
    If (keyword_set(mag)) Then Begin
       me = data(7,is:ie)
       bmax = data(8,is:ie)
    EndIf

    If (keyword_set(psf)) Then Begin
        xsize=8.5*2.54
        ysize=4.25d0*2.54

        !P.FONT=3

        filename = 'scalar.eps'
        Set_Plot, "PS"
        choose_color, /ash_color, /quiet
        Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize
        !P.CHARSIZE=1
    Endif Else Begin
        window, 1, xsize=1000, ysize=500
        !P.CHARSIZE=1.75
    EndElse
    day = 24d0*3600d0
    docolor
    window, 1
    If (keyword_set(mag)) Then Begin
       yrng = [min([min(cke),min(me)]),max(ke)]
    Endif Else Begin
       yrng = [min(cke),max(ke)]
    EndElse
    plot, times/day, ke, /xs, title='Volume Averaged KE', xtitle='time (days)', ytitle='<KE>', ylog=log, thick=2, yrange=yrng
    oplot, times/day, cke, color=120, thick=2
    oplot, times/day, mcke, color=150, thick=2
    oplot, times/day, drke, color=100, thick=2
    If (keyword_set(mag)) Then Begin
       oplot, times/day, me, color=200, thick=2
       oplot, times/day, me+ke, linestyle=2, thick=2
    Endif

    window, 2
    plot, times/day, mach, /xs, title='Maximum Mach #', xtitle='time (days)', ytitle='<M>', ylog=log,thick=2

    If (keyword_set(mag)) Then Begin
       window, 3
       plot, times/day, bmax, yrange=minmax(bmax), /xs, title='Maximum B Field', xtitle='time (days)', ytitle='B!Dmax!N', ylog=ylog, thick=2
    Endif

    If (keyword_set(psf)) Then Begin
        DEVICE, /close
        SET_PLOT, 'x'
        print, 'Image ready.'
        ;SPAWN, 'kghostview '+filename+' &'
    EndIf

    If (keyword_set(stop)) Then Begin
        stop
    EndIf

end
