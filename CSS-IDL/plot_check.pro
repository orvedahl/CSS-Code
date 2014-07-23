pro plot_check, filename, nr

    openr, file_unit, filename, /get_lun
    dint = 0L
    data = dblarr(3)
    vmode = dblarr(nr)
    rho = dblarr(nr)
    mfr = dblarr(nr)
    for ir=0,nr-1 do begin
        readf, file_unit, dint, data, format='(i4,3e19.10)'
        vmode(ir) = data(0)
        mfr(ir) = data(1)
        rho(ir) = data(2)
    endfor
    close, file_unit
    free_lun, file_unit

    window, 0, xsize=1200, ysize=500
    !P.CHARSIZE=2
    !P.MULTI=[0,3,1]
    plot, vmode, psym=2
    plot, mfr
    plot, rho
    !P.MULTI=0

end
