pro simpson_int, nr

    onethird = 1d0/3d0
    fout = dblarr(nr)

    r = 2d0*dindgen(nr)/double(nr-1)-1d0
    dr = 2d0/double(nr-1)

    fin = dblarr(nr)
    fin(*) = !dpi*cos(!dpi*r)

    fout(0L) = 0d0
    for ir=1L,nr-1L do begin
        fout(ir) = fout(ir-1)+dr*fin(ir)
    endfor

    fout = fout

    fexc = sin(!dpi*r)

    plot, r, fexc
    oplot, r, fout, psym=2

    stop

end

