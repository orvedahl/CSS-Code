pro load_fields, file, nth, nphi, nt
common cblk, fftcube, ffttot, kphi, radii, theta, phi, quantities
    root_dir = '/albireo/augustso/CSSDE/'
    filename = root_dir + file
    openr, lun, filename, /get_lun
    a = assoc(lun, fltarr(nth,nphi,nt))
    data = a[0]
    close, lun
    free_lun, lun

    dth = 2d0*!dpi/double(nth)/3d0
    theta = dth*dindgen(nth) + !dpi/6d0

    apod = (1d0-cos(!dpi*(theta-!dpi/6d0)/(max(theta)-min(theta)))^64)

    fftcube = dblarr(nth/2,nphi/2, nt)
    for it=0,nt-1 do begin
;        for ith=nth/6-1,5*nth/6-1 do begin
        tdat = reform(data[*,*,it])
        for ip=0,nphi-1 do begin
            tdat[*,ip] = tdat[*,ip]*apod
        endfor
        tmp = fft(tdat,-1,/double)
        fftcube[0,0,it] = fftcube[0,0,it]+ tmp[0,0]
        fftcube[1:nth/2-1,1:nphi/2-1,it] = fftcube[1:nth/2-1,1:nphi/2-1,it] + (tmp[1:nth/2-1,1:nphi/2-1] + reverse(reverse(tmp[nth/2:*,nphi/2:*],1),2))


        ;    fftcube[0,it] = fftcube[0,it]+ tmp[0,0]
        ;    fftcube[1:nphi/2-1,it] = fftcube[1:nphi/2-1,it] + (tmp[1:nphi/2-1] + reverse(tmp[nphi/2:*]))
;        endfor
    endfor

    ffttot = total(total(fftcube,3),1)/double(nt*nth/2d0)
    kphi = dindgen(nphi/2)*2d0*!dpi/6.9d10/(!dpi/4d0)

end

pro spectral_analysis
common cblk, fftcube, ffttot, kphi, radii, theta, phi, quantities

    set_plot, 'ps'
    Device, filename='spectra_sld_low.eps', /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=15, YSIZE=15, SET_FONT='Helvetica Italic', /TT_FONT
    choose_color, /rain
    !P.FONT=1
    !P.CHARSIZE=1.5
    !P.THICK=4
    !X.THICK=4
    !Y.THICK=4
    !P.POSITION=[0.12,0.12,0.95,0.92]
    plot, kphi*1d8, smooth(abs(ffttot),5), /xlog, /ylog, xrange=[1d-2,1d2], yrange=[1d-3,1d4],  xtitle='k (Mm!U-1!N)', ytitle='<KE>', title='Convection Spectra'
    oplot, 5d0*kphi*1d8, 5d2*(5d8*kphi)^(-3.5d0), thick=4, color=50
    oplot, 5d0*kphi*1d8, 6d1*(5d8*kphi)^(-5d0/3d0), thick=4, color=200

    xyouts, /data, 2d1, 0.8d0, 'k!U-5/3'
    xyouts, /data, 2d1, 3d-2, 'k!U-7/2'

    device,/close
    spawn, 'evince spectra_sld_low.eps &'

stop
end
