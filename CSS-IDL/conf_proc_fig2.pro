pro conf_proc_fig2, filename, nth, nph,r,ph2,ph1

    uftavg = dcomplexarr(nph,4)
    files =['case1.ss','case2.ss','case3.ss','case4.ss']
    for i=0,3 do begin
       openr,file_unit,files(i),/get_lun
       af = assoc(file_unit,dcomplexarr(nth-8,nph))
       urfft = af(0)
       close,file_unit
       free_lun,file_unit
       uftavg(*,i) = dcomplex(0d0)
       for j=0,nth-9 do begin
          uftavg(*,i) = uftavg(*,i)+urfft(j,*)
       endfor
       uftavg(*,i) = uftavg(*,i)/dcomplex(nth-8)
    endfor
    uftavg = abs(uftavg)
    
    xsize=10d0
    ysize=10d0

    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=1, XSIZE=xsize, YSIZE=ysize

    x0 = 0.14
    y0 = 0.14
    x1 = 0.96
    y1 = 0.96
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]

    maxft=max(uftavg)
    minft=min(uftavg)
    uftavg = uftavg/maxft
    linestyles = lindgen(4)
    dph = r*(ph2-ph1)*!dpi/180d0/double(nph-1)
    rx = r*(ph2-ph1)*!dpi/180d0/(dindgen(nph/2)+1d0)/1d8
    plot, rx, uftavg(0:nph/2,0), xtitle=textoidl('r (Mm)'), ytitle='Power', yrange=[0d0,1d0], $
          linestyle=linestyles[3], xs=1, thick=3, /xlog, xthick=3, ythick=3, charthick=3
    oplot, rx, uftavg(0:nph/2,1), linestyle=linestyles[2], thick=3
    oplot, rx, uftavg(0:nph/2,2), linestyle=linestyles[1], thick=3
    oplot, rx, uftavg(0:nph/2,3), linestyle=linestyles[0], thick=3
    names = ['1','2','3','4']
    linestyles=reverse(linestyles)
    legend,names,LIN = [linestyles],/top,/left,thick=3,charthick=3

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

end
