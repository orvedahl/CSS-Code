pro load_fields
    common sohofig, fig1, fig2, fig3, fig4, fen, ke, phi, theta, r

    ;Load newtonian case
    file = '~/cssworkdir/solar_wvh_bio/Shell_Slices/shell_slice_2309000'
    read_shell_slice, file, phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities
    print, 'file loaded'
    print, 'size of data', size(data)
    num_records = n_elements(data(0,0,0,0,*))
    N_R = n_elements(r)
    N_Theta = n_elements(theta)
    N_Phi = n_elements(phi)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)
    fen = dblarr(N_Phi,N_Theta,4)
    ke = dblarr(N_Phi,N_Theta,4)
    ;pick out radial velocities
    fig1 = reform(data(*,*,3,4,0))
    ur = fig1
    ut = reform(data(*,*,3,5,0))
    up = reform(data(*,*,3,6,0))
    rho = reform(data(*,*,3,0,0))
    ke(*,*,0) = 0.5d0*rho*(ur^2+ut^2+up^2)
    T = reform(data(*,*,3,1,0))
    Tavg = mean(T)
    Cp = 3.4d8
    lsun = 3.86d33
    fen(*,*,0) = 4d0*!dpi*r(3)^2*Cp*rho*ur*(T-Tavg)/lsun
    fig2 = reform(data(*,*,N_R-1,4,0))
    ur = fig2
    ut = reform(data(*,*,N_R-1,5,0))
    up = reform(data(*,*,N_R-1,6,0))
    rho = reform(data(*,*,N_R-1,0,0))
    ke(*,*,1) = 0.5d0*rho*(ur^2+ut^2+up^2)
    T = reform(data(*,*,N_R-1,1,0))
    Tavg = mean(T)
    fen(*,*,1) = 4d0*!dpi*r(N_R-1)^2*Cp*rho*ur*(T-Tavg)/lsun

    ;Load sld case
    file = '~/cssworkdir/hyper_test4/Shell_Slices/shell_slice_2549000'
    read_shell_slice, file, phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities
    print, 'file loaded'
    print, 'size of data', size(data)
    num_records = n_elements(data(0,0,0,0,*))
    N_R = n_elements(r)
    N_Theta = n_elements(theta)
    N_Phi = n_elements(phi)
    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

    ;pick out radial velocities
    fig3 = reform(data(*,*,3,4,0))
    fig4 = reform(data(*,*,N_R-1,4,0))

    ur = fig3
    ut = reform(data(*,*,3,5,0))
    up = reform(data(*,*,3,6,0))
    rho = reform(data(*,*,3,0,0))
    ke(*,*,2) = 0.5d0*rho*(ur^2+ut^2+up^2)
    T = reform(data(*,*,3,1,0))
    Tavg = mean(T)
    fen(*,*,2) = 4d0*!dpi*r(3)^2*Cp*rho*ur*(T-Tavg)/lsun
    ur = fig4
    ut = reform(data(*,*,N_R-1,5,0))
    up = reform(data(*,*,N_R-1,6,0))
    rho = reform(data(*,*,N_R-1,0,0))
    ke(*,*,3) = 0.5d0*rho*(ur^2+ut^2+up^2)
    T = reform(data(*,*,N_R-1,1,0))
    Tavg = mean(T)
    fen(*,*,3) = 4d0*!dpi*r(N_R-1)^2*Cp*rho*ur*(T-Tavg)/lsun

end

pro soho2010_fig1
    common sohofig, fig1, fig2, fig3, fig4, fen, ke, phi, theta, r
    space=0.01d0
    rstar = 6.96d10
    xsize = 8d0*2.54d0
    ysize = xsize/4d0
    csize = 1.25
    cs2 = 1.5

    !P.MULTI=0
    !P.FONT=1
    Set_Plot, "PS"
    filename = "soho2010_fig1.eps"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    choose_color, /ash_color, /quiet
    css_lon = phi*180d0/!dpi
    css_lat = reverse(90d0-theta*!radeg)

    limits=[min(css_lat),min(css_lon),max(css_lat),max(css_lon)]

    mnv1 = 0.8*min(fig1)
    mxv1 = max(fig1)

    mnv2 = 0.6*min(fig2)
    mxv2 = 0.9*max(fig2)

    mnv3 = 0.8*min(fig1)
    mxv3 = max(fig1)

    mnv4 = 0.6*min(fig2)
    mxv4 = 0.9*max(fig2)
    
    print, mnv1, mxv1
    print, mnv2, mxv2

    mnv = mnv4
    mxv = mxv4
    draw = transpose(fig4)

 ;Figure out bounds from size of domain
    x0 = 0.01d0
    y0 = 0.08d0
    x1 = 0.24d0
    y1 = 0.99d0
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox

    xyouts,/NORMAL,x0-0.5d0*space,0.01,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,x1-3.5d0*space,0.01,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    mnv = mnv2
    mxv = mxv2
    draw = transpose(fig2)

 ;Figure out bounds from size of domain
    x0 = 0.26
    x1 = 0.49
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox

    xyouts,/NORMAL,x0-0.5d0*space,0.01,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,x1-3.5d0*space,0.01,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    mnv = mnv3
    mxv = mxv3
    draw = transpose(fig3)

 ;Figure out bounds from size of domain
    x0 = 0.51
    x1 = 0.74
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox
    
    xyouts,/NORMAL,x0-0.5d0*space,0.01,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,x1-3.5d0*space,0.01,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    mnv = mnv1
    mxv = mxv1
    draw = transpose(fig1)

 ;Figure out bounds from size of domain
    x0 = 0.76
    x1 = 0.99
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    highlightbox = limits

    display_css_shell_slice, css_lat, css_lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
      orthographic=[0.5d0*(min(css_lon)+max(css_lon)),0.5d0*(min(css_lat)+max(css_lat))], latlon_window=limits, highlightbox=highlightbox
    
    xyouts,/NORMAL,x0-0.5d0*space,0.01,textoidl('-10\circ'),charsize=csize,charthick=1.0,color=0
    xyouts,/NORMAL,x1-3.5d0*space,0.01,textoidl('+10\circ'),charsize=csize,charthick=1.0,color=0

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

end

pro soho2010_fig2
    common sohofig, fig1, fig2, fig3, fig4, fen, ke, phi, theta, r

    fig1m = fig1-mean(fig1)
    fig2m = fig2-mean(fig2)
    fig3m = fig3-mean(fig3)
    fig4m = fig4-mean(fig4)

    ph1 = phi(0)
    ph2 = max(phi)
    nph = n_elements(phi)
    dph = (ph2-ph1)/double(nph-1)
    rx = (ph2-ph1)/(dindgen(nph/2)+1d0)/1d8

    nth = n_elements(theta)
    ;Avoid theta boundaries
    fig1a = fig1m(*,8:nth-8)
    fig2a = fig2m(*,8:nth-8)
    fig3a = fig3m(*,8:nth-8)
    fig4a = fig4m(*,8:nth-8)

    ;Compute spectra
    fig1fft = fft(fig1a,/double)
    fig2fft = fft(fig2a,/double)
    fig3fft = fft(fig3a,/double)
    fig4fft = fft(fig4a,/double)

    ;Average in theta
    fig1fftavg = total(fig1fft,2)
    fig2fftavg = total(fig2fft,2)
    fig3fftavg = total(fig3fft,2)
    fig4fftavg = total(fig4fft,2)

    fig1fftavg = fig1fftavg(0:nph/2L-1L)+reverse(fig1fftavg(nph/2L:*))
    fig2fftavg = fig2fftavg(0:nph/2L-1L)+reverse(fig2fftavg(nph/2L:*))
    fig3fftavg = fig3fftavg(0:nph/2L-1L)+reverse(fig3fftavg(nph/2L:*))
    fig4fftavg = fig4fftavg(0:nph/2L-1L)+reverse(fig4fftavg(nph/2L:*))

    fig1pwr = abs(fig1fftavg)
    fig2pwr = abs(fig2fftavg)
    fig3pwr = abs(fig3fftavg)
    fig4pwr = abs(fig4fftavg)

    fig1pwr = fig1pwr/max(fig1pwr)
    fig2pwr = fig2pwr/max(fig2pwr)
    fig3pwr = fig3pwr/max(fig3pwr)
    fig4pwr = fig4pwr/max(fig4pwr)

    rx1 = rx*r(3)
    nr = n_elements(r)
    rx2 = rx*r(nr-1)

    ;Pdfs
    ntot = double(nth)*double(nph)
    mnv1 = min([min(fig1m),min(fig3m)])
    mxv1 = max([max(fig1m),max(fig3m)])
    mnv2 = min([min(fig2m),min(fig4m)])
    mxv2 = max([max(fig2m),max(fig4m)])
    fig1h = histogram(fig1m,nbins=200,loc=vs1,min=mnv1,max=mxv1)/ntot+1d-6
    fig2h = histogram(fig2m,nbins=200,loc=vs2,min=mnv2,max=mxv2)/ntot+1d-6
    fig3h = histogram(fig3m,nbins=200,loc=vs1,min=mnv1,max=mxv1)/ntot+1d-6
    fig4h = histogram(fig4m,nbins=200,loc=vs2,min=mnv2,max=mxv2)/ntot+1d-6

    fen1 = fen(*,*,0)
    fen2 = fen(*,*,1)
    fen3 = fen(*,*,2)
    fen4 = fen(*,*,3)

    mnv1 = min([min(fen1),min(fen3)])
    mxv1 = max([max(fen1),max(fen3)])
    mnv2 = min([min(fen2),min(fen4)])
    mxv2 = max([max(fen2),max(fen4)])
    fen1h = histogram(fen1,nbins=200,loc=fs1,min=mnv1,max=mxv1)/ntot+1d-6
    fen2h = histogram(fen2,nbins=200,loc=fs2,min=mnv2,max=mxv2)/ntot+1d-6
    fen3h = histogram(fen3,nbins=200,loc=fs1,min=mnv1,max=mxv1)/ntot+1d-6
    fen4h = histogram(fen4,nbins=200,loc=fs2,min=mnv2,max=mxv2)/ntot+1d-6

    ke1 = ke(*,*,0)
    ke2 = ke(*,*,1)
    ke3 = ke(*,*,2)
    ke4 = ke(*,*,3)

    mnv1 = min([min(ke1),min(ke3)])
    mxv1 = max([max(ke1),max(ke3)])
    mnv2 = min([min(ke2),min(ke4)])
    mxv2 = max([max(ke2),max(ke4)])
    ke1h = histogram(ke1,nbins=200,loc=ks1,min=mnv1,max=mxv1)/ntot+1d-6
    ke2h = histogram(ke2,nbins=200,loc=ks2,min=mnv2,max=mxv2)/ntot+1d-6
    ke3h = histogram(ke3,nbins=200,loc=ks1,min=mnv1,max=mxv1)/ntot+1d-6
    ke4h = histogram(ke4,nbins=200,loc=ks2,min=mnv2,max=mxv2)/ntot+1d-6

    vs1 = vs1/1d5

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

    plot, rx2, 0d0*fig4pwr, /xlog, /xs, xtitle='k (Mm)', ytitle='Power', thick=1, xthick=3, ythick=3
    oplot, rx2, fig2pwr, color=125, thick=3
    oplot, rx2, fig4pwr, thick=3

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

    xsize=16d0
    ysize=16d0
    !P.MULTI=0
    !P.FONT=1
    !P.CHARSIZE=2
    Set_Plot, "PS"
    filename = 'soho2010_pdf.eps'
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
    !P.POSITION=[0.15,0.15,0.95,0.95]
    choose_color, /ash_color, /quiet
    
    plot, vs1, 0d0*fig3h, /ylog, xstyle=8, yrange=[1d-5,1d-1], xtitle=textoidl('v_r (km s^{-1})'), ytitle='PDF', $
      thick=1, xthick=3, ythick=3
    oplot, vs1, fig1h, color=125, thick=3
    oplot, vs1, fig3h, thick=3

    axis, xaxis=1, xrange=[-20d0,300d0]
    oplot, vs1, fen1h, color=60, thick=3
    oplot, vs1, fen3h, color=170, thick=3

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'



    ;plot, vs2, fig2h, /ylog, /xs, yrange=[1d-5,1d0]
    ;oplot, vs2, fig4h, color=125

    stop

end
