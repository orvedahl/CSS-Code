function average_quantity_over_records, az_avgs, n_r, n_theta, iq, num_records
    result = dblarr(n_theta, n_r)
 ;Average                                                                                                                                                         
    for i=0,num_records-1 do begin
        result(*,*) = result(*,*) + az_avgs(*,*,iq,i)
    endfor
    result = result/double(num_records)
    return, result
end

function average_over_files, files, n_files, n_r, n_theta, n_q, num_records
    result = dblarr(n_theta,n_r,n_q)
    print, 'Averaging'
    if(n_files gt 1) then begin
 ;Average                                                                                                                                                        
        print, 'Percent complete:'
        for i=0,n_files-1 do begin
            read_az_avg, files(i),data = data,num_recs=num_records
            for j=0,n_q-1 do begin
                result(*,*,j) = result(*,*,j) + average_quantity_over_records(data, n_r, n_theta, j, num_records)
            endfor
            print, (100*(i+1))/n_files
        endfor
        result = result/double(n_files)
        print, 'Finished averaging and reading'
    endif else begin
        read_az_avg, files(0),data = data
        for j=0,n_q-1 do begin
            result(*,*,j) = average_quantity_over_records(data, n_r, n_theta, j, num_records)
        endfor
    endelse

    return, result
end


pro conf_proc_fig, filename, omega0, nr, dir=dir

    ;obar = dblarr(nr,4)
    ;for i=0,3 do begin
    ;   root_directory = '/freyr3/augustso/CSSDE/'
    ;   files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Azimuthal Average Files to Read')
    ;   n_files = n_elements(files)
    ;
    ;   if(n_files eq 1) then begin
    ;      if(files eq '') then return
    ;   endif

    ;   read_az_avg, files(0), theta = theta, radius = r, data = values, QUANTITIES = quantities
    ;   costheta = cos(theta)
    ;   num_records = n_elements(values(0,0,0,*))
    ;   NR = n_elements(r)
    ;   NTheta = n_elements(theta)
    ;   NQ = n_elements(quantities)
    ;   data = average_over_files(files,n_files,NR,NTheta,NQ,num_records)

    ;   nlevs=10                   ;  number of contour levels
    ;   qloc = where(quantities eq 6) ;get vphi
    ;   vphi = data(*,*,qloc)
    ;   omega = dblarr(NTheta,NR)
    ;   conang = dblarr(NTheta,NR)
       
    ;   for ir=0,NR-1 do begin
    ;      for it=0,NTheta-1 do begin
    ;         omega(it,ir) = (omega0+vphi(it,ir)/r(ir)/sin(theta(it)))/2.0d0/!PI*1.0d9
    ;         obar(ir,i) = obar(ir,i)+omega(it,ir)
    ;      endfor
    ;   endfor
    ;   obar(*,i) = obar(*,i)/double(NTheta)
    ;endfor

    obar = dblarr(nr)
    root_directory = '/freyr3/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Azimuthal Average Files to Read')
    n_files = n_elements(files)
    
    if(n_files eq 1) then begin
        if(files eq '') then return
    endif

    read_az_avg, files(0), theta = theta, radius = r, data = values, QUANTITIES = quantities
    costheta = cos(theta)
    num_records = n_elements(values(0,0,0,*))
    NR = n_elements(r)
    NTheta = n_elements(theta)
    NQ = n_elements(quantities)
    data = average_over_files(files,n_files,NR,NTheta,NQ,num_records)

    nlevs=10 ;  number of contour levels
    qloc = where(quantities eq 6) ;get vphi
    vphi = data(*,*,qloc)
    qloc = where(quantities eq 1) ;get rho
    rho = data(*,*,qloc)
    omega = dblarr(NTheta,NR)
    conang = dblarr(NTheta,NR)
    cabar = dblarr(NR)
    rhobar = dblarr(NR)
    for ir=0,NR-1 do begin
        for it=0,NTheta-1 do begin
            omega(it,ir) = (omega0+vphi(it,ir)/r(ir)/sin(theta(it)))/2.0d0/!PI*1.0d9
            rhobar(ir) = rhobar(ir)+rho(it,ir)
            ;conang(it,ir) = rho(it,NR-1)*r(NR-1)^2*omega0/rho(it,ir)/r(ir)^2/2.0d0/!PI*1.0d9
            ;cabar(ir) = cabar(ir)+conang(it,ir)
            obar(ir) = obar(ir)+omega(it,ir)
        endfor
    endfor
    obar = obar/double(NTheta)
    rhobar = rhobar/double(NTheta)
    mass = 1.989d33 ;g
    mr = dblarr(NR)
    mr(0) = mass
    for ir=1,NR-1 do begin
        mr(ir) = mr(0)+4d0*!dpi*int_tabulated(r(0:ir),r(0:ir)^2*rhobar(0:ir),/DOUBLE)
    endfor
    cabar = r(NR/3)^2*mr(Nr/3)*omega0/r^2/mr/2d0/!dpi*1d9-20d0

    xsize=10d0
    ysize=10d0
    !P.FONT=1
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Helvetica Italic', /TT_FONT
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
 ;Create a plot of omega at slices
 ;degree angles (from equator = 0)                                       
    max_omega=max(obar)
    min_omega=min(obar)
    linestyles = lindgen(4)
    sunsym = sunsymbol()
    r2 = 6.96d10
    plot, r/r2, omega0*1.0d9*make_array(NR,value=1.0d0)/2.0d0/!PI, xtitle=textoidl('r/R_')+sunsym, ytitle=textoidl('\Omega (nHz)'), $
          yrange=[min_omega,max_omega], xrange=[0.91d0,1d0], linestyle=0, xs=1, thick=2, xthick=3, ythick=3, charthick=3, color=0
    oplot, r/r2, cabar, color=200, linestyle=0, thick=3
    oplot, r/r2, obar, color=45, linestyle=0, thick=3

    xyouts,/NORMAL,0.1,0.925,'d',charsize=1.5, charthick=1.0,color=0

    ;plot, r/r2, obar(*,0), xtitle=textoidl('r/R_')+sunsym, ytitle=textoidl('\Omega (nHz)'), $
    ;      yrange=[min_omega,max_omega], linestyle=linestyles[3], xs=1, thick=3, xthick=3, ythick=3, charthick=3
    ;oplot, r/r2, obar(*,1), linestyle=linestyles[2], thick=3
    ;oplot, r/r2, obar(*,2), linestyle=linestyles[1], thick=3
    ;oplot, r/r2, obar(*,3), linestyle=linestyles[0], thick=3
    ;oplot, r/r2, omega0*1.0d9*make_array(NR,value=1.0d0)/2.0d0/!PI, color=0, linestyle=0, thick=2
;    names = ['1','2','3','4']
;    linestyles = reverse(linestyles)
;    legend,names,LIN = [linestyles],/bottom,/left,charthick=2,thick=3

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'
stop
end
