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


pro diff_rot_fig, filename, omega0, dir=dir

    root_directory = '/freyr3/augustso/CSSDE/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Azimuthal Average Files to Read')
    
 ;Determine the number of files.                                                                                                                                     
    n_files = n_elements(files)
    
 ;Handle the event that the user cancels.                                                                                                                            
    if(n_files eq 1) then begin
       if(files eq '') then return
   endif

    ;We assume the files read in are from the same run.
    read_az_avg, files(0), theta = theta, radius = r, data = values, QUANTITIES = quantities
    costheta = cos(theta)
    num_records = n_elements(values(0,0,0,*))
 ;Determine the number of radial points.                                                                       
    NR = n_elements(r)
 ;Determine the number of meridional points.                                                                                                                         
    NTheta = n_elements(theta)
 ;Determine the number of quantities.                                                                                                                                
    NQ = n_elements(quantities)
    
    data = average_over_files(files,n_files,NR,NTheta,NQ,num_records)

    nlevs=10                    ;  number of contour levels
    qloc = where(quantities eq 6) ;get vphi
    vphi = data(*,*,qloc)
    omega = dblarr(NTheta,NR)
    obar = dblarr(NR)
    for ir=0,NR-1 do begin
        for it=0,NTheta-1 do begin
            omega(it,ir) = (omega0+vphi(it,ir)/r(ir)/sin(theta(it)))/2.0d0/!PI*1.0d9
            obar(ir) = obar(ir)+omega(it,ir)
        endfor
    endfor
    obar = obar/double(NTheta)

    xsize=18d0
    ysize=14d0
    choose_color, /rainbow, /quiet
    Set_Plot, "PS"
    Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize

    x0 = 0.1
    y0 = 0.1
    x1 = 0.695
    y1 = 0.9
    !P.NOERASE=1
    !P.MULTI=0
    !P.POSITION=[x0,y0,x1,y1]
    r1 = min(r)
    r2 = max(r)
 ;Create a plot of omega at slices
 ;degree angles (from equator = 0)                                       
    max_omega=max(omega(NTheta/12:NTheta/2,*))
    min_omega=min(omega(NTheta/12:NTheta/2,*))
    linestyles = lindgen(4)+1
    sunsym = sunsymbol()
    plot, r/r2, omega(NTheta/2,*), title='Rotation Profiles', xtitle=textoidl('r/R_')+sunsym, ytitle='Omega (nHz)', color=0, $
      yrange=[min_omega,max_omega], linestyle=linestyles[0], xs=1, thick=3
    oplot, r/r2, (omega(NTheta/3,*)+omega(2*NTheta/3,*))/2., color=30, linestyle=linestyles[1], thick=3
    oplot, r/r2, (omega(NTheta/4,*)+omega(3*NTheta/4,*))/2., color=50, linestyle=linestyles[2], thick=3
    oplot, r/r2, (omega(NTheta/6,*)+omega(5*NTheta/6,*))/2., color=70, linestyle=linestyles[3], thick=3
    oplot, r/r2, obar, color=0, linestyle=0, thick=2
    oplot, r/r2, omega0*1.0d9*make_array(NR,value=1.0d0)/2.0d0/!PI, color=0, linestyle=0, thick=1
    theta_p = [theta(NTheta/2),theta(NTheta/3),theta(NTheta/4),theta(NTheta/6)]*180d0/!dpi-90d0
    theta_p = -theta_p
    names = strarr(5)
    for i=0,3 do begin
        names(i) = strtrim(string(round(theta_p(i)),'(i3)'),2)
    endfor
    names(4) = 'Avg.'
    legend,names,LIN = [linestyles,0]

   !P.NOERASE=1
    !P.MULTI=0
    x0 = 0.725d0
    y0 = 0.05d0
    x1 = 1.1d0
    y1 = 0.95d0

    margin = 0.1d0
    th2 = max(theta)
    bound=[r1*cos(th2)*(1d0-margin),-r2*sin(th2)*(1d0+margin),r2*(1d0+margin),r2*sin(th2)*(1d0+margin)]
    bound=bound/r2
    aspect = (bound(2)-bound(0))/(bound(3)-bound(1))

    dx = x1-x0
    dy = y1-y0
 ;if (aspect lt 1d0) then begin                                                                                                                                      
 ;    x0 = x0+dx*(1d0-aspect)*0.5d0                                                                                                                                  
 ;    x1 = x0+dx*aspect                                                                                                                                              
 ;endif else begin                                                                                                                                                   
 ;    y0 = y0+dy*(1d0-1d0/aspect)*0.5d0                                                                                                                              
 ;    y1 = y0+dy/aspect                                                                                                                                              
 ;endelse                                                                                                                                                            

    !P.POSITION=[x0,y0,x1,y1]
    costheta = cos(theta)
    mnv = min(omega)
    mxv = max(omega)
    image_polar, transpose(omega), r, costheta, mini=mnv, max=mxv, spacing=0.3d0*(r2-r1)/double(nr-1), /nocontour
    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'

end
