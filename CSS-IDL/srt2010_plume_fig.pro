function average_quantity_over_records, az_avgs, n_r, n_theta, iq, num_records
    result = dblarr(n_theta, n_r)
    result[*,*] = 0d0
    for i=0,num_records-1 do begin
        result(*,*) = result(*,*) + az_avgs(*,*,iq,i)
    endfor
    result = result/double(num_records)
    return, result
end

function average_over_files, files, n_files, n_r, n_theta, n_q, num_records
    result = dblarr(n_theta,n_r,n_q)
    print, 'Averaging'
    result[*,*,*] = 0d0
    if(n_files gt 1) then begin
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


pro srt2010_plume_fig, filename, omega0, nr, dir=dir
    obar = dblarr(nr,4)
    for i=0,3 do begin
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

       nlevs=10                   ;  number of contour levels
       If ((i ge 0)and(i le 2)) Then Begin
           qloc = where(quantities eq 6) ;get vphi
       Endif Else Begin
           qloc = where(quantities eq 7) ;get vphi
       EndElse
       vphi = data(*,*,qloc)
       omega = dblarr(NTheta,NR)
       conang = dblarr(NTheta,NR)
       omega[*,*] = 0d0
       obar[*,i] = 0d0
       for ir=0,NR-1 do begin
          for it=0,NTheta-1 do begin
             omega[it,ir] = (omega0+vphi[it,ir]/r[ir]/sin(theta[it]))/2.0d0/!PI*1.0d9
          endfor
          obar[ir,i] = Mean(omega[*,ir])
       endfor
   endfor
   obar = obar+30d0
   obar(*,3) = (obar(*,3)-470d0)*50d0/(max(obar(*,3))-min(obar(*,3)))+475d0
   conang = dblarr(nr)
   conang = 470d0*(min(r)/r)^2
   rsun = 6.96d10
   rmid = 0.94d0*rsun
   func = 1d0/(1d0+exp(-50d0*(r-rmid)/(max(r)-min(r))))
   solarang = 470d0-7d3*(abs(1d0-r/rmid)^2)*func-3d0*(abs(r-rmid)/(max(r)-min(r)))^1.5*(1d0-func)

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
    linestyles = lindgen(5)*0
    linestyles = [2,3,0,0,0]
    colors = [0,0,40,190,0]
    thickn = [2,3,3,3,5]
    sunsym = sunsymbol()
    rsun = 6.96d10
    xyouts,/NORMAL,0.1,0.925,'d',charsize=1.5, charthick=1.0,color=0

    plot, r/rsun, solarang, xtitle=textoidl('r/R_')+sunsym, ytitle=textoidl('\Omega (nHz)'), $
          yrange=[min_omega,max_omega], linestyle=linestyles[0], xs=1, thick=2, xthick=3, ythick=3, charthick=3, color=0
    oplot, r/rsun, obar(*,0), linestyle=linestyles[1], thick=thickn[1], color=colors[1]
    oplot, r/rsun, obar(*,1), linestyle=linestyles[2], thick=thickn[2], color=colors[2]
    oplot, r/rsun, obar(*,2), linestyle=linestyles[3], thick=thickn[3], color=colors[3]
    oplot, r/rsun, obar(*,3), linestyle=linestyles[4], thick=thickn[4], color=colors[4]
    ;oplot, r/rsun, obar(*,4), linestyle=linestyles[4], thick=5, color=colors[4]
    names = ['Sun','1','2','3','4']
    legend,names,COLORS=[colors], LIN=[linestyles],/bottom,/left,charthick=2,thick=thickn

    DEVICE, /close
    SET_PLOT, 'x'
    print, 'Image ready.'
    SPAWN, 'kghostview '+filename+' &'
stop
end
