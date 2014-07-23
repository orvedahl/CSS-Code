;----------------------------------------------------------------
;
;plot_azimuthal_averages.pro
;
;Kyle Augustson @ University of Colorado at Boulder July 2006
;
;This is the top level program that reads in azimuthal average
;files using a file selection dialog and then either averages
;the files and plots the result.
;
;----------------------------------------------------------------
;Import some necessary procedures.
@read_shell_avg.pro
@get_quantity_name.pro
@UniformDerivative6.pro

pro get_state, state
    drawing_widget = 0L
    quantity_widget = 0L
    options_id = 0L
    N_Q = 1
    n_f = 0
    x = 600
    y = 600
    window_size = {x:x,y:y}
    x = 700
    y = 700
    base_size = {x:x,y:y}
    n_r = 0
    num_records = 0
    perturb = 0
    lumin = 0
    current_quantity_selection = 0
    log10 = 0
    stopnow = 0
    root_dir = 'tmp'
    state = {drawing_widget:drawing_widget,quantity_widget:quantity_widget,window_size:window_size, $
             base_size:base_size,options_id:options_id,nq:N_Q,nf:n_f,nr:n_r,log10:log10,num_records:num_records,$
             current_quantity_selection:current_quantity_selection, perturb:perturb, lumin:lumin, stopnow:stopnow, root_dir:root_dir}
end

pro base_event, ev
    WIDGET_CONTROL, ev.id, GET_UVALUE = uvalue
    case uvalue of
        'exit' : WIDGET_CONTROL, ev.top, /DESTROY
        'animate' : handle_animation, ev
        'load_file' : handle_file_load_event, ev
        'quantity_selection' : handle_quantity_change, ev
        'print' : handle_print_event, ev
        'perturb' : handle_perturb_event, ev
        'luminosity' : handle_luminosity_event, ev
        'diffrot' : handle_diffrot_event,ev
        'logval' : handle_log10_event, ev
        'stopnow' : handle_stop_event, ev
        'average' : handle_average_event, ev
        else:                   ;do nothing
    endcase
end

pro handle_average_event, ev
  common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
  WIDGET_CONTROL, ev.top, GET_UVALUE=state
  nelem = n_elements(data[0,0,*])
  iq = state.current_quantity_selection
  nq = state.nq
  If (iq eq nq) Then Begin
     avg = total(domega,2)/double(nelem)
  Endif Else Begin
     avg = total(reform(data[*,iq,*]),2)/double(nelem)
  EndElse
  rmax = 6.96d10
  plot, r/rmax, avg, title=quantity_names[iq], ytitle=quantity_names[iq], xtitle='r/R', yrange=minmax(avg),/xs  
end


pro handle_stop_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    state.stopnow = 1-state.stopnow
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
end

pro handle_log10_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    state.log10 = 1-state.log10
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
end

pro handle_diffrot_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    nr = state.num_records
    domega = dblarr(n_elements(data(*,0,0)),nr)
    for ir=0,nr-1 do begin
        domega(*,ir) = omega0 + data(*,6,ir)/r
    endfor
    domega = domega/omega0 ;1d9*domega/2d0/!dpi
    nq = state.nq
    state.current_quantity_selection = nq
    quantity_names = [quantity_names,'Omega']
    quantities = [quantities,quantities(nq-1)+1]
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    handle_drawing_event,state.num_records, state.current_quantity_selection, state.drawing_widget, state.perturb, state.lumin, state.log10, state.stopnow
end

pro handle_animation, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    if (n_elements(data) eq 1) then return
    
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    window_id = state.drawing_widget
    nr = state.num_records
    nq = state.nq
    iq = state.current_quantity_selection
    xsize=state.window_size.x
    ysize=state.window_size.y
    iqr = quantities(iq)
    If (iqr le quantities(nq-1)) Then Begin
        If ((state.lumin eq 1) and ((iqr ge 15) and (iqr le 19))) Then Begin
            lum = dblarr(nr,n_elements(r))
            lum(*,*) = 0d0
            minv = 0d0
            maxv = 0d0
            for ir=0L,nr-1L do begin
                lum(ir,*) = 4d0*!dpi*r^2*data(*,iq,ir)
                tempmax = max(lum(ir,*))
                tempmin = min(lum(ir,*))
                if (tempmin lt minv) then minv=tempmin
                if (tempmax gt maxv) then maxv=tempmax
            endfor
            
            for i=0,nr-1 do begin
                WSET, window_id
                plot, r, lum(i,*), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r (cm)', thick=2, yrange=[minv,maxv]
                wait, 0.025
            endfor
        Endif Else Begin
            if ((state.perturb eq 1)) then begin
                minv = 0d0
                maxv = 0d0
                for ir=0L,nr-1L do begin
                    tempmax = max(data(*,iq,ir)-data(*,iq,0))
                    tempmin = min(data(*,iq,ir)-data(*,iq,0))
                    if (tempmin lt minv) then minv=tempmin
                    if (tempmax gt maxv) then maxv=tempmax
                endfor
                for ir=0L,nr-1L do begin
                    WSET, window_id
                    colors = 205*indgen(nr-1)/(nr-2)+50
                    plot, r, reform(data(*,iq,ir)-data(*,iq,0)), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r (cm)',thick=2, yrange=[minv,maxv]
                    wait, 0.025
                endfor
            endif else begin
                for ir=0L,nr-1L do begin
                    WSET, window_id
                    colors = 205*indgen(nr-1)/(nr-2)+50
                    plot, r, reform(data(*,iq,ir)), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r (cm)', yrange=[min(data(*,iq,*)),max(data(*,iq,*))], thick=2
                    wait, 0.025
                endfor
            endelse
        Endelse
    Endif Else Begin
       for ir=0L,nr-1L do begin
          WSET, window_id
          colors = 205*indgen(nr-1)/(nr-2)+50
          plot, r, reform(domega(*,ir))-domega[*,0]*state.perturb, title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r (cm)', yrange=[min(domega),max(domega)]-minmax(domega[*,0]*state.perturb), thick=2
          wait, 0.025
       endfor
    EndElse
end

pro handle_perturb_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    state.perturb = 1-state.perturb
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    handle_drawing_event,state.num_records, state.current_quantity_selection, state.drawing_widget, state.perturb, state.lumin, state.log10, state.stopnow
end

pro handle_luminosity_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    state.lumin = 1-state.lumin
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    handle_drawing_event,state.num_records, state.current_quantity_selection, state.drawing_widget, state.perturb, state.lumin, state.log10, state.stopnow
end

pro handle_quantity_change, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    ;Handle the event that there are no items in the combobox
    if((WIDGET_INFO(state.quantity_widget,/COMBOBOX_NUMBER) eq 0) or (ev.index eq 0) ) then return
    state.current_quantity_selection = ev.index-1
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    handle_drawing_event,state.num_records, state.current_quantity_selection, state.drawing_widget, state.perturb, state.lumin, state.log10, state.stopnow
end

pro handle_drawing_event, nr, iq, window_id, prtb, lumin, log10, stopatend
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    if(n_elements(data) eq 1 ) then return
    WSET, window_id
    !P.MULTI=0
    !P.POSITION=[0.125,0.1,0.975,0.95]
    !P.NOERASE=0
    colors = long(205d0*indgen(nr-1)/double(nr-2)+50)
    rmax = 6.96d10
    iqr = quantities(iq)
    nq = n_elements(data(0,*,0))
    If (stopatend) Then stop
    If (iqr le quantities(nq-1)) Then Begin
        If ((lumin eq 1) and ((iqr ge 15) and (iqr le 20))) Then Begin
            lum = dblarr(nr,n_elements(r))
            lum(*,*) = 0d0
            minv = 0d0
            maxv = 0d0
            for ir=0L,nr-1L do begin
                lum(ir,*) = 4d0*!dpi*r^2*data(*,iq,ir)
                tempmax = max(lum(ir,*))
                tempmin = min(lum(ir,*))
                if (tempmin lt minv) then minv=tempmin
                if (tempmax gt maxv) then maxv=tempmax
            endfor
            plot, r/rmax, lum(0,*), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r/R', yrange=[minv,maxv],/xs
            for i=1,nr-1 do begin
                oplot, r/rmax, lum(i,*), color=colors(i-1)
            endfor
        Endif Else Begin
            If (prtb eq 1) Then Begin
                minv = 0d0
                maxv = 0d0
                for ir=0L,nr-1L do begin
                    tempmax = max(data(*,iq,ir)-data(*,iq,0))
                    tempmin = min(data(*,iq,ir)-data(*,iq,0))
                    if (tempmin lt minv) then minv=tempmin
                    if (tempmax gt maxv) then maxv=tempmax
                endfor
                plot, r/rmax, data(*,iq,nr-1)-data(*,iq,0), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r/R', yrange=[minv,maxv],/xs
                for i=1,nr-1 do begin
                    oplot, r/rmax, data(*,iq,i)-data(*,iq,0), color=colors(i-1)
                endfor
            EndIf Else Begin
                If (log10 eq 1) Then Begin
                    plot, r/rmax, alog10(abs(data(*,iq,0))), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r/R', yrange=[min(alog10(abs(data(*,iq,*)))),max(alog10(abs(data(*,iq,*))))],/xs
                    for i=1,nr-1 do begin
                        oplot, r/rmax, alog10(abs(data(*,iq,i))), color=colors(i-1)
                    endfor
                Endif Else Begin
                    plot, r/rmax, data(*,iq,0), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r/R',/xs,/ys , yrange=[min(data(*,iq,*)),max(data(*,iq,*))]
                    for i=1,nr-1 do begin
                        oplot, r/rmax, data(*,iq,i), color=colors(i-1)
                    endfor
                EndElse
            Endelse
        Endelse
    Endif Else Begin
       temp = domega
       If (prtb eq 1) Then Begin
          for ir=1,nr-1 do begin
             temp[*,ir] = temp[*,ir]-temp[*,0]
          endfor
          temp[*,0] = 0d0
       Endif
       mnmx = minmax(temp)
        plot, r/rmax, temp[*,0], title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r/R', yrange=mnmx,/xs
        for i=1,nr-1 do begin
            oplot, r/rmax, temp[*,i], color=colors(i-1)
        endfor
    EndElse

end

pro handle_print_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    SPAWN, 'pwd', root_directory
    psfile = DIALOG_PICKFILE(PATH=root_directory,FILTER='*.eps',/FIX_FILTER)
    nf = n_elements(psfile)
    if(nf gt 1) then begin
       print, 'Too many PostScript file names'
    endif
    if(nf le 1) then begin
       if(psfile eq '') then return
    endif
    if(n_elements(data) le 1) then return
    
    aspect = float(state.window_size.y)/state.window_size.x
    xsize = 6.0*2.54
    ysize = xsize
    xoffset = (8.5*2.54 - xsize)/2.0
    yoffset = (11.0*2.54 - ysize)/2.0
    print, 'Output size (cm)',xsize, ysize, aspect
    SET_PLOT, 'ps'
    DEVICE, filename = psfile, BITS = 8,/COLOR,/TT_Font,Set_Font='Helvetica' ;,/ENCAPSULATED
    DEVICE, XSIZE = xsize, YSIZE = ysize, XOFFSET=xoffset, YOFFSET=yoffset
    !P.FONT=1
    !P.CHARSIZE=2
    !P.THICK=2
    !X.THICK=4
    !Y.THICK=4
    pos = !P.POSITION
    !P.POSITION=[0.22,0.15,0.98,0.98]
    nr = state.num_records
    iq = state.current_quantity_selection
    colors = 205*indgen(nr-1)/(nr-2)+50
    rmax = 6.96d10
    iqr = quantities[iq]
    nq = n_elements(data[0,*,0])
    If (iqr le quantities(nq-1)) Then Begin
       plot, r/rmax, data(*,iq,0), title=quantity_names(iq), ytitle=quantity_names(iq), xtitle='r/R',/xs,/ys , yrange=[min(data(*,iq,*)),max(data(*,iq,*))]
       for i=1,nr-1 do begin
          oplot, r/rmax, data(*,iq,i), color=colors(i-1)
       endfor
    Endif Else Begin
       temp = domega
       mnmx = minmax(temp)
       plot, r/rmax, temp[*,0], ytitle=textoidl('\Omega/\Omega_0'), xtitle='r/R', yrange=mnmx,/xs
       for i=1,nr-1 do begin
          oplot, r/rmax, temp[*,i], color=colors(i-1)
       endfor
    EndElse

    DEVICE, /close
    SET_PLOT, 'x'
    !P.FONT=0
    !P.CHARSIZE=1
    !P.POSITION=pos
    command = 'evince ' + psfile + ' &'
    SPAWN, command
end

pro handle_file_load_event, ev
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    WIDGET_CONTROL,ev.top,GET_UVALUE=state

    ;Open a file dialog box to select (a) file(s).
    files = DIALOG_PICKFILE(PATH=state.root_dir,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)
    state.nf = n_files

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;Remove the current list of names in the combobox.
    if(state.nq ge 1) then begin
        for i=state.nq,1,-1 do begin
            WIDGET_CONTROL, state.quantity_widget, COMBOBOX_DELETEITEM=i
            WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=i-1
        endfor
    endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius = r, data = values, QUANTITIES = quantities

    num_records = n_elements(values(0,0,*))
    state.num_records = num_records

    ;Determine the number of radial points.
    N_R = n_elements(r)
    state.nr = N_R

    ;Determine the number of quantities.
    N_Q = n_elements(quantities)
    state.nq = N_Q

    data = dblarr(N_R,N_Q,num_records*n_files)
    data(*,0:N_Q-1,0:num_records-1) = values

    if (n_files gt 1) then begin
        for i=1,n_files-1 do begin
            read_shell_avg, files(i), data = values
            num_rec = n_elements(values(0,0,*))
            data(*,0:N_Q-1,num_records:num_records+num_rec-1) = values
            num_records = num_records + num_rec
        endfor
    endif

    state.num_records = num_records

    ;Create a string array of quantity names.
    ;Add the current name list to the combobox.
    WIDGET_CONTROL, state.quantity_widget, GET_VALUE=temp
    quantity_names = strarr(N_Q)
    for i=0,N_Q-1 do begin
        quantity_names(i) = get_quantity_name(quantities(i))
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM=quantity_names(i)
        WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=i+1
    endfor

    ;Default to the first quantity for initialization
    WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=1
    WIDGET_CONTROL,ev.top,SET_UVALUE = state

    state.current_quantity_selection=0
    handle_drawing_event,state.num_records, state.current_quantity_selection, state.drawing_widget, state.perturb, state.lumin, state.log10, state.stopnow   
end

pro plot_shell_avg_gui, omega=omega, root=root
    common common_block, data, r, files, quantity_names, quantities, animated, domega, omega0
    data=0
    domega = 0d0
    If (keyword_set(omega)) Then Begin
        omega0 = omega
    Endif Else Begin
        omega0 = 2.6662237d-6
    Endelse



    r=0
    files=''    
    ;for now use docolor, eventually replace with customizable color table
    docolor
    !P.CHARSIZE=1.5
    ;Create the parent widget.
    base = WIDGET_BASE(/ROW,TITLE = 'Shell Average Viewer',/TLB_SIZE_EVENTS)

    get_state, state

    If (keyword_set(root)) Then Begin
        state.root_dir = root
    Endif Else Begin
        state.root_dir = '/home8/ryor5023/Programs/CSS_Code/Runs/'
    EndElse

    ;Create two columns, one for drawing, and one for buttons and options
    drawing_column = WIDGET_BASE(base,/COLUMN)
    options_column = WIDGET_BASE(base,/COLUMN)
    ;Create a combobox that displays the current quantity
    quantity_widget = WIDGET_COMBOBOX(drawing_column,/DYNAMIC_RESIZE,UVALUE='quantity_selection',VALUE='Select a Quantity')
    ;A drawing widget attached to the drawing column
    drawing_widget = WIDGET_DRAW(drawing_column, XSIZE=state.window_size.x, YSIZE=state.window_size.y, FRAME = 2,UVALUE='drawing')
    ;A button to print the current plot
    print_button = WIDGET_BUTTON(options_column,VALUE='Print',UVALUE='print')
    ;Create a button widget to load files.
    load_button = WIDGET_BUTTON(options_column,VALUE='Load File(s)',UVALUE='load_file')
    ;Create a button widget to animate plots.
    animate_button = WIDGET_BUTTON(options_column,VALUE='Animate',UVALUE='animate')
    ;Create a button widget show perturbation plots.
    perturb_button = WIDGET_BUTTON(options_column,VALUE='Perturb',UVALUE='perturb')
    ;Create a button widget show luminosity plots.
    luminosity_button = WIDGET_BUTTON(options_column,VALUE='Luminosity',UVALUE='luminosity')
    ;Create a button widget show omega plots.
    omega_button = WIDGET_BUTTON(options_column,VALUE='Omega',UVALUE='diffrot')
    ;Create a button widget show log plots.
    log_button = WIDGET_BUTTON(options_column,VALUE='Log10',UVALUE='logval')
    ;Create a button widget show log plots.
    avg_button = WIDGET_BUTTON(options_column,VALUE='Average',UVALUE='average')
    ;Create a button widget show stop allowing data access.
    stop_button = WIDGET_BUTTON(options_column,VALUE='Stop',UVALUE='stopnow')
    ;Create a button widget to exit the IDL program.
    exit_button = WIDGET_BUTTON(options_column,VALUE='Exit',UVALUE='exit')
    ;Create widget event controller
    WIDGET_CONTROL,base,/REALIZE
    ;Get drawing_widget id
    WIDGET_CONTROL, drawing_widget, GET_VALUE=drawing_widget_id

    ;Save widget ids
    state.drawing_widget = drawing_widget_id
    state.quantity_widget = quantity_widget
    state.options_id = options_column
    WSET, drawing_widget_id
    XMANAGER, 'base', base, /NO_BLOCK

    WIDGET_CONTROL,base,TLB_GET_SIZE=base_size

    state.base_size.x = base_size(0)
    state.base_size.y = base_size(1)
    WIDGET_CONTROL,base,SET_UVALUE=state
end
