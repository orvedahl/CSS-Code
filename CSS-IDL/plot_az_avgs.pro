
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
@read_az_avg.pro
@get_quantity_name.pro
@ctmaker_cw2.pro
@TextBox.pro
@streamfunction.pro
@image_polar.pro
@UniformDerivative6.pro
pro get_state, state
    drawing_widget = 0L
    quantity_widget = 0L
    max_slider = 0L
    min_slider = 0L
    options_id = 0L
    ctmaker_id = 0L
    N_Q = 1
    n_f = 0
    x = 500
    y = 500
    window_size = {x:x,y:y}
    x = 500
    y = 500
    base_size = {x:x,y:y}
    n_r = 0
    n_theta = 0
    num_records = 0
    current_quantity_selection = 0
    aspect = 0.0
    twind=0
    do_omega=0
    state = {drawing_widget:drawing_widget,quantity_widget:quantity_widget,window_size:window_size, $
             base_size:base_size,options_id:options_id,nq:N_Q,nf:n_f,nr:n_r,ntheta:n_theta, $
             num_records:num_records,current_quantity_selection:current_quantity_selection,$ 
             max_slider:max_slider, min_slider:min_slider,ctmaker_id:ctmaker_id,aspect:aspect,twind:twind,do_omega:do_omega}
end

function average_quantity_over_records, az_avgs, n_r, n_theta, iq, num_records, do_wind=do_wind,r=r,theta=theta
    If (keyword_set(do_wind)) Then Begin
        result = dblarr(n_theta, n_r, 2)
        dr = r(1)-r(0)
        dt = theta(1)-theta(0)
        for i=0,num_records-1 do begin
            vphi = reform(az_avgs(*,*,5,i))
            S = reform(az_avgs(*,*,2,i))
            diff = result
            diff(*,*,*) =0d0
            for ir=0,n_r-1 do begin
                diff(*,ir,0) = -sin(theta)*ud6(dt,vphi(*,ir),1)/r(ir)
                diff(*,ir,1) = ud6(dt,S(*,ir),1)/r(ir)
            endfor
            for it=0,n_theta-1 do begin
                diff(it,*,0) = diff(it,*,0)+cos(theta)*ud6(dr,vphi(it,*),1)
            endfor
            result = result + diff
        endfor
    Endif Else Begin
        result = dblarr(n_theta, n_r)
       ;Average
        for i=0,num_records-1 do begin
            result(*,*) = result(*,*) + az_avgs(*,*,iq,i)
        endfor
    EndElse
    result = result/double(num_records)
    return, result

end

function average_over_files, files, n_files, n_r, n_theta, n_q, num_records,wind=wind,r=r,theta=theta
    result = dblarr(n_theta,n_r,n_q)
    If (keyword_set(wind)) Then Begin 
        wind = dblarr(n_theta,n_r,2)
        wind(*,*,*) = 0d0
    EndIf
    print, 'Averaging'
    if(n_files gt 1) then begin
        ;Average
        print, 'Percent complete:'
        for i=0,n_files-1 do begin
            read_az_avg, files(i),data = data,num_recs=num_records
            for j=0,n_q-1 do begin
                result(*,*,j) = result(*,*,j) + average_quantity_over_records(data, n_r, n_theta, j, num_records)
                If (n_elements(wind) gt 1) Then Begin
                    wind = wind + average_quantity_over_records(data, n_r, n_theta, j, num_records,/do_wind,r=r,theta=theta)
                EndIf
            endfor
            print, (100*(i+1))/n_files
        endfor
        If (n_elements(wind) gt 1) Then Begin
            wind = wind/double(n_files)
        Endif
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

pro base_event, ev
    if(tag_names(ev, /STRUCTURE_NAME) eq 'WIDGET_BASE') then begin
        handle_resize_event, ev
    endif else begin
       WIDGET_CONTROL, ev.id, GET_UVALUE = uvalue
       case uvalue of
            'exit' : WIDGET_CONTROL, ev.top, /DESTROY
            'load_file' : handle_file_load_event, ev
            'load_ct' : handle_load_color_table_event, ev
            'quantity_selection' : handle_quantity_change, ev
            'max_slider' : handle_max_event, ev
            'min_slider' : handle_min_event, ev
            'print' : handle_print_event, ev
            'ctmaker' : handle_ctmaker_event, ev
            'textbox' : differential_rotation, ev
            'thermal' : thermal_wind, ev
            'streamlines' : streamlines, ev
            else: ;do nothing
       endcase
    endelse
end

pro streamlines, event
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    WIDGET_CONTROL, event.top, GET_UVALUE=state
    iqr = where(quantities eq 1)
    iqvr = where(quantities eq 5)
    Ntheta = n_elements(costheta)
    Nr = n_elements(r)
    rhovr = dblarr(Ntheta, Nr)
    rhovr(*,*) = data(*,*,iqvr)
    ;rhovr(*,*) = data(*,*,iqr)*data(*,*,iqvr)
    iqvt = where(quantities eq 6)
    rhovtheta = dblarr(Ntheta, Nr)
    rhovtheta = data(*,*,iqvt)
    ;rhovtheta = data(*,*,iqr)*data(*,*,iqvt)
    psi = dblarr(Ntheta, Nr)
    psi = streamfunction(rhovr, rhovtheta, r, costheta, /double)

    WSET, state.drawing_widget
    min_value=min(psi)
    max_value=max(psi)
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    state.current_quantity_selection = 31
    WIDGET_CONTROL, event.top, SET_UVALUE=state
    handle_drawing_event,state.nr, state.current_quantity_selection, state.drawing_widget,state.do_omega,state.twind,state.nq
end

pro handle_ctmaker_event, ev
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    red = reform(ev.value.table[1,*])
    green = reform(ev.value.table[2,*])
    blue = reform(ev.value.table[3,*])
    tvlct, red, green, blue
    handle_drawing_event,state.nr, state.current_quantity_selection, state.drawing_widget,state.do_omega,state.twind,state.nq
end

pro handle_max_event, ev
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    if(n_elements(data) le 1) then return
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    WIDGET_CONTROL, ev.id, GET_VALUE=max_value
    handle_drawing_event,state.nr, state.current_quantity_selection, state.drawing_widget,state.do_omega,state.twind,state.nq
end

pro handle_min_event, ev
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    if(n_elements(data) le 1) then return
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    WIDGET_CONTROL, ev.id, GET_VALUE=min_value
    handle_drawing_event,state.nr, state.current_quantity_selection, state.drawing_widget,state.do_omega,state.twind,state.nq
end

pro handle_quantity_change, ev
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    ;Handle the event that there are no items in the combobox
    if((WIDGET_INFO(state.quantity_widget,/COMBOBOX_NUMBER) eq 0) or (ev.index eq 0) ) then return
    state.current_quantity_selection = ev.index-1
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    iq = state.current_quantity_selection
    case (iq) of 
        29: begin
            max_value = max(wind(*,*,0))
            min_value = min(wind(*,*,0))
            end
        30: begin
            max_value = max(wind(*,*,1))
            min_value = min(wind(*,*,1))
        end
        31: begin
            max_value = max(psi(*,*))
            min_value = min(psi(*,*))
            end
        Else: Begin    
            max_value = max(data(*,*,iq))
            min_value = min(data(*,*,iq))
            End
    EndCase
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    handle_drawing_event,state.nr, state.current_quantity_selection, state.drawing_widget,state.do_omega,state.twind,state.nq
end

pro handle_drawing_event, nr, iq, window_id,domg,dowind, nq
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    if(n_elements(data) eq 1 ) then return
    WSET, window_id
    cbplace = [0.85,0.4,0.88,0.6]
    nlevs=10                    ;  number of contour levels
    if(domg eq 1) then begin
        image_polar, transpose(reform(omega)), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(nr-1.0d0), /nocontour
        colourbar2,cbplace,mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)
    endif else begin
        case (iq) of
            29: begin 
                image_polar, transpose(reform(wind(*,*,0))), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(nr-1.0d0), /nocontour
                colourbar2,cbplace,mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)
                end
            30: begin 
                  image_polar, transpose(reform(wind(*,*,1))), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(nr-1.0d0), /nocontour
                  colourbar2,cbplace,mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)
                end
            31: begin 
                  image_polar, transpose(reform(psi)), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(nr-1.0d0), /nocontour
                  colourbar2,cbplace,mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)
                end
            else: begin
                  image_polar, transpose(reform(data(*,*,iq))), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(nr-1.0d0), /nocontour
                  colourbar2,cbplace,mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)
                  end
        endcase
    EndElse
end

pro handle_print_event, ev
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    SPAWN, 'pwd', root_directory
    psfile = DIALOG_PICKFILE(PATH=root_directory,FILTER='*.ps',/FIX_FILTER)
    nf = n_elements(psfile)
    if(nf gt 1) then begin
       print, 'Too many PostScript file names'
    endif
    if(nf le 1) then begin
       if(psfile eq '') then return
    endif
    if(n_elements(data) le 1) then return
    
    aspect = float(state.window_size.y)/state.window_size.x
    xsize = 5.5*2.54
    ysize = xsize*aspect
    xoffset = (8.5*2.54 - xsize)/2.0
    yoffset = (11.0*2.54 - ysize)/2.0
    print, 'Output size (cm)',xsize, ysize, aspect
    SET_PLOT, 'ps'
    DEVICE, filename = psfile, BITS = 8,/COLOR,/HELVETICA
    DEVICE, XSIZE = xsize, YSIZE = ysize, XOFFSET=xoffset, YOFFSET=yoffset
    image_polar, data(*,*,state.current_quantity_selection), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(2.0d0*(state.nr-1.0d0))
    colourbar2,[0.1,0.4,0.13,0.6],mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value),plttype='ps',text_char_size=0.75

    DEVICE, /close
    SET_PLOT, 'x'
    command = 'ggv ' + psfile
    SPAWN, command
end

pro thermal_wind, event
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi

    WIDGET_CONTROL, event.top, GET_UVALUE=state
    WSET, state.drawing_widget
    state.twind = 2-state.twind
    ;nlevs = 10
    ;levs=get_contour_levels(omega,nlevs,mn=min_value,mx=max_value,pct=0.04)
    min_value=min(wind)
    max_value=max(wind)
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    image_polar, transpose(reform(wind(*,*,0))), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(2.0d0*(state.nr-1.0d0))
    colourbar2,[0.1,0.4,0.13,0.6],mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)
    WIDGET_CONTROL, event.top, SET_UVALUE=state

end

pro differential_rotation, event
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi

    WIDGET_CONTROL, event.top, GET_UVALUE=state
    qloc = where(quantities eq 7) ;get vphi
    ;qloc = where(quantities eq 6) ;get vphi
    vphi = data(*,*,qloc)
    omega0 = 0D
    ty = size(omega0, /TYPE)
    state.do_omega = 1-state.do_omega
    omega0 = TextBox(Title='Enter \Omega_{0}', Group_Leader=event.top, $
    Label='\Omega_{0}: ', Cancel=cancelled, XSize=200, Value='2.66E-6')

    if not cancelled then begin
        ;stop
        omega0 = fix(omega0, type=ty)
        NR = state.nr
        NTheta = state.ntheta
        omega = dblarr(NTheta,NR)
        theta = acos(costheta)
        for ir=0,NR-1 do begin
            for it=0,NTheta-1 do begin
                omega(it,ir) = (omega0+vphi(it,ir)/r(ir)/sin(theta(it)))/2.0d0/!PI*1.0d9
            endfor
        endfor

        WSET, state.drawing_widget
        ;nlevs = 10
        ;levs=get_contour_levels(omega,nlevs,mn=min_value,mx=max_value,pct=0.04)
        min_value=min(omega)
        max_value=max(omega)
        WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
        WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
        image_polar, transpose(reform(omega)), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(2.0d0*(state.nr-1.0d0))
        colourbar2,[0.1,0.4,0.13,0.6],mn=min_value,mx=max_value,textmx=string(max_value),textmn=string(min_value)

        ;Create a plot of omega at slices
        ;degree angles (from equator = 0)
        max_omega=max(omega(NTheta/12:NTheta/2,*))
        min_omega=min(omega(NTheta/12:NTheta/2,*))
        window, 0,xsize=500,ysize=500
        linestyles = lindgen(4)+1
        plot, r, omega(NTheta/2,*), title='Rotation Profiles', xtitle='r (cm)', ytitle='Omega (nHz)', color=255, yrange=[min_omega,max_omega], linestyle=linestyles[0], xs=1
        oplot, r, (omega(NTheta/3,*)+omega(2*NTheta/3,*))/2., color=230, linestyle=linestyles[1]
        oplot, r, (omega(NTheta/4,*)+omega(3*NTheta/4,*))/2., color=200, linestyle=linestyles[2]
        oplot, r, (omega(NTheta/6,*)+omega(5*NTheta/6,*))/2., color=170, linestyle=linestyles[3]
        oplot, r, omega0*1.0d9*make_array(state.nr,value=1.0d0)/2.0d0/!PI, color=255, linestyle=0
        theta_p = [theta(NTheta/2),theta(NTheta/3),theta(NTheta/4),theta(NTheta/6)]*180d0/!dpi-90d0
        names = strarr(4)
        for i=0,3 do begin
            names(i) = strtrim(string(round(theta_p(i)),'(i3)'),2)
        endfor

        legend,names,LIN = linestyles
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM='Omega'
        WIDGET_CONTROL,state.quantity_widget,SET_COMBOBOX_SELECT=state.nq+1
        WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=state.current_quantity_selection
        
    endif
    WIDGET_CONTROL, event.top, SET_UVALUE=state
end

pro handle_resize_event, ev, init=init
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    WIDGET_CONTROL, ev.top, GET_UVALUE=state

    wsize=0
    If (keyword_set(init)) Then Begin
        WIDGET_CONTROL, ev.top, TLB_GET_SIZE=base_size
        WIDGET_CONTROL, state.ctmaker_id, GET_VALUE=ct_value    
        WIDGET_CONTROL, ev.top, /DESTROY
        state.base_size.x = base_size(0)
        state.base_size.y = base_size(1)
        dy = state.window_size.x/state.aspect-state.window_size.y
        dx = 0.3d0*dy
        print, 'resizing with dy=', dy        
        state.window_size.x = state.window_size.x + dx
        state.window_size.y = min([1100,fix(state.window_size.y + dy)])
        state.base_size.x = state.base_size.x + dx
        state.base_size.y = min([1110,fix(state.base_size.y + dy)])
 ;Create the parent widget.
        base = WIDGET_BASE(/ROW,TITLE = 'Azimuthal Average Viewer',/TLB_SIZE_EVENTS,SCR_XSIZE=state.base_size.x,SCR_YSIZE=state.base_size.y)
    Endif Else Begin
        if (state.aspect eq 0) then begin
            image_polar, transpose(reform(data(*,*,0))), r, costheta, mini=min_value, max=max_value, spacing=(max(r)-min(r))/(state.nr-1.0d0), bound=bound
            bound=bound/max(r)
            aspect = (bound(2)-bound(0))/(bound(3)-bound(1))
            state.aspect=aspect
            nx = state.window_size.x
            if (aspect lt 1d0) then begin
                wsize=[nx,fix(nx/aspect)]
            endif else begin
                wsize=[fix(nx*aspect),nx]
            endelse
        endif
        
        if (n_elements(wsize) gt 1) then begin
            WIDGET_CONTROL, state.ctmaker_id, GET_VALUE=ct_value    
            WIDGET_CONTROL, ev.id, /DESTROY
            dx = wsize(0)-state.window_size.x
            dy = wsize(1)-state.window_size.y
            state.window_size.x = wsize(0)
            state.window_size.y = wsize(1)
            base_size = intarr(2)
            base_size(0) = state.base_size.x + dx
            base_size(1) = state.base_size.y + dy
 ;Create the parent widget.
            base = WIDGET_BASE(/ROW,TITLE = 'Azimuthal Average Viewer',/TLB_SIZE_EVENTS,SCR_XSIZE=base_size(0),SCR_YSIZE=base_size(1))
        endif else begin
            WIDGET_CONTROL, ev.id, TLB_GET_SIZE=base_size
            WIDGET_CONTROL, state.ctmaker_id, GET_VALUE=ct_value    
            WIDGET_CONTROL, ev.id, /DESTROY
            dx = base_size(0) - state.base_size.x
            dy = base_size(1) - state.base_size.y
            state.window_size.x = state.window_size.x + dx
            state.window_size.y = state.window_size.y + dy
 ;Create the parent widget.
            base = WIDGET_BASE(/ROW,TITLE = 'Azimuthal Average Viewer',/TLB_SIZE_EVENTS)
        endelse
    Endelse

    ;Create two columns, one for drawing, and one for buttons and options
    drawing_column = WIDGET_BASE(base,/COLUMN)
    options_column = WIDGET_BASE(base,/COLUMN)
    ;Create a combobox that displays the current quantity
    quantity_widget = WIDGET_COMBOBOX(drawing_column,/DYNAMIC_RESIZE,UVALUE='quantity_selection',VALUE='Select a Quantity')
    WIDGET_CONTROL, quantity_widget, /UPDATE
    ;A drawing widget attached to the drawing column
    
    drawing_widget = WIDGET_DRAW(drawing_column, XSIZE=state.window_size.x, YSIZE=state.window_size.y, FRAME = 2,UVALUE='drawing')
    ;A button to print the current plot
    print_button = WIDGET_BUTTON(options_column,VALUE='Print',UVALUE='print')
    ;Create a button widget to load files.
    load_button = WIDGET_BUTTON(options_column,VALUE='Load File(s)',UVALUE='load_file')

    ;Create a button widget to exit the IDL program
    exit_button = WIDGET_BUTTON(options_column,VALUE='Exit',UVALUE='exit')
    ctmaker_widget = ctmaker(BASE_WIDGET=options_column,UVALUE='ctmaker')

    ;Create a button widget to compute omega.
    omega_button = WIDGET_BUTTON(options_column,VALUE='Omega',UVALUE='textbox')
    therm_button = WIDGET_BUTTON(options_column,VALUE='T. Wind',UVALUE='thermal')
    v_stream_button = WIDGET_BUTTON(options_column, VALUE='V Streamlines', UVALUE='streamlines')
    ;Min, max sliders
    max_slider = CW_FSLIDER(options_column,MAXIMUM=max_value,MINIMUM=min_value,VALUE=max_value,FORMAT='(G13.6)',UVALUE='max_slider')
    min_slider = CW_FSLIDER(options_column,MAXIMUM=max_value,MINIMUM=min_value,VALUE=min_value,FORMAT='(G13.6)',UVALUE='min_slider')

    ;Create widget event controller
    WIDGET_CONTROL,base,/REALIZE

    ;Get drawing_widget id
    WIDGET_CONTROL, drawing_widget, GET_VALUE=drawing_widget_id

    ;Save widget ids
    state.drawing_widget = drawing_widget_id
    state.quantity_widget = quantity_widget
    state.max_slider = max_slider
    state.min_slider = min_slider
    state.ctmaker_id = ctmaker_widget
    WSET, drawing_widget_id
    XMANAGER, 'base', base, /NO_BLOCK

    if state.nq gt 1 then begin
       for i=0,state.nq-1 do begin
           WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM=quantity_names(i)
           WIDGET_CONTROL,state.quantity_widget,SET_COMBOBOX_SELECT=i+1
       endfor
    endif
    If (state.twind eq 1) Then Begin
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM='dOmega/dz'
        WIDGET_CONTROL,state.quantity_widget,SET_COMBOBOX_SELECT=state.nq+1
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM='dS/dtheta'
        WIDGET_CONTROL,state.quantity_widget,SET_COMBOBOX_SELECT=state.nq+2
    EndIf
    WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=state.current_quantity_selection
    WIDGET_CONTROL, state.ctmaker_id, SET_VALUE=ct_value
    state.base_size.x = base_size(0)
    state.base_size.y = base_size(1)
    WIDGET_CONTROL,base,SET_UVALUE=state
    handle_drawing_event,state.nr, state.current_quantity_selection, drawing_widget_id,state.do_omega,state.twind,state.nq
end

pro handle_file_load_event, ev
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    WIDGET_CONTROL,ev.top,GET_UVALUE=state

    ;Open a file dialog box to select (a) file(s).
    root_directory = '/home8/ryor5023/Programs/CSS_Code/Runs/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Azimuthal Average Files to Read')
 
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
    read_az_avg, files(0), theta = theta, radius = r, data = values, QUANTITIES = quantities
    costheta = cos(theta)
    num_records = n_elements(values(0,0,0,*))
    state.num_records = num_records
    tp = 0.5d0*!dpi-theta
    rmin = min(r)
    rmax = max(r)
    bound=[rmin*cos(max(tp)),rmax*sin(min(tp)),rmax,rmax*sin(max(tp))]/rmax
    state.aspect = (bound(2)-bound(0))/(bound(3)-bound(1))
    print, 'aspect=', state.aspect
    ;Determine the number of radial points.
    N_R = n_elements(r)
    state.nr = N_R

    ;Determine the number of meridional points.
    N_Theta = n_elements(theta)
    state.ntheta = N_Theta

    ;Determine the number of quantities.
    N_Q = n_elements(quantities)
    state.nq = N_Q

    ;Create a string array of quantity names.
    ;Add the current name list to the combobox.
    WIDGET_CONTROL, state.quantity_widget, GET_VALUE=temp
    quantity_names = strarr(N_Q)
    for i=0,N_Q-1 do begin
        quantity_names(i) = get_quantity_name(quantities(i))
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM=quantity_names(i)
        WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=i+1
    endfor

    If (state.twind eq 1) Then Begin
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM='dOmega/dz'
        WIDGET_CONTROL,state.quantity_widget,SET_COMBOBOX_SELECT=state.nq+1
        WIDGET_CONTROL,state.quantity_widget,COMBOBOX_ADDITEM='dS/dtheta'
        WIDGET_CONTROL,state.quantity_widget,SET_COMBOBOX_SELECT=state.nq+2
        WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=state.current_quantity_selection
    EndIf

    ;Default to the first quantity for initialization
    WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=1
    WIDGET_CONTROL,ev.top,SET_UVALUE = state
    If (state.twind eq 1) Then Begin
        wind = dblarr(N_Theta,N_R,2)
        omega0 = 0D
        ty = size(omega0, /TYPE)
        omega0 = TextBox(Title='Enter Omega0', Group_Leader=ev.top, $
                         Label='Omega0: ', Cancel=cancelled, XSize=200, Value='2.6E-6')
        Cp = 0D
        ty = size(Cp, /TYPE)
        Cp = TextBox(Title='Enter Cp', Group_Leader=ev.top, $
                         Label='Cp: ', Cancel=cancelled, XSize=200, Value='3.5E8')

        ;Load background file
        file = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select Background File')
        openr, file_unit, file(0), /get_lun
        snr = 0
        grav = dblarr(N_R)
        readf,file_unit,snr,format='(i5)'
        readf,file_unit,grav,format='(e16.9)'
        close, file_unit
        free_lun, file_unit
        data = average_over_files(files,n_files,N_R,N_Theta,N_Q,num_records,wind=wind,r=r,theta=theta)

        wind(*,*,0) = 2d0*omega0*wind(*,*,0)
        for ir=0,N_R-1 do begin
            wind(*,ir,1) = grav(ir)*wind(*,ir,1)/Cp
        endfor
        
    Endif Else Begin
        data = average_over_files(files,n_files,N_R,N_Theta,N_Q,num_records)
    EndElse

    nlevs=10                    ;  number of contour levels
    levs=get_contour_levels(data[*,*,0],nlevs,mn=min_value,mx=max_value,pct=0.04)

    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    state.current_quantity_selection=0
    handle_drawing_event,state.nr, state.current_quantity_selection, state.drawing_widget,state.do_omega,state.twind,state.nq
    handle_resize_event,ev, /init
end

pro plot_az_avgs, thermal_wind=thermal_wind
    common common_block, data, r, costheta, files, quantity_names, min_value, max_value, quantities, omega, wind, psi
    data=0
    r=0
    theta=0
    omega=0        

    files=''
    ;for now use docolor, eventually replace with customizable color table
    docolor
    ;Create the parent widget.
    base = WIDGET_BASE(/ROW,TITLE = 'Azimuthal Average Viewer',/TLB_SIZE_EVENTS)

    get_state, state
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

    ;Create a button widget to exit the IDL program.
    exit_button = WIDGET_BUTTON(options_column,VALUE='Exit',UVALUE='exit')
    ;Incorporate the colortable maker
    ctmaker_widget = ctmaker(BASE_WIDGET=options_column,UVALUE='ctmaker')

    ;Create a button widget to compute omega.
    omega_button = WIDGET_BUTTON(options_column,VALUE='Omega',UVALUE='textbox')
    therm_button = WIDGET_BUTTON(options_column,VALUE='T. Wind',UVALUE='thermal')
    v_stream_button = WIDGET_BUTTON(options_column, VALUE='V Streamlines', UVALUE='streamlines')
    ;Min, max sliders
    max_slider = CW_FSLIDER(options_column,VALUE=0.,FORMAT='(G13.6)',UVALUE='max_slider')
    min_slider = CW_FSLIDER(options_column,VALUE=0.,FORMAT='(G13.6)',UVALUE='min_slider')
    
    ;Create widget event controller
    WIDGET_CONTROL,base,/REALIZE
    ;Get drawing_widget id
    WIDGET_CONTROL, drawing_widget, GET_VALUE=drawing_widget_id

    ;Save widget ids
    If (keyword_set(thermal_wind)) Then Begin
        state.twind=1
    EndIf
    state.drawing_widget = drawing_widget_id
    state.quantity_widget = quantity_widget
    state.options_id = options_column
    state.max_slider = max_slider
    state.min_slider = min_slider
    state.ctmaker_id = ctmaker_widget
    WSET, drawing_widget_id
    XMANAGER, 'base', base, /NO_BLOCK

    WIDGET_CONTROL,base,TLB_GET_SIZE=base_size

    state.base_size.x = base_size(0)
    state.base_size.y = base_size(1)
    WIDGET_CONTROL,base,SET_UVALUE=state
end
