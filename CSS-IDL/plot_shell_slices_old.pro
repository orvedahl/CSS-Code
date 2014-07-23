
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
@read_shell_slice_old.pro
@display_shell_slice.pro
@get_quantity_name.pro
@ctmaker_cw2.pro
@compact_fd6.pro

pro get_shell_state, state
    drawing_widget = 0L
    quantity_widget = 0L
    radius_widget = 0L
    rec_slider = 0L
    max_slider = 0L
    min_slider = 0L
    options_id = 0L
    ctmaker_id = 0L
    N_Q = 1
    n_f = 0
    x = 500
    y = 500
    window_size = {x:x,y:y}
    x = 600
    y = 600
    base_size = {x:x,y:y}
    n_r = 0
    n_theta = 0
    n_phi = 0
    num_records = 0
    current_radius_selection = 0
    current_quantity_selection = 0
    current_record_selection = 0
    state = {drawing_widget:drawing_widget,quantity_widget:quantity_widget,window_size:window_size, radius_widget:radius_widget, $
             base_size:base_size,options_id:options_id,nq:N_Q,nf:n_f,nr:n_r,ntheta:n_theta, nphi:n_phi, $
             num_records:num_records,current_quantity_selection:current_quantity_selection,$ 
             max_slider:max_slider, min_slider:min_slider,ctmaker_id:ctmaker_id,rec_slider:rec_slider, $
             current_record_selection:current_record_selection, current_radius_selection:current_radius_selection}
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
            'radius_selection' : handle_radius_change, ev
            'rec_slider' : handle_record_event, ev
            'max_slider' : handle_max_event, ev
            'min_slider' : handle_min_event, ev
            'print' : handle_print_event, ev
            'ctmaker' : handle_ctmaker_event, ev
            'vort' : handle_vort_event, ev
            else: ;do nothing
       endcase
    endelse
end

pro handle_vort_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    if(n_elements(data) le 1) then return
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    state.current_quantity_selection=16
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    wr = dblarr(state.ntheta,state.nphi,state.nr)
    sint = sin(theta)
    up = reform(data(*,*,*,5,state.current_record_selection))
    ut = reform(data(*,*,*,4,state.current_record_selection))
    dpi = 1d0/(max(phi)-min(phi))
    dti = 1d0/(max(theta)-min(theta))
    for ir=0,state.nr-1 do begin
        for ip=0,state.nphi-1 do begin
            wr(*,ip,ir) = compact_fd6(dti,sint*up(*,ip,ir),0,0,1,0)/r(ir)/sint
        endfor
        for it=0,state.ntheta-1 do begin
            wr(it,*,ir) = wr(it,*,ir)-compact_fd6(dpi,ut(it,*,ir),0,0,1,0)/r(ir)/sint
        endfor
    endfor
    max_value = max(wr(*,*,state.current_radius_selection))
    min_value = min(wr(*,*,state.current_radius_selection))
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_record_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    if(n_elements(data) le 1) then return
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    rec_value=0
    WIDGET_CONTROL,state.rec_slider, GET_VALUE=rec_value
    state.current_record_selection = rec_value
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_ctmaker_event, ev
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    red = reform(ev.value.table[1,*])
    green = reform(ev.value.table[2,*])
    blue = reform(ev.value.table[3,*])
    tvlct, red, green, blue
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_max_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    if(n_elements(data) le 1) then return
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    WIDGET_CONTROL, ev.id, GET_VALUE=max_value
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_min_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    if(n_elements(data) le 1) then return
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    WIDGET_CONTROL, ev.id, GET_VALUE=min_value
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_quantity_change, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    ;Handle the event that there are no items in the combobox
    if((WIDGET_INFO(state.quantity_widget,/COMBOBOX_NUMBER) eq 0) or (ev.index eq 0) ) then return
    state.current_quantity_selection = ev.index-1
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    max_value = max(data(*,*,state.current_radius_selection,state.current_quantity_selection))
    min_value = min(data(*,*,state.current_radius_selection,state.current_quantity_selection))
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_radius_change, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    ;Handle the event that there are no items in the combobox
    if((WIDGET_INFO(state.radius_widget,/COMBOBOX_NUMBER) eq 0) or (ev.index eq 0) ) then return
    state.current_radius_selection = ev.index-1
    WIDGET_CONTROL, ev.top, SET_UVALUE=state
    If (state.current_quantity_selection eq 16) Then Begin
        max_value = max(wr(*,*,state.current_radius_selection))
        min_value = min(wr(*,*,state.current_radius_selection))
    Endif Else Begin
        max_value = max(data(*,*,state.current_radius_selection,state.current_quantity_selection,state.current_record_selection))
        min_value = min(data(*,*,state.current_radius_selection,state.current_quantity_selection,state.current_record_selection))
    EndElse
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro handle_drawing_event,ev, ir, iq, irec, window_id
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    WIDGET_CONTROL,ev.top,GET_UVALUE=state
    if(n_elements(data) eq 1 ) then return
    WSET, window_id

    lon = phi*180d0/!dpi
    lat = reverse(90d0-theta*!radeg)
    
    limits=[min(lat),min(lon),max(lat),max(lon)]
    
    mnv = min_value ;min(data(*,*,ir,iq,irec))
    mxv = max_value ;max(data(*,*,ir,iq,irec))
    !P.NOERASE=0
    !P.MULTI=0
    !P.POSITION=[0.05,0.05,0.95,0.95]
    If (iq eq 16) Then Begin
        draw = Reverse(Transpose(Reform(wr(*,*,ir))),2)
    Endif Else Begin
        draw = Reverse(Transpose(Reform(data(*,*,ir,iq,irec))),2)
    Endelse
    If (mnv ne mxv) Then Begin
        display_shell_slice, lat, lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
          orthographic=[0.5d0*(min(lon)+max(lon)),0.5d0*(min(lat)+max(lat))], latlon_window=limits
    Endif
    ;stop
end

pro handle_print_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
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
    xsize = 6.0*2.54
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
    command = 'kghostview ' + psfile
    SPAWN, command
end

pro handle_resize_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    WIDGET_CONTROL, ev.top, GET_UVALUE=state
    WIDGET_CONTROL, state.ctmaker_id, GET_VALUE=ct_value
    WIDGET_CONTROL, ev.id, TLB_GET_SIZE=base_size
    WIDGET_CONTROL, ev.id, /DESTROY
    ;Create the parent widget.
    base = WIDGET_BASE(/ROW,TITLE = 'Shell Slice Viewer',/TLB_SIZE_EVENTS)
    ;Create two columns, one for drawing, and one for buttons and options
    drawing_column = WIDGET_BASE(base,/COLUMN)
    options_column = WIDGET_BASE(base,/COLUMN)
    ;Create a combobox that displays the current quantity
    quantity_widget = WIDGET_COMBOBOX(drawing_column,/DYNAMIC_RESIZE,UVALUE='quantity_selection',VALUE='Select a Quantity')
    radius_widget = WIDGET_COMBOBOX(drawing_column,/DYNAMIC_RESIZE,UVALUE='radius_selection',VALUE='Select a Radius')
    WIDGET_CONTROL, quantity_widget, /UPDATE
    ;A drawing widget attached to the drawing column
    
    dx = base_size(0) - state.base_size.x
    dy = base_size(1) - state.base_size.y
    state.window_size.x = state.window_size.x + dx
    state.window_size.y = state.window_size.y + dy


    drawing_widget = WIDGET_DRAW(drawing_column, XSIZE=state.window_size.x, YSIZE=state.window_size.y, FRAME = 2,UVALUE='drawing')
    ;A button to print the current plot
    print_button = WIDGET_BUTTON(options_column,VALUE='Print',UVALUE='print')
    ;Create a button widget to load files.
    load_button = WIDGET_BUTTON(options_column,VALUE='Load File(s)',UVALUE='load_file')
    ;Create a button widget to compute vorticity.
    vort_button = WIDGET_BUTTON(options_column,VALUE='Compute w_r',UVALUE='vort')
    ;Create a button widget to exit the IDL program
    exit_button = WIDGET_BUTTON(options_column,VALUE='Exit',UVALUE='exit')
    ctmaker_widget = ctmaker(BASE_WIDGET=options_column,UVALUE='ctmaker')
    rec_slider = WIDGET_SLIDER(options_column,MAXIMUM=state.num_records-1,MINIMUM=0,VALUE=state.current_record_selection,UVALUE='rec_slider')
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
    state.radius_widget = radius_widget
    state.rec_slider = rec_slider
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
    WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=state.current_quantity_selection
    if state.nr gt 1 then begin
       for i=0,state.nr-1 do begin
           WIDGET_CONTROL,state.radius_widget,COMBOBOX_ADDITEM=string(r(i))
           WIDGET_CONTROL,state.radius_widget,SET_COMBOBOX_SELECT=i+1
       endfor
    endif
    WIDGET_CONTROL, state.radius_widget, SET_COMBOBOX_SELECT=state.current_radius_selection
    WIDGET_CONTROL, state.ctmaker_id, SET_VALUE=ct_value
    state.base_size.x = base_size(0)
    state.base_size.y = base_size(1)
    WIDGET_CONTROL,base,SET_UVALUE=state
end

pro handle_file_load_event, ev
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    WIDGET_CONTROL,ev.top,GET_UVALUE=state

    ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr3/augustso/CSSDE/'
    file = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select Shell Slice File to Read')
 
    ;Determine the number of files.
    n_files = n_elements(file)
    state.nf = n_files

    print, 'loading file:',file
    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(file eq '') then return
     endif

    ;Remove the current list of names in the combobox.
    if(state.nq ge 1) then begin
        for i=state.nq,1,-1 do begin
            WIDGET_CONTROL, state.quantity_widget, COMBOBOX_DELETEITEM=i
            WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=i-1
        endfor
    endif

    ;Remove the current list of names in the combobox.
    if(state.nr ge 1) then begin
        for i=state.nr,1,-1 do begin
            WIDGET_CONTROL, state.radius_widget, COMBOBOX_DELETEITEM=i
            WIDGET_CONTROL, state.radius_widget, SET_COMBOBOX_SELECT=i-1
        endfor
    endif

    ;We assume the files read in are from the same run.
    read_shell_slice_old, file, phi=phi, theta = theta, radii = r, data = data, QUANTITIES = quantities
    print, 'file loaded'
    print, 'size of data', size(data)
    ;phi = 180.0d0*phi/!PI
    ;theta = 180.0d0*theta/!PI-90.0d0    

    num_records = n_elements(data(0,0,0,0,*))
    state.num_records = num_records

    ;Determine the number of radial points.
    N_R = n_elements(r)
    state.nr = N_R

    ;Determine the number of meridional points.
    N_Theta = n_elements(theta)
    state.ntheta = N_Theta

    N_Phi = n_elements(phi)
    state.nphi = N_Phi

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

    ;Default to the first quantity for initialization
    WIDGET_CONTROL, state.quantity_widget, SET_COMBOBOX_SELECT=1

    ;Create a string array of radii names.
    ;Add the current name list to the combobox.
    WIDGET_CONTROL, state.radius_widget, GET_VALUE=temp
    for i=0,N_R-1 do begin
        WIDGET_CONTROL,state.radius_widget,COMBOBOX_ADDITEM=string(r(i))
        WIDGET_CONTROL, state.radius_widget, SET_COMBOBOX_SELECT=i+1
    endfor

    ;Default to the first quantity for initialization
    WIDGET_CONTROL, state.radius_widget, SET_COMBOBOX_SELECT=1

    nlevs=10                    ;  number of contour levels
    min_value = min(data(*,*,0))
    max_value = max(data(*,*,0))
    Print, 'min=',min_value,' max=',max_value
 ;levs=get_contour_levels(data[*,*,0],nlevs,mn=min_value,mx=max_value,pct=0.04)
    If (max_value eq min_value) Then Begin
        max_value = 1.1d0*max_value
        min_value = 0.9d0*max_value
    Endif
    WIDGET_CONTROL, state.max_slider, SET_VALUE=[max_value,min_value,max_value]
    WIDGET_CONTROL, state.min_slider, SET_VALUE=[min_value,min_value,max_value]
    WIDGET_CONTROL, state.rec_slider, SET_SLIDER_MAX=num_records-1, SET_SLIDER_MIN=0
    state.current_quantity_selection=0
    state.current_record_selection=0
    state.current_radius_selection=N_R-1
    WIDGET_CONTROL,ev.top,SET_UVALUE = state
    handle_drawing_event,ev,state.current_radius_selection, state.current_quantity_selection,state.current_record_selection, state.drawing_widget
end

pro plot_shell_slices_old
    common shell_slice_block, data, r, phi, theta, file, quantity_names, min_value, max_value, quantities, wr
    data=0
    r=0
    theta=0
    phi=0
    num_records=10
    file=''
    ;for now use docolor, eventually replace with customizable color table
    docolor
    ;Create the parent widget.
    base = WIDGET_BASE(/ROW,TITLE = 'Shell Slice Viewer',/TLB_SIZE_EVENTS)

    get_shell_state, state
    ;Create two columns, one for drawing, and one for buttons and options
    drawing_column = WIDGET_BASE(base,/COLUMN)
    options_column = WIDGET_BASE(base,/COLUMN)
    ;Create a combobox that displays the current quantity
    quantity_widget = WIDGET_COMBOBOX(drawing_column,/DYNAMIC_RESIZE,UVALUE='quantity_selection',VALUE='Select a Quantity')
    radius_widget = WIDGET_COMBOBOX(drawing_column,/DYNAMIC_RESIZE,UVALUE='radius_selection',VALUE='Select a Radius')
    ;A drawing widget attached to the drawing column
    drawing_widget = WIDGET_DRAW(drawing_column, XSIZE=state.window_size.x, YSIZE=state.window_size.y, FRAME = 2,UVALUE='drawing')
    ;A button to print the current plot
    print_button = WIDGET_BUTTON(options_column,VALUE='Print',UVALUE='print')
    ;Create a button widget to load files.
    load_button = WIDGET_BUTTON(options_column,VALUE='Load File(s)',UVALUE='load_file')
    ;Create a button widget to compute vorticity.
    vort_button = WIDGET_BUTTON(options_column,VALUE='Compute w_r',UVALUE='vort')
    ;Create a button widget to exit the IDL program.
    exit_button = WIDGET_BUTTON(options_column,VALUE='Exit',UVALUE='exit')
    ;Incorporate the colortable maker
    ctmaker_widget = ctmaker(BASE_WIDGET=options_column,UVALUE='ctmaker')
    record_slider = WIDGET_SLIDER(options_column,MINIMUM=0,MAXIMUM=num_records,VALUE=0,UVALUE='rec_slider')
    ;Min, max sliders
    max_slider = CW_FSLIDER(options_column,VALUE=0.,FORMAT='(G13.6)',UVALUE='max_slider')
    min_slider = CW_FSLIDER(options_column,VALUE=0.,FORMAT='(G13.6)',UVALUE='min_slider')
    
    ;Create widget event controller
    WIDGET_CONTROL,base,/REALIZE
    ;Get drawing_widget id
    WIDGET_CONTROL, drawing_widget, GET_VALUE=drawing_widget_id

    ;Save widget ids
    state.drawing_widget = drawing_widget_id
    state.quantity_widget = quantity_widget
    state.radius_widget = radius_widget
    state.options_id = options_column
    state.rec_slider = record_slider
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
