;----------------------------------------------------------------
;
;plot_azimuthal_averages.pro
;
;Kyle Augustson @ University of Colorado at Boulder April 2009
;
;This is the top level program that reads in azimuthal average
;files using a file selection dialog and then either averages
;the files and plots the result.
;
;----------------------------------------------------------------
;Import some necessary procedures.
@read_flux_bal.pro
@legend.pro

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
    current_quantity_selection = 0
    state = {drawing_widget:drawing_widget,window_size:window_size,base_size:base_size,$
             options_id:options_id,nq:N_Q,nf:n_f,nr:n_r, num_records:num_records}
end

pro base_event, ev
    if(tag_names(ev, /STRUCTURE_NAME) eq 'WIDGET_BASE') then begin
        handle_resize_event, ev
    endif else begin
       WIDGET_CONTROL, ev.id, GET_UVALUE = uvalue
       case uvalue of
            'exit' : WIDGET_CONTROL, ev.top, /DESTROY
            'load_file' : handle_file_load_event, ev
            'print' : handle_print_event, ev
            else: ;do nothing
       endcase
    endelse
end

pro handle_drawing_event, window_id
    common common_block, data, r, files, quantity_names, times, avgs
    if(n_elements(data) eq 1 ) then return
    WSET, window_id

    nq = n_elements(quantity_names)
    nr = n_elements(r)

    colors = 205*indgen(nq)/(nq-1)+50
    lsun = 3.839d33
    Lnorm = dindgen(nr)
    Lnorm(*) = 1d0
    minv = 1d5
    maxv = -1d5
    for iq=0,nq-1 do begin
        temp = min(4d0*!PI*r^2*avgs(*,iq)/lsun)
        if (temp lt minv) then minv=temp
        temp = max(4d0*!PI*r^2*avgs(*,iq)/lsun)
        if (temp gt maxv) then maxv=temp
    endfor
    plot, r, Lnorm, title='Fluxes', ytitle='Fluxes', xtitle='r (cm)', yrange=[minv,maxv], /xs
    for iq=1,nq-1 do begin
        oplot, r, 4d0*!PI*r^2*avgs(*,iq)/lsun, color=colors(iq)
    endfor
    lines = intarr(nq)
    lines(*) = 0
    legend,quantity_names,linestyle=lines,colors=colors

    stop
end

pro handle_print_event, ev
    common common_block, data, r, files, quantity_names, times, avgs
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
    DEVICE, filename = psfile, BITS = 8,/COLOR,/HELVETICA,/ENCAPSULATED
    DEVICE, XSIZE = xsize, YSIZE = ysize, XOFFSET=xoffset, YOFFSET=yoffset

    nq = n_elements(quantity_names)
    colors = 205*indgen(nq-1)/(nq-2)+50

    colors = 205*indgen(nq-1)/(nq-2)+50
    plot, r, avgs(*,0), title='Fluxes', ytitle='Fluxes', xtitle='r (cm)', yrange=[min(avgs(*,*)),max(avgs(*,*))]
    for iq=1,nq-1 do begin
        oplot, r, avgs(*,iq), color=colors(iq)
    endfor
    lines = intarr(nq)
    lines(*) = 0
    legend,quantity_names,linestyle=lines,colors=colors

    DEVICE, /close
    SET_PLOT, 'x'
    command = 'kghostview ' + psfile + ' &'
    SPAWN, command
end

pro handle_file_load_event, ev
    common common_block, data, r, files, quantity_names, times, avgs
    WIDGET_CONTROL,ev.top,GET_UVALUE=state

    ;Open a file dialog box to select (a) file(s).
    root_directory = '/home8/ryor5023/Programs/CSS_Code/Runs/'
    files = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select Azimuthal Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)
    state.nf = n_files

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_flux_bal, files(0), radius = r, data = values, QUANTITIES = quantities, q_names=quantity_names, times=times

    num_records = n_elements(values(0,0,*))
    state.num_records = num_records

    ;Determine the number of radial points.
    N_R = n_elements(r)
    state.nr = N_R

    ;Determine the number of quantities.
    N_Q = n_elements(quantities)
    state.nq = N_Q

    WIDGET_CONTROL,ev.top,SET_UVALUE = state
    data = values

    ;Build averages
    tmax = max(times)
    if (tmax eq 0d0) then begin
       dt = dblarr(state.num_records)
       dt(*) = 1d0
       deltat = state.num_records
    endif else begin
       dt = times-shift(times,1)
       dt(0)=times(1)-times(0)
       deltat=times(state.num_records-1)-times(0)
    endelse
    avgs = dblarr(state.nr,state.nq)
    for iq=0,state.nq-1 do begin
        avgs(*,iq) = 0d0
        for ir=0,state.num_records-1 do begin
            avgs(*,iq)=avgs(*,iq)+dt(ir)*data(*,iq,ir)
        endfor
        avgs(*,iq) = avgs(*,iq)/deltat
    endfor

    handle_drawing_event,state.drawing_widget
    
end

pro plot_flux_bal_gui
    common common_block, data, r, files, quantity_names, times, avgs
    data=0
    r=0
    files=''
    ;for now use docolor, eventually replace with customizable color table
    docolor
    !P.CHARSIZE=1.5
    ;Create the parent widget.
    base = WIDGET_BASE(/ROW,TITLE = 'Flux Balance Viewer',/TLB_SIZE_EVENTS)

    get_state, state
    ;Create two columns, one for drawing, and one for buttons and options
    drawing_column = WIDGET_BASE(base,/COLUMN)
    options_column = WIDGET_BASE(base,/COLUMN)
    ;A drawing widget attached to the drawing column
    drawing_widget = WIDGET_DRAW(drawing_column, XSIZE=state.window_size.x, YSIZE=state.window_size.y, FRAME = 2,UVALUE='drawing')
    ;A button to print the current plot
    print_button = WIDGET_BUTTON(options_column,VALUE='Print',UVALUE='print')
    ;Create a button widget to load files.
    load_button = WIDGET_BUTTON(options_column,VALUE='Load File(s)',UVALUE='load_file')
    ;Create a button widget to exit the IDL program.
    exit_button = WIDGET_BUTTON(options_column,VALUE='Exit',UVALUE='exit')
    ;Create widget event controller
    WIDGET_CONTROL,base,/REALIZE
    ;Get drawing_widget id
    WIDGET_CONTROL, drawing_widget, GET_VALUE=drawing_widget_id

    ;Save widget ids
    state.drawing_widget = drawing_widget_id
    state.options_id = options_column
    WSET, drawing_widget_id
    XMANAGER, 'base', base, /NO_BLOCK

    WIDGET_CONTROL,base,TLB_GET_SIZE=base_size

    state.base_size.x = base_size(0)
    state.base_size.y = base_size(1)
    WIDGET_CONTROL,base,SET_UVALUE=state
end
