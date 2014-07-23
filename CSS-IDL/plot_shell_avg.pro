@docolor.pro
@read_shell_avg.pro
pro plot_shell_avg, qoi=qoi

    If (not keyword_set(qoi)) Then Begin
       qoi=2 
    EndIf
    ;Open a file dialog box to select (a) file(s).
    root_directory = '/alphacen/augustso/CSS_volume/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius=r, data=values, q_names=q_names, quantities=quantities
    num_records = n_elements(values(0,0,*))

    ;Determine the number of radial points.
    nr = n_elements(r)

    ;Determine the number of quantities.
    nq = n_elements(quantities)

    ;average = dblarr(nr,nq)
    ;average = total(values,2)/(1.0d0*num_records)

    iq = where(quantities eq qoi)

    docolor
    colors = 205*indgen(num_records-1)/(num_records-2)+50
    plot, r, values(*,iq,0), ytitle=q_names(iq), xtitle='r (cm)'
    for i=1,num_records-1 do begin
        oplot, r, values(*,iq,i), color=colors(i-1)
    endfor
    stop
end
