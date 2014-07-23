pro background_compare, css_case, gamma=gamma, Cp=Cp

    ;Open a file dialog box to select (a) file(s).
    root_directory = '/freyr3/augustso/CSSDE/'+css_case+'/Shell_Avgs/'
    files = DIALOG_PICKFILE(PATH=root_directory,/MULTIPLE_FILES,TITLE='Select Shell Average Files to Read')
 
    ;Determine the number of files.
    n_files = n_elements(files)

    ;Handle the event that the user cancels.
    if(n_files eq 1) then begin
       if(files eq '') then return
     endif

    ;We assume the files read in are from the same run.
    read_shell_avg, files(0), radius = r, data = values, QUANTITIES = quantities

    num_records = n_elements(values(0,0,*))

    ;Determine the number of radial points.
    N_R = n_elements(r)

    ;Determine the number of quantities.
    N_Q = n_elements(quantities)

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

    root_directory = '/freyr3/augustso/CSSDE/'+css_case
;Read in background file
    openr, funit, root_directory+'/background.dat', /get_lun
    readf, funit, format='(i5)', n
    gravn = dblarr(n)
    rhon = dblarr(n)
    Sn = dblarr(n)
    dsdrn = dblarr(n)
    mur = dblarr(n)
    dmudr = dblarr(n)
    opacity = dblarr(n)
    dopdlnd = dblarr(n)
    dopdlnT = dblarr(n)
    kappas = dblarr(n)
    dkappasdr = dblarr(n)
    kappa0 = dblarr(n)
    dkappa0dr = dblarr(n)
    epsn = dblarr(n)
    readf, funit, format='(E16.9)', gravn
    readf, funit, format='(E16.9)', rhon
    readf, funit, format='(E16.9)', Sn
    readf, funit, format='(E16.9)', dsdrn
    readf, funit, format='(E16.9)', mur
    readf, funit, format='(E16.9)', dmudr
    readf, funit, format='(E16.9)', opacity
    readf, funit, format='(E16.9)', dopdlnd
    readf, funit, format='(E16.9)', dopdlnT
    readf, funit, format='(E16.9)', kappas
    readf, funit, format='(E16.9)', dkappasdr
    readf, funit, format='(E16.9)', kappa0
    readf, funit, format='(E16.9)', dkappa0dr
    readf, funit, format='(E16.9)', epsn
    readf, funit, format='(E16.9)', Lradn
    close, funit
    free_lun, funit 

    Cv = Cp/gamma
    gam1 = gamma-1d0

    ;get density
    rho1 = data[*,0,num_records-1L]

    ;get temperature
    Pn = rhon^gamma*exp(Sn/Cv)
    Tn = Pn/gam1/Cv/rhon
    T1 = data[*,1,num_records-1L]

    ;get entropy
    S1 = data[*,2,num_records-1L]

    ;get pressure
    P1 = data[*,3,num_records-1L]

    opac1 = opacity+(alog(rho1)-alog(rhon))*dopdlnd + (alog(T1)-alog(Tn))*dopdlnT

    docolor
    ;window, 1, xsize=1280, ysize=960
    xsize = 10.0*2.54
    ysize = xsize/1.5
    xoffset = (8.5*2.54 - xsize)/2.0
    yoffset = (11.0*2.54 - ysize)/2.0
    print, 'Output size (cm)',xsize, ysize
    SET_PLOT, 'ps'
    DEVICE, filename = 'bkgd_compare.eps', BITS = 8,/COLOR,/HELVETICA,/ENCAPSULATED
    DEVICE, XSIZE = xsize, YSIZE = ysize, XOFFSET=xoffset, YOFFSET=yoffset
    !P.MULTI=[0,2,1]
    !P.CHARSIZE=1
    !P.THICK=4
    colors = [0,80,110,150,190]
    plot, r, (rho1-rhon)/rhon, /xs, yrange=[-0.1,0.1], title='Relative Thermodynamic Changes', xtitle='R (cm)', ytitle='Relative Change'
    oplot, r, (T1-Tn)/Tn, color=colors[1]
    oplot, r, (S1-Sn)/Sn, color=colors[2]
    oplot, r, (P1-Pn)/Pn, color=colors[3]
    oplot, r, (opac1-opacity)/opacity, color=colors[4]
    oplot, r, 0d0*r, linestyle=2
    names = ['Density','Temperature','Entropy','Pressure','Opacity']
    legend, names, linestyle=colors*0, color=colors, /bottom, /left

    plot, r, dsdrn,/xs, yrange=[min(deriv(r,S1)),max(deriv(r,S1))], title='Entropy Gradient', xtitle='R (cm)', ytitle='dS/dr (erg/g K cm)'
    oplot, r, deriv(r,S1), color=150
    oplot, r, 0d0*r, linestyle=2
    names = ['dsdr0','dsdr1']
    legend, names, linestyle=colors[0:1]*0, color=[0,150], /bottom, /left
    DEVICE, /close
    SET_PLOT, 'x'
    command = 'kghostview bkgd_compare.eps &'
    SPAWN, command

    stop
end
