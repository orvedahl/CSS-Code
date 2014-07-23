@write_bulk_cartesian.pro
@get_checkpoint_quantity_name.pro

Pro make_vapor_movie, iq, istart, iend, vdfname, iskip=iskip, read_dir=read_dir, write_dir=write_dir, triangle_file=triangle_file

    If (not keyword_set(iskip)) Then iskip=500L
    If (not keyword_set(read_dir)) Then read_dir='/alphacen/augustso/CSSD/entropy_base_case/Checkpoints/'
    If (not keyword_set(write_dir)) Then write_dir='/alphacen/augustso/CSSD/Vapor/entropy_base_case/movie/'

    exists = file_test(read_dir, /directory)
    If (exists eq 0) Then Begin 
        Print, 'Read directory does not exist.' 
        stop 
    Endif

    exists = file_test(write_dir, /directory)
    If (exists eq 0) Then Begin 
        Print, 'Write directory does not exist.' 
        Print, 'Making write directory...'
        file_mkdir, write_dir
        file_chmod, write_dir, '777'o
    EndIf

    Print, 'Reading in files... assuming all from the same run'

    k=100000L
    If (istart lt k) Then Begin
        zero_char = '0'
        i=5L
        pretemp = ' '
        stemp = ' '
        While ((istart lt k) and (i gt 0L)) Do Begin
            short_temp = string(i,format='(i1)')
            format = '(i'+strtrim(short_temp,2)+')'
            fmt = strtrim(format,2)
            stemp = string(istart,format=fmt)
            pretemp = strtrim(pretemp+zero_char,2)
            i=i-1L
            k=k/10L 
        EndWhile
    EndIf Else Begin
        stemp = string(istart, format='i6')
    EndElse

    fname = strtrim(read_dir+'checkpoint_'+strtrim(pretemp+stemp,2)+'_header',2)

    If (byte(1,0)) Then Begin
        openr, file_unit, fname, /get_lun, /swap_endian
    Endif Else Begin
        openr, file_unit, fname, /get_lun
    EndElse

    temp = intarr(5)
    readf, file_unit, temp,format='(5i6)'

    nth  = fix(temp(0))
    nph  = fix(temp(1))
    nr   = fix(temp(2))
    nq   = fix(temp(3))
    iter = fix(temp(4))

    print, temp

    temp = dblarr(8)
    readf, file_unit, temp, format='(8e14.7)'
    dt   = temp(0)
    time = temp(1)
    r1   = temp(2)
    r2   = temp(3)
    th1  = temp(4)
    th2  = temp(5)
    ph1  = temp(6)
    ph2  = temp(7)

    print, temp

    close, file_unit
    free_lun, file_unit
    
    bounds = [nr,nth,nph,r1,r2,th1,th2,ph1,ph2]

    print, 'bounds = [', bounds, ']'

    ;Read in quantity from files
    qname = get_checkpoint_quantity_name(iq)
    iq = iq-1L
    nf = (iend-istart)/iskip+1L
    clearline = fifteenb()
    For ifile=istart,iend,iskip Do Begin

        k=100000L
        If (ifile lt k) Then Begin
            zero_char = '0'
            i=5L
            pretemp = ' '
            stemp = ' '
            While ((ifile lt k) and (i gt 0L)) Do Begin
                short_temp = string(i,format='(i1)')
                format = '(i'+strtrim(short_temp,2)+')'
                fmt = strtrim(format,2)
                stemp = string(ifile,format=fmt)
                pretemp = strtrim(pretemp+zero_char,2)
                i=i-1L
                k=k/10L
            EndWhile
        EndIf Else Begin
            stemp = string(ifile, format='i6')
        EndElse

        fname = strtrim(read_dir+'checkpoint_'+strtrim(pretemp+stemp,2),2)
        If (byte(1,0)) Then Begin
            openr, file_unit, fname, /get_lun
        EndIf Else Begin
            openr, file_unit, fname, /get_lun, /swap_endian
        EndElse

        volumes = assoc(file_unit,dblarr(nth,nph,nr,nq))
        data = reform(volumes(*,*,*,iq,0))

        close, file_unit
        free_lun, file_unit

        filename = strtrim(qname+'_'+strtrim(pretemp+stemp,2)+'.dat',2)

        If (ifile eq istart) Then Begin
            write_bulk_cartesian, data, bounds, filename, triangle_file=triangle_file, dir=write_dir, cube_size=cube_size, extents=extents
            dim = strtrim(strtrim(string(cube_size(0),format='(i5)'),2)+'x'+strtrim(string(cube_size(1),format='(i5)'),2)+'x'+strtrim(string(cube_size(2),format='(i5)'),2),2)
            cmd = 'vdfcreate -coordsystem "cartesian" -dimension ' + dim + ' -level 1 -order 0:1:2 '+extents+'-numts '+ $
              strtrim(string(nf,format='(i5)'),2) + ' -vars3d ' + strtrim(qname,2) + ' ' + write_dir + vdfname
            print, cmd
            spawn, cmd
        Endif Else Begin
            write_bulk_cartesian, data, bounds, filename, /triangle_file, dir=write_dir, cube_size=cube_size, extents=extents
        EndElse

        If (byte(1,0)) Then Begin
            cmd = 'raw2vdf -dbl -swapbyte -varname ' + strtrim(qname,2) + ' -ts ' + strtrim(string((ifile-istart)/iskip,format='(i5)'),2) + $
              ' ' + write_dir + vdfname + ' ' + write_dir + filename
        Endif Else Begin
            cmd = 'raw2vdf -dbl -varname ' + strtrim(qname,2) + ' -ts ' + strtrim(string((ifile-istart)/iskip,format='(i5)'),2) + $
              ' ' + write_dir + vdfname + ' ' + write_dir + filename
        EndElse

        print, cmd
        spawn, cmd

        per = strtrim(round((1d0*(ifile-istart)/iskip)/(1d0*(nf-1d0))*100.0d0),2)
        lenp = strtrim(strlen(strtrim(per,2)),2)
        form="($,a"+lenp+",' % Completed',a,a)"
        print, form=form, per, '         ', clearline
    EndFor

End

function fifteenb

   return, string("15b)

end

