pro build_cartesian_chunk, iter, dir1=dir1, dir2=dir2

    if (not keyword_set(dir1)) then begin
       dir1 = '/alphacen/augustso/CSSDE/Vapor/new_ebc_small/Checkpoints/'
    endif

    if (not keyword_set(dir2)) then begin
       dir2 = '/alphacen/augustso/CSSDE/Vapor/new_ebc_small/Checkpoints/'
    endif

    k=100000L
    If (iter lt k) Then Begin
        zero_char = '0'
        i=5L
        pretemp = ' '
        stemp = ' '
        While ((iter lt k) and (i gt 0L)) Do Begin
            short_temp = string(i,format='(i1)')
            format = '(i'+strtrim(short_temp,2)+')'
            fmt = strtrim(format,2)
            stemp = string(iter,format=fmt)
            pretemp = strtrim(pretemp+zero_char,2)
            i=i-1L
            k=k/10L
        EndWhile
    EndIf Else Begin
        pretemp = ' '
        stemp = string(iter, format='(i6)')
    EndElse
    
    dir1 = strtrim(dir1,2)+strtrim(pretemp+stemp+'/',2)
    dir2 = strtrim(dir2,2)+strtrim(pretemp+stemp+'/',2)+'Cartesian/'

    fname = strtrim(dir1+'checkpoint_'+strtrim(pretemp+stemp,2)+'_header',2)

    openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format

    temp = intarr(5)
    readf, file_unit, temp,format='(5i6)'

    nr   = fix(temp(0))
    nth  = fix(temp(1))
    nph  = fix(temp(2))
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

    fname = strtrim(dir1+'checkpoint_'+strtrim(pretemp+stemp,2),2)
    
    openr, file_unit, fname, /get_lun, /swap_if_big_endian
    volumes = assoc(file_unit,dblarr(nr,nth,nph,nq))

    names  = ['rho_cart.dat','ut_cart.dat','up_cart.dat','ur_cart.dat','S_cart.dat','rho_fluc_cart.dat','S_fluc_cart.dat']
    qnames = ['rho','ut','up','ur','S','rho_fluc','S_fluc']

    data = dblarr(nr,nth,nph,nq+2)
    data(*,*,*,0:nq-1) = reform(volumes(*,*,*,*,0))
    
    close, file_unit
    free_lun, file_unit

    rho_avg = dblarr(nr)
    s_avg = dblarr(nr)

    for ir=0,nr-1 do begin
        rho_avg(ir) = mean(data(ir,*,*,0))
        s_avg(ir) = mean(data(ir,*,*,4))
        data(ir,*,*,nq) = data(ir,*,*,0)-rho_avg(ir)
        data(ir,*,*,nq+1) = data(ir,*,*,4)-s_avg(ir)
    endfor

    ;Set up r, theta, phi arrays, adjusting coordinate system so that ph1=0
    dr = (r2-r1)/nr
    dth = (th2-th1)/nth
    dph = (ph2-ph1)/nph
    r = dr*dindgen(nr)+r1
    theta = dth*dindgen(nth)+th1
    phi = dph*dindgen(nph)

    sint = sin(theta)
    cost = cos(theta)
    sinp = sin(phi)
    cosp = cos(phi)

    ;Determine nx,ny,nz from nr,nth,nphi

    deltax = max(r2*sint-r1*sin(th1)*cos(ph2-ph1))
    nx = fix(double(nr)*deltax/(r2-r1))
    
    it_max = where(r2*sint-r1*sin(th1)*cos(ph2-ph1) eq deltax)

    dy = min(r1*min(sint)*abs(sinp-shift(sinp,1)))
    ny = fix(r2*sint(it_max)*max(sinp)/dy)

    dz = min(r1*abs(shift(cost,1)-cost))
    nz = fix(r2*(cos(th1)-cos(th2))/dz)

    Print, 'nx,ny,nz=',nx,ny,nz

    ;Set up r,theta,phi

end
