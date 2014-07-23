pro write_bulk_cartesian, data, bounds, filename, triangle_file=triangle_file, resize=resize, dir=dir, stop=stop, add_inner=add_inner,cube_size=cube_size,extents=extents

    nr = round(bounds(0))
    nth = round(bounds(1))
    nph = round(bounds(2))
    r1 = bounds(3)
    r2 = bounds(4)
    th1 = bounds(5)
    th2 = bounds(6)
    ph1 = bounds(7)
    ph2 = bounds(8)
    

    If not keyword_set(resize) Then Begin
       resize=1
    Endif
 
    nr_old = nr
    nth_old = nth
    nph_old = nph
    nr = LONG(nr)
    nth = LONG(nth/resize)
    nph = LONG(nph/resize)

    If not keyword_set(dir) Then Begin
       dir='/alphacen/augustso/CSSD/base_case/Checkpoints/'
    EndIf
 
    ;Set up grid (assumes regular r, theta, phi grid (so pre-interpolate!))
    radius = (r2-r1)*dindgen(nr)/double(nr-1)+r1
    theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
    phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1

    radius = radius/r2 ;normalize radius

    ;Set up cartesian grid to interpolate to
    n_tot = LONG(nr*nth*nph)
    ;print, 'Generating grid:'
    x = dblarr(n_tot)
    y = dblarr(n_tot)
    z = dblarr(n_tot)
    sint = sin(theta)
    cost = cos(theta)
    sinp = sin(phi)
    cosp = cos(phi)

    i=0L
    for ir = 0, nr-1 do begin
        for ip = 0, nph-1 do begin
            x(i:i+1L*(nth-1)) = radius(ir)*sint*cosp(ip)
            y(i:i+1L*(nth-1)) = radius(ir)*sint*sinp(ip)
            z(i:i+1L*(nth-1)) = radius(ir)*cost
            i=i+1L*(nth)
         endfor
        ;print, 'Radius = ', radius(ir)
    endfor

    x(n_tot-1L) = radius(nr-1)*sint(nth-1)*cosp(nph-1)
    y(n_tot-1L) = radius(nr-1)*sint(nth-1)*sinp(nph-1)
    z(n_tot-1L) = radius(nr-1)*cost(nth-1)

    If (not keyword_set(triangle_file)) Then Begin
       ;Construct convex hulls of the volume given the points above
       print, 'Creating a convex hull of tetrahedra of the grid'

       qhull, x, y, z, tetrahedra, /DELAUNAY

       ;write triangles
       openw, file_unit, dir+'tetrahedra.head', /get_lun
       temp = size(tetrahedra)
       writeu, file_unit, temp(2)
       close, file_unit
       free_lun, file_unit

       openw, file_unit, dir+'tetrahedra.dat', /get_lun
       a = assoc(file_unit, lonarr(4,temp(2)))
       a(0) = tetrahedra
       close, file_unit
       free_lun, file_unit
    EndIf Else Begin
       openr, file_unit, dir+'tetrahedra.head', /get_lun
       num_triangles = 1L*intarr(1)
       readu, file_unit, num_triangles
       close, file_unit
       free_lun, file_unit

       openr, file_unit, dir+'tetrahedra.dat', /get_lun
       temp = assoc(file_unit, lonarr(4,num_triangles))
       tetrahedra = reform(temp(0))
       close, file_unit
       free_lun, file_unit
    EndElse

    ;print, 'Reading in the data'

    long_data = dblarr(n_tot)
    i=0L
    If (resize gt 1) Then Begin
       new_data = congrid(data,nth,nph,nr)
    EndIf Else Begin
       new_data = data
    EndElse
    for ir = 0L, nr-1L do begin
        for ip = 0L, nph-1L do begin
            long_data(i:i+1L*(nth-1)) = reform(new_data(*,ip,ir))
            i=i+1L*(nth)
        endfor
    endfor

     ;Volume parameters
     mxr = max(radius)
     mnr = min(radius)
     mxth = max(theta)
     mnth = min(theta)
     mxph = max(phi)
     mnph = min(phi)

     xmax = max(x)
     xmin = min(x)
     ymax = max(y)
     ymin = min(y)
     zmax = max(z)
     zmin = min(z)

     nx = min([round(nr*(xmax-xmin)/(mxr-mnr)),512L])
     ny = round(nx*(ymax-ymin)/(xmax-xmin))
     nz = round(nx*(zmax-zmin)/(xmax-xmin))
     cube_size = [nx,ny,nz]
     print, 'Cube_Size = [', cube_size, ']'
     extents = '-extents ' + strtrim(string(min(x)),2) + ':' + strtrim(string(min(y)),2) + ':' + strtrim(string(min(z)),2) + ':' + $
            strtrim(string(max(x)),2) + ':' + strtrim(string(max(y)),2) + ':' + strtrim(string(max(z)),2) + ' '
     print, extents

     mind = min(long_data)
     missing_value = max([mind,0.0d0])
     interp_data = qgrid3(x,y,z,long_data,tetrahedra,DIMENSION=cube_size,MISSING=missing_value)

     ;Make sure that anything with r less than r1 is minimum value!
     x = (xmax-xmin)*dindgen(nx)/(nx-1.0)+xmin
     y = (ymax-ymin)*dindgen(ny)/(ny-1.0)+ymin
     z = (zmax-zmin)*dindgen(nz)/(nz-1.0)+zmin
     for iz=0,nz-1 do begin
         for iy=0,ny-1 do begin
             rtemp = sqrt(x*x+y(iy)^2+z(iz)^2)            
             indices = where(rtemp lt radius(0))
             If (n_elements(indices) gt 1) Then Begin
                interp_data(indices,iy,iz)=missing_value
             Endif
         endfor
      endfor

     ;Make sure that anything with r gt than r2 is minimum value!
     for iz=0,nz-1 do begin
         for iy=0,ny-1 do begin
             rtemp = sqrt(x*x+y(iy)^2+z(iz)^2)            
             indices = where(rtemp gt radius(nr-1))
             If (n_elements(indices) gt 1) Then Begin
                interp_data(indices,iy,iz)=missing_value
             Endif
         endfor
     endfor
 
                                ;Write out data cube
     If (byte(1,0)) Then Begin
         openw, file_unit, dir+filename, /get_lun, /swap_endian
     EndIf Else Begin
         openw, file_unit, dir+filename, /get_lun
     EndElse

     a = assoc(file_unit, dblarr(nx,ny,nz))
     a(0) = interp_data
     close,file_unit
     free_lun, file_unit 

     If (keyword_set(stop)) Then Begin
        stop
     EndIf
end
