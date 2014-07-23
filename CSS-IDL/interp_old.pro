pro test2d

  arr = dblarr(128,128)
  ox = 4d0*dindgen(128)/127d0-2d0
  oy = 4d0*dindgen(128)/127d0-2d0
  for ix=0,127 do begin
     arr[ix,*] = sin(!dpi*ox[ix])*sin(!dpi*oy)
  endfor

  nx = 2d0*dindgen(512)/511d0-1d0
  ny = 2d0*dindgen(512)/511d0-1d0
  interp2d, arr, ox, oy, nx, ny, newarr

stop

end

pro interp2d, arr, ox, oy, nx, ny, newarr
  nnx = n_elements(nx)
  nny = n_elements(ny)
  newarr = dblarr(nnx,nny)
  const = [nx[0],nx[nnx-1],ny[0],ny[nny-1]]

  onx = n_elements(ox)
  ony = n_elements(oy)
  
  x1 = const[0]
  x2 = const[1]
  y1 = const[2]
  y2 = const[3]
  odx = ox[1]-ox[0]
  ody = oy[1]-oy[0]

  print, odx, ody, x1, x2, y1, y2
  ix = intarr(4)
  iy = intarr(4)
  ry = dblarr(4)

  for i2=0,nny-1 do begin
     yc = ny[i2]
     for i1=0,nnx-1 do begin
        xc = nx[i1]        
 ;For bicubic convolution interpolation, need nearest 2 points from old grid
 ;  in each direction. Find them in O(1) time.
        ix[1] = ((Floor((xc-ox[0])/odx)) mod onx)
        iy[1] = Max([Floor((yc-oy[0])/ody), 0])
        ix[2] = ((ix[1]+1) mod onx)
        iy[2] = Min([iy[1]+1,ony-1])
        ix[0] = ix[1]-1
        If (ix[0] lt 1) Then ix[0] = onx+ix[0]-1
        iy[0] = Max([iy[1]-1,0])
        ix[3] = ((ix[1]+2) mod onx)
        iy[3] = Min([iy[1]+2,ony-1])

        ry[*] = 0d0
        for ii=0,3 do begin
           weights = 0d0
           for ij=0,3 do begin
              temp = (yc-oy[iy[1]])/ody ;(oy[iy[2]]-oy[iy[1]])
              temp = kernelo3(temp)
              weights = weights + temp
              ry[ii] = ry[ii] + temp*arr[ix[ii],iy[ij]]
           endfor
           ry[ii] = ry[ii]/weights
        endfor

        ux = 0d0
        weights = 0d0
        for ii=0,3 do begin
           temp = (xc-ox[ix[1]])/odx ;(ox[ix[2]]-ox[ix[1]])
           temp = kernelo3(temp)
           weights = weights + temp
           ux = ux + temp*ry[ii]
        endfor
        newarr[i1,i2] = ux/weights
     endfor
  endfor
  stop
end

function kernelo3, ky
  kx = abs(ky)
  If (kx le 1d0) Then Begin
     res = 1d0+kx*kx*(-2.5d0+1.5d0*kx)
  Endif Else Begin
     If (kx lt 2d0) Then Begin
        res = 2d0+kx*(-4d0+2.5d0*kx-0.5d0*kx*kx)
     Endif Else Begin
        res = 0d0
     EndElse
  Endelse
  return, res
end
