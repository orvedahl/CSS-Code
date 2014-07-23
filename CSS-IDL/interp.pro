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

        mu = (xc-ox[ix[1]])/(ox[ix[2]]-ox[ix[1]])
        mu2 = mu*mu
        mu3 = mu*mu2
        ry[*] = 0d0
        for ij=0,3 do begin
           akm1 = arr[ix[0],iy[ij]]
           ak   = arr[ix[1],iy[ij]]
           akp1 = arr[ix[2],iy[ij]]
           akp2 = arr[ix[3],iy[ij]]
        
           a3 = ak 
           a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
           a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
           a2 = -0.5d0*akm1 + 0.5d0*akp1
           ry[ij] = a0*mu3+a1*mu2+a2*mu+a3           
        endfor

        mu = (yc-oy[iy[1]])/(oy[iy[2]]-oy[iy[1]])
        mu2 = mu*mu
        mu3 = mu*mu2

        akm1 = ry[0]
        ak   = ry[1]
        akp1 = ry[2]
        akp2 = ry[3]
        
        a3 = ak 
        a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
        a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
        a2 = -0.5d0*akm1 + 0.5d0*akp1
        newarr[i1,i2] = a0*mu3+a1*mu2+a2*mu+a3     
     endfor
  endfor
  
end

pro interp3d, arr, ox, oy, oz, nx, ny, nz, newarr
  nnx = n_elements(nx)
  nny = n_elements(ny)
  nnz = n_elements(nz)
  newarr = dblarr(nnx,nny,nnz)
  const = [nx[0],nx[nnx-1],ny[0],ny[nny-1],nz[0],nz[nnz-1]]

  onx = n_elements(ox)
  ony = n_elements(oy)
  onz = n_elements(oz)
  
  x1 = const[0]
  x2 = const[1]
  y1 = const[2]
  y2 = const[3]
  z1 = const[4]
  z2 = const[5]

  odx = ox[1]-ox[0]
  ody = oy[1]-oy[0]
  odz = oz[1]-oz[0]

  print, odx, ody, odz, x1, x2, y1, y2, z1, z2
  ix = intarr(4)
  iy = intarr(4)
  iz = intarr(4)
  ry = dblarr(4,4)
  rz = dblarr(4)

  for i3=0,nnz-1 do begin
     zc = nz[i3]
     for i2=0,nny-1 do begin
        yc = ny[i2]
        for i1=0,nnx-1 do begin
           xc = nx[i1]        
 ;For tricubic convolution interpolation, need nearest 2 points from old grid
 ;  in each direction. Find them in O(1) time.
           ix[1] = ((Floor((xc-ox[0])/odx)) mod onx)
           iy[1] = Max([Floor((yc-oy[0])/ody), 0])
           iz[1] = Max([Floor((zc-oz[0])/odz), 0])
           ix[2] = ((ix[1]+1) mod onx)
           iy[2] = Min([iy[1]+1,ony-1])
           iz[2] = Min([iz[1]+1,onz-1])
           ix[0] = ix[1]-1
           If (ix[0] lt 1) Then ix[0] = onx+ix[0]-1
           iy[0] = Max([iy[1]-1,0])
           iz[0] = Max([iz[1]-1,0])
           ix[3] = ((ix[1]+2) mod onx)
           iy[3] = Min([iy[1]+2,ony-1])
           iz[3] = Min([iz[1]+2,onz-1])
           
           mu = (xc-ox[ix[1]])/(ox[ix[2]]-ox[ix[1]])
           mu2 = mu*mu
           mu3 = mu*mu2
           ry[*] = 0d0
           for ik=0,3 do begin
              for ij=0,3 do begin
                 akm1 = arr[ix[0],iy[ij],iz[ik]]
                 ak   = arr[ix[1],iy[ij],iz[ik]]
                 akp1 = arr[ix[2],iy[ij],iz[ik]]
                 akp2 = arr[ix[3],iy[ij],iz[ik]]
              
                 a3 = ak 
                 a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                 a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                 a2 = -0.5d0*akm1 + 0.5d0*akp1
                 ry[ij,ik] = a0*mu3+a1*mu2+a2*mu+a3           
              endfor
           endfor

           mu = (yc-oy[iy[1]])/(oy[iy[2]]-oy[iy[1]])
           mu2 = mu*mu
           mu3 = mu*mu2           
           for ik=0,3 do begin       
              akm1 = ry[0,ik]
              ak   = ry[1,ik]
              akp1 = ry[2,ik]
              akp2 = ry[3,ik]
           
              a3 = ak 
              a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
              a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
              a2 = -0.5d0*akm1 + 0.5d0*akp1
              rz[ik] = a0*mu3+a1*mu2+a2*mu+a3
           endfor

           mu = (zc-oz[iz[1]])/(oz[iz[2]]-oz[iz[1]])
           mu2 = mu*mu
           mu3 = mu*mu2

           akm1 = rz[0]
           ak   = rz[1]
           akp1 = rz[2]
           akp2 = rz[3]
           
           a3 = ak 
           a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
           a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
           a2 = -0.5d0*akm1 + 0.5d0*akp1
           newarr[i1,i2,i3] = a0*mu3+a1*mu2+a2*mu+a3
        endfor
     endfor
     Print, 'At z level', i3
     Print, fix(1d2*double(i3)/double(nnz-1)), '% complete'
  endfor
end
