pro triangle_test, nx, nz, nvert

  vert=dblarr(2,nvert) ;CCW list of vertices
  vert[0,*] = randomu(systime(/seconds),nvert,/DOUBLE,/UNIFORM) ;[0.7d0,0.9d0,0.25d0] ;x values
  vert[1,*] = randomu(2*systime(/seconds),nvert,/DOUBLE,/UNIFORM) ;[0.2d0,0.7d0,0.8d0] ;y values
  ;Sort the vertices by x value
  sinds = sort(reform(sqrt((vert[0,*]-0.5d0)^2+(vert[1,*]-0.5d0)^2)))
  vert[0,*] = vert[0,sinds]
  vert[1,*] = vert[1,sinds]

  ;Define plane
  ny = nx
  dx = 1d0/double(nx-1)
  dy = 1d0/double(ny-1)     
  xarr = dx*dindgen(nx)
  yarr = dy*dindgen(ny)

  plane = dblarr(nx,ny)
  plane[*,*] = 0d0

 ;loop over edges
  for iv=0,nvert-1 do begin
 ;iterate along edge
 ;Determine number of points crossed     
     If (iv ne nvert-1) Then Begin
        dvx = (vert[0,iv+1]-vert[0,iv]) 
        dvy = (vert[1,iv+1]-vert[1,iv]) 
     Endif Else Begin
        dvx = (vert[0,0]-vert[0,iv])
        dvy = (vert[1,0]-vert[1,iv])
     Endelse
     mx = dvx/dx
     my = dvy/dy
     my = fix(ceil(sqrt(mx^2+my^2)))/30
     for ip=0L,my-1L do begin
        tmp = double(ip)/double(my)
        xv = vert[0,iv]+dvx*tmp
        yv = vert[1,iv]+dvy*tmp

        ix = fix(xv/dx)
        iy = fix(yv/dy)
        for ity=max([iy-nz,0]),min([iy+nz,ny-1]) do begin
           z1 = (yarr[ity]-yv)/(2d0*nz*dy)
           ;z1 = min([z1^2,1d0])
           ;z1 = z1^2-2d0*z1+1d0
           for itx=max([ix-nz,0]),min([ix+nz,nx-1]) do begin
              z2 = (xarr[itx]-xv)/(2d0*nz*dx)
              z2 = min([z2^2+z1^2,1d0])
              z2 = z2^2-2d0*z2+1d0
              plane[itx,ity] = plane[itx,ity]+z2
           endfor
        endfor
     endfor
  endfor
  stop

;        xinds = ix + zs
;        xti = where((xinds lt nx)and(xinds ge 0))
;        xinds = xinds[xti]
;        ;zx = xinds-ix+nz/2

;        yinds = iy + zs
;        yti = where((yinds lt ny)and(yinds ge 0))
;        yinds = yinds[yti]
;        nty = n_elements(yinds)
        ;zy = yinds-iy+nz/2

;        for ity=0L,nty-1L do begin
           ;plane[xinds,yinds[ity]] = plane[xinds,yinds[ity]] + kernel[zx,zy[ity]]
;           plane[xinds,yinds[ity]] = plane[xinds,yinds[ity]] + exp(-1d4*((xarr[xinds]-xv)^2+(yarr[yinds[ity]]-yv)^2))
;        endfor
;     endfor
;  endfor
;  stop
end
