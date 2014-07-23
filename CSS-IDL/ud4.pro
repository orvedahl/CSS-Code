function ud4, dx, arr
  nx = n_elements(arr)
  darr = dblarr(nx)

  darr[0] = (2d0*arr[3]-9d0*arr[2]+18d0*arr[1]-11d0*arr[0])/6d0
  darr[1] = (2d0*arr[4]-9d0*arr[3]+18d0*arr[2]-11d0*arr[1])/6d0
  for ix=2,nx-3 do begin
     darr[ix] = arr[ix+2]+2d0*arr[ix+1]-2d0*arr[ix-1]-arr[ix-2]
  endfor
  darr[2:nx-3] = darr[2:nx-3]*0.25d0
  ix = nx-2
  darr[ix] = -(2d0*arr[ix-3]-9d0*arr[ix-2]+18d0*arr[ix-1]-11d0*arr[ix])/6d0
  ix = nx-1
  darr[ix] = -(2d0*arr[ix-3]-9d0*arr[ix-2]+18d0*arr[ix-1]-11d0*arr[ix])/6d0

  darr = darr/dx
  return, darr
end

