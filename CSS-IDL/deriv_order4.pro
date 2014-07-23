function do4, dr, x, b1, b2, bdry, dtype

  mynr = n_elements(x)
  dxdy = 0d0*x
  coefs = dblarr(5)
  dri = 1d0/dr

  If ((dtype eq 2)or(dtype eq 0)) Then Begin
      If ((bdry mod 2) eq 0) Then Begin
          dxdy[0] = 0d0
          coefs = [-19d0,6d0,15d0,-2d0,0d0]
          coefs = coefs/24d0
          dxdy[1] = coefs[0]*x[0] + coefs[1]*x[1] + coefs[2]*x[2] + coefs[3]*x[3]
          If (bdry gt 0) Then Begin
              dxdy[0] = dr*b1[1]
              dxdy[1] = dxdy[1] - 0.25d0*dr*b1[1]
          EndIf 
      EndIf Else Begin
          coefs = [-25d0,48d0,-36d0,16d0,-3d0]
          coefs = coefs/12d0
          dxdy[0] = coefs[0]*x[0] + coefs[1]*x[1] + coefs[2]*x[2] + coefs[3]*x[3] + coefs[4]*x[4]
          coefs = [-3d0,-10d0,18d0,-6d0,1d0]
          coefs = coefs/12d0
          dxdy[1] = coefs[0]*x[0] + coefs[1]*x[1] + coefs[2]*x[2] + coefs[3]*x[3] + coefs[4]*x[4]
      EndElse
          
      coefs = [1d0,-8d0,8d0,-1d0,0d0]
      coefs = coefs/12d0
  Endif Else Begin
      coefs = [1d0,-8d0,8d0,-1d0,0d0]
      coefs = coefs/12d0 
      dxdy[0] = coefs[0]*b1[0] + coefs[1]*b1[1] + coefs[2]*x[1] + coefs[3]*x[2]
      dxdy[1] = coefs[0]*b1[1] + coefs[1]*x[0]  + coefs[2]*x[2] + coefs[3]*x[3]
  EndElse
  
  for j=2,mynr-3 do begin
     dxdy[j] = coefs[0]*x[j-2] + coefs[1]*x[j-1] + coefs[2]*x[j+1] + coefs[3]*x[j+2]
  Endfor
  
  If ((dtype eq 1)or(dtype eq 0)) Then Begin
      If ((bdry eq 0)or(bdry eq 3)or(bdry eq 4)) Then Begin
          coefs = [2d0,-15d0,-6d0,19d0,0d0]
          coefs = coefs/24d0
          dxdy[mynr-2] = coefs[0]*x[mynr-4] + coefs[1]*x[mynr-3] + coefs[2]*x[mynr-2] + coefs[3]*x[mynr-1]
          dxdy[mynr-1] = 0d0
          If ((bdry eq 3)or(bdry eq 4)) Then Begin
              dxdy[mynr-2] = dxdy[mynr-2] - 0.25d0*dr*b2[0]
              dxdy[mynr-1] = dr*b2[0]
          Endif
      EndIf Else Begin
          coefs = [-1d0,6d0,-18d0,10d0,3d0]
          coefs = coefs/12d0
          dxdy[mynr-2] = coefs[0]*x[mynr-5] + coefs[1]*x[mynr-4] + coefs[2]*x[mynr-3] + coefs[3]*x[mynr-2] + coefs[4]*x[mynr-1]
          coefs = [3d0,-16d0,36d0,-48d0,25d0]
          coefs = coefs/12d0
          dxdy[mynr-1]   = coefs[0]*x[mynr-5] + coefs[1]*x[mynr-4] + coefs[2]*x[mynr-3] + coefs[3]*x[mynr-2] + coefs[4]*x[mynr-1]
      EndElse
  Endif Else Begin
      dxdy[mynr-2] = coefs[0]*x[mynr-4] + coefs[1]*x[mynr-3] + coefs[2]*x[mynr-1] + coefs[3]*b2[0]
      dxdy[mynr-1]   = coefs[0]*x[mynr-3] + coefs[1]*x[mynr-2] + coefs[2]*b2[0]   + coefs[3]*b2[1]
  EndElse
  
  dxdy = dri*dxdy
  return, dxdy
end
