@read_checkpoint.pro

function compute_diffusion_speed, css_case=css_case, iter, Temperature_Top, gamma, Cp, ss_factor, mag=mag

  read_checkpoint, iter, css_case=css_case, data=data, rad=rad, theta=theta, phi=phi, time=time

  nr  = n_elements(rad)
  nth = n_elements(theta)
  nph = n_elements(phi)

  ;Compute averages
  meanV = dblarr(nth,nr,3)
  meanB = dblarr(nth,nr,3)

  for iv=0,2 do begin
     meanV[*,*,iv] = total(data[*,*,*,iv+1],2)/double(nph)
  endfor

  If (keyword_set(mag)) Then Begin
     for iv=0,2 do begin
        meanB[*,*,iv] = total(data[*,*,*,iv+5],2)/double(nph)
     endfor
  Endif

  cst = ss_factor*(gamma-1d0)*Cp*Temperature_Top
  diffusion_speed = dblarr(nth,nph,nr)
  for ip=0,nph-1 do begin
     diffusion_speed[*,ip,*] = cst + (data[*,ip,*,1]-meanV[*,*,0])^2 + (data[*,ip,*,2]-meanV[*,*,1])^2 + (data[*,ip,*,3]-meanV[*,*,2])^2
  endfor

  If (keyword_set(mag)) Then Begin
     diffusion_speed[*,ip,*] = diffusion_speed[*,ip,*] + (data[*,ip,*,5]-meanB[*,*,0])^2 + (data[*,ip,*,6]-meanB[*,*,1])^2 + (data[*,ip,*,7]-meanB[*,*,2])^2
  EndIf

  diffusion_speed = sqrt(diffusion_speed)

  return, diffusion_speed

end
