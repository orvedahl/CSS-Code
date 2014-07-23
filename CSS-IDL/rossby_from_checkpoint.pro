@vectorcurl.pro

pro rossbyfc, iter, css_case=css_case, omega=omega

  vectorcurl, iter, css_case=css_case, curlB=curl, rad=rn, theta=theta, phi=phi

  ens = sqrt(curl[*,*,*,0]^2+curl[*,*,*,1]^2+curl[*,*,*,2]^2)
  dr = rn[1]-rn[0]
  dth = theta[1]-theta[0]
  dphi = phi[1]-phi[0]

  nr = n_elements(rn)
  rossby = dblarr(nr)
  for ir=0,nr-1 do begin
     rossby[ir] = mean(ens[*,*,ir])
  endfor

  rossby = rossby/omega

  stop
end
