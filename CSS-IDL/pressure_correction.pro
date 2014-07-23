pro pressure_correction, nth, th1, th2
  omega_0 = 2.666d-6
  r1 = 6.195d10
  rhor1 = 2.8996d-2
  th1 = !dpi*th1/180d0
  th2 = !dpi*th2/180d0

  theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
  sines = sin(theta)
  cosines = cos(theta)

  cv = 467.425d0
  av = -51.25d0*cv/455.425d0
  bv = -81.525d0*cv/455.425d0
  cv = 2d0*!dpi*cv/1d9-omega_0
  av = 2d0*!dpi*av/1d9
  bv = 2d0*!dpi*bv/1d9
  vphi_prof = cv+av*cosines^2+bv*cosines^4 ;at 0.89
  vphi_prof = r1*vphi_prof

  pprime = -0.5d0*omega_0^2*cosines^2 + (0.375d0*av^2 + 29d0*av*bv/32d0 + 65d0*bv^2/128d0 + 0.5d0*av*cv + 0.75d0*bv*cv)*cos(2d0*theta)
  pprime = pprime + (av^2/32d0 + 0.125d0*av*bv + 23d0*bv^2/256d0 + bv*cv/16d0)*cos(4d0*theta) + (av*bv/96d0 + 5d0*bv^2/384d0)*cos(6d0*theta)
  pprime = pprime + bv^2*cos(8d0*theta)/1024d0 + (av+bv+cv)^2*alog(sines) + omega_0*(1.5d0*av + 1.25d0*bv + 2d0*cv)*sines 
  pprime = pprime + omega_0*(av/6d0 + 5d0*bv/24d0)*sin(3d0*theta) + omega_0*bv*sin(5d0*theta)/40d0
  pprime = r1^2*rhor1*pprime

  gamma = 1.6595d0
  Cp = 3.3d8
  sr1 = 6.853897056d9
  pr1 = rhor1^gamma*exp(sr1*gamma/Cp)

  ;pprime = pr1 pprime

  G = 6.67d-8
  g1 = G*1.98d33/r1^2

  eprime = (cv^2+omega_0^2)*cosines^2+av*(cv+omega_0)*cosines^4 + (av^2 + 2d0*bv*(cv+omega_0))*cosines^6/3d0
  eprime = eprime + 0.5d0*av*bv*cosines^8 + 0.2d0*bv^2*cosines^10 + cv*omega_0*cos(2d0*theta)
  eprime = eprime + sines^2*(cv+omega_0+av*cosines^2 + bv*cosines^4)^2
  eprime = -r1*Cp*eprime/g1
  
  ;eprime = sr1+eprime-eprime[3764]



stop
end
