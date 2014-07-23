function kappa, k, a, b, c, alp, bet

  kap = (a*sin(k) + (b/2d0)*sin(2*k) + (c/3d0)*sin(3*k))/(1+2d0*alp*cos(k)+2d0*bet*cos(2*k))/!dpi

  return, kap

end

pro dispersion_fig

  filename='dispersion_fig.eps'

  xsize = 3.5d0
  ysize = xsize

  Set_Plot, "PS"
  !P.FONT=1
  Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, /TT_FONT, SET_FONT='Helvetica', /CMYK, /INCHES
  !P.FONT=1
  !P.MULTI=0
  !P.NOERASE=1
  !P.THICK=3
  !X.THICK=3
  !Y.THICK=3
  !P.CHARSIZE=1
  docolor
  choose_color_ben, /rain

  !P.POSITION=[0.18d0,0.15d0,0.96d0,0.96d0]

  k = dindgen(128)/128d0+1d0/256d0
  invk = 1d0/k
  plot, invk, k, /xs, /ys, xtitle='N!Dpts', ytitle=textoidl('|k-\kappa|k^{-1}'), /nodata, /xlog, /ylog, yrange=[1d-8,1d0]
  oplot, invk, abs(k-kappa(!dpi*k,4d0/3d0,-1d0/3d0,0d0,0d0,0d0))/k, color=40, linestyle=1 ;Fourth order FD
  oplot, invk, abs(k-kappa(!dpi*k,3d0/2d0,-3d0/5d0,1d0/10d0,0d0,0d0))/k, color=40 ;Sixth order FD
  oplot, invk, abs(k-kappa(!dpi*k,14d0/9d0,1d0/9d0,0d0,1d0/3d0,0d0))/k, color=220 ;Sixth order CFD
  oplot, invk, abs(k-kappa(!dpi*k,17d0/12d0,101d0/150d0,1d-2,1d0/2d0,1d0/20d0))/k, color=220, linestyle=1 ;Tenth order CFD

  DEVICE, /close
  print, 'Image ready.'
  SPAWN, 'evince '+filename+' &'
  close, /all
  stop
end
