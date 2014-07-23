 pro css_perf

   pleiades256_ncpu = [96,192,384,384,768,768,1536,1536,3072]
   pleiades512_ncpu = [96,192,384,384,768,768,1536,1536,3072,3072,6144]
   pleiades1024_ncpu = [192,384,384,768,768,1536,1536,3072,3072,6144,3072,6144,6144,12288]

   pleiades256_ithr = [4649.56,9000.94,17277.01,16861.97,32763.38,32695.17,59887.61,61771.81,101408.73]
   pleiades512_ithr = [1146.02,2132.37,4425.83,4065.86,8707.15,8147.76,16627.66,15771.54,29209.99,28688.47,56272.31]
   pleiades1024_ithr = [424.58,984.54,838.3,2075.39,1861.03,4169.37,3943.14,7443.47,7572.34,13736.6,7733.82,14221.18,12789.64,24578.53]

   kraken256_ncpu = [96,192,384,384,768,768,1536,1536,3072]
   kraken512_ncpu = [96,192,384,768,1152,2304,4608,9216]
   kraken1024_ncpu = [192,384,768,1536,3072,6144]

   kraken256_ithr = [2300,4911.43,8831.87,8655,14370.07,13672.83,25376.07,28178.49,41516.4]
   kraken512_ithr = [885.18,1603.96,3464.22,6504.46,8728.42,16100.59,28841.67,32909.5]
   kraken1024_ithr = [415.74,848.12,1560.52,3039.06,5259.65,7182.28]

   Set_Plot, 'PS'
   xsize=13.33
   ysize=10.00
   filename = 'css_perf.eps'
   Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, /INCHES, XSIZE=xsize, YSIZE=ysize
   !P.FONT=1
   !P.MULTI=0
   !P.NOERASE=0
   !P.POSITION=[0.135,0.17,0.92,0.92]
   !P.CHARSIZE=3
   !P.CHARTHICK=1.5
   docolor
   choose_color, /rainbow, /quiet

   plot, pleiades1024_ncpu, pleiades1024_ithr, xtitle='N (Cores)', ytitle='Time Steps per Wall Hour', $
     thick=1, psym=0, /xlog, /ylog, symsize=2, yrange=[2d2,2d5], xrange=[50,20000], /xs, /ys, xthick=8, ythick=8, $
     ytickname=['1,000','10,000','100,000'], xtickname=['100','1,000','10,000']
   oplot, pleiades256_ncpu, pleiades256_ithr, psym=-1, symsize=2, thick=8, color=20
   oplot, pleiades512_ncpu, pleiades512_ithr, psym=-2, symsize=2, thick=8, color=60
   oplot, pleiades1024_ncpu, pleiades1024_ithr, psym=-4, symsize=2, thick=8, color=100

   oplot, kraken256_ncpu, kraken256_ithr, psym=-5, symsize=2, thick=8, color=160
   oplot, kraken512_ncpu, kraken512_ithr, psym=-6, symsize=2, thick=8, color=200
   oplot, kraken1024_ncpu, kraken1024_ithr, psym=-7, symsize=2, thick=8, color=0

   names = ['Pleiades 256', 'Pleiades 512', 'Pleiades 1024', 'Kraken 256', 'Kraken 512', 'Kraken 1024']
   syms  = -[1,2,4,5,6,7]
   colors = [20,60,100,160,200,0]
   legend, names, color=colors, psym=syms, thick=8*[1,1,1,1,1,1], charsize=2.5

   xyouts, 0.01, 0.94, 'B', /NORMAL, charsize=5
   xyouts, 0.1175, 0.1075, '50', /NORMAL
   xyouts, 0.345, 0.1075, '300', /NORMAL
   xyouts, 0.64, 0.1075, '3,000', /NORMAL
   xyouts, 0.89, 0.1075, '20,000', /NORMAL

   DEVICE, /close
   SET_PLOT, 'x'
   print, 'Image ready.'
   SPAWN, 'kghostview '+filename+' &'

end
