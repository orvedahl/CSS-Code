pro css2_scaling
;ASH lmax=340, nr=200
ash340_cpusk = [128,256,512,1024,2040,3072]*1d0
ash340_timesk = [864.000,494.000,295.160,192.480,137.600,198.400]
ash340_itersk = 0d0*ash340_cpusk+400d0
ash340_iterhrk = ash340_itersk/ash340_timesk*3600d0
ash340_gridptsk = 340d0*340d0*200d0/ash340_cpusk

ash340_cpusr = [130,258,514,770,1026,1538,2050,3074,4096]
ash340_iterhrr = [6.2d2,1.1d3,2.1d3,2.8d3,3.5d3,4.8d3,6d3,7.6d3,7d3]*257d0/200d0
ash340_gridptsr = 340d0*340d0*257d0/ash340_cpusr

ash340_cpusp  = 1d0*[64,128,256, 512, 1024,2048,3048,4088]
ash340_itersp = [50.0,100.0,250.0,250.0,500.0,61389.0,40595.0,73828.0]
nt64   = 161.0
nt128  = 188.0
nt256  = 180.0+46.0
nt512  = 137.0
nt1024 = 132.0
nt2048 = 3600.0*2.0
nt3048 = 3600.0-60.0
nt4088 = 7200.0
ash340_timesp = [nt64, nt128,nt256,nt512,nt1024,nt2048,nt3048,nt4088]
ash340_iterhrp = ash340_itersp/ash340_timesp*3600d0
ash340_gridptsp = 340d0*340d0*200d0/ash340_cpusp

;ASH, lmax=680, nr=400

ash680_cpusk = [1024.00,2048.00,3072.00,4096.00,8192.00, 10240.0]
ash680_iterhrk = [0.294074,0.473934,0.582072,0.674309,0.989609,0.958313]*3600d0
ash680_gridptsk = 680d0*680d0*400d0/ash680_cpusk

ash680_cpusr = [1538,2050,3074,4096]
ash680_iterhrr = [1200d0,1200d0,1600d0,2300d0]
ash680_gridptsr = 680d0*680d0*257d0/ash680_cpusr

ash680_cpusp = 1d0*[1024,2048,4096,8192,10000]
ash680_timesp = [1154.8300,660.96997,372.25272,214.43140,188.22000]
ash680_itersp = 0d0*ash680_cpusp+500d0
ash680_iterhrp = ash680_itersp/ash680_timesp*3600d0
ash680_gridptsp = 680d0*680d0*400d0/ash680_cpusp

;nr = 192
cpus256k = 1d0*[12,24,48,96,192,384,768,1536,3072]
cpus512k = 1d0*[96,192,384,768,1536,3072,6144,12288]
cpus1024k = 1d0*[384,768,1536,3072,6144,12288,24576,49152]

iterhr256k = [599.01,1240.73,2478.68,4457.16,8420.64,17667.85,33379.22,65615.57,110584.01]*2d0
iterhr512k = [1202.45,2398.27,4322.97,8628.88,17037.16,32319.99,52641.83,77282.22]*2d0
iterhr1024k = [1201.31,2171.47,4262.59,8333.33,15492.77,31249.78,49720.10,68872.54]*2d0

gridpts256k = 1d0*[1048576,524288,262144,131072,65536,32768,16384,8192,4096]
gridpts512k = 1d0*[524288.00,262144.00,131072.00,65536.00,32768.00,16384.00,8192.00,4096.00]
gridpts1024k = 1d0*[524288.00,262144.00,131072.00,65536.00,32768.00,16384.00,8192.00,4096.00]

;nr = 576
cpus512p  = [768.00,1536.00,3072.00,6144.00,12288.00,18432.00,24576.00]
cpus1024p = [384.00,768.00,1536.00,3072.00,6144.00,12288.00,18432.00,24576.00,36864.00]
cpus2048p = [1536.00,3072.00,6144.00,12288.00,18432.00,24576.00,36864.00]

iterhr512p  = [5498.53,10900.88,21906.30,41982.33,87305.90,125257.49,142435.60]*2d0
iterhr1024p = [628.42,1328.89,2739.29,5480.31,10668.40,21526.41,31512.84,41688.19,63691.79]*2d0
iterhr2048p = [587.05,1138.75,2356.89,4937.56,7845.31,9882.60,14637.04]*2d0

gridpts512p  = [196608.00,98304.00,49152.00,24576.00,12288.00,8192.00,6144.00]
gridpts1024p = [1572864.00,786432.00,393216.00,196608.00,98304.00,49152.00,32768.00,24576.00,16384.00]
gridpts2048p = [1572864.00,786432.00,393216.00,196608.00,131072.00,98304.00,65536.00]

space=0.01d0    
xsize = 7.5d0*2.54d0
ysize = xsize/2.5d0
csize = 1.0
Set_Plot, "PS"
!P.MULTI=0
!P.FONT=0
!P.CHARSIZE=csize
Device, filename='ashfdncss2_scaling.eps', /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Times-Roman', /CMYK

choose_color, /rainbow

!P.NOERASE=1
!P.THICK=4
!X.THICK=2
!Y.THICK=2
dx = 0.39
x0 = 0.09
x1 = x0+dx
y0 = 0.20
y1 = 0.95
!P.POSITION=[x0,y0,x1,y1]


colors=[0,0,50,50,200,200]
names=['340P','680P','340K','680K','340R','680R']
syms=-[2,1,2,1,2,1]
stys = [0,0,0,0,0,0]
plot, ash340_cpusp, ash340_iterhrp, yrange=[1d3,1d5], xrange=[1d1,5d4], /xlog, /ylog, psym=syms[0], $
  ytitle='Timesteps Hr!U-1!N', xtitle='N (Cores)', color=colors[0], thick=4, charsize=csize, /ys, /xs, linestyle=stys[0]
oplot, ash680_cpusp, ash680_iterhrp, psym=syms[1], color=colors[1], thick=4, linestyle=stys[1]
oplot, ash340_cpusk, ash340_iterhrk, psym=syms[2], color=colors[2], thick=4, linestyle=stys[2]
oplot, ash680_cpusk, ash680_iterhrk, psym=syms[3], color=colors[3], thick=4, linestyle=stys[3]
oplot, ash340_cpusr, ash340_iterhrr, psym=syms[4], color=colors[4], thick=4, linestyle=stys[4]
oplot, ash680_cpusr, ash680_iterhrr, psym=syms[5], color=colors[5], thick=4, linestyle=stys[5]
legend,names,psym=syms,colors=colors, linestyle=stys, /top, /left, charsize=csize, symsize=0.75+0d0*stys
xyouts, /NORMAL, x0-0.07, 0.95, '(a)', charsize=csize,charthick=1,color=0

;x0 = x1+0.06d0
;x1 = x0 + dx
;!P.POSITION=[x0,y0,x1,y1]

;plot, ash340_gridptsp, ash340_iterhrp, yrange=[1d3,1d5], xrange=[1d3,1d6], /xlog, /ylog, psym=syms[0], $
;  ytitle='Iterations Hr!U-1!N', xtitle='Gridpoints Core!U-1!N', color=colors[0], thick=4, charsize=csize
;oplot, ash680_gridptsp, ash680_iterhrp, psym=syms[1], color=colors[1], thick=4
;oplot, ash340_gridptsk, ash340_iterhrk, psym=syms[2], color=colors[2], thick=4, linestyle=2
;oplot, ash680_gridptsk, ash680_iterhrk, psym=syms[3], color=colors[3], thick=4, linestyle=2
;legend,names,psym=syms,colors=colors, linestyle=stys, /top, /right
;xyouts, /NORMAL, x0-0.04, 0.95, 'B', charsize=csize2,charthick=1,color=0

x0 = x1+0.1d0
x1 = x0 + dx
!P.POSITION=[x0,y0,x1,y1]

colors=[50,50,50,200,200,200]
names=['512P','1024P','2048P','256K','512K','1024K']
syms=-[2,1,4,5,2,1]
stys = [0,0,0,0,0,0]
plot, cpus512p, iterhr512p, yrange=[1d3,1d6], xrange=[1d1,1d5], /xlog, /ylog, psym=syms[0], $
  ytitle='Timesteps Hr!U-1!N', xtitle='N (Cores)', color=0, thick=4, charsize=csize
oplot, cpus512p, iterhr512p, psym=syms[0], color=colors[0], thick=4, linestyle=stys[0]
oplot, cpus1024p, iterhr1024p, psym=syms[1], color=colors[1], thick=4, linestyle=stys[1]
oplot, cpus2048p, iterhr2048p, psym=syms[2], color=colors[2], thick=4, linestyle=stys[2]
oplot, cpus256k, iterhr256k, psym=syms[3], color=colors[3], thick=4, linestyle=stys[3]
oplot, cpus512k, iterhr512k, psym=syms[4], color=colors[4], thick=4, linestyle=stys[4]
oplot, cpus1024k, iterhr1024k, psym=syms[5], color=colors[5], thick=4, linestyle=stys[5]
legend,names,psym=syms,colors=colors, linestyle=stys, /top, /left, charsize=csize, symsize=0.75+0d0*stys
xyouts, /NORMAL, x0-0.07, 0.95, '(b)', charsize=csize,charthick=1,color=0

;x0 = x1+0.06d0
;x1 = x0 + dx
;!P.POSITION=[x0,y0,x1,y1]

;plot, gridpts512p, iterhr512p, yrange=[1d3,1d6], xrange=[1d3,1d7], /xlog, /ylog, psym=syms[0], $
;  ytitle='Iterations Hr!U-1!N', xtitle='Gridpoints Core!U-1!N', color=colors[0], thick=4, charsize=csize
;oplot, gridpts1024p, iterhr1024p, psym=syms[1], color=colors[1], thick=4
;oplot, gridpts2048p, iterhr2048p, psym=syms[2], color=colors[2], thick=4
;oplot, gridpts256k, iterhr256k, psym=syms[3], color=colors[3], thick=4, linestyle=2
;oplot, gridpts512k, iterhr512k, psym=syms[4], color=colors[4], thick=4, linestyle=2
;oplot, gridpts1024k, iterhr1024k, psym=syms[5], color=colors[5], thick=4, linestyle=2
;legend,names,psym=syms,colors=colors, linestyle=stys, /top, /right
;xyouts, /NORMAL, x0-0.04, 0.95, 'D', charsize=csize2,charthick=1,color=0

DEVICE, /close
SET_PLOT, 'x'
print, 'Image ready.'
SPAWN, 'evince ashfdncss2_scaling.eps &'

stop

end
