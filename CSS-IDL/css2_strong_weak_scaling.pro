pro css2_strong_weak_scaling

;nr = 192
;cpus256k = 1d0*[12,24,48,96,192,384,768,1536,3072]
;cpus512k = 1d0*[96,192,384,768,1536,3072,6144,12288]
;cpus1024k = 1d0*[384,768,1536,3072,6144,12288,24576,49152]

;iterhr256k = [599.01,1240.73,2478.68,4457.16,8420.64,17667.85,33379.22,65615.57,110584.01]*2d0
;iterhr512k = [1202.45,2398.27,4322.97,8628.88,17037.16,32319.99,52641.83,77282.22]*2d0
;iterhr1024k = [1201.31,2171.47,4262.59,8333.33,15492.77,31249.78,49720.10,68872.54]*2d0

;gridpts256k = 1d0*[1048576,524288,262144,131072,65536,32768,16384,8192,4096]
;gridpts512k = 1d0*[524288.00,262144.00,131072.00,65536.00,32768.00,16384.00,8192.00,4096.00]
;gridpts1024k = 1d0*[524288.00,262144.00,131072.00,65536.00,32768.00,16384.00,8192.00,4096.00]

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
Device, filename='css_scaling.eps', /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize, SET_FONT='Times-Roman', /CMYK

choose_color, /rainbow

!P.NOERASE=1
!P.THICK=4
!X.THICK=4
!Y.THICK=4
dx = 0.39
x0 = 0.10
x1 = x0+dx
y0 = 0.16
y1 = 0.95
!P.POSITION=[x0,y0,x1,y1]

colors = [50,90,200]
names  = ['512','1024','2048']

syms = -[4,2,1]
stys =  [0,0,0]

ncpus = 1d2 + (1d5-1d2)*dindgen(128)/127d0
ideal512 = 3600d0*cpus512p[0]/iterhr512p[0]/ncpus
ideal1024 = 3600d0*cpus1024p[0]/iterhr1024p[0]/ncpus
ideal2048 = 3600d0*cpus2048p[0]/iterhr2048p[0]/ncpus

plot, cpus512p, 3600d0/iterhr512p, yrange=[1d-2,10d0], xrange=[1d2,1d5], /xlog, /ylog, psym=syms[0], $
  ytitle='Elapsed Time (s!E !N/!E !Nstep)', xtitle='N (cores)', color=0, thick=4, charsize=csize
oplot, cpus512p, 3600d0/iterhr512p, psym=syms[0], color=colors[0], thick=4, linestyle=stys[0]
oplot, ncpus, ideal512, linestyle=2, thick=2
oplot, cpus1024p, 3600d0/iterhr1024p, psym=syms[1], color=colors[1], thick=4, linestyle=stys[1]
oplot, ncpus, ideal1024, linestyle=2, thick=2
oplot, cpus2048p, 3600d0/iterhr2048p, psym=syms[2], color=colors[2], thick=4, linestyle=stys[2]
oplot, ncpus, ideal2048, linestyle=2, thick=2

oplot, [1d2,1d5], 1d0+[0d0,0d0], linestyle=3, thick=1
oplot, [1d2,1d5], 0.1d0+[0d0,0d0], linestyle=3, thick=1

oplot, 1d3+[0d0,0d0], [0.01d0,10d0], linestyle=3, thick=1
oplot, 1d4+[0d0,0d0], [0.01d0,10d0], linestyle=3, thick=1

legend,names,psym=syms,colors=colors, linestyle=stys, /bottom, /left, charsize=0.8*csize, symsize=0.8+0d0*stys
xyouts, /NORMAL, x0-0.09, 0.95, '(a)', charsize=csize,charthick=1,color=0

x0 = x1+0.1d0
x1 = x0 + dx
!P.POSITION=[x0,y0,x1,y1]

int512  = 3600d0*cpus512p[0]/iterhr512p[0]/cpus512p
int1024 = 3600d0*cpus1024p[0]/iterhr1024p[0]/cpus1024p
int2048 = 3600d0*cpus2048p[0]/iterhr2048p[0]/cpus2048p

plot, gridpts512p, (3600d0/iterhr512p-int512)/int512, yrange=[-0.15d0,0.25d0], xrange=[1d3,1d7], /xlog, psym=syms[0], $
  ytitle='Relative Difference', xtitle='Gridpoints!E !N/!E !NCore', color=0, thick=4, charsize=csize, /ys
oplot, gridpts512p, (3600d0/iterhr512p-int512)/int512, psym=syms[0], color=colors[0], thick=4
oplot, gridpts1024p, (3600d0/iterhr1024p-int1024)/int1024, psym=syms[1], color=colors[1], thick=4
oplot, gridpts2048p, (3600d0/iterhr2048p-int2048)/int2048, psym=syms[2], color=colors[2], thick=4

oplot, [1d3,1d7], [0d0,0d0], linestyle=3, thick=1
oplot, 1d4+[0d0,0d0], [-1d0,1d0], linestyle=3, thick=1
oplot, 1d5+[0d0,0d0], [-1d0,1d0], linestyle=3, thick=1
oplot, 1d6+[0d0,0d0], [-1d0,1d0], linestyle=3, thick=1

legend,names,psym=syms,colors=colors, linestyle=stys, /top, /right, charsize=0.8*csize, symsize=0.8+0d0*stys
xyouts, /NORMAL, x0-0.075, 0.95, '(b)', charsize=csize,charthick=1,color=0

DEVICE, /close
SET_PLOT, 'x'
print, 'Image ready.'
SPAWN, 'evince css_scaling.eps &'

stop

end
