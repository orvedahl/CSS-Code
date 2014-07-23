@display_shell_slice.pro
@choose_color.pro
pro spd_fig3

;Read in the two cooling_function files

nr  = 48
nth = 512
nph = 1024
th1 = 60d0*!dpi/180d0
th2 = 120d0*!dpi/180d0
ph1 = 0d0*!dpi/180d0
ph2 = 120d0*!dpi/180d0
r1  = 9d0
r2  = 10d0

rsun = 6.96d10
lsun = 3.846d33

r = (r2-r1)*dindgen(nr)/double(nr-1)+r1
r = r*rsun/r2
theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1

t1 = 60d0
t2 = 96d0
sigmar = 0.03*(max(r)-r(0))

cft1 = dblarr(nr,nth,nph)
cft2 = dblarr(nr,nth,nph)

openr,file_unit,'cooling_function_t1.dat',/get_lun
a = assoc(file_unit,dblarr(nr,nth,nph))
cft1 = a(0)
close,file_unit
free_lun, file_unit

openr,file_unit,'cooling_function_t2.dat',/get_lun
a = assoc(file_unit,dblarr(nr,nth,nph))
cft2 = a(0)
close,file_unit
free_lun, file_unit

;Average in theta and phi
cft1_avg = dblarr(nr)
cft2_avg = dblarr(nr)

for ir=0L,nr-1L do begin
    cft1_avg(ir) = mean(reform(cft1(ir,*,*)))
    cft2_avg(ir) = mean(reform(cft2(ir,*,*)))
endfor

;True time averaged cf_avg
nr=256
r1 = 9.8
rn = (r2-r1)*dindgen(nr)/double(nr-1)+r1
rn = rn*rsun/r2

cft1_avg = interpol(cft1_avg,r,rn,/SPLINE)
cft2_avg = interpol(cft2_avg,r,rn,/SPLINE)
cf_avg = lsun/5d0/!dpi/rn^3*exp(-(rn-rsun)^2/sigmar^2/2d0)/sqrt(2d0*!dpi)/sigmar

xsize=8.5*2.54
ysize=4.25d0*2.54

!P.FONT=3

filename = 'spd_fig3.eps'
Set_Plot, "PS"
choose_color, /ash_color, /quiet
Device, filename=filename, /ENCAPSULATED, BITS_PER_PIXEL=8, /COLOR, XSIZE=xsize, YSIZE=ysize

!P.MULTI=0
!P.NOERASE=1
!P.POSITION=[0.01,0.1,0.57,0.9]
draw = bytscl(transpose(alog10(reform(cft1(47,*,*)))))
    
lon = phi*180d0/!dpi
lat = reverse(90d0-theta*!radeg)

limits=[-30d0,0d0,30d0,120d0]

mnv = 0.0d0
mxv = 0.0d0
nlevs=15                        ;  number of contour levels
levs=get_shell_slice_levels(alog10(cft1),lat,nlevs,mn=mnv,mx=mxv,pct=0.04,/uneven_range)
display_shell_slice, lat, lon, draw, /publication, /grid, minv=mnv, maxv=mxv, res_scale=0.3, $
         /white_bg, mollweide=[60d0,0], latlon_window=limits, /psf


!P.NOERASE=1
!P.POSITION=[0.65,0.20,0.95,0.85]
plot, rn/rsun, 1d9*cf_avg, xtitle='R/R'+sunsymbol(), ytitle=textoidl(' \epsilon 10^{-9} erg cm^{-3} s^{-1}'),thick=3
oplot, rn/rsun, 1d9*cft1_avg, color=50, thick=3
oplot, rn/rsun, 1d9*cft2_avg, color=120, thick=3

DEVICE, /close
SET_PLOT, 'x'
print, 'Image ready.'
SPAWN, 'kghostview '+filename+' &'

stop

end
