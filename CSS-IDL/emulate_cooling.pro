@view_shell2.pro
pro emulate_cooling
common cooling, phiseed, thetaseed, fluxseed, timeseed, phi, theta, r
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

r = (r2-r1)*dindgen(nr)/double(nr-1)+r1
r = r*rsun/r2
theta = (th2-th1)*dindgen(nth)/double(nth-1)+th1
phi = (ph2-ph1)*dindgen(nph)/double(nph-1)+ph1

costheta = cos(theta)

tmax = 120d0 ;hours
ntime  = 10000L
time = tmax*dindgen(ntime)/double(ntime-1L)
dt = 2d0

sgs = initialize_cooling(tmax)
nsgs = n_elements(sgs(*,0))
sigmar = 0.03*(max(r)-r(0))

eps = dblarr(nth,nph)
it = ntime/2L
eps(*,*) = 0d0
for is=0L,nsgs-1L do begin
    th0=sgs(is,1)
    ph0=sgs(is,0)
    t0=sgs(is,5)
    deltat=sgs(is,4)
    sigmas = sgs(is,2)
    eps = eps + sgs(is,3)*time_envelope(time(it),t0,deltat,dt)*spatial_envelope(phi,ph0,theta,th0,sigmas)
    print, 100d0*is/double(nsgs),'% done'
endfor

epstot = dblarr(nr,nth,nph)
for ir=0L,nr-1L do begin
    epstot(ir,*,*) = 2d0*eps(*,*)/r(ir)/sqrt(2d0*!dpi)/sigmar*exp(-(r(ir)-rsun)^2/sigmar^2/2d0)
endfor

;Save into file
openw, file_unit, 'cooling_function_t1.dat', /get_lun
a = assoc(file_unit, dblarr(nr,nth,nph))
a(0) = epstot
close,file_unit
free_lun, file_unit

eps = dblarr(nth,nph)
it = 4L*ntime/5L
eps(*,*) = 0d0
for is=0L,nsgs-1L do begin
    th0=sgs(is,1)
    ph0=sgs(is,0)
    t0=sgs(is,5)
    deltat=sgs(is,4)
    sigmas = sgs(is,2)
    eps = eps + sgs(is,3)*time_envelope(time(it),t0,deltat,dt)*spatial_envelope(phi,ph0,theta,th0,sigmas)
    print, 100d0*is/double(nsgs),'% done'
endfor

epstot = dblarr(nr,nth,nph)
for ir=0L,nr-1L do begin
    epstot(ir,*,*) = 2d0*eps(*,*)/r(ir)/sqrt(2d0*!dpi)/sigmar*exp(-(r(ir)-rsun)^2/sigmar^2/2d0)
endfor

;Save into file
openw, file_unit, 'cooling_function_t2.dat', /get_lun
a = assoc(file_unit, dblarr(nr,nth,nph))
a(0) = epstot
close,file_unit
free_lun, file_unit

end

function initialize_cooling, tmax
     common cooling, phiseed, thetaseed, fluxseed, timeseed, phi, theta, r
     ;Generate the random sequences for placement (uniformally distributed)
     phiseed   = 100101001L
     thetaseed = 201030401L
     fluxseed  = 567839121L
     timeseed  = 874359101L

     ph1 = phi(0)
     ph2 = max(phi)
     th1 = theta(0)
     th2 = max(theta)
     rsun = 6.96d10
     Lsun = 3.86d33
     nth = n_elements(theta)
     nph = n_elements(phi)
     area0 = (ph2-ph1)*(cos(th1)-cos(th2))*rsun^2
     area0 = area0/double(nph/5L)
     nsg = 4d0*!dpi*rsun^2/area0

     ;nsg = (nth/10L)*(nph/10L) ;nsg*(ph2-ph1)*(cos(th1)-cos(th2))/4d0/!dpi
     
     Q0 = Lsun/4d0/!dpi/rsun^2/nsg
     deltat = 20d0 ;mean lifetime
     size0   = sqrt(area0)/2d0/!dpi/rsun ;mean angular size
     ntot = nsg*tmax/deltat
     dt = 2d0

     phiseq = (ph2-ph1)*randomn(phiseed,ntot,/UNIFORM,/DOUBLE)+ph1
     thetaseq = (th2-th1)*randomn(phiseed,ntot,/UNIFORM,/DOUBLE)+th1

     ;Generate the random sequence for flux and duration (normally distributed)
     fluxes = 2d0*Q0*randomn(fluxseed,ntot,/UNIFORM,/DOUBLE)
     sizes  = 2d0*size0*randomn(fluxseed,ntot,/UNIFORM,/DOUBLE)
     lifetimes = 2d0*deltat*randomn(timeseed,ntot,/UNIFORM,/DOUBLE)
     starttimes = tmax*randomu(timeseed,ntot,/UNIFORM,/DOUBLE)
     
     sgs = dblarr(ntot,6)
     for is=0L,ntot-1L do begin
         sgs(is,0) = phiseq(is)
         sgs(is,1) = thetaseq(is)
         sgs(is,2) = sizes(is)
         sgs(is,3) = fluxes(is)
         sgs(is,4) = lifetimes(is)
         sgs(is,5) = starttimes(is)
     endfor
     return, sgs
end

function spatial_envelope, phi, phi0, theta, theta0, sigma

    np = n_elements(phi)
    nt = n_elements(theta)
    func = dblarr(nt,np)
    for ip=0L,np-1L do begin
        func(*,ip) = exp(-((phi(ip)-phi0)^2+(theta-theta0)^2)/sigma^2/2d0)/!dpi/2d0/sigma^2
    endfor

    return, func

end

function time_envelope, t, t0, deltat, dt

     env = 0d0
     If ((t ge t0) and (t le  t0+dt)) Then Begin
         env = sin(!dpi*(t-t0)/dt/2d0)^2
     EndIf Else If ((t ge t0+deltat-dt) and (t le t0+deltat)) Then Begin
         env = sin(!dpi*(t-t0-deltat)/dt/2d0)^2
     EndIf Else If ((t ge t0+dt) and (t le t0+deltat-dt)) Then Begin
         env = 1d0
     EndIf Else Begin
         env = 0d0
     EndElse
     return, env
end
