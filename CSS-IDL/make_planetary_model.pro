pro make_planetary_model, bounds, constg=constg, make_plot=make_plot

nr = fix(bounds(0))
nth = fix(bounds(1))
nph = fix(bounds(2))
th1 = fix(bounds(3))
th2 = fix(bounds(4))
ph1 = fix(bounds(5))
ph2 = fix(bounds(6))

;Create pressure array
lnP0 = alog10(5d2) ;bar
lnPs = alog10(1d-5) ;bar

nlnP = 2048L
lnP = (lnPs-lnP0)*dindgen(nlnP)/double(nlnP-1L)+lnP0

P = 10d0^(lnP)

gs = 9.42d2
Rp = 9.44d9
G = 6.67d-8
Mp = Rp^2*gs/G

gamma = 1.473d0
Cp = 1.431d8
Cv = Cp/gamma
gam1 = gamma-1d0

;Dayside at theta=90, phi=180 maximum T (presumably)
aday = [2149.9581,4.1395571,-186.24851,135.52524,106.20433,-35.851966,-50.022826,-18.462489,-3.3319965,-0.30295925,-0.011122316]
ndayf = n_elements(aday)
Tday = dblarr(nlnP)
Tday(*) = aday(0)
for in=1,ndayf-1 do begin
    Tday = Tday+aday(in)*lnP^in
endfor

anight = [1388.2145,267.66586,-215.53357,61.814807,135.68661,2.0149044,-40.907246,-19.015628,-3.8771634,-0.38413901,-0.015089084]
nnf = n_elements(anight)
Tnight = dblarr(nlnP)
Tnight(*) = anight(0)
for in=1,nnf-1 do begin
    Tnight = Tnight+anight(in)*lnP^in
endfor

;ain = [5529.7168,-6869.6504,4142.7231,-936.23053,87.120975]
ain = [4314.4632,-14966.974,34564.817,-41540.983,28526.979,-11542.497,2734.4161,-351.28564,18.915705]
ninf = n_elements(ain)
Tin = dblarr(nlnP)
Tin(*) = ain(0)
for in=1,ninf-1 do begin
    Tin = Tin+ain(in)*lnP^in
endfor

Td = dblarr(nlnP)
Tn = dblarr(nlnP)

ind1 = where((P le 10d0),complement=ind2)
Td(ind1) = Tday(ind1)
Td(ind2) = Tin(ind2)
Tn(ind1) = Tnight(ind1)
Tn(ind2) = Tin(ind2)

smoo = 40L
Tds = smooth(Td,smoo)
Tds = smooth(Tds,smoo)
Tds = smooth(Tds,smoo)

Tns = smooth(Tn,smoo)
Tns = smooth(Tns,smoo)
Tns = smooth(Tns,smoo)


intrd = dblarr(nlnP)
intrd(0) = 0d0
lnPr = reverse(lnP)
Tdsr = reverse(Tds)
for in=1L,nlnP-1L do begin
    intrd(in) = int_tabulated(lnPr(0L:in),Tdsr(0L:in),/DOUBLE)
endfor

intrd = reverse(intrd)

intrn = dblarr(nlnP)
intrn(0) = 0d0
lnPr = reverse(lnP)
Tnsr = reverse(Tns)
for in=1L,nlnP-1L do begin
    intrn(in) = int_tabulated(lnPr(0L:in),Tnsr(0L:in),/DOUBLE)
endfor

intrn = reverse(intrn)

If (keyword_set(constg)) Then Begin
    rd = Rp-gam1*Cv*alog(10d0)*intrd/gs
    rn = Rp-gam1*Cv*alog(10d0)*intrn/gs
    rn = rn-rn(0)+rd(0)
    mr = (rn(nlnP-1)+Rp)/2d0
    drd = Rp-mr
    rd = rd+drd
    rn = rn+drd
Endif Else Begin
    rd = Rp/(1d0+gam1*Cv*alog(10d0)*Rp*intrd/G/Mp)
    rn = Rp/(1d0+gam1*Cv*alog(10d0)*Rp*intrn/G/Mp)
    rn = rn-rn(0)+rd(0)
    mr = (rn(nlnP-1)+Rp)/2d0
    drd = Rp-mr
    rd = rd+drd
    rn = rn+drd
EndElse

P_cgs = P*1d6

rhod = P_cgs/gam1/Cv/Tds
rhon = P_cgs/gam1/Cv/Tns

Sd = Cv*alog(gam1*Cv*Tds/rhod^gam1)
Sn = Cv*alog(gam1*Cv*Tns/rhon^gam1)

If (keyword_set(make_plot)) Then Begin
    set_plot,'ps'
    device,filename='hd_model_T_logP.eps',/encapsul
    plot, Tns, lnP, ytitle=textoidl('log_{10} P (bar)'), xtitle=textoidl('T (K)'), yrange=[lnP(0),lnP(nlnP-1)], /xs
    oplot, Tds, lnP
    device,/close 
    set_plot,'x'
    spawn, 'kghostview hd_model_T_logP.eps &'

    set_plot,'ps'
    device,filename='hd_model_r_S.eps',/encapsul
    plot, rd, Sd, xtitle=textoidl('r (cm)'), ytitle=textoidl('S (erg/K/cm^3)'), yrange=[min(Sn),max(Sd)], /xs
    oplot, rn, Sn
    device,/close 
    set_plot,'x'
    spawn, 'kghostview hd_model_r_S.eps &'

    set_plot,'ps'
    device,filename='hd_model_r_rho.eps',/encapsul
    plot, rd, alog10(rhod), xtitle=textoidl('r (cm)'), ytitle=textoidl('log_{10} \rho (g/cm^3)'), /xs
    oplot, rn, alog10(rhon)
    device,/close 
    set_plot,'x'
    spawn, 'kghostview hd_model_r_rho.eps &'
EndIf

;Build initial CSS model

stop

end
