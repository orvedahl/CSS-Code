@read_checkpoint.pro
@write_checkpoint.pro

;Replicate a checkpoint to expand phi domain
pro stitch_checkpoint, iter, css_case

    read_checkpoint, iter, css_case=css_case, data=data, rad=rad, theta=theta, phi=phi, time=time

    ;Stitch
    nv = n_elements(data(0,0,0,*))
    nr = n_elements(rad)
    nth = n_elements(theta)
    nph = n_elements(phi)

    np2 = 2*nph
    data2 = dblarr(nr,nth,np2,nv)
    for iv=0,nv-1 do begin
        data2(*,*,0:nph-1,iv) = data(*,*,*,iv)
        data2(*,*,nph:np2-1,iv) = data(*,*,*,iv)
    endfor

    ph1 = min(phi)
    ph2 = max(phi)
    phi2 = 2d0*(ph2-ph1)*dindgen(np2)/double(np2-1)+ph1

    ;Write out new checkpoint and header
    write_checkpoint,iter,css_case=css_case,data=data2,rad=rad,theta=theta,phi=phi2,time=time,tag='_stitch'

end
