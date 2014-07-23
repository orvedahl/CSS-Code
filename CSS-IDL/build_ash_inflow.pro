@read_shell_ash.pro
@cheby_trans.pro
pro build_ash_inflow, outfile=outfile

  read_shell_ash, '/alphacen/augustso/solar_p4_170/Shell_Slices/2474784', [1,2,3,4,28], 600, [0], data, costheta, times=times, iters=iters, radii=radii, quan=quan, title=title
  
  nth = n_elements(costheta)
  nph = n_elements(data[0,*,0,0,0])
  
  ashtheta = acos(costheta)
  phi=shift((dindgen(nphi+1))*90d0/double(nphi),1)  &  phi[0]=0d0
  phi=phi+.5*(phi[0]+phi[1])-45d0
  phi=[phi, phi[0]+90d0] ;  replicate first index (for smooth map wrapping)
  ashphi=phi

 ;Put the data in CSS form
  reddata = fltarr(nth,nph,5,590)
  reddata[*,*,0,*] = data[*,*,3,9:*] ;Get rid of first 10 due to ringing left over from lmax170 case
  reddata[*,*,1:3,*] = data[*,*,0:2,9:*]
  reddata[*,*,4,*] = data[*,*,4,9:*]
  
 ;Free memory
  data = 0
  
 ;Apply CSS theta boundary conditions, save for future implementation?
  
 ;Write out header
  ;nth
  ;nphi
  niter=590
  ;ashtheta
  ;ashphi

  ;Write out data
  
end


