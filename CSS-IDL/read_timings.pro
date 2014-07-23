pro read_timings, css_case=css_case

  dir = '/nobackupp6/kaugusts/'+css_case

  openr, fun, dir+'/timings.dat', /get_lun

  ;Read first line for number of processors
  temp = 0d0
  readu, fun, temp

  nproc = fix(temp(0))
  
  print, 'nproc=', nproc
  close, fun
  free_lun, fun

  openr, fun, dir+'/timings.dat', /get_lun
  a = assoc(fun,dblarr(17,nproc+1))

  timings = a[0]
  timings = timings[*,1:nproc]
  close, fun
  free_lun, fun

stop
end
