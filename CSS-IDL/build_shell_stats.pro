pro build_shell_stats, statsfile

  root_directory = '/freyr3/augustso/CSSDE/'

  ;Open a file dialog box to select (a) file(s).
  
  files = DIALOG_PICKFILE(PATH=root_directory,TITLE='Select Shell Slice File to Read', /MULTI)
 
    ;Determine the number of files.
  n_files = n_elements(files)

  print, 'loading file:',files
                                ;Handle the event that the user cancels.
  if(n_files eq 1) then begin
     if(files eq '') then return
  endif

                                ;We assume the files read in are from the same run.
  read_shell_slice, files(0), phi=phi, theta = theta, radii = r, data = temp, QUANTITIES = quantities
  numrec = n_elements(temp(0,0,0,0,*))
  nr = n_elements(r)
  nth = n_elements(theta)
  nph = n_elements(phi)
  nq = n_elements(quantities)
  
  quantity_names = strarr(nq)
  for i=0,nq-1 do begin
     quantity_names(i) = get_quantity_name(quantities(i))
  endfor
  ntot = long(numrec*n_files)
  ur = dblarr(nth,nph,ntot)
  urfft = dcomplexarr(nth-8,nph)
  urfft(*,*) = dcomplex(0d0)
  ur(*,*,0:numrec-1) = temp(*,*,nr-1,3,*)
  for i=1L,long(n_files)-1L do begin
     read_shell_slice, files(i), data = temp
     ur(*,*,i*numrec:(i+1L)*numrec-1L) = temp(*,*,nr-1,3,*)
     print, i, ' of', n_files
  endfor

  ;average in time avoiding theta boundaries
  for i=0L,ntot-1L do begin
     urplane = reform(ur(4:nth-5,*,i))
     ;inds = where(urplane lt 0d0)
     ;mneg = mean(urplane(inds))
     ;inds = where(urplane lt mneg,complement=cinds)
     ;urplane(inds) = 1d0
     ;urplane(cinds) = 0d0
     urfft = urfft+fft(urplane,/double)
  endfor
  urfft=urfft/dcomplex(ntot)

  ;average in theta and phi, and average the two
  uftavg = dcomplexarr(nph)
  uftavg(*) = dcomplex(0d0)
  for i=0,nth-9 do begin
     uftavg = uftavg+urfft(i,*)
  endfor
  uftavg = uftavg/dcomplex(nth-8)
  
  ;Save urfft to disk
  openw,file_unit,statsfile,/get_lun
  afile = assoc(file_unit,dcomplexarr(nth-8,nph))
  afile(0) = urfft
  close, file_unit
  free_lun, file_unit
  
stop
end
