@compact_fd6.pro

pro vectorcurl, iter, css_case=css_case, curlB=curlB, mag=mag, rad=rn, theta=theta, phi=phi

  dir1 = '/freyr1/augustso/CSSDE/'+css_case+'/Checkpoints/'

  k=1000000L
  If (iter lt k) Then Begin
     zero_char = '0'
     i=6L
     pretemp = ' '
     stemp = ' '
     While ((iter lt k) and (i gt 0L)) Do Begin
        short_temp = string(i,format='(i1)')
        format = '(i'+strtrim(short_temp,2)+')'
        fmt = strtrim(format,2)
        stemp = string(iter,format=fmt)
        pretemp = strtrim(pretemp+zero_char,2)
        i=i-1L
        k=k/10L
     Endwhile
  Endif Else Begin
     pretemp = ' '
     stemp = string(iter, format='(i7)')
  EndElse

  dir1 = strtrim(dir1,2)+strtrim(pretemp+stemp+'/',2)

  fname = strtrim(dir1+'header',2)

  openr, file_unit, fname, /get_lun, /swap_if_big_endian ;open the file in read-only mode for big endian format                                                 

  temp = lonarr(5)
  readf, file_unit, temp,format='(5i9)'
  
  nr   = fix(temp(0))
  nt  = fix(temp(1))
  np  = fix(temp(2))
  nq   = fix(temp(3))
  iter = long(temp(4))
  
  print, temp
  
  temp = dblarr(8)
  readf, file_unit, temp, format='(8e14.7)'
  dt   = temp(0)
  time = temp(1)
  r1   = temp(2)
  r2   = temp(3)
  th1  = temp(4)
  th2  = temp(5)
  phi1  = temp(6)
  phi2  = temp(7)
  
  print, temp
  
  close, file_unit
  free_lun, file_unit

  rn = (r2-r1)*dindgen(nr)/double(nr-1)+r1
  theta = (th2-th1)*dindgen(nt)/double(nt-1) + th1
  phi = (phi2-phi1)*dindgen(np)/double(np-1) + phi1

  dri = 1d0/(r2-r1)
  dti = 1d0/(th2-th1)
  dpi = 1d0/(phi2-phi1)
  If (keyword_set(mag)) Then begin
     qnames = ['bt','bp','br']
  Endif Else Begin
     qnames = ['w','v','u']
  EndElse
  data = dblarr(nt,np,nr,3)

  for iq=0,2 do begin
     tempname = strtrim(dir1+'checkpoint.'+qnames(iq),2)
     openr, file_unit, tempname, /get_lun, /swap_if_big_endian
     a = assoc(file_unit,dblarr(nt,np,nr))
     data[*,*,*,iq] = a[0]
     close, file_unit
     free_lun, file_unit
  endfor
  
  clearline = fifteenb()

  sines = sin(theta)
  cosines = cos(theta)
  cosec = 1d0/sines
  cotan = cosines/sines
  rinv = 1d0/rn
  curlB = dblarr(nt,np,nr,3)
  for ir=0,nr-1 do begin
     for ip=0,np-1 do begin
        tmp = sines*data[*,ip,ir,1]
        bc1 = -cosines(0)*data[0,ip,ir,2]/dti
        bc2 = -cosines(nt-1)*data[nt-1,ip,ir,2]/dti
        curlB[*,ip,ir,2] = rinv[ir]*cosec*compact_fd6(dti,tmp,bc1,bc2,4)
        tmp = data[*,ip,ir,2]
        bc1 = 0d0
        bc2 = 0d0
        curlB[*,ip,ir,1] = -rinv[ir]*compact_fd6(dti,tmp,bc1,bc2,0)
     endfor
     per = strtrim(round((1.0d0*ir)/(1.0d0*(nr-1.0d0))*100.0d0),2)
     lenp = strtrim(strlen(strtrim(per,2)),2)
     form="($,a"+lenp+",' % Completed',a,a)"
     print, form=form, per, '         ', clearline
  endfor

  print, '          '
  print, 'Theta done'

  for ir=0,nr-1 do begin
     for it=0,nt-1 do begin
        tmp = data[it,*,ir,2]
        bc1 = data[it,np-4:np-1,ir,2]
        bc2 = data[it,0:3,ir,2]
        curlB[it,*,ir,0] = rinv[ir]*cosec[it]*compact_fd6(dpi,tmp,bc1,bc2,4,dtype=3)
        tmp = data[it,*,ir,0]
        bc1 = data[it,np-4:np-1,ir,0]
        bc2 = data[it,0:3,ir,0]
        curlB[it,*,ir,2] = curlB[it,*,ir,2]-rinv[ir]*cosec[it]*compact_fd6(dpi,tmp,bc1,bc2,4,dtype=3)
     endfor
     per = strtrim(round((1.0d0*ir)/(1.0d0*(nr-1.0d0))*100.0d0),2)
     lenp = strtrim(strlen(strtrim(per,2)),2)
     form="($,a"+lenp+",' % Completed',a,a)"
     print, form=form, per, '         ', clearline
  endfor

  print, '          '
  print, 'Phi done'

  for ip=0,np-1 do begin
     for it=0,nt-1 do begin
        tmp = rn*data[it,ip,*,2]
        bc1 = -data[it,ip,0,2]/dri
        bc2 = -data[it,ip,nr-1,2]/dri
        curlB[it,ip,*,0] = curlB[it,ip,*,0]-rinv*compact_fd6(dri,tmp,bc1,bc2,4)
        tmp = rn*data[it,ip,*,0]
        bc1 = -data[it,ip,0,0]/dri
        bc2 = -data[it,ip,nr-1,0]/dri
        curlB[it,ip,*,1] = curlB[it,ip,*,1]+rinv*compact_fd6(dri,tmp,bc1,bc2,0)
     endfor
     per = strtrim(round((1.0d0*ip)/(1.0d0*(np-1.0d0))*100.0d0),2)
     lenp = strtrim(strlen(strtrim(per,2)),2)
     form="($,a"+lenp+",' % Completed',a,a)"
     print, form=form, per, '         ', clearline
  endfor

  print, '          '
  print, 'Radial done'

end

function fifteenb

   return, string("15b)

end
