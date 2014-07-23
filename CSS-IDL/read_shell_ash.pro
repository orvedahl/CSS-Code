@get_quantity_name_ash.pro
pro read_shell_ash, file, qs, recnums,shell_levels,big_array, costheta, TIMES = times, TITLES = titles, $
                ITERS = iters, Radii = radii, n_shells = n_shells, n_q = n_q, QUANTITIES = quantities, $
                all_qs = all_qs, all_qs_rs=all_qs_rs, all_rs = all_rs
  
; This routine reads the shell data from "file."
; INPUTS:
;	file - the file containing the shell data
;	qs   - The indice array you want from the file (staring with 1)
;	recnums - the record numbers you want to extract (starting with 0)
;	shell_levels - the shell indices you want to extract (starting with 0)
;
; OUTPUTS:
;	  big_array    - an array of dimenions [recnums,shell_levels,qs,ntheta, nphi]
;			this is the primary output - contains all your shell info
;	  costheta - the cosine theta array read from the file	  
;	  times  - the times corresponding to each record number
;	  titles - the titles for your quantities in case you care
;	  iters - the iteration number cooresponding to each record number
;	  radii = the radii corresponding to each shell_level indice
  nqs   = n_elements(qs)
  n_shell_levels = n_elements(shell_levels)
  nrecs = recnums
;---------------------------------------------
;----------    First Open The File    --------
;---------------------------------------------
  
  OPENR,file_lun,file,byteorder,/swap_if_big_endian,/get_lun
  READU,file_lun,temp  &  n_phi=fix(temp)
  READU,file_lun,temp  &  n_theta=fix(temp)
  READU,file_lun,temp  &  n_shells=fix(temp)
  READU,file_lun,temp  &  n_q=fix(temp)
  
  Print, n_phi, n_theta, n_shells, n_q
  
  If (n_q eq 0) Then Begin
     close, file_lun
     free_lun, file_lun
     openr, file_lun, file, /swap_endian, /get_lun
     READU,file_lun,temp  &  n_phi=fix(temp)
     READU,file_lun,temp  &  n_theta=fix(temp)
     READU,file_lun,temp  &  n_shells=fix(temp)
     READU,file_lun,temp  &  n_q=fix(temp)
     
     Print, n_phi, n_theta, n_shells, n_q
  EndIf
  
  IF (KEYWORD_SET(all_qs)) THEN BEGIN
     nqs   = n_q
  ENDIF
  IF (KEYWORD_SET(all_rs)) THEN BEGIN
     shell_levels = indgen(n_shells)
     n_shell_levels = n_elements(shell_levels)
  ENDIF
  IF (KEYWORD_SET(all_qs_rs)) THEN BEGIN
     shell_levels = indgen(n_shells)
     nqs   = n_q
     n_shell_levels = n_elements(shell_levels)
  ENDIF
  
  big_array    = fltarr(n_theta,n_phi,n_shell_levels,nqs,nrecs)
  shell_slices =assoc(file_lun,fltarr(n_theta,n_phi+1,n_shells,n_q))
  costheta=shell_slices[0:n_theta-1,0,0,1,0]
  quantities=fix(shell_slices[0:n_q-1,0,0,3,0])
  IF (KEYWORD_SET(all_qs) or KEYWORD_SET(all_qs_rs)) THEN BEGIN
     qs = quantities
  ENDIF
  radii = fltarr(n_shell_levels)
  titles = strarr(nqs)
  iters = lindgen(nrecs)
  times = fltarr(nrecs)
  FOR m = 0, nqs -1 DO BEGIN
     iq = where(quantities eq qs[m])
     titles[m]=get_quantity_name(quantities[iq])
  ENDFOR
  FOR m = 0, n_shell_levels-1 do begin
     radii[m]=shell_slices[8,0,shell_levels[m],0,0]
  ENDFOR
;------------------------------------------------------------
;----------   Now Read Everything In
;------------------------------------------------------------
  FOR i = 0, nrecs -1 do begin
     iters[i] = shell_slices[4,0,0,0,i]
     times[i]  = shell_slices[5,0,0,0,i]
     
     FOR k = 0, nqs -1 do begin
        FOR j = 0, n_shell_levels -1 DO BEGIN
           iq = where(quantities eq qs[k])
           big_array[*,*,j,k,i] = shell_slices[*,1:n_phi,shell_levels[j],iq,i]
        ENDFOR
     ENDFOR
  ENDFOR
;--------------------------------------------------------------
;       Now close the file, and thank me later for this nice interface!!
;--------------------------------------------------------------
  close, file_lun
  free_lun, file_lun
END
