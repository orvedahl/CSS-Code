function get_quantity_name, quantity, tex=tex
  nq = n_elements(quantity)
  names = strarr(nq)
  for iq=0,nq-1 do begin
   case quantity[iq] of
     1    : name='Rho'
     2    : name='Temperature'
     3    : name='Entropy'
     4    : name='Pressure'
     5    : name='Vr'
     6    : name='Vtheta'
     7    : name='Vphi'
     8    : name='KE'
     9    : name='ds/dr'
     10   : name='Br'
     11   : name='Btheta'
     12   : name='Bphi'
     13   : name='ME'
     14   : name='Poloidal Mag'
     15   : name='F_ks'
     16   : name='F_en'
     17   : name='F_ke'
     18   : name='F_ac'
     19   : name='F_sl'
     else : name='Unknown'

   endcase
   names[iq] = name
endfor
   return, names

end






