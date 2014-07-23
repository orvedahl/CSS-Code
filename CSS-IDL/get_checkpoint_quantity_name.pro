function get_checkpoint_quantity_name, qs

   nq = n_elements(qs)
   string_name=strarr(nq)
   for i=0,nq-1 do begin
       case qs(i) of
            1 : string_name(i) = 'rho'
            2 : string_name(i) = 'u_theta'
            3 : string_name(i) = 'u_phi'
            4 : string_name(i) = 'u_r'
            5 : string_name(i) = 'S'
            6 : string_name(i) = 'B_r'
            7 : string_name(i) = 'B_theta'
            8 : string_name(i) = 'B_phi'
           else : string_name(i) = 'Unknown'
       endcase
   endfor

   return, string_name

end
