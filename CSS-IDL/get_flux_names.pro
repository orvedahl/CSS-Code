function get_flux_names, qs

   nq = n_elements(qs)
   string_name=strarr(nq)
   for i=0,nq-1 do begin
       case qs(i) of
            1 : string_name(i) = 'KE'
            2 : string_name(i) = 'Enth'
            3 : string_name(i) = 'Visc'
            4 : string_name(i) = 'Rad'
            5 : string_name(i) = 'Unr'
            6 : string_name(i) = 'Poyn'
           else : string_name(i) = 'Unknown'
       endcase
   endfor

   return, string_name

end
