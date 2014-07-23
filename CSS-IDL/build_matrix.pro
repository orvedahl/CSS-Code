pro build_matrix, order, acc, bc

;Solve for all points, and build a Fortran code snippet written to file
;Build matrices for each point up to where centered difference can be used
;Do inner boundary
nd = acc+order-1
if (bc eq 'dirichlet') then begin
   print, 'inner'
   inner = dblarr(nd,nd/2)
   for k=0,nd/2-1 do begin
       A = dblarr(nd,nd)
       b = dblarr(nd)
       pts  = dblarr(nd)
       for i=0,nd-1 do begin
           m=i+2        
           if (k ge 1) then begin
              n=-k
           endif else begin
              n=1
           endelse
           for j=0,nd-1 do begin
               A(i,j) =  (1.0d0*n)^m/double(factorial(m))
               pts(j) = n
               n=n+1
               if (n eq 0) then begin
                  n=1 
               endif
           endfor        
           b(i) = 0.0d0
       endfor
       b(nd-1) = -4.0d0/7.0d0
       A = transpose(A)
       ;Solve for a_n
       LUDC, A, index

       temp = LUSOL(A,index,b)
       inner(*,k) = temp(*)
       print, 'x='
       print, reform(inner(*,k))
       print, pts
       print, -total(inner(*,k))
       sum = 0.0d0
       if (k ge 1) then begin
          n=-k
       endif else begin
          n=1
       endelse
       for i=0,nd-1 do begin
           sum = sum + n*inner(i,k)
           n=n+1
           if (n eq 0) then begin
              n=1 
           endif
       endfor
       print, sum 
   endfor

   ;Build matrix for interior (centered difference)
   A = dblarr(nd,nd)
   b = dblarr(nd)
   pts = dblarr(nd)
   for i=0,nd-1 do begin

       if (i eq 0) then begin
          if (order eq 1) then begin
             m=2
          endif else begin
             m=1
          endelse
       endif else
          m = i+2
          if (m eq order) then begin
             m=order+1
          endif
       endelse

       n = -nd/2
       for j=0,nd-1 do begin       
           A(i,j) = (1.0d0*n)^m/double(factorial(m))
           pts(j) = n
           n=n+1
           if (n eq 0) then begin
              n=1 
           endif
       endfor
       b(i) = 0.0d0
   endfor 
   b(nd-1) = -2.0d0/5.0d0
   A = transpose(A)
   ;Solve for a_n
   LUDC, A, index
   interior = LUSOL(A,index,b)
   print, 'interior'
   print, 'x='
   print, interior
   print, pts
   print, -total(interior)
   sum = 0.0d0
   n=-nd/2
   for i=0,nd-1 do begin
       sum = sum + n*interior(i)
       n=n+1
       if (n eq 0) then begin
          n=1 
       endif
   endfor
   print, sum 

   print, 'outer'
   outer = dblarr(nd,nd/2)
   for k=nd/2-1,0,-1 do begin
       A = dblarr(nd,nd)
       b = dblarr(nd)
       pts = dblarr(nd)
       for i=0,nd-1 do begin
           m=i+2        
           n=-nd+k
           for j=0,nd-1 do begin
               A(i,j) = (1.0d0*n)^m/double(factorial(m))
               pts(j) = n
               n=n+1
               if (n eq 0) then begin
                  n=1 
               endif
           endfor
           b(i) = 0.0d0
       endfor
       b(nd-1) = -60.0d0/7.0d0
       A = transpose(A)
       ;Solve for a_n
       LUDC, A, index
       temp = LUSOL(A,index,b)
       outer(*,k) = temp(*)
       print, 'x='
       print, reform(outer(*,k))
       print, pts
       print, -total(outer(*,k))
       sum = 0.0d0
       n=-nd+k
       for i=0,nd-1 do begin
           sum = sum + n*outer(i,k)
           n=n+1
           if (n eq 0) then begin
              n=1 
           endif
       endfor
       print, sum 
   endfor
endif else if (bc eq 'neumann') then begin
   print, 'inner'
   inner = dblarr(nd,nd/2)
   for k=0,nd/2-1 do begin
       A = dblarr(nd,nd)
       b = dblarr(nd)
       for i=0,nd-1 do begin
           m=i+2        
           if (k ge 1) then begin
              n=-k
           endif else begin
              n=1
           endelse
           for j=0,nd-1 do begin
               A(i,j) =  (1.0d0*n)^m/double(factorial(m))
               n=n+1
               if (n eq 0) then begin
                  n=1 
               endif
           endfor        
           b(i) = 0.0d0
       endfor
       b(nd-1) = 60.0d0
       A = transpose(A)
       ;Solve for a_n
       LUDC, A, index
       temp = LUSOL(A,index,b)
       inner(*,k) = temp(*)
       print, 'x='
       print, reform(inner(*,k))
       print, total(inner(*,k))
       sum = 0.0d0
       if (k ge 1) then begin
          n=-k
       endif else begin
          n=1
       endelse
       for i=0,nd-1 do begin
           sum = sum + n*inner(i,k)
           n=n+1
           if (n eq 0) then begin
              n=1 
           endif
       endfor
       print, sum 
   endfor

   ;Build matrix for interior (centered difference)
   A = dblarr(nd,nd)
   b = dblarr(nd)
   for i=0,nd-1 do begin
       m=i+2
       n = -nd/2
       for j=0,nd-1 do begin       
           A(i,j) = (1.0d0*n)^m/double(factorial(m))
           n=n+1
           if (n eq 0) then begin
              n=1 
           endif
       endfor
       b(i) = 0.0d0
   endfor 
   b(nd-1) = 3.0d0
   A = transpose(A)
   ;Solve for a_n
   LUDC, A, index
   interior = LUSOL(A,index,b)
   print, 'interior'
   print, 'x='
   print, interior
   print, total(interior)
   sum = 0.0d0
   n=-nd/2
   for i=0,nd-1 do begin
       sum = sum + n*interior(i)
       n=n+1
       if (n eq 0) then begin
          n=1 
       endif
   endfor
   print, sum 

   print, 'outer'
   outer = dblarr(nd,nd/2)
   for k=nd/2-1,0,-1 do begin
       A = dblarr(nd,nd)
       b = dblarr(nd)
    for i=0,nd-1 do begin
        m=i+2        
        n=-nd+k
        for j=0,nd-1 do begin
            A(i,j) =  (1.0d0*n)^m/double(factorial(m))
            n=n+1
            if (n eq 0) then begin
               n=1 
            endif
        endfor
        b(i) = 0.0d0
    endfor
    b(nd-1) = 60.0d0/7.0d0
    ;print, 'A='
    A = transpose(A)
    ;print, A
    ;Solve for a_n
    LUDC, A, index
    temp = LUSOL(A,index,b)
    outer(*,k) = temp(*)
    print, 'x='
    print, reform(outer(*,k))
    print, total(outer(*,k))
    sum = 0.0d0
    n=-nd+k
    for i=0,nd-1 do begin
        sum = sum + n*outer(i,k)
         n=n+1
         if (n eq 0) then begin
            n=1 
         endif
    endfor
    print, sum 
    endfor
endif else begin
    print, 'Pick a boundary type'
endelse

end

