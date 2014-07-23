@NonUniformDerivative6.pro

function cheby_trans, x, inverse=inverse

   sx = size(x)
   ;Case X is a column vector.
   If ((sx[0] eq 2) and (sx[1] eq 1)) Then Begin
       n1=sx[2]
       sel='col'
   Endif Else Begin
       If ((sx[0] eq 2) and (sx[1] ne sx[2])) Then Message,'Array must be a square array or a vector.'
       n1=sx[1]
       sel='array'
   EndElse
   
   ;Case X is a row vector.
   If (sx[0] eq 1) Then Begin
       n1=n_elements(x)
       sel='row'
   End

   n = n1-1
   
   If (not keyword_set(inverse)) Then Begin
       cj = dindgen(n1)*0d0+1d0
       cj(0) = 1d0/2d0
       cj(n) = 1d0/2d0

       cjk = cj#replicate(1,n1)
       cjk = cjk*transpose(cjk)
   EndIf

   theta = !dpi/double(n)

   thetajk = dindgen(n1)#replicate(1,n1)
   thetajk = thetajk*transpose(thetajk)

   Ajk = cos(theta*thetajk)

   If (not keyword_set(inverse)) Then Begin
       Ajk = cjk*Ajk*2d0/double(n)
   EndIf

   Case (sel) Of
       'col': Begin
           ;For an inverse transform compute Transpose(Ajk)##x
           If (keyword_set(inverse)) Then return,Transpose(Ajk)##x

   	   ;For a forward transform compute C##X
           Return,Ajk##x
	End

        'row': Begin
	    ;For an inverse transform compute x##Ajk
            If (keyword_set(inverse)) Then Return,x##Ajk

  	    ;For a forward transform compute x##Transpose(Ajk)
            Return,x##Transpose(Ajk)
        End

        'array': Begin
            ;For an inverse transform compute Transpose(Ajk)##x##Ajk
            If (keyword_set(inverse)) Then Return,Transpose(Ajk)##x##Ajk

            ;For a forward transform compute C##X##Transpose(Ajk)
            Return,Ajk##x##Transpose(Ajk)
	End

        Else: Message,'Input must be a square array or a vector.'

    EndCase

end

function d_by_dr, x, ro, ri, nr2=nr2, ri1=ri1
 
  If (not keyword_set(nr2)) Then Begin
      ri1 = ri
      nr1 = n_elements(x)
  Endif Else Begin
      nr1 = n_elements(x)-nr2
  EndElse
  
  x = 1d0*x
  dx = dblarr(n_elements(x))
  dx(*) = 0d0

  nr = nr1
  nm = nr1

  scaling = 2d0/(ro-ri1)
  c1 = 2d0*double(nr-1)*scaling
  dx(nm-1) = 0d0
  dx(nm-2) = c1*x(nm-1)/2d0

  for r=nm-3, nm-nr, -1 do begin
      c1 = c1-2d0*scaling
      dx(r) = dx(r+2)+c1*x(r+1)
  endfor

  dx(0) = dx(0)/2d0

  If (keyword_set(nr2)) Then Begin
     nr = nr2
     nm = nr1+nr2

     scaling = 2d0/(ro-ri1)
     c1 = 2d0*double(nr-1)*scaling

     dx(nm-1) = 0d0
     dx(nm-2) = c1*x(nm-1)/2d0

     for r=nm-3, nm-nr, -1 do begin
         c1 = c1 - 2d0*scaling
         dx(r) = dx(r+2)+c1*x(r+1)
     endfor

     dx(nm-nr) = dx(nm-nr)/2d0

  EndIf

  return, dx

end

function compute_density, r_m, grav_m, rho_m, gamma_m, Cp_m, S_m, dsdr_m, N_R_1, N_R_1_aliased, Anchor_Index, $
                          N_R_2 = N_R_2, N_R_2_aliased = N_R_2_aliased, InnerRadius_1=InnerRadius_1, ncgcp=ncgcp

   If (keyword_set(N_R_2)) Then Begin
       N_R = N_R_1+N_R_2
       Dual_Chebyshev = 1
   Endif Else Begin
       N_R = N_R_1
       InnerRadius_1 = r_m(N_R-1)
       Dual_Chebyshev = 0
   EndElse

   Matrix = dblarr(N_R,N_R)
   rhs = dblarr(N_R)
   rho = rhs
   pivot = intarr(N_R)
   drho = rho
   gam = rho
   c1 = rho
   c2 = rho
   T = 0d0
   n2 = 0d0
   TC1 = dblarr(N_R_1)

   If (keyword_set(N_R_2)) Then TC2 = dblar(N_R_2)

   dT1 = dblarr(N_R_1,N_R_1)

   If (keyword_set(N_R_2)) Then dT2 = dblarr(N_R_2,N_R_2)

   scaling = 0d0
   d1 = 0d0
   z = rho
   Pi_over_N_1 = !dpi/double(N_R_1-1)

   If (keyword_set(N_R_2)) Then Pi_over_N_2 = !dpi/double(N_R_2-1)

   iteration = 0
   np = 0
   rp = 0
   M = 0
   zero = 0d0
   one = 1d0
   two = 2d0
   half = 1d0/2d0
   error_for_variations=one
   error_for_functions=one
   
        
   If (Anchor_Index < 0) Then Anchor_Index = 1

   ; use the solar model interpolated onto our grid as an initial guess
   z = r_m
   OuterRadius = z(0)
   InnerRadius = z(N_R-1)
   rho = rho_m

   ; get rid of higher order modes
   drho = cheby_trans(rho(0:N_R_1-1))
   drho(N_R_1_aliased-1:N_R_1-1) = zero
   rho(0:N_R_1-1) = cheby_trans(drho,/inverse)

   If (Dual_Chebyshev) Then Begin
       drho = cheby_trans(rho(N_R_1:*))
       drho(N_R_2_aliased-1:N_R_2-1) = zero
       rho(N_R_1:*) = cheby_trans(drho,/inverse)
   EndIf

   ;set up the chebyshev derivative matrices
   TC1(*) = one 
   TC1(0)=half
   TC1(N_R_1-1)=half

   M = N_R_1-1
   dT1(*,*) = zero
   scaling = two/(OuterRadius-InnerRadius_1)

   dT1(0,*) = zero
   for n = 1, M-1, 2 do begin
       d1 = two*n*scaling
       for r=1,n,2 do begin
           dT1(n-1,r-1) = d1
       endfor
   endfor
   
   for n = 2, M, 2 do begin
       d1 = two*n*scaling
       for r=2,n,2 do begin
           dT1(n-1,r-1) = d1
       endfor
   endfor

   for r=0,M do begin
       dT1(*,r) = cheby_trans(reform(dT1(*,r)))
   endfor

   ;dT1 = transpose(dT1)

   If (Dual_Chebyshev) Then Begin
       TC2(*) = one
       TC2(0)=half
       TC2(N_R_2-1)=half
       M = N_R_2-1
       dT2(*,*) = zero
       scaling = two/(InnerRadius_1-InnerRadius)
       dT2(0,*) = zero
       for n = 1, M-1, 2 do begin
           d1 = two*n*scaling
           for r=0,n-1,2 do begin
               dT2(n,r) = d1
           endfor
       endfor

       for n = 2, M, 2 do begin
           d1 = two*n*scaling
           for r=1,n-1,2 do begin
               dT2(n,r) = d1
           endfor
       endfor
       dT2 = cheby_trans(dT2,/inverse)
   EndIf

   ;initialize a few constants
   gam = two - gamma_m
   If (keyword_set(ncgcp)) Then Begin
       dlngammadr = deriv(r_m,gamma_m)/gamma_m
       dlnCpdr = deriv(r_m,Cp_m)/Cp_m
       c1 = dsdr_m/Cp_m + (alog(rho) + S_m/Cp_m + one)*dlngammadr - S_m/Cp_m*dlnCpdr
       c2 = (grav_m/gamma_m)*exp(-gamma_m*S_m/Cp_m)
   Endif Else Begin
       c1 = dsdr_m/Cp_m
       c2 = (grav_m/gamma_m)*exp(-gamma_m*S_m/Cp_m)
   EndElse

   ; now iterate until convergence
   iteration = 1
   solar_max_iter = 10
   convergence_for_variations = 1d-7
   convergence_for_functions = 1d-13
   While ((iteration le solar_max_iter) and $
             ((error_for_variations gt Convergence_for_variations)  $
             and (error_for_functions  gt Convergence_for_functions) )) Do Begin

       ; first set up the coefficient matrix
       matrix(*,*) = zero
       ; the first eqn...
       for r = 0,N_R_1-1 do begin
           for n = 0,N_R_1-1 do begin
               T = cos(n*r*Pi_over_N_1) * TC1(n)
               matrix(r,n) = dT1(n,r) + (c1(r) + gam(r)*c2(r)*rho(r)^(gam(r)-one) ) * T
           endfor
       endfor

       If (Dual_Chebyshev) Then Begin
           for r = N_R_1,N_R-1 do begin
               rp = r-N_R_1
               for n = 0,N_R_2-1 do begin
                   T = cos(n*rp*Pi_over_N_2) * TC2(n)
                   np = n+N_R_1
                   matrix(r,np) = dT2(n,rp+1) + (c1(r) + gam(r)*c2(r)*rho(r)^(gam(r)-one) ) * T
               endfor
           endfor
       EndIf

       ; now set up the right-hand side
       rhoct = cheby_trans(rho)
       If (keyword_set(N_R_2)) Then Begin
           drho = d_by_dr(rhoct,OuterRadius,InnerRadius,nr2=N_R_2,ri1=InnerRadius_1)
       Endif Else Begin
           drho = d_by_dr(rhoct,OuterRadius,InnerRadius)
       EndElse
       drho = cheby_trans(drho,/inverse)
       drho = nud6(r_m,rho)
       for r = 0,N_R-1 do begin
           rhs(r) = -(drho(r) + c1(r)*rho(r) + c2(r)*rho(r)^(gam(r)))
       endfor
       ; This is the Boundary Condition
       matrix(anchor_index,*) = zero
       If (anchor_index gt N_R_1) Then Begin
           rp = anchor_index-1-N_R_1
           for n=0,N_R_2_aliased-1 do begin
               T = cos(n*rp*Pi_over_N_2) * TC2(n)     
               matrix(anchor_index,n+N_R_1) = T
           endfor
       Endif Else Begin
           rp = anchor_index-1
           for n=0,N_R_1_aliased-1 do begin
               T = cos(n*rp*Pi_over_N_1) * TC1(n)                 
               matrix(anchor_index,n) = T
           EndFor
       EndElse
       rhs(anchor_index) = zero
        
       ; and this is the continuity condition
       If (Dual_Chebyshev) Then Begin
           matrix(N_R_1-1,*) = zero 
           For n=0,N_R_1_aliased-1 Do Begin
               T = TC1(n) * (-one)^n
               n2 = n*n
               matrix(N_R_1,n) = -T
           EndFor
           For n=0,N_R_2_aliased-1 Do Begin
               T = TC2(n)
               n2 = n*n
               np = n+N_R_1
               matrix(N_R_1,np) = T
           EndFor
           rhs(N_R_1) = rho(N_R_1-1)-rho(N_R_1)
       EndIf

       error_for_functions = total(abs(rhs)) 
       If (Dual_Chebyshev) Then error_for_functions = error_for_functions-abs(rhs(N_R_1+1))
       error_for_functions = error_for_functions / total(abs(drho))

       ; now solve for the variations,...
       ;matrix = transpose(matrix)
       ;rhs = transpose(rhs)
       ludc,matrix,pivot
       rhs = lusol(matrix, pivot, rhs)
       ;remove the highest few modes
       If (N_R_1_aliased lt N_R_1) Then rhs(N_R_1_aliased-1:N_R_1-1) = zero
       If (Dual_Chebyshev) Then Begin
           If (N_R_2_aliased lt N_R_2) Then rhs(N_R_1+N_R_2_aliased-1:N_R-1) = zero
       EndIf

       ; ...move them into physical space...      
       drho = dblarr(N_R)
       drho(0:N_R_1-1) = cheby_trans(rhs(0:N_R_1-1),/INVERSE)

       If (Dual_Chebyshev) Then Begin
           drho(N_R_1:*) = cheby_trans(rhs(N_R_1:*),/INVERSE)
       EndIf

       ; add in the variations
       rho = rho + drho           
stop
       ; check for sanity
       If (min(rho) le zero) Then Begin
           print, 'ERROR!! Density went negative in compute_density'
           for r=0,N_R-1 do begin
               print, r,r_m(r),rho(r)
           endfor
           stop
       EndIf

       ; check for convergence
       error_for_variations = total(abs(drho))/(total(abs(rho)))
           
       print, 'compute_density,(vars,funcs) : ', iteration, error_for_variations, error_for_functions

       iteration = iteration + 1
           
   EndWhile
        
   return, rho

end
