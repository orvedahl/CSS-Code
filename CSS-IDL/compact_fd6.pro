;===============================================================
;dri = 1/(x2-x1) where x2 and x1 are the endpoints of the domain
;arr = array to be differenced
;b1  = boundary value on x1 (specify if ibc = 3 or 4)
;b2  = boundary value on x2 (specify if ibc = 2 or 4) 
;
;ibc = boundary condition type:  
; 0 - darr/dr = 0
; 1 - Dirichlet on x1 & x2
; 2 - Neumann on x1 & Dirichlet on x2
; 3 - Vice Versa
; 4 - Neumann on x1 & x2
;
;dtype = domain type (0 for full domain, 1 or 2 for boundary 3 for internal)
;(typically use dtype = 0 
;=====================================================

function compact_fd6,dri,arr,b1,b2,ibc,dtype=dtype,darr=darr

    If (not keyword_set(dytpe)) Then dtype = 0

    n1 = n_elements(arr)
    If (not keyword_set(darr)) Then Begin
        darr = dblarr(n1)
    EndIf

    If (ibc eq 5) Then Begin
        d1i = double(n1)
    EndIf Else Begin
        d1i = double(n1-1)
    EndElse

    c1 = d1i*7d0/9d0
    c2 = d1i/36d0

    If ((dtype eq 1) or (dtype eq 3)) Then Begin
        darr(0) = c1*(arr(1)-b1(3))+c2*(arr(2)-b1(2))
        darr(1) = c1*(arr(2)-arr(0))+c2*(arr(3)-b1(3))
    EndIf
    If ((dtype eq 2) or (dtype eq 3)) Then Begin
        darr(n1-2) = c1*(arr(n1-1)-arr(n1-3))+c2*(b2(0)-arr(n1-4))
        darr(n1-1) = c1*(b2(0)-arr(n1-2))+c2*(b2(1)-arr(n1-3))
    EndIf

    For i=2,n1-3 Do Begin
        darr(i) = c1*(arr(i+1)-arr(i-1))+c2*(arr(i+2)-arr(i-2))
    EndFor

    alpha = 1d0/3d0
    uppr = dblarr(n1)
    lowr = dblarr(n1)
    diag = dblarr(n1)
    diag(*) = 1d0
    uppr(*) = alpha
    lowr(*) = alpha   

    If (dtype lt 3) Then Begin
        If ((dtype eq 1) or (dtype eq 0)) Then Begin
            c1 = -43d0/96d0*d1i
            c2 = -5d0/6d0*d1i
            c3 = 9d0/8d0*d1i
            c4 = 1d0/6d0*d1i
            c5 = -1d0/96d0*d1i
            darr(n1-2) =-(c1*arr(n1-1) + c2*arr(n1-2) + c3*arr(n1-3) + c4*arr(n1-4) + c5*arr(n1-5))
            
                                ;fifth order coeff.
            c1 = -10d0/3d0*d1i
            c2 = -3d0*d1i
            c3 = 6d0*d1i
            c4 = 1d0/3d0*d1i
            If (ibc eq 0) Then Begin
                darr(n1-1)=0d0
                                ;if ibc = 3,4 Then dxdy(*,n2) must be set by the calling routine
            EndIf
            If ((ibc eq 1) or (ibc eq 2)) Then Begin
                darr(n1-1)=-(c1*arr(n1-1)+c2*arr(n1-2)+c3*arr(n1-3)+c4*arr(n1-4))
            Endif

            If (dtype eq 1) Then Begin
                darr(0) = darr(0)-alpha*(-b1(0)+9d0*b1(1)-45d0*b1(2)+45d0*arr(0)-9d0*arr(1)+arr(2))*d1i/60d0
            Endif

                                ;--------------------------------------------------------------
                                ;       Now we set up the matrix
                                ;--------------------------------------------------------------
                                ; here is the sixth order interior value     
                                ;here are the pentadiagonal and fifth order values for the boundary.
            alpha2=3d0/4d0
            gamma2=1d0/8d0
            alpha1=6d0
            beta1=3d0
                                ;precondition the matrix to make it tridiagonal
            const=1d0/(alpha2-beta1*gamma2)
            up1=(alpha1*alpha2-beta1)*const

            If ((ibc eq 1) or (ibc eq 2)) Then Begin
                darr(n1-1)=(alpha2*darr(n1-1)-beta1*darr(n1-2))*const
            Endif
            lowr(n1-2)=alpha2
            uppr(n1-2)=gamma2

            If (ibc ge 0) Then Begin
                If ((ibc ne 1) and (ibc ne 2)) Then Begin
                    lowr(n1-1) = 0d0
                EndIf Else Begin
                    lowr(n1-1) = up1
                EndElse
            EndIf

        EndIf
            
        If ((dtype eq 2) or (dtype eq 0)) Then Begin
            c1 = -43d0/96d0*d1i
            c2 = -5d0/6d0*d1i
            c3 = 9d0/8d0*d1i
            c4 = 1d0/6d0*d1i
            c5 = -1d0/96d0*d1i
            darr(1) = c1*arr(0) + c2*arr(1) + c3*arr(2) + c4*arr(3) + c5*arr(4) 
                
                                ;fifth order coeff.
            c1 = -10d0/3d0*d1i
            c2 = -3d0*d1i
            c3 = 6d0*d1i
            c4 = 1d0/3d0*d1i
            If (ibc eq 0) Then Begin
                darr(0)=0d0
            EndIf
                                ;if ibc = 2,4 Then dxdy(*,1) must be set by the calling routine
            If ((ibc eq 1) or (ibc eq 3)) Then Begin
                darr(0)=c1*arr(0)+c2*arr(1)+c3*arr(2)+c4*arr(3)
            EndIf
            If (dtype eq 2) Then Begin
                darr(n1-1) = darr(n1-1)-alpha*(-arr(n1-3)+9d0*arr(n1-2)-45d0*arr(n1-1)+45d0*b2(1)-9d0*b2(2)+b2(3))*d1i/60d0
            Endif

                                ;--------------------------------------------------------------
                                ;       Now we set up the matrix
                                ;--------------------------------------------------------------
                                ; here is the sixth order interior value     
                                ;here are the pentadiagonal and fifth order values for the boundary.
            alpha2=3d0/4d0
            gamma2=1d0/8d0
            alpha1=6d0
            beta1=3d0

                                ;precondition the matrix to make it tridiagonal
            const=1d0/(alpha2-beta1*gamma2)
            up1=(alpha1*alpha2-beta1)*const
            If ((ibc mod 2) ne 0) Then Begin
                darr(0)=(alpha2*darr(0)-beta1*darr(1))*const
            EndIf
            If ((ibc mod 2) eq 0) Then Begin
                uppr(0) = 0d0
            EndIf Else Begin
                                ; fifth order bc.
                uppr(0)=up1
            EndElse
            uppr(1)=alpha2
            lowr(1)=gamma2
        EndIf

        lowr(0) = 0d0
        uppr(n1-1) = 0d0

        darr = tri_solver(lowr,diag,uppr,darr,n1,is=0)
    EndIf Else Begin
         lowr(0) = 0d0
         uppr(n1-1) = 0d0
                                ;Dtype = 3 (internal boundaries only)        
         darr(0) = darr(0)-alpha*(-b1(0)+9d0*b1(1)-45d0*b1(2)+45d0*arr(0)-9d0*arr(1)+arr(2))*d1i/60d0
         darr(n1-1) = darr(n1-1)-alpha*(-arr(n1-3)+9d0*arr(n1-2)-45d0*arr(n1-1)+45d0*b2(1)-9d0*b2(2)+b2(3))*d1i/60d0

         darr = tri_solver(lowr,diag,uppr,darr,n1)
    EndElse

    ;stop

    return, darr*dri

end
