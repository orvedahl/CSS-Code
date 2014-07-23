function compact_fd6_2,dri,arr,b1,b2,ibc,dtype, darr=darr

    n1 = n_elements(arr)
    If (not keyword_set(darr)) Then Begin
        darr = dblarr(n1)
    EndIf

    If (ibc eq 5) Then Begin
        d1i = double(n1)
    EndIf Else Begin
        d1i = double(n1-1)
    EndElse

    d1i2 = d1i^2
    c1 = d1i2*12d0/11d0
    c2 = d1i2*3d0/44d0

    If ((dtype eq 1) or (dtype eq 3)) Then Begin
        darr(0) = c1*(arr(1)-2d0*arr(0)+b1(3))+c2*(arr(2)-2d0*arr(0)+b1(2))
        darr(1) = c1*(arr(2)-2d0*arr(1)+arr(0))+c2*(arr(3)-2d0*arr(1)+b1(3))
    EndIf
    If ((dtype eq 2) or (dtype eq 3)) Then Begin
        darr(n1-2) = c1*(arr(n1-1)-2d0*arr(n1-2)+arr(n1-3))+c2*(b2(0)-2d0*arr(n1-2)+arr(n1-4))
        darr(n1-1) = c1*(b2(0)-2d0*arr(n1-1)+arr(n1-2))+c2*(b2(1)-2d0*arr(n1-1)+arr(n1-3))
    EndIf

    For i=2,n1-3 Do Begin
        darr(i) = c1*(arr(i+1)-2d0*arr(i)+arr(i-1))+c2*(arr(i+2)-2d0*arr(i)+arr(i-2))
    EndFor

    alpha = 2d0/11d0
    uppr = dblarr(n1)
    lowr = dblarr(n1)
    diag = dblarr(n1)
    diag(*) = 1d0
    uppr(*) = alpha
    lowr(*) = alpha   

    If (dtype lt 3) Then Begin
        If ((dtype eq 1) or (dtype eq 0)) Then Begin
            
            c1=99d0/80d0*d1i2
            c2=-3d0*d1i2
            c3=93d0/40d0*d1i2
            c4=-0.6d0*d1i2
            c5=3d0/80d0*d1i2            

            darr(n1-2) = c1*arr(n1-1) + c2*arr(n1-2) + c3*arr(n1-3) + c4*arr(n1-4) + c5*arr(n1-5)
     
            gamma2=0.1d0
            alpha2=-7d0/20d0

            uppr(n1-2)=gamma2
            lowr(n1-2)=alpha2
       

            If ((ibc ne 1) and (ibc ne 2)) Then Begin
                lowr(n1-1) = 6d0
                c1=89d0/18d0*d1i2
                c2=-12d0*d1i2
                c3=15d0/2d0*d1i2
                c4=-4d0/9d0*d1i2
                c5=0d0
                
                If (ibc eq 0) Then Begin
                    c6=0d0
                Endif Else Begin
                    c6=5d0/3d0*d1i
                EndElse

                darr(n1-1)=c1*arr(n1-1)+c2*arr(n1-2)+c3*arr(n1-3)+c4*arr(n1-4)+c5*arr(n1-5)+c6*darr(n1-1)
            Endif Else Begin
                ; fourth order boundary
                alpha1=10d0
                lowr(n1-1)= alpha1
                c1 = 145d0/12d0*d1i2
                c2 = -76d0/3d0*d1i2
                c3 = 29d0/2d0*d1i2
                c4 = -4d0/3d0*d1i2
                c5 = 1d0/12d0*d1i2
                
                darr(n1-1)=c1*arr(n1-1)+c2*arr(n1-2)+c3*arr(n1-3)+c4*arr(n1-4)+c5*arr(n1-5)
            EndElse
            
            If (dtype eq 1) Then Begin
                darr(0) = darr(0)-alpha*(2d0*b1(0)-27d0*b1(1)+270d0*b1(2)-490d0*b1(3)+270d0*arr(0)-27d0*arr(1)+2d0*arr(2))*d1i2/180d0
            EndIf

        EndIf
            
        If ((dtype eq 2) or (dtype eq 0)) Then Begin

            c1=99d0/80d0*d1i2
            c2=-3d0*d1i2
            c3=93d0/40d0*d1i2
            c4=-0.6d0*d1i2
            c5=3d0/80d0*d1i2

            darr(1) = c1*arr(0) + c2*arr(1) + c3*arr(2) + c4*arr(3) + c5*arr(4) 
             
            gamma2=0.1d0
            alpha2=-7d0/20d0

            uppr(1)=alpha2
            lowr(1)=gamma2

            If ((ibc mod 2) eq 0) Then Begin
                uppr(0) = 6d0
                c1=89d0/18d0*d1i2
                c2=-12d0*d1i2
                c3=15d0/2d0*d1i2
                c4=-4d0/9d0*d1i2
                c5=0d0

                If (ibc eq 0) Then Begin
                    c6=0d0
                Endif Else Begin
                    c6=-5d0/3d0*d1i
                EndElse

                darr(0)=c1*arr(0)+c2*arr(1)+c3*arr(2)+c4*arr(3)+c5*arr(4)+c6*darr(0)
            Endif Else Begin
                ; fourth order boundary
                alpha1=10d0
                uppr(0) = alpha1
                c1 = 145d0/12d0*d1i2
                c2 = -76d0/3d0*d1i2
                c3 = 29d0/2d0*d1i2
                c4 = -4d0/3d0*d1i2
                c5 = 1d0/12d0*d1i2
        
                darr(0)=c1*arr(0)+c2*arr(1)+c3*arr(2)+c4*arr(3)+c5*arr(4)
            EndElse

            If (dtype eq 2) Then Begin
                darr(n1-1) = darr(n1-1)-alpha*(2d0*arr(n1-3)-27d0*arr(n1-2)+270d0*arr(n1-1)-490d0*b2(0)+270d0*b2(1)-27d0*b2(2)+2d0*b2(3))*d1i2/180d0
            EndIf

        EndIf

        lowr(0) = 0d0
        uppr(n1-1) = 0d0

        ;stop

        darr = tri_solver(lowr,diag,uppr,darr,n1,is=2)

    EndIf Else Begin
         lowr(0) = 0d0
         uppr(n1-1) = 0d0
                                ;Dtype = 3 (internal boundaries only)        
         darr(0) = darr(0)-alpha*(2d0*b1(0)-27d0*b1(1)+270d0*b1(2)-490d0*b1(3)+270d0*arr(0)-27d0*arr(1)+2d0*arr(2))*d1i2/180d0
         darr(n1-1) = darr(n1-1)-alpha*(2d0*arr(n1-3)-27d0*arr(n1-2)+270d0*arr(n1-1)-490d0*b2(0)+270d0*b2(1)-27d0*b2(2)+2d0*b2(3))*d1i2/180d0

         darr = tri_solver(lowr,diag,uppr,darr,n1,is=0)
    EndElse

    ;stop

    return, darr*dri^2

end
