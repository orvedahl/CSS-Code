Function Compute_Quadratic_Smoothing, quan, q1, q2, r_ind
    n = n_elements(quan)
    ssmoo = dblarr(n)
    qsmoo = dblarr(n)
    rsmoo = dblarr(n)
    rsmoo = r_ind
    ssmoo = quan
    ssmoo(0) = q1
    ssmoo(n-1) = q2
    qsmoo(0) = q1
    qsmoo(n-1) = q2

    For ii=1,n-2 Do Begin
       ;Left
        a1 = 0d0
        b1 = 0d0
        c1 = 0d0
        If (ii gt 1) Then Begin
            y1 = ssmoo(ii-2)
            y2 = ssmoo(ii-1)
            y3 = ssmoo(ii+1)
            x1 = rsmoo(ii-2)
            x2 = rsmoo(ii-1)
            x3 = rsmoo(ii+1)
            d = (x1-x2)*(x1-x3)*(x2-x3)
            a1 = x3*x2*(x2-x3)*y1+x1*x3*(x3-x1)*y2+x1*x2*(x1-x2)*y3
            a1 = a1/d
            b1 = x3^2*(y1-y2)+x1^2*(y2-y3)+x2^2*(y3-y1)
            b1 = b1/d
            c1 = x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2)
            c1 = c1/d
        EndIf
        ;Right
        a2 = 0d0
        b2 = 0d0
        c2 = 0d0
        If (ii lt n-2) Then Begin
            y1 = ssmoo(ii-1)
            y2 = ssmoo(ii+1)
            y3 = ssmoo(ii+2)
            x1 = rsmoo(ii-1)
            x2 = rsmoo(ii+1)
            x3 = rsmoo(ii+2)
            d = (x1-x2)*(x1-x3)*(x2-x3)
            a2 = x3*x2*(x2-x3)*y1+x1*x3*(x3-x1)*y2+x1*x2*(x1-x2)*y3
            a2 = a2/d
            b2 = x3^2*(y1-y2)+x1^2*(y2-y3)+x2^2*(y3-y1)
            b2 = b2/d
            c2 = x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2)
            c2 = c2/d
        EndIf
        a1 = a1+a2
        b1 = b1+b2
        c1 = c1+c2
        qsmoo(ii) = 0.5d0*(a1+b1*rsmoo(ii)+c1*rsmoo(ii)^2)
    EndFor
    qsmoo(1) = 2d0*qsmoo(1)
    qsmoo(n-2) = 2d0*qsmoo(n-2)

    return, qsmoo
End

Pro Test_Smoothing

    nr = 128
    r1 = 6.1950362d10
    r2 = 6.8882930d10
    r_ind = (r2-r1)*dindgen(nr)/double(nr-1)+r1

    sawtooth = cos(nr*!dpi*(r_ind-r1)/(r2-r1))
    func = (r1/r_ind)^3+1d-2*sawtooth

    fsmoo = compute_quadratic_smoothing(func,1d0,(r1/r2)^3,r_ind)
    fsmoo = compute_quadratic_smoothing(fsmoo,1d0,(r1/r2)^3,r_ind)
    fsmoo = compute_quadratic_smoothing(fsmoo,1d0,(r1/r2)^3,r_ind)
    fsmoo = compute_quadratic_smoothing(fsmoo,1d0,(r1/r2)^3,r_ind)
    fsmoo = compute_quadratic_smoothing(fsmoo,1d0,(r1/r2)^3,r_ind)
    stop
End
