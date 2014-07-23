@tri_solver.pro
@compact_fd6.pro
@compact_fd6_2.pro

pro test_compact_diff, nr, ibc, dtype

    r1 = 0.75d0
    r2 = 1d0
    dr = (r2-r1)/double(nr-1)
    dri = 1d0/(r2-r1)
    r = dr*dindgen(nr)+r1

    ;f = r

    f = sin(2d0*!dpi*(r-r1)/(r2-r1))

    b1 = r1-dr*(dindgen(4)+1d0)
    b1 = reverse(b1)
    b1 = sin(2d0*!dpi*(b1-r1)/(r2-r1))
    b2 = r2+dr*(dindgen(4)+1d0)
    b2 = sin(2d0*!dpi*(b2-r1)/(r2-r1))

    If (dtype lt 3) Then Begin
        If (ibc eq 0) Then Begin
            b1 = 0d0
            b2 = 0d0
        Endif

        If (ibc eq 2) Then Begin
            b1 = 0d0
        Endif

        If (ibc eq 3) Then Begin
            b2 = 0d0
        Endif

        If (ibc eq 4) Then Begin
            b1 = 2d0*!dpi/dri/(r2-r1) ; 1d0/dri
            b2 = 2d0*!dpi/dri/(r2-r1) ; 1d0/dri
            df = f
            df(*) = 0d0
            df(0) = b1
            df(nr-1) = b2
        Endif
    EndIf

    If (dtype eq 2) Then Begin
        b2 = r2+dr*(dindgen(4)+1d0)
    Endif

    If (dtype eq 1) Then Begin
        b1 = r1-dr*(dindgen(4)+1d0)
        b1 = reverse(b1)
    EndIf

    df1 = compact_fd6(dri,f,b1,b2,ibc,dtype=dtype,darr=df)

    ;exact = 0d0*dblarr(nr)+1d0
    exact1 = 2d0*!dpi*cos(2d0*!dpi*(r-r1)/(r2-r1))/(r2-r1)
    local_error = df1-exact1
    global_error = sqrt(total(local_error^2))

    print, 'global error first derivative =', global_error

    ;f = r^2
    ;b1 = r1-dr*(dindgen(4)+1d0)
    ;b1 = reverse(b1)
    ;b1 = b1^2
    ;b2 = r2+dr*(dindgen(4)+1d0)
    ;b2 = b2^2

    b1 = r1-dr*(dindgen(4)+1d0)
    b1 = reverse(b1)
    b1 = sin(2d0*!dpi*(b1-r1)/(r2-r1))
    b2 = r2+dr*(dindgen(4)+1d0)
    b2 = sin(2d0*!dpi*(b2-r1)/(r2-r1))

    If (dtype lt 3) Then Begin
        If (ibc eq 0) Then Begin
            b1 = 0d0
            b2 = 0d0
        Endif

        If (ibc eq 2) Then Begin
            b1 = 0d0
        Endif

        If (ibc eq 3) Then Begin
            b2 = 0d0
        Endif

        If (ibc eq 4) Then Begin
            b1 = 2d0*!dpi/dri/(r2-r1)
            b2 = 2d0*!dpi/dri/(r2-r1)
            ;b1 = 2d0*r1/dri
            ;b2 = 2d0*r2/dri
            df = f
            df(*) = 0d0
            df(0) = b1
            df(nr-1) = b2
        Endif
    EndIf

    If (dtype eq 2) Then Begin
        b2 = r2+dr*(dindgen(4)+1d0)
        b2 = b2^2
    Endif

    If (dtype eq 1) Then Begin
        b1 = r1-dr*(dindgen(4)+1d0)
        b1 = reverse(b1)
        b1 = b1^2
    EndIf

    df2 = compact_fd6_2(dri,f,b1,b2,ibc,dtype,darr=df)

    ;exact = 0d0*dblarr(nr)+2d0
    exact2 = -((2d0*!dpi)/(r2-r1))^2*sin(2d0*!dpi*(r-r1)/(r2-r1))
    local_error = df2-exact2
    global_error = sqrt(total(local_error^2))

    print, 'global error second derivative =', global_error

    ;stop

end
