pro test_part, nr, nnr

    mynr = nr/nnr
    myd1 = intarr(nnr)
    myd2 = intarr(nnr)
    mynd = intarr(nnr)

    irem = (nr Mod nnr)

    myd1(0) = 1
    myd2(0) = mynr+1
    mynd(0) = mynr+1
    for ii=1,nnr-1 do begin
        If (ii le (irem-1)) Then Begin
            mynd(ii) = mynr+1
        Endif Else Begin
            mynd(ii) = mynr
        Endelse
        myd1(ii) = myd2(ii-1)+1
        myd2(ii) = myd1(ii)+mynd(ii)-1
    Endfor
stop
end
