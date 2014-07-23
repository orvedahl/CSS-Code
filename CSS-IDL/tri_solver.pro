function tri_solver,lowr,diag,uppr,rhs,n,is=is

    If (not keyword_set(is)) Then is = 0

    If (is eq 3) Then Begin
        bp = dblarr(n)
        cp = dblarr(n-1)
        dp = dblarr(n)
        result = dblarr(n)
        
        bp(0) = diag(0)
        cp(0) = uppr(0)
        dp(0) = rhs(0)

        for i=1,n-2 do begin
            cp(i) = uppr(i)*bp(i-1)
            bp(i) = diag(i)*bp(i-1)-cp(i-1)*lowr(i)
            dp(i) = rhs(i)*diag(i-1)-dp(i-1)*lowr(i)
        endfor
        
        bp(n-1) = diag(n-1)*bp(n-2)-cp(n-2)*lowr(n-1)
        dp(n-1) = rhs(n-1)*diag(n-2)-dp(n-2)*lowr(n-1)

        result(n-1) = dp(n-1)/bp(n-1)
        for i=n-2,2,-1 do begin
            result(i) = (dp(i)-cp(i)*result(i+1))/bp(i)
        endfor

        ;zero pivot
        result(1) = 11d0/2d0*rhs(2)-11d0/2d0*result(2)-result(3)
        result(0) = rhs(0)-10d0*result(1)
        stop

    Endif Else If (is eq 2) Then Begin
        bp = dblarr(n)
        cp = dblarr(n)
        dp = dblarr(n)
        result = dblarr(n)

        dp(0)=rhs(0)/diag(0)
        bp(0)=1d0

        bet = diag(1)
        dp(1)=rhs(1)/bet
        cp(1)=uppr(0)/bet
        bp(1)=lowr(1)/bet
        for i=2,n-1 do begin
            cp(i)=uppr(i-1)/bet
            bet=diag(i)-lowr(i)*cp(i)
            dp(i)=(rhs(i)-lowr(i)*dp(i-1))/bet
            bp(i)=-bp(i-1)*lowr(i)/bet
        endfor
        
        for i=n-2,0,-1 do begin
            dp(i)=dp(i)-cp(i+1)*dp(i+1)
            bp(i)=bp(i)-cp(i+1)*bp(i+1)
        endfor

        result(0)=dp(0)/bp(0)
        for i=1,n-1 do begin
            result(i)=dp(i)-result(0)*bp(i)
        endfor

    Endif Else If (is eq 1) Then Begin
        bp = dblarr(n)
        cp = dblarr(n)
        dp = dblarr(n)
        result = dblarr(n)

        dp(n-1)=rhs(n-1)/diag(n-1)
        bp(n-1)=1d0
        bet = diag(n-2)
        dp(n-2)=rhs(n-2)/bet
        cp(n-2)=uppr(n-1)/bet
        bp(n-2)=lowr(n-2)/bet
        for i=n-3,0,-1 do begin
            cp(i)=uppr(i+1)/bet
            bet=diag(i)-lowr(i)*cp(i)
            dp(i)=(rhs(i)-lowr(i)*dp(i+1))/bet
            bp(i)=-bp(i+1)*lowr(i)/bet
        endfor
        
        for i=1,n-1 do begin
            dp(i)=dp(i)-cp(i-1)*dp(i-1)
            bp(i)=bp(i)-cp(i-1)*bp(i-1)
        endfor

        result(n-1)=dp(n-1)/bp(n-1)
        for i=n-2,0,-1 do begin
            result(i)=dp(i)-result(n-1)*bp(i)
        endfor

    Endif Else Begin
        cp = dblarr(n-1)
        dp = dblarr(n)
        result = dblarr(n)
        
        cp(0) = uppr(0)/diag(0)
        dp(0) = rhs(0)/diag(0)

        for i=1,n-2 do begin
            cp(i) = uppr(i)/(diag(i)-cp(i-1)*lowr(i))
            dp(i) = (rhs(i)-dp(i-1)*lowr(i))/(diag(i)-cp(i-1)*lowr(i))
        endfor

        dp(n-1) = (rhs(n-1)-dp(n-2)*lowr(n-1))/(diag(n-1)-cp(n-2)*lowr(n-1))
        
        result(n-1) = dp(n-1)
        for i=n-2,0,-1 do begin
            result(i) = dp(i)-cp(i)*result(i+1)
        endfor

    EndElse

    return, result

end
