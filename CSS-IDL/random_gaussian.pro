pro random_gaussian, ngr, lower, upper
    values = dblarr(ngr)
    jj = 0L
    ii = 0L
    while (ii le ngr-1) do begin
       x = randomu(long(1234123489d0*randomu(1234218L+jj,1,/DOUBLE,/UNIFORM)),1,/DOUBLE,/UNIFORM)
       y = randomu(long(12343489d0*randomu(123217L+jj,1,/DOUBLE,/UNIFORM)),1,/DOUBLE,/UNIFORM)
       x = 2d0*(x-0.5d0)
       y = 2d0*(y-0.5d0)
       s = x*x+y*y
       If (s lt 1d0) Then Begin
          tmp = x*sqrt(-2d0*alog(s)/s)
          If ((tmp ge lower)and(tmp le upper)) Then Begin
             values[ii] = tmp
             ii=ii+1L
          EndIf
       Endif
       jj = jj+1L
       ;Endif
    endwhile
stop
end

pro random_exp, ngr, lower, upper
    values = dblarr(ngr)
    jj = 0L
    ii = 0L
    while (ii le ngr-1) do begin
       x = randomu(long(1234123489d0*randomu(1234218L+jj,1,/DOUBLE,/UNIFORM)),1,/DOUBLE,/UNIFORM)
       tmp = -alog(x)
       If ((tmp ge lower)and(tmp le upper)) Then Begin
          values[ii] = tmp
          ii=ii+1L
       EndIf
       jj = jj+1L
       ;Endif
    endwhile
stop
end
