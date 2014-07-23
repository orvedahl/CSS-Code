function build_cd_mat, n
   mat = dblarr(n,n)

   for ii=0L,n-2L do begin
       for jj=0L,(n-ii-2L)/2L do begin
           mat[ii,ii+2L*jj+1L] = 2d0*double(ii+2L*jj+1L)
       endfor
       mat[ii,*] = cheby_trans(mat[ii,*])
   endfor

   return, mat
end

function int_dr, x, ro, ri
  nr = n_elements(x)
  intx = dblarr(nr)
  intx(*) = 0d0
  scaling = 0.5d0*(ro-ri)
  for ir=2L,nr-2L do begin
      intx[ir] = 0.5d0*(x[ir-1L]-x[ir+1L])/double(ir)
  endfor
  intx[nr-1L] = 0.5d0*x[nr-2L]/double(nr-1L)
  ;intx[2L] = 0.5d0*x[1L]-0.25d0*x[3L]
  intx[1L] = x[0L]-0.5d0*x[2L]
  intx[0L] = 0.5d0*x[1L]
  intx = scaling*intx
  return, intx
end

function nearest, x,sign=sign,eps=eps
  If (not keyword_set(eps)) Then Begin
     ;Find machine epsilon for doubles
      meps=(machar(/double)).machep
      eps=2d0^(meps)
  Endif

  If (keyword_set(sign)) Then Begin 
      res=x-eps
  Endif Else Begin
      res=x+eps
  EndElse
  return,res 
end

function find_n0, x, eps=eps
  nc = n_elements(x)
  mxv = max(abs(x))
  n0  = nc-1L
  nloop = 0L
  ii = n0
  while ((ii gt 0)and(nloop lt nc)) do begin
      tmp=abs(x[ii])/mxv
      tmp1=nearest(tmp,/sign, eps=eps)
      tmp2=nearest(tmp,eps=eps)
      iprev = ii
      If ((tmp1 lt 0d0)and(tmp2 gt 0d0)) Then Begin
          n0 = ii
          ii = ii-1L
      Endif
      If (ii eq iprev) Then nloop = nc+1
  endwhile
  ;print, n0

  return, n0
end

function d_by_dr, x, ro, ri
 
  ri1 = ri
  nr1 = n_elements(x)
  
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

  return, dx

end


;function d_by_dr, x, dx, ro, ri, eps=eps
;  nx = n_elements(x)
;  nd = n_elements(dx(0,*));;

;  x = 1d0*x
;  dx(*,*) = 0d0

;  scaling = 2d0/(ro-ri)
;  n0=find_n0(x,eps=eps)
  ;n0 = min([nx-2L,n0])
;  n0 = nx-1L
;  cnt = 2d0*double(n0)
  ;dx(n0:nx-1L,0) = 0d0
;  dx(n0,0) = 0d0
;  dx(n0-1L,0) = cnt*x(n0)
;  for ii=n0-2L,0L,-1L do begin
;      cnt = cnt-2d0
;      dx(ii,0) = dx(ii+2L,0) + cnt*x(ii+1L)
;  endfor

;  If (nd gt 2) Then Begin
;      for jj=1,nd-1 do begin
;          n0 = find_n0(dx(*,jj-1),eps=eps)
;          n2 = min([n0,nx-jj-1])
;          cnt = 2d0*double(n2)
;          dx(n2:nx-1L,jj) = 0d0
;          dx(n2-1L,jj) = cnt*dx(n2,jj-1)
;          for ii=n2-2L,0L,-1L do begin
;              cnt = cnt-2d0
;              dx(ii,jj) = dx(ii+2L,jj) + cnt*dx(ii+1L,jj-1)
;          endfor
;      endfor
;  EndIf;

;  for jj=0,nd-1 do begin
;      dx(*,jj) = dx(*,jj)*scaling^(jj+1)
;      dx(0,jj) = dx(0,jj)/2d0
;  endfor

;  return, 1
;end

function dual_intdr, x, nr2, ro, ri1, ri
    n = n_elements(x)
    nr1 = n-nr2
    y = dblarr(n)
    y(0:nr1-1) = int_dr(x(0:nr1-1),ro,ri1)
    y(nr1:*) = int_dr(x(nr1:*),ri1,ri)
    return, y
end

function dual_dbydr, x, nr2, ro, ri1, ri
    n = n_elements(x)
    nr1 = n-nr2
    y = dblarr(n)
    res = d_by_dr(x(0:nr1-1),ro,ri1)
    y(0:nr1-1) = res
    res = d_by_dr(x(nr1:*),ri1,ri)
    y(nr1:*) =res
    return, y
end

function dual_ct, x, nr2, inverse=inverse
    n = n_elements(x)
    nr1 = n-nr2
    y = dblarr(n)
    y(0:nr1-1) = cheby_trans(x(0:nr1-1),inverse=inverse)
    y(nr1:*) = cheby_trans(x(nr1:*),inverse=inverse)
    return, y
end

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

   	   ;For a forward transform compute Ajk##X
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

            ;For a forward transform compute Ajk##X##Transpose(Ajk)
            Return,Ajk##x##Transpose(Ajk)
	End

        Else: Message,'Input must be a square array or a vector.'

    EndCase
end

