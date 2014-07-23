pro read_background, filename, nonuni=nonuni, r1=r1, r2=r2

   openr, fun, filename, /get_lun

   nr = 0
   readf, fun, nr, format='(i5)'

   nq = 15
   If (keyword_set(nonuni)) Then nq = 18

   data = dblarr(nr,nq)

   for iq=0,nq-1 do begin
       for ir=0,nr-1 do begin
           tmp=0d0
           readf, fun, tmp, format='(e16.9)'
           data(ir,iq) = tmp
       endfor
   endfor

   gravity = data(*,0)
   rho    = data(*,1)
   s      = data(*,2)
   dsdr   = data(*,3)
   mu     = data(*,4)
   dmudr  = data(*,5)
   opac   = data(*,6)
   dodlnd = data(*,7)
   dodlnT = data(*,8)
   ks = data(*,9)
   dksdr = data(*,10)
   k0 = data(*,11)
   dk0dr = data(*,12)
   eps = data(*,13)
   Lrad = data(*,14)

   If (keyword_set(nonuni)) Then Begin
       r = data(*,15)
       dxdr = data(*,16)
       d2xdr2 = data(*,17)
   EndIf
stop
end
