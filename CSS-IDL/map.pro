;****************************************************
   function map,a,nx,ny
;----------------------------------------------------
;
; Maps array a onto nx by ny box
;----------------------------------------------------
   s=size(a)
   nx0=s(1)
   ny0=s(2)
;
   bx=[0,1,1,0]
   by=[0,0,1,1]
;
   polywarp,nx0*bx,ny0*by,nx*bx,ny*by,1,kx,ky
;
   aa=poly_2d(a,kx,ky,1,nx,ny)
;
   return,aa
   end

