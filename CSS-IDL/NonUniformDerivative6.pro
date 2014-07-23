function nud6, rad, arr

   np = n_elements(rad)
   if (n_elements(arr) ne np) then begin
      return, -1
   endif

   r0 = max(rad)
   arr0 = max(arr)   
   r = rad(*)/r0   
   array = arr(*)/arr0

   der = dblarr(np)
   der(*) = 0.0d0

   ;Lower boundary
   d3  = r(6)-r(0)
   d2  = r(5)-r(0)
   d1  = r(4)-r(0)
   dm1 = r(3)-r(0)
   dm2 = r(2)-r(0)
   dm3 = r(1)-r(0)

   a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
   a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
   a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
   am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
   am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
   am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

   a0  = a3+a2+a1+am1+am2+am3
   a0  = -a0
   der(0) = a3*array(6)+a2*array(5)+a1*array(4)+a0*array(0)+am1*array(3)+am2*array(2)+am3*array(1)
   
   d3  = r(6)-r(1)
   d2  = r(5)-r(1)
   d1  = r(4)-r(1)
   dm1 = r(3)-r(1)
   dm2 = r(2)-r(1)
   dm3 = r(0)-r(1)

   a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
   a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
   a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
   am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
   am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
   am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

   a0  = a3+a2+a1+am1+am2+am3
   a0  = -a0
   der(1) = a3*array(6)+a2*array(5)+a1*array(4)+a0*array(1)+am1*array(3)+am2*array(2)+am3*array(0)

   d3  = r(6)-r(2)
   d2  = r(5)-r(2)
   d1  = r(4)-r(2)
   dm1 = r(3)-r(2)
   dm2 = r(1)-r(2)
   dm3 = r(0)-r(2)

   a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
   a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
   a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
   am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
   am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
   am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

   a0  = a3+a2+a1+am1+am2+am3
   a0  = -a0
   der(2) = a3*array(6)+a2*array(5)+a1*array(4)+a0*array(2)+am1*array(3)+am2*array(1)+am3*array(0)
   
   ;Interior
   for i=3L,np-4L do begin
       d3  = r(i+3)-r(i)
       d2  = r(i+2)-r(i)
       d1  = r(i+1)-r(i)
       dm1 = r(i-1)-r(i)
       dm2 = r(i-2)-r(i)
       dm3 = r(i-3)-r(i)

       a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
       a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
       a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
       am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
       am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
       am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

       a0  = a3+a2+a1+am1+am2+am3
       a0  = -a0
       ;c1  = -d3*d2*d1*dm1*dm2*dm3

       der(i) = a3*array(i+3)+a2*array(i+2)+a1*array(i+1)+a0*array(i)+am1*array(i-1)+am2*array(i-2)+am3*array(i-3)
   endfor

   ;Upper boundary
   d3  = r(np-7L)-r(np-3L)
   d2  = r(np-6L)-r(np-3L)
   d1  = r(np-5L)-r(np-3L)
   dm1 = r(np-4L)-r(np-3L)
   dm2 = r(np-2L)-r(np-3L)
   dm3 = r(np-1L)-r(np-3L)

   a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
   a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
   a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
   am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
   am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
   am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

   a0  = a3+a2+a1+am1+am2+am3
   a0  = -a0
   der(np-3L) = a3*array(np-7L)+a2*array(np-6L)+a1*array(np-5L)+a0*array(np-3L)+am1*array(np-4L)+am2*array(np-2L)+am3*array(np-1L)

   d3  = r(np-7L)-r(np-2L)
   d2  = r(np-6L)-r(np-2L)
   d1  = r(np-5L)-r(np-2L)
   dm1 = r(np-4L)-r(np-2L)
   dm2 = r(np-3L)-r(np-2L)
   dm3 = r(np-1L)-r(np-2L)

   a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
   a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
   a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
   am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
   am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
   am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

   a0  = a3+a2+a1+am1+am2+am3
   a0  = -a0
   der(np-2L) = a3*array(np-7L)+a2*array(np-6L)+a1*array(np-5L)+a0*array(np-2L)+am1*array(np-4L)+am2*array(np-3L)+am3*array(np-1L)

   d3  = r(np-7L)-r(np-1L)
   d2  = r(np-6L)-r(np-1L)
   d1  = r(np-5L)-r(np-1L)
   dm1 = r(np-4L)-r(np-1L)
   dm2 = r(np-3L)-r(np-1L)
   dm3 = r(np-2L)-r(np-1L)

   a3  = -d2*d1*dm1*dm2*dm3/(d3*(d1-d3)*(d2-d3)*(d3-dm1)*(d3-dm2)*(d3-dm3))
   a2  = -d3*d1*dm1*dm2*dm3/(d2*(d2-d3)*(d2-d1)*(d2-dm1)*(d2-dm2)*(d2-dm3))
   a1  = -d3*d2*dm1*dm2*dm3/(d1*(d1-d3)*(d1-d2)*(d1-dm1)*(d1-dm2)*(d1-dm3))
   am1 = -d3*d2*d1*dm2*dm3/(dm1*(dm1-d3)*(dm1-d2)*(dm1-d1)*(dm1-dm2)*(dm1-dm3))
   am2 = -d3*d2*d1*dm1*dm3/(dm2*(dm2-d3)*(dm2-d2)*(dm2-d1)*(dm2-dm1)*(dm2-dm3))
   am3 = -d3*d2*d1*dm1*dm2/(dm3*(dm3-d3)*(dm3-d2)*(dm3-d1)*(dm3-dm1)*(dm3-dm2))

   a0  = a3+a2+a1+am1+am2+am3
   a0  = -a0
   der(np-1L) = a3*array(np-7L)+a2*array(np-6L)+a1*array(np-5L)+a0*array(np-1L)+am1*array(np-4L)+am2*array(np-3L)+am3*array(np-2L)

   der = arr0*der/r0

   return, der
end
