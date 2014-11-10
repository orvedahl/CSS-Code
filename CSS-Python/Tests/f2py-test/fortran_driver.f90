!
! fortran version of python_driver.py
!
! R. Orvedahl 11-7-2014

program integrate

   use functionf90

   implicit none

   integer :: nx
   double precision :: xlo, xhi
   character(len=8) :: method

   double precision :: numerical, exact

   nx = 127

   xlo = 0.0d0
   xhi = 1.0d0

   method = "simp"

   call def_func_and_integrate(nx, xlo, xhi, method, numerical, exact)

   write(*,*)
   write(*,*) "Pure Fortran:"
   write(*,*) nx, numerical, exact, abs(exact-numerical)
   write(*,*)

end program integrate

