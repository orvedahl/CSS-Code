!
! functions to perform integration
!
! R. Orvedahl 11-6-2014

module functionf90

   use integratef90

   implicit none

   contains

   !=========================================================================
   ! main interface to 1D integration routines
   !=========================================================================
   subroutine def_func_and_integrate(nx, xlo, xhi, method, integral, exact)

      ! nx     -->  number of grid points
      ! xlo    -->  lower boundary
      ! xhi    -->  higher boundary
      ! method -->  what method to use when integrating
      !
      ! integral --> output: numerical value of the integral
      ! exact    --> output: exact value of the integral

      integer, intent(in) :: nx
      double precision, intent(in) :: xlo, xhi
      character(len=*), intent(in) :: method

      double precision, intent(out) :: integral, exact

   !f2py intent(in) :: nx, xlo, xhi, method
   !f2py intent(out) :: integral, exact

      double precision :: dx, pi=dacos(-1.0d0)
      double precision :: fx(0:nx-1), x(0:nx-1)
      integer :: i

      ! stepsize
      dx = (xhi - xlo) / dble(nx-1)

      ! define the function
      do i=0,nx-1
         x(i) = xlo + dble(i)*dx
         fx(i) = dsin(pi*x(i))
      enddo

      ! exact integral
      exact = 1.0d0/pi*(-dcos(pi*xhi) + dcos(pi*xlo))

      ! numerically integrate
      call integrate1D(integral, nx, fx, x, dx, method)

      return

   end subroutine def_func_and_integrate

end module functionf90


