!
! functions to perform integration
!
! R. Orvedahl 11-6-2014

!=========================================================================
! main interface to 1D integration routines
!=========================================================================
subroutine integrate1D(integral, npts, func, x, dx, method)

   ! func   -->  function to integrate (1d array)
   ! x      -->  x values (1d array)
   ! npts   -->  size of the function and x value arrays
   ! dx     -->  specify the uniform spacing of the x-array
   ! method -->  what method to use when integrating
   !
   ! integral --> output: value of the integral

   double precision, intent(in) :: dx
   integer, intent(in) :: npts
   character(len=*), intent(in), optional :: method
   double precision, dimension(0:npts-1), intent(in) :: func, x

   double precision, intent(out) :: integral

!f2py intent(in) :: dx, npts, func, x
!f2py intent(in), optional :: method
!f2py depends(npts) :: func, x
!f2py intent(out) :: integral

   character(len=32) :: int_method

   if (present(method) .and. (method /= "simp")) then
      ! method was passed and holds value other then default
      int_method = method
   else
      ! method was not passed or it has the default value
      int_method = "simp"
   endif

   ! call the appropriate routine
   if (int_method == "simp") then
      call simpson_integrate(integral, npts, func, dx)

   elseif (int_method == "trap") then
      call trap_integrate(integral, npts, func, x)

   else
      write(*,*)
      write(*,*) "---ERROR: unknown integration method (f90): ",int_method
      write(*,*)
      stop
   endif

   return

end subroutine integrate1D


!=========================================================================
! simpson integration 1D
!=========================================================================
subroutine simpson_integrate(integral, npts, func, dx)

   integer, intent(in) :: npts
   double precision, intent(in) :: dx
   double precision, dimension(0:npts-1), intent(in) :: func
   double precision, intent(out) :: integral

!f2py intent(in) :: npts, func, dx
!f2py depends(npts) :: func
!f2py intent(out) :: integral

   integer :: nslabs, i
   logical :: odd
   double precision :: third

   third = 1.d0/3.d0

   nslabs = npts - 1

   integral = 0.0d0

   if (mod(nslabs,2)==0) then
      odd = .false.
   else
      odd = .true.
   endif

   i = 0
   do while (i <= nslabs - 2)
      ! simpson integration over even number of slabs
      integral = integral + third*dx*(func(i) + 4.0d0*func(i+1) + func(i+2))
      i = i + 2
   enddo

   if (odd) then
      integral = integral + dx/12.0d0*(&
                 -func(nslabs-2) + 8.0d0*func(nslabs-1) + 5.0d0*func(nslabs))
   endif

   return

end subroutine simpson_integrate


!=========================================================================
! trapezoidal rule integration 1D
!=========================================================================
subroutine trap_integrate(integral, npts, func, x)

   integer, intent(in) :: npts
   double precision, dimension(0:npts-1), intent(in) :: func, x
   double precision, intent(out) :: integral

!f2py intent(in) :: npts, func, x
!f2py depends(npts) :: func, x
!f2py intent(out) :: integral

   integer :: i, nslabs

   integral = 0.0d0

   nslabs = npts - 1

   i = 0
   do while (i < nslabs)
      !                     <----width----> <----height in middle---->
      integral = integral + (x(i+1) - x(i))*0.5d0*(func(i) + func(i+1))
      i = i + 1
   enddo

   return

end subroutine trap_integrate


!=========================================================================
! integrate 1./v(r) over r
!=========================================================================
subroutine integrate_inv_vr(integral, npts, vr, rad, method)

   integer, intent(in) :: npts
   character(len=*), intent(in), optional :: method
   double precision, dimension(0:npts-1), intent(in) :: vr, rad
   double precision, intent(out) :: integral

!f2py intent(in) :: npts, vr, rad
!f2py intent(in), optional :: method
!f2py depends(npts) :: vr, rad
!f2py intent(out) :: integral

   double precision, dimension(0:npts-1) :: intgrnd

   intgrnd = vr

   ! avoid division by zero
   where(intgrnd == 0.0d0)
      intgrnd = huge(1.0d0)
   endwhere

   intgrnd = 1.0d0 / intgrnd

   dr = rad(1) - rad(0)
   call integrate1D(integral, npts, vr, rad, dr, method)

   return

end subroutine integrate_inv_vr


