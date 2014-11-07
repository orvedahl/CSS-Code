!
! functions to perform shell averages and angular averages
!
! R. Orvedahl 11-6-2014

module shell_avgf90

   use integratef90

   implicit none

   contains

   !=========================================================================
   ! shell average of values: integrate over theta & phi
   !=========================================================================
   subroutine full_shell_avg(avgvalues, nr, nth, nphi, nq, method, &
                                  values, theta, phi)

      ! avgvalues --> output: shell avg over theta and phi of values
      ! values    --> data to average
      ! theta     --> 1D array of theta values
      ! phi       --> 1D array of phi values
      ! nr, nth, nphi, nq --> defines size/shape of values
      ! method -->  what method to use when integrating
      !
      ! idea:               integral( data sin(th) dth dphi)
      !       avgdata(r) = ----------------------------------
      !                        integral( sin(th) dth dphi)
      !
      ! the r**2 in the area element cancels, since both int are over angles
      !
      ! int( sin(th) dth dphi) = (phi2 - phi1)*(cos(th1) - cos(th2)

      integer, intent(in) :: nr, nth, nphi, nq
      character(len=*), intent(in), optional :: method
      double precision, intent(in) :: values(0:nth-1, 0:nphi-1, 0:nr-1, 0:nq-1)
      double precision, intent(in) :: theta(0:nth-1), phi(0:nphi-1)

      double precision, intent(out) :: avgvalues(0:nr-1, 0:nq-1)

   !f2py intent(in) :: nr, nth, nphi, nq, values, theta, phi
   !f2py intent(in), optional :: method
   !f2py depends(nth, nphi, nr, nq) :: values
   !f2py depends(nth) :: theta
   !f2py depends(nphi) :: phi
   !f2py depends(nr, nq) :: avgvalues
   !f2py intent(out) :: avgvalues

      character(len=32) :: int_method
      double precision :: phihi, philo, thhi, thlo, dphi, dth
      double precision :: area, invarea
      double precision :: stheta(0:nth-1)
      double precision :: tmp(0:nth-1, 0:nr-1, 0:nq-1), val
      double precision :: vals(0:nth-1, 0:nphi-1, 0:nr-1, 0:nq-1)
      integer :: i,j,k

      if (present(method) .and. (method /= "")) then
         ! method was passed and holds value other then default
         int_method = method
      else
         ! method was not passed or it has the default value
         int_method = "simp"
      endif

      ! get integral bounds
      phihi = phi(nphi-1)
      philo = phi(0)
      thhi  = theta(nth-1)
      thlo  = theta(0)
      dphi  = phi(1) - phi(0)
      dth   = theta(1) - theta(0)

      stheta(:) = dsin(theta)

      ! denominator of average
      area = (phihi - philo)*(dcos(thlo) - dcos(thhi))
      invarea = 1.d0/area

      ! multiply by sin(th)
      do k=0,nq-1
         do j=0,nr-1
            do i=0,nphi-1
               vals(:,i,j,k) = stheta(:)*values(:,i,j,k)
            enddo
         enddo
      enddo

      ! integrate over phi
      do k=0,nq-1
         do j=0,nr-1
            do i=0,nth-1
               call integrate1D(val, nphi, vals(i,:,j,k), phi, dphi, &
                                int_method)
               tmp(i,j,k) = val
            enddo
         enddo
      enddo

      ! integrate over theta
      do j=0,nq-1
         do i=0,nr-1
            call integrate1D(val, nth, tmp(:,i,j), theta, dth, int_method)
            avgvalues(i,j) = val
         enddo
      enddo

      ! normalize by area
      avgvalues(:,:) = invarea*avgvalues(:,:)

      return

   end subroutine full_shell_avg


end module shell_avgf90


