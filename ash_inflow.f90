Module ASH_Inflow

  Use Constants
  Use Derivatives
  Use Interpolation

  Implicit None

  Integer, Parameter :: TrueSingle    = Selected_real_kind(6,15)
  Integer, Public :: ashnth, ashnphi, ashrecl !Store iashtime during checkpointing for restart
  Real*8 :: ashstartt
  Real*8, Dimension(:), Allocatable :: ashtimes, ashtheta, ashphi
  Real*8, Dimension(:,:,:), Allocatable, Public :: ashvars  
  Real*8, Dimension(:,:,:,:), Allocatable :: ashslices

Contains

  Subroutine Initialize_ASH_Inflow()
    Implicit None
    Real*8, Dimension(:,:,:), Allocatable :: fullslice
    Real(TrueSingle), Dimension(:,:,:,:), Allocatable :: readdata
    Real*8, Dimension(:), Allocatable :: costheta, buf1d
    Real*8 :: tmp
    Integer :: iq, ip, it, is, unit=11

    If (.not. Start_From_Checkpoint) Then
       iashtime = 3
    EndIf
    ashrecl = ashnth*(ashnphi+1)
    Allocate(ashtimes(4),ashtheta(ashnth),ashphi(ashnphi))

    odphi = (phi2-phi1)/Dble(ashnphi)
    ashphi = odphi*(/(Dble(ip),ip=0,ashnphi-1)/)+phi1

    Allocate(oldtheta(ashnth),oldphi(ashnphi))
    oldphi = ashphi
    onth = ashnth
    onphi = ashnphi

    !Read first 4 records for cubic interpolation
    Allocate(ashslices(4,mynth,mynphi,5),fullslice(nth,nphi,5),costheta(ashnth))
    If (mybot_rank .eq. 0) Then
       Allocate(readdata(ashnth,ashnphi+1,1,5),olddata(ashnth,ashnphi,5))
       Open(unit,file='ash_inflow_data',form='unformatted',access='direct',status='old',recl=ashrecl)
       Read(unit,rec=6) readdata(:,:,1,1) 
       ashstartt = readdata(6,1,1,1)
       Close(unit)
    EndIf

    Do is=1,4
       fullslice = 0d0
       If (mybot_rank .eq. 0) Then
          Open(unit,file='ash_inflow_data',form='unformatted',access='direct',status='old',recl=ashrecl)
          readdata = 0.0
          do iq=1,5
             Read(unit,rec=((iashtime-3)+is-1)*5+iq) readdata(:,:,1,iq)
          enddo
          olddata = readdata(:,2:,1,:)
          costheta = readdata(:,1,1,2)          
          ashtimes(is) = readdata(6,1,1,1)
          Close(unit)

          ashtheta = acos(costheta)
          odth = ashtheta(ashnth/2)-ashtheta(ashnth/2-1)
          oldtheta = ashtheta

          !Do spatial interpolation
          Do iq=1,5
             Call Interp2d(fullslice(:,:,iq),olddata(:,:,iq),ashphi,ashtheta,phi,theta)
          EndDo

          !Apodize v_theta, recall that ash data is lain out as vr,vtheta,vphi,entropy,density
          !Account for changes in the dynamic pressure in entropy (gives higher entropy at higher lats)
          Do it=1,16
             tmp = sin(0.5d0*pi*(theta(it)-th1)/(theta(16)-th1))**2
             fullslice(it,:,4) = fullslice(it,:,4) + (Cv*(1d0-tmp)*rhor(1)/pr(1))*fullslice(it,:,2)**2
             fullslice(nth-it+1,:,4) = fullslice(nth-it+1,:,4) + (Cv*(1d0-tmp)*rhor(1)/pr(1))*fullslice(nth-it+1,:,2)**2
             fullslice(it,:,2) = tmp*fullslice(it,:,2)
             fullslice(nth-it+1,:,2) = tmp*fullslice(nth-it+1,:,2)
          EndDo
       EndIf
       Call Barrier(bot_comm)
       
       Call MPI_BCAST(fullslice(1,1,1),nth*nphi*5,MPI_REAL8,0,bot_comm,ierr)
       ashslices(is,:,:,:) = fullslice(myth1:myth2,myphi1:myphi2,:)
    EndDo

    Deallocate(fullslice)
    If (mybot_rank .eq. 0) Then
       Deallocate(readdata,olddata)
    EndIf

    Allocate(buf1d(ashnth+4))
    If (mybot_rank.eq.0) Then
       ashtimes = ashtimes-ashstartt
       buf1d(1:4) = ashtimes
       buf1d(5:) = costheta
    Else
       buf1d = 0d0
    EndIf
    Call MPI_Bcast(buf1d,ashnth+4,MPI_REAL8,0,bot_comm,ierr)
    costheta = buf1d(5:)
    ashtimes = buf1d(1:4)
    ashtheta = acos(costheta)
    odphi = (phi2-phi1)/Dble(ashnphi)
    ashphi = odphi*(/(Dble(ip),ip=0,ashnphi-1)/)+phi1
    Deallocate(buf1d)
    Deallocate(costheta)

    odth = ashtheta(ashnth/2)-ashtheta(ashnth/2-1)
    oldtheta = ashtheta
    oldphi = ashphi
    onth = ashnth
    onphi = ashnphi

    Allocate(ashvars(mynth,mynphi,5))
  End Subroutine Initialize_ASH_Inflow

  Subroutine Finalize_ASH_Inflow()
    Implicit None
    Deallocate(ashtimes,ashtheta,ashphi,ashslices,ashvars)
  End Subroutine Finalize_ASH_Inflow

  Subroutine Read_ASH_Inflow(tc)
    Implicit None
    Real*8, Intent(In) :: tc !Current sim time
    Real*8, Dimension(:,:,:), Allocatable :: fullslice
    Real*4, Dimension(:,:,:,:), Allocatable :: readdata
    Real*8 :: tmp, sndtmp(1)
    Integer :: iq, it, unit=11

    !Check if we need a new slice
    If (tc .ge. ashtimes(3)) Then
       !Shift slices, if necessary
       ashtimes(1) = ashtimes(2)
       ashtimes(2) = ashtimes(3)
       ashtimes(3) = ashtimes(4)
       ashslices(1,:,:,:) = ashslices(2,:,:,:)
       ashslices(2,:,:,:) = ashslices(3,:,:,:)
       ashslices(3,:,:,:) = ashslices(4,:,:,:)

       !Step forward in ASH time
       iashtime = iashtime + 1

       !Read next time slice
       Allocate(fullslice(nth,nphi,5))
       fullslice = 0d0
       If (mybot_rank .eq. 0) Then
          Allocate(readdata(ashnth,ashnphi+1,1,5),olddata(ashnth,ashnphi,5))
          
          Open(unit,file='ash_inflow_data',form='unformatted',access='direct',status='old',recl=ashrecl)
          Do iq=1,5
             Read(unit,rec=(iashtime)*5+iq) readdata(:,:,1,iq)
          EndDo
          olddata = readdata(:,2:,1,:)
          ashtimes(4) = readdata(6,1,1,1)-ashstartt
          Deallocate(readdata)
          Close(unit)

          !Do spatial interpolation
          Do iq=1,5
             Call Interp2d(fullslice(:,:,iq),olddata(:,:,iq),ashphi,ashtheta,phi,theta)
          EndDo
          Deallocate(olddata)

          !Apodize v_theta, recall that ash data is lain out as vr,vtheta,vphi,entropy,density
          !Account for changes in the dynamic pressure in entropy (gives higher entropy at higher lats)
          Do it=1,16
             tmp = sin(0.5d0*pi*(theta(it)-th1)/(theta(16)-th1))**2
             fullslice(it,:,4) = fullslice(it,:,4) + (Cv*(1d0-tmp)*rhor(1)/pr(1))*fullslice(it,:,2)**2
             fullslice(it,:,2) = tmp*fullslice(it,:,2)
             fullslice(nth-it+1,:,4) = fullslice(nth-it+1,:,4) + (Cv*(1d0-tmp)*rhor(1)/pr(1))*fullslice(nth-it+1,:,2)**2
             fullslice(nth-it+1,:,2) = tmp*fullslice(nth-it+1,:,2)
          EndDo
       EndIf
       Call Barrier(bot_comm)

       Call MPI_BCAST(fullslice(1,1,1),nth*nphi*5,MPI_REAL8,0,bot_comm,ierr)
       ashslices(4,:,:,:) = fullslice(myth1:myth2,myphi1:myphi2,:)
       Deallocate(fullslice)

       If (mybot_rank .eq. 0) Then
          sndtmp(1) = ashtimes(4)
       Else
          sndtmp(1) = 0d0
       EndIf

       Call MPI_Bcast(sndtmp,1,MPI_REAL8,0,bot_comm,ierr)
       ashtimes(4) = sndtmp(1)
    EndIf

    !Interpolate in time
    Call Cubic_Time_Interp(tc)
  End Subroutine Read_ASH_Inflow

  Subroutine Cubic_Time_Interp(tc)
    Implicit None
    Real*8, Intent(In) :: tc
    Integer :: iq, it, ip
    Real*8 :: a0, a1, a2, a3, mu, mu2, mu3, tmp(4)

    mu = (tc-ashtimes(2))/(ashtimes(3)-ashtimes(2))
    mu2 = mu**2
    mu3 = mu*mu2
    ashvars = 0d0
    Do iq=1,5
       Do ip=1,mynphi
          Do it=1,mynth
             tmp = ashslices(:,it,ip,iq)
             a3 = tmp(2)
             a0 = -0.5d0*tmp(1) + 1.5d0*a3 - 1.5d0*tmp(3) + 0.5d0*tmp(4)
             a1 = tmp(1) - 2.5d0*a3 + 2d0*tmp(3) - 0.5d0*tmp(4)
             a2 = -0.5d0*tmp(1) + 0.5d0*tmp(3)
             ashvars(it,ip,iq) = a0*mu3+a1*mu2+a2*mu+a3
          EndDo
       EndDo
    EndDo

  End Subroutine Cubic_Time_Interp

End Module ASH_Inflow
