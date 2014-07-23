!***********************************************************************
!
!       Copyright 2009, Kyle Augustson
!	JILA -- University of Colorado
!
!***********************************************************************
Module Time_Step_RK4

  Use Physics
  Use Diffusion

  Implicit None

  Logical :: RK3TVD = .False., RK2TVD=.False.
  Real*8, Allocatable :: vars_save(:,:,:,:)
  !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: vars_save
Contains

  Subroutine Initialize_Step_RK4()
    Implicit None
    Allocate(vars_save(mynth,mynphi,mynr,nv))
  End Subroutine Initialize_Step_RK4

  Subroutine Finalize_Step_RK4()
    Implicit None
    Deallocate(vars_save)
  End Subroutine Finalize_Step_RK4

  !------------------------------------------------------------------------

  Subroutine Step_RK2_TVD()
    Implicit None
    Integer :: ir, it, ip, iv

    If (Do_Timings) tstart = MPI_WTIME()

    If ((Bottom_Specified).and.(myr1 .eq. 1)) Call Read_ASH_Inflow(sim_time)

    !  This is the fourth order Runge-Kutta method
    dtcur = deltat
    dvars1 = 0d0    
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf
    vars1 = vars0
    
    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(1)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    If (SLDiffusion) Then
       !Explicit Euler
       vars2 = dvars1
       dvars1 = 0d0

       If (Do_Timings) Then
          tstop = MPI_WTIME()
          timings(16) = timings(16) + tstop-tstart
       EndIf

       Call Compute_Diffusion(1)

       If (Do_Timings) Then
          tstart = MPI_WTIME()
       EndIf

       vars_save = dvars1
       dvars1 = vars2
    EndIf

    !#DIR$ VECTOR ALIGNED
    vars1 = vars0 + dtcur*dvars1

    Call Bound()

    !k2
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf

    dtcur = 0.5d0*deltat

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(2)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    !#DIR$ VECTOR ALIGNED
    Do iv=1,nv
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          Do ip=1,mynphi
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do it=1,mynth
                vars1(it,ip,ir,iv) = 0.5d0*vars0(it,ip,ir,iv) + 0.5d0*vars1(it,ip,ir,iv) + dtcur*dvars1(it,ip,ir,iv)
             EndDo
          EndDo
       EndDo
    EndDo

    If (SLDiffusion) Then
       dtcur = deltat
       
       !#DIR$ VECTOR ALIGNED
       Do iv=1,nv
          !#DIR$ VECTOR ALIGNED
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do it=1,mynth
                   vars1(it,ip,ir,iv) = vars1(it,ip,ir,iv) + dtcur*vars_save(it,ip,ir,iv)
                EndDo
             EndDo
          EndDo
       EndDo
    Endif

    Call Bound()

    vars0 = vars1

    If ((Bottom_Specified) .and. (mod(istep,ufreq).eq.0)) Then
       If (myr1 .eq. 1) Then 
          Call ASH_Enthalpy_Flux()
       EndIf
       Call MPI_Bcast(ash_fen_mod,1,MPI_REAL8,0,node_comm,ierr)
    EndIf

    If ((do_entropy_rain).and.(myr2.eq.nr).and.(mod(istep,ufreq/10).eq.0)) Then
       Call Update_Rain_Squalls()
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

  End Subroutine Step_RK2_TVD

  !------------------------------------------------------------------------

  Subroutine Step_RK3_TVD()
    Implicit None
    Integer :: ir, it, ip, iv

    If (Do_Timings) tstart = MPI_WTIME()

    If ((Bottom_Specified).and.(myr1 .eq. 1)) Call Read_ASH_Inflow(sim_time)

    !  This is the fourth order Runge-Kutta method
    dtcur = deltat
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf
    vars1 = vars0
    
    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(1)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    If (SLDiffusion) Then
       !Explicit Euler
       vars2 = dvars1
       dvars1 = 0d0

       If (Do_Timings) Then
          tstop = MPI_WTIME()
          timings(16) = timings(16) + tstop-tstart
       EndIf

       Call Compute_Diffusion(1)

       If (Do_Timings) Then
          tstart = MPI_WTIME()
       EndIf

       vars_save = dvars1
       dvars1 = vars2
    EndIf

    !#DIR$ VECTOR ALIGNED
    vars1 = vars0 + dtcur*dvars1

    Call Bound()

    !k2
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf

    dtcur = deltat/4d0

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(2)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    !#DIR$ VECTOR ALIGNED
    Do iv=1,nv
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          Do ip=1,mynphi
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do it=1,mynth
                vars1(it,ip,ir,iv) = 0.75d0*vars0(it,ip,ir,iv) + 0.25d0*vars1(it,ip,ir,iv) + dtcur*dvars1(it,ip,ir,iv)
             EndDo
          EndDo
       EndDo
    EndDo

    Call Bound()

    !k3
    dtcur = 2d0*deltat/3d0
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(3)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    !#DIR$ VECTOR ALIGNED
    Do iv=1,nv
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          Do ip=1,mynphi
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do it=1,mynth
                vars1(it,ip,ir,iv) = onethird*vars0(it,ip,ir,iv) + twothirds*vars1(it,ip,ir,iv) + dtcur*dvars1(it,ip,ir,iv)
             EndDo
          EndDo
       EndDo
    EndDo

    If (SLDiffusion) Then
       dtcur = deltat
       
       !#DIR$ VECTOR ALIGNED
       Do iv=1,nv
          !#DIR$ VECTOR ALIGNED
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do it=1,mynth
                   vars1(it,ip,ir,iv) = vars1(it,ip,ir,iv) + dtcur*vars_save(it,ip,ir,iv)
                EndDo
             EndDo
          EndDo
       EndDo
    Endif

    Call Bound()

    vars0 = vars1

    If ((Bottom_Specified) .and. (mod(istep,ufreq).eq.0)) Then
       If (myr1 .eq. 1) Then 
          Call ASH_Enthalpy_Flux()
       EndIf
       Call MPI_Bcast(ash_fen_mod,1,MPI_REAL8,0,node_comm,ierr)
    EndIf

    If ((do_entropy_rain).and.(myr2.eq.nr).and.(mod(istep,ufreq/10).eq.0)) Then
       Call Update_Rain_Squalls()
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

  End Subroutine Step_RK3_TVD

  !------------------------------------------------------------------------

  Subroutine Step_RK4()
    Implicit None
    Integer :: ir, it, ip, iv

    If (Do_Timings) tstart = MPI_WTIME()

    If ((Bottom_Specified).and.(myr1 .eq. 1)) Call Read_ASH_Inflow(sim_time)

    !  This is the fourth order Runge-Kutta method
    dtcur = deltat/2d0
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf
    vars1 = vars0
    
    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(1)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    If (SLDiffusion) Then
       !Explicit Euler
       vars2 = dvars1
       dvars1 = 0d0

       If (Do_Timings) Then
          tstop = MPI_WTIME()
          timings(16) = timings(16) + tstop-tstart
       EndIf

       Call Compute_Diffusion(1)

       If (Do_Timings) Then
          tstart = MPI_WTIME()
       EndIf

       vars_save = dvars1 + vars2/6d0
       dvars1 = vars2
    Else
       vars_save = dvars1   
    EndIf

    !#DIR$ VECTOR ALIGNED
    vars1 = vars0 + dtcur*dvars1

    Call Bound()

    !k2
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf

    dtcur = deltat/2d0

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(2)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf
    
    vars_save = vars_save + dvars1/3d0
    vars1 = vars0 + dtcur*dvars1

    Call Bound()

    !k3
    dtcur = deltat
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(3)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    vars_save = vars_save + dvars1/3d0
    vars1 = vars0 + dtcur*dvars1

    Call Bound()

    !k4
    dtcur = deltat
    dvars1 = 0d0
    If (Laplacian) Then
       vars2 = 0d0
       dvars2 = 0d0
    ElseIf(do_radheat)Then
       dvars2 = 0d0
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

    Call RHS(4)

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    vars_save = vars_save + dvars1/6d0
    vars1 = vars0 + dtcur*vars_save

    Call Bound()

    vars0 = vars1

    If ((Bottom_Specified) .and. (mod(istep,ufreq).eq.0)) Then
       If (myr1 .eq. 1) Then 
          Call ASH_Enthalpy_Flux()
       EndIf
       Call MPI_Bcast(ash_fen_mod,1,MPI_REAL8,0,node_comm,ierr)
    EndIf

    If ((do_entropy_rain).and.(myr2.eq.nr).and.(mod(istep,ufreq/10).eq.0)) Then
       Call Update_Rain_Squalls()
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(16) = timings(16) + tstop-tstart
    EndIf

  End Subroutine Step_RK4

  Subroutine ASH_Enthalpy_Flux()
    Implicit None
    Real*8 :: tconst, uravg, ntot, Len, snd(2), rcv(2)

    ntot = Dble(nphi*nth)
    snd(1) = Sum(vars1(:,:,1,4))
    snd(2) = Sum(temperature(:,:,1))
    rcv = 0d0
    Call Global_AllReduce(snd,rcv,bot_comm,'sum')
    uravg = rcv(1)/ntot
    T2avg = rcv(2)/ntot
    
    Len = Cp*Sum(vars1(:,:,1,1)*(vars1(:,:,1,4)-uravg)*(temperature(:,:,1)-T2avg))
    snd(1) = Len
    rcv = 0d0
    Call Global_AllReduce(snd,rcv,bot_comm,'sum')
    Len = rcv(1)/ntot
    Len = 4d0*pi*r1**2*Len/lumr1
    Len = 1.3d0*Min(1d0/Min(2d0,Max(0d0,abs(Len))),2d0)
    ash_fen_mod = 0.5d0*(Len+ash_fen_mod)
    tconst = Min(abs(sim_time)/1d5,1d0)
    If (tconst .lt. 1d0) ash_fen_mod = 1d0
    If (myrank .eq. 0) Print*, 'ash_fen_mod=', ash_fen_mod
  End Subroutine ASH_Enthalpy_Flux

End Module Time_Step_RK4
