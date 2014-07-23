Module Physics

  Use Constants !Flags and constants
  Use Derivatives !Derivatives and derived functions
  Use Entropy_Rain
  Use ASH_Inflow

  Implicit None

  Real*8 :: mass, m0, mtot, mrsint, Ltot, L0, E0, vol, xi_cs, divBtot, divBmax, k0expon, damping_time
  Real*8 :: cfl_safety, dtmaxb, dtmaxc, dtmaxd, deltat_old, delta_s
  Real*8, Allocatable, Dimension(:) :: Lrad_new
  Logical :: ell0_heating=.False., lrad_update=.False.
Contains

  !------------------------------------------------------------------------

  Subroutine Initialize_Physics()
    Implicit None
    !  initialize time
    sim_time=0d0
    !  allocate main data variables: temps and vars
    Allocate(temperature(mynth,mynphi,mynr),invrho(mynth,mynphi,mynr))
    Allocate(vars0(mynth,mynphi,mynr,nv),vars1(mynth,mynphi,mynr,nv),vars2(mynth,mynphi,mynr,nv),dvars1(mynth,mynphi,mynr,nv))

    If (Magnetic) Then
       Allocate(maga(mynth,mynphi,mynr,3))
    EndIf

    If (do_radheat) Then
       Allocate(dvars2(mynth,mynphi,mynr,2))
    EndIf

    !Initialize background state variables
    Allocate(gravity(nr),rhor(nr),sr(nr),tr(nr),pr(nr),dsdrm(nr),drhodrm(nr),inv_drag_time(mynr))
    Allocate(mur(nr),dlnmu(nr),ks(nr),dlnks(nr),etar(nr),eps(nr),k0(nr))
    Allocate(Lrad(nr),opac0(nr),opac(nr),dopdlnd(nr),dopdlnT(nr),kr(nr),dlnkr(nr))
    Allocate(sr1(mynth),sr2(mynth),pr1(mynth),pr2(mynth),tr1(mynth),tr2(mynth),rhor1(mynth),rhor2(mynth),dsdr1(mynth),dsdr2(mynth))
    Allocate(rhot1(mynr),tt1(mynr),pt1(mynr),st1(mynr))

    If (lrad_update) Allocate(Lrad_new(nr))
    sr=1d0
    dsdrm=1d0
    Lrad=1d0
    tr=1d0
    rhor=1d0
    pr = 1d0
    gravity=0d0
    kr=0d0
    dlnkr=0d0
    ks=0d0
    dlnks=0d0
    k0=0d0
    mur=0d0
    dlnmu=0d0
    etar=0d0
    eps=0d0
    opac=0d0
    dopdlnd=0d0
    dopdlnT=0d0
    dxdr=1d0
    d2xdr2=1d0
    Cv = Cp/gamma
    invCp = 1d0/Cp
    invCv = 1d0/Cv
    gamm1 = gamma-1d0
    If (Unstratified) Then
       pconst = g0/T0
       invpconst = T0/g0
       invtconst = 1d0
    Else
       pconst = gamm1*Cp
       invpconst = invCp/gamm1
       invtconst = invCv/gamm1
    EndIf
    divBtot = 0d0
    divBmax = 0d0
  End Subroutine Initialize_Physics

  !Clean up arrays at shutdown
  Subroutine Finalize_Physics()
    Implicit None
    Deallocate(vars0,vars1,vars2,dvars1,temperature,invrho)
    Deallocate(gravity,kr,dlnkr,mur,dlnmu,etar,tr,rhor,ks,dlnks,eps,sr,k0)
    Deallocate(dsdrm,opac,dopdlnd,dopdlnT,opac0,sr1,sr2,tr1,tr2,pr1,pr2,rhor1,rhor2,dsdr1,dsdr2)
    If (Allocated(vphi_prof)) Deallocate(vphi_prof) 
    If (do_radheat) Deallocate(dvars2)
  End Subroutine Finalize_Physics

  !========================================================================

  Subroutine Compute_Step_Size(dtmax)
    Implicit None
    Real*8, Intent(Out) :: dtmax
    Real*8 :: Etot, vmax, vamax, dxmin, rbuff(6),sbuff(6), drad, radius, sinth
    Real*8, Dimension(mynphi) :: ut, vt, wt, brt, btt, bpt, tt, rhot
    Real*8 :: nanarr(nv), tmp(nv), mins(nv), maxs(nv), gmins(nv), gmaxs(nv)
    Logical :: nan
    Integer :: ir, it, ip, iv
    dtmax=1d4
    dtmaxb=1d4
    dtmaxc=1d4

    !  Courant limits (sound speed + speed)
    If (Unstratified) Then
       temperature = T0
    Else
       temperature = invtconst*(vars0(:,:,:,1)**gamm1)*dexp(vars0(:,:,:,5)*invCv)
    EndIf
    If (magnetic) Call Compute_B(0)

    If (Unstratified) Then
       !Check for energy and angular momentum conservation
       Etot = 0d0
       Ltot = 0d0
       mtot = 0d0
       mrsint = 0d0
       !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
       Do ir=1,mynr
          radius = r_ind(myr1+ir-1)
          drad = dr_ind(myr1+ir-1)
          !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=64
          Do it=1,mynth
             sinth = sines(myth1+it-1)
             dxmin = Min(drad,radius*Min(dth,sinth*dphi)) !4
             rhot = vars0(it,:,ir,1)
             mtot = mtot + radius**2*sinth*Sum(rhot)*drad*dth*dphi !nops = 7*mynphi
             !mrsint = mrsint + radius**3*sinth**2*Sum(rhot)*drad*dth*dphi
             mrsint = mrsint + radius**3*sinth**3*Sum(rhot)*drad*dth*dphi
             wt = vars0(it,:,ir,2)**2 !2*mynphi
             vt = vars0(it,:,ir,3)
             ut = vars0(it,:,ir,4)**2 !2*mynphi
             tt = temperature(it,:,ir)
             !Ltot = Ltot + radius**3*sinth**2*Sum(rhot*vt)*drad*dth*dphi
             Ltot = Ltot + radius**3*sinth*Sum(rhot*Sqrt(wt+sinth**2*(vt+radius*sinth*omega_0)**2))*drad*dth*dphi !8+29*mynphi
             vt = vt*vt !2*mynphi
             Etot = Etot + radius**2*sinth*Sum(rhot*(Cv*tt+0.5d0*(ut+wt+vt)))*drad*dth*dphi !6+6*mynphi
             tt = g0
             ut = ut+vt+wt+tt !3*mynphi
             vmax  = Minval(dxmin/Sqrt(ut)) !Max sound speed !16*mynphi+2
             dtmaxc= Min(dtmaxc,vmax) !1
             If (magnetic) Then
                btt = vars0(it,:,ir,6)**2
                bpt = vars0(it,:,ir,7)**2
                brt = vars0(it,:,ir,8)**2
                Etot = Etot + 0.5d0*radius**2*sinth*Sum(brt+btt+bpt)*drad*dth*dphi*one_over_4pi
                vamax = Minval(dxmin/Sqrt(ut+(brt+btt+bpt)*one_over_4pi/rhot))
                dtmaxb = Min(dtmaxb,vamax)
             EndIf
          EndDo !71*mynphi+21 (hydro)
       EndDo !Total ops = mynr*mynth*(71*mynphi+21)
    Else
       !Check for energy and angular momentum conservation
       Etot = 0d0
       Ltot = 0d0
       mtot = 0d0
       mrsint = 0d0
       !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
       Do ir=1,mynr
          radius = r_ind(myr1+ir-1)
          drad = dr_ind(myr1+ir-1)
          !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=64
          Do it=1,mynth
             sinth = sines(myth1+it-1)
             dxmin = Min(drad,radius*Min(dth,sinth*dphi)) !4
             rhot = vars0(it,:,ir,1)
             mtot = mtot + radius**2*sinth*Sum(rhot)*drad*dth*dphi !nops = 7*mynphi
             !mrsint = mrsint + radius**3*sinth**2*Sum(rhot)*drad*dth*dphi
             mrsint = mrsint + radius**3*sinth**3*Sum(rhot)*drad*dth*dphi
             wt = vars0(it,:,ir,2)**2 !2*mynphi
             vt = vars0(it,:,ir,3)
             ut = vars0(it,:,ir,4)**2 !2*mynphi
             tt = temperature(it,:,ir)
             !Ltot = Ltot + radius**3*sinth**2*Sum(rhot*vt)*drad*dth*dphi
             Ltot = Ltot + radius**3*sinth*Sum(rhot*Sqrt(wt+sinth**2*(vt+radius*sinth*omega_0)**2))*drad*dth*dphi !8+29*mynphi
             vt = vt*vt !2*mynphi
             Etot = Etot + radius**2*sinth*Sum(rhot*(Cv*tt+0.5d0*(ut+wt+vt)))*drad*dth*dphi !6+6*mynphi
             tt = gamm1*Cp*xi_cs*xi_cs*tt !4*mynphi
             ut = ut+vt+wt+tt !3*mynphi
             vmax  = Minval(dxmin/Sqrt(ut)) !Max sound speed !16*mynphi+2
             dtmaxc= Min(dtmaxc,vmax) !1
             If (magnetic) Then
                btt = vars0(it,:,ir,6)**2
                bpt = vars0(it,:,ir,7)**2
                brt = vars0(it,:,ir,8)**2
                Etot = Etot + 0.5d0*radius**2*sinth*Sum(brt+btt+bpt)*drad*dth*dphi*one_over_4pi
                vamax = Minval(dxmin/Sqrt(ut+(brt+btt+bpt)*one_over_4pi/rhot))
                dtmaxb = Min(dtmaxb,vamax)
             EndIf
          EndDo !71*mynphi+21 (hydro)
       EndDo !Total ops = mynr*mynth*(71*mynphi+21)
    EndIf

    If (Laplacian) Then
       dtmaxd = Min(Minval(dr_ind),r1*Min(Min(sines(1),sines(nth))*dphi,dth))**2/Max(Maxval(ks),Maxval(mur/rhor),Maxval(etar))
    Else
       dtmaxd = 1d10
    EndIf

    !All Reduce constants
    sbuff = (/Etot,Ltot,mtot,mrsint,divBtot,1d0/)
    rbuff = 0d0
    Call Global_AllReduce(sbuff,rbuff,node_comm,'sum')
    Etot = rbuff(1)
    Ltot = rbuff(2)
    mtot = rbuff(3)
    mrsint = rbuff(4)
    divBtot = rbuff(5)

    !Print conservation properties
    If (myrank .eq. 0) Then
       Print*, 'E_tot = ', Etot, 'dE_tot = ', abs(Etot-E0)/E0
       Print*, 'L_tot = ', Ltot, 'dL_tot = ', abs(Ltot-L0)/L0
       Print*, 'm_tot = ', mtot, 'dm_tot = ', abs(mtot-m0)/m0       
       If (magnetic) Then
          Print*, '||div.B||_L1 = ', divBtot*dth*dphi/vol
       EndIf
    EndIf

    !  find global limit (store in dtmax) and apply the safety factor
    If (magnetic) Then
       sbuff = (/dtmaxb,dtmaxc,dtmaxd,1d0/divBmax,1d0,1d0/)
    Else
       sbuff = (/dtmaxb,dtmaxc,dtmaxd,1d0,1d0,1d0/)
    EndIf
    rbuff = 0d0
    Call Global_AllReduce(sbuff,rbuff,node_comm,'min')
    dtmaxb=rbuff(1)
    dtmaxc=rbuff(2)
    dtmaxd=rbuff(3)
    dtmax=Min(dtmaxb,dtmaxc,dtmaxd)*cfl_safety
    If (magnetic) Then
       divBmax = rbuff(4)
       If (myrank .eq. 0) Then
          Print*, '||div.B||_Linf = ', 1d0/divBmax
       EndIf
    EndIf

    If (debug) Then
       nanarr = 0d0
       Do iv=1,nv
          mins(iv) = Minval(vars0(:,:,:,iv))
          maxs(iv) = Maxval(vars0(:,:,:,iv))
          nan = ANY(ISNAN(vars0(:,:,:,iv)))
          If (nan) nanarr(iv) = myrank
       EndDo

       Call Global_AllReduce(nanarr,tmp,node_comm,'sum')
       nanarr=tmp
       tmp = 0d0
       Call Global_AllReduce(mins,tmp,node_comm,'min')
       gmins = tmp
       tmp = 0d0
       Call Global_AllReduce(maxs,tmp,node_comm,'max')
       gmaxs = tmp

       If (myrank .eq. 0) Then
          Print*, 'ISNAN Vars0', nanarr
          Print*, 'Min Vars0', gmins
          Print*, 'Max Vars0', gmaxs
          Print*, '            '
       EndIf
    EndIf

    !Restore vector potential
    If (magnetic) vars0(:,:,:,6:8) = maga

    !  check to see if time step is too small
    If ((sim_time.Gt.0).And.(dtmax .Lt. 1e-3)) Then
       If (myrank.Eq.0) Print*,'ERROR in Compute_Step_Size: time step too small'
       If (myrank.Eq.0) Print*,'      sim_time = ',sim_time,'  time step = ',dtmax
       Call Graceful_Exit()
    EndIf
    Call Barrier(node_comm)
  End Subroutine Compute_Step_Size

  !=======================================================================

  Subroutine RHS(call_num)
    Implicit None
    Integer :: call_num
    Real*8 :: nanarr(nv), tmp(nv), mins(nv), maxs(nv), gmins(nv), gmaxs(nv), sndr(3), rcvr(3), rtmp1, rtmp2, rstmp(mynth,mynphi,2)
    Logical :: nan
    Integer :: ir, it, ip, iv

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rstmp

    If (debug) Then
       nanarr = 0d0
       Do iv=1,nv
          mins(iv) = Minval(vars1(:,:,:,iv))
          maxs(iv) = Maxval(vars1(:,:,:,iv))
          nan = ANY(ISNAN(vars1(:,:,:,iv)))
          If (nan) Then
             nanarr(iv) = myrank
             Print*, 'NaN detected on rank', myrank, 'for variable iv=', iv
          Endif
       EndDo

       Call Global_AllReduce(nanarr,tmp,node_comm,'sum')
       nanarr=tmp
       tmp = 0d0
       Call Global_AllReduce(mins,tmp,node_comm,'min')
       gmins = tmp
       tmp = 0d0
       Call Global_AllReduce(maxs,tmp,node_comm,'max')
       gmaxs = tmp

       Call Barrier(node_comm)

       If (myrank .eq. 0) Then
          Print*, 'ISNAN Vars1', nanarr
          Print*, 'Min Vars1', gmins
          Print*, 'Max Vars1', gmaxs
          Print*, '            '
          If (Sum(nanarr) .gt. 0d0) Then
             Print*, 'NaNaNaNa NaNaNaNa Hey Hey Goodbye...'
          EndIf
       EndIf

       If (Sum(nanarr) .gt. 0d0) Call Graceful_Exit()
    EndIf

    If (Do_Timings) tstart =  MPI_WTIME()
    If (magnetic) Call Compute_B(1)

    If ((mod(istep,ufreq) .eq. 0).and.(call_num.eq.1)) Then
       mtot = 0d0
       Ltot = 0d0
       mrsint = 0d0
       !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          rtmp1 = dr_ind(myr1+ir-1)*r_ind(myr1+ir-1)**2
          rtmp2 = rtmp1*r_ind(myr1+ir-1)
          rstmp(:,:,1) = vars1(:,:,ir,1)
          rstmp(:,:,2) = vars1(:,:,ir,3)

          !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
          !#DIR$ VECTOR ALIGNED
          Do ip=1,mynphi
             mtot = mtot + Sum(sines(myth1:myth2)*rstmp(:,ip,1))*rtmp1
             mrsint = mrsint + Sum(rstmp(:,ip,1)*sines(myth1:myth2)**3)*rtmp2
             Ltot = Ltot + Sum(rstmp(:,ip,1)*rstmp(:,ip,2)*sines(myth1:myth2)**2)*rtmp2
          EndDo
       EndDo
       sndr(1) = mtot*dth*dphi
       sndr(2) = Ltot*dth*dphi
       sndr(3) = mrsint*dth*dphi
       rcvr = 0d0
       Call Global_AllReduce(sndr,rcvr,node_comm,'sum')
       mtot = rcvr(1)
       Ltot = rcvr(2)
       mrsint = rcvr(3)
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(5) = timings(5) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    !Exchange boundary info (theta)
    Call Boundary_Exchange_Init(2)

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(3) = timings(3) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf
    If (Unstratified) Then
       temperature = T0
       invrho = 1d0/vars1(:,:,:,1)
    Else
       !Compute temperature
       !#DIR$ VECTOR ALIGNED
       invrho = vars1(:,:,:,1)**gamm1
       !#DIR$ VECTOR ALIGNED
       temperature = invtconst*invrho
       !#DIR$ VECTOR ALIGNED
       invrho = vars1(:,:,:,5)*invCv
       !#DIR$ VECTOR ALIGNED
       invrho = dexp(invrho)
       !#DIR$ VECTOR ALIGNED
       temperature = temperature*invrho
       !#DIR$ VECTOR ALIGNED
       invrho = 1d0/vars1(:,:,:,1)
       If (do_radheat) Then
          Do ir=1,mynr
             dvars2(:,:,ir,1) = (4d0*a_const*c_light*kr(myr1+ir-1)/3d0)*(temperature(:,:,ir)**3)*(invrho(:,:,ir)**2)
          EndDo
       EndIf
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(5) = timings(5) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    Call Boundary_Exchange_Finish(2)

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(3) = timings(3) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    !Set boundary info
    If (myth1 .eq. 1) Then
       If (Unstratified) Then
          binfo(2,1)%bdata(1,:,:,3)=0d0
       Else
          binfo(2,1)%bdata(1,:,:,3)=cosines(1)/sines(1)*vars1(1,:,:,3)
       EndIf
       binfo(2,1)%bdata(1,:,:,4)=0d0
       binfo(2,1)%bdata(1,:,:,5)=0d0

       If ((magnetic).and.(Theta_Normal_Field)) Then
          binfo(2,1)%bdata(1,:,:,6) = 0d0
          binfo(2,1)%bdata(1,:,:,7)=-cosines(1)/sines(1)*vars1(1,:,:,7)
          binfo(2,1)%bdata(1,:,:,8) = 0d0
       EndIf
    EndIf

    If (myth2 .eq. nth) Then
       If (Theta_Symmetric) Then
          Do it=1,4
             binfo(2,2)%bdata(it,:,:,1) = vars1(mynth-it,:,:,1)
             binfo(2,2)%bdata(it,:,:,2) = -vars1(mynth-it,:,:,2)
             binfo(2,2)%bdata(it,:,:,3) = vars1(mynth-it,:,:,3)
             binfo(2,2)%bdata(it,:,:,4) = vars1(mynth-it,:,:,4)
             binfo(2,2)%bdata(it,:,:,5) = vars1(mynth-it,:,:,5)
             If (magnetic) Then
                binfo(2,2)%bdata(it,:,:,6) = -vars1(mynth-it,:,:,6) !Is this correct, field totally anti-symmtric about equator?
                binfo(2,2)%bdata(it,:,:,7) = -vars1(mynth-it,:,:,7)
                binfo(2,2)%bdata(it,:,:,8) = -vars1(mynth-it,:,:,8)
             EndIf
          EndDo
       Else
          If (Unstratified) Then
             binfo(2,2)%bdata(1,:,:,3)=0d0
          Else
             binfo(2,2)%bdata(1,:,:,3)=cosines(nth)/sines(nth)*vars1(mynth,:,:,3)
          EndIf
          binfo(2,2)%bdata(1,:,:,4)=0d0
          binfo(2,2)%bdata(1,:,:,5)=0d0

          If ((magnetic).and.(Theta_Normal_Field)) Then
             binfo(2,2)%bdata(1,:,:,6) = 0d0
             binfo(2,2)%bdata(1,:,:,7)=-cosines(nth)/sines(nth)*vars1(mynth,:,:,7)
             binfo(2,2)%bdata(1,:,:,8) = 0d0
          EndIf
       EndIf
    EndIf

    !Set temperature boundary info
    If (do_radheat) Then
       If (myth1 .eq. 1) Then
          tbcinfo(2,1)%bdata(1,:,:,1) = temperature(1,:,:)
       Else
          tbcinfo(2,1)%bdata(:,:,:,1) = (binfo(2,1)%bdata(:,:,:,1)**gamm1)*dexp(binfo(2,1)%bdata(:,:,:,5)*invCv)*invtconst
       EndIf

       If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
          tbcinfo(2,2)%bdata(1,:,:,1) = temperature(mynth,:,:)
       Else
          tbcinfo(2,2)%bdata(:,:,:,1) = (binfo(2,2)%bdata(:,:,:,1)**gamm1)*dexp(binfo(2,2)%bdata(:,:,:,5)*invCv)*invtconst
       EndIf
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(5) = timings(5) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    Call Boundary_Exchange_Init(3)

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(4) = timings(4) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    !Compute source terms and theta derivatives and corresponding rhs terms
    Call Source_Theta_Terms()

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(6) = timings(6) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. nnodes-1) Then
          Print*, 'Source'
          Print*, 'minmax vars1'
          Do iv=1,nv
             Print*, Minval(vars1(:,:,:,iv)), Maxval(vars1(:,:,:,iv))
          EndDo
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    !Exchange Phi boundary info
    Call Boundary_Exchange_Finish(3)

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(4) = timings(4) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (do_radheat) Then
       tbcinfo(3,1)%bdata(:,:,:,1) = (binfo(3,1)%bdata(:,:,:,1)**gamm1)*dexp(binfo(3,1)%bdata(:,:,:,5)*invCv)*invtconst
       tbcinfo(3,2)%bdata(:,:,:,1) = (binfo(3,2)%bdata(:,:,:,1)**gamm1)*dexp(binfo(3,2)%bdata(:,:,:,5)*invCv)*invtconst
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(5) = timings(5) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    Call Boundary_Exchange_Init(1)

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(2) = timings(2) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    !Compute phi derivatives and corresponding rhs terms
    Call Phi_Terms()

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(7) = timings(7) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. nnodes-1) Then
          Print*, 'Phi'
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    Call Boundary_Exchange_Finish(1)

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(2) = timings(2) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (do_radheat) Then
       If (myr1 .eq. 1) Then
          tbcinfo(1,1)%bdata(:,:,1,1) = temperature(:,:,1)
       Else
          tbcinfo(1,1)%bdata(:,:,:,1) = (binfo(1,1)%bdata(:,:,:,1)**gamm1)*dexp(binfo(1,1)%bdata(:,:,:,5)*invCv)*invtconst
       EndIf

       If (myr2 .eq. nr) Then
          tbcinfo(1,2)%bdata(:,:,1,1) = temperature(:,:,mynr)
       Else
          tbcinfo(1,2)%bdata(:,:,:,1) = (binfo(1,2)%bdata(:,:,:,1)**gamm1)*dexp(binfo(1,2)%bdata(:,:,:,5)*invCv)*invtconst
       EndIf
    EndIf

    If (myr1 .eq. 1) Then
       If (Bottom_Stress_Free) Then
          binfo(1,1)%bdata(:,:,1,2) = vars1(:,:,1,2)/r1/dxdr(1)
          binfo(1,1)%bdata(:,:,1,3) = vars1(:,:,1,3)/r1/dxdr(1)
       EndIf

       If (Bottom_Constant_Flux) Then
          Do ip=1,mynphi
             If (Laplacian .or. SLaplacian) Then                
                binfo(1,1)%bdata(:,ip,1,5) = -lumr1*one_over_4pi/(r1**2*ks(1)*dxdr(1))
             Else
                binfo(1,1)%bdata(:,ip,1,5) = dsdr1/dxdr(1)
             EndIf
          EndDo
       EndIf

       If (bot_inout) Then
          binfo(1,1)%bdata(:,:,1,2:5) = 0d0 !Everything but density
       EndIf

       If (magnetic) Then
          If ((Bottom_Radial_Field).and.(.not.bot_inout)) Then
             binfo(1,1)%bdata(:,:,1,6) = 0d0
             binfo(1,1)%bdata(:,:,1,7) = 0d0
             binfo(1,1)%bdata(:,:,1,8) = -2*vars1(:,:,1,8)/r1/dxdr(1)
          EndIf
       EndIf
    EndIf

    If (myr2 .eq. nr) Then
       If (Top_Stress_Free) Then
          binfo(1,2)%bdata(:,:,1,2) = vars1(:,:,mynr,2)/r2/dxdr(nr)
          binfo(1,2)%bdata(:,:,1,3) = vars1(:,:,mynr,3)/r2/dxdr(nr)
       Else
          binfo(1,2)%bdata(:,:,1,2) = vars1(:,:,mynr,2)
          binfo(1,2)%bdata(:,:,1,3) = vars1(:,:,mynr,3)
       EndIf

       If (Top_Constant_Flux) Then
          Do ip=1,mynphi
             If (Laplacian .or. SLaplacian) Then
                binfo(1,2)%bdata(:,ip,1,5) = -lumr2*one_over_4pi/(r2**2*ks(nr)*dxdr(nr))
             Else
                binfo(1,2)%bdata(:,ip,1,5) = dsdr2/dxdr(nr)
             EndIf
          EndDo
       EndIf

       If (top_inout) Then
          binfo(1,2)%bdata(:,:,1,2:5) = 0d0 !Everything but density
       EndIf

       If ((magnetic).and.(Top_Radial_Field)) Then
          If (.not. top_inout) Then
             binfo(1,2)%bdata(:,:,1,6) = 0d0
             binfo(1,2)%bdata(:,:,1,7) = 0d0
             binfo(1,2)%bdata(:,:,1,8) = -2*vars1(:,:,mynr,8)/r2/dxdr(nr)
          EndIf
       EndIf
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(5) = timings(5) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf


    Call Radial_Terms()

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(8) = timings(8) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. nnodes-1) Then
          Print*, 'Radial'
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    If ((solar_prof).and.(myr1.eq.1).and.(.not. Bottom_Specified)) Then
       Call Drag_Toward_Solar()
    EndIf

    If (Laplacian) Then
       !Do cross-derivatives
       Call Laplacian_Cross_Terms_Radial()

       If (Debug) Then
          If (myrank .eq. nnodes-1) Then
             Print*, 'Laplacian Radial'
             Print*, 'minmax dvars1'
             Do iv=1,nv
                Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
             EndDo
          EndIf
       EndIf

       Call Laplacian_Cross_Terms_Theta()

       If (Debug) Then
          If (myrank .eq. nnodes-1) Then
             Print*, 'Laplacian Theta'
             Print*, 'minmax dvars1'
             Do iv=1,nv
                Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
             EndDo
          EndIf
       EndIf
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(5) = timings(5) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    !If (((myth1 .eq. 1).or.(myth2.eq.nth)).and.(Theta_Open)) Then
    !   Call Theta_Inflow_Outflow()
    !Else
    Call Density_Terms_Theta()
    !EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(9) = timings(9) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. nnodes-1) Then
          Print*, 'Density Theta'
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    If ((bot_inout).and.(myr1 .eq. 1)) Then
       Call Bottom_Inflow_Outflow()
    ElseIf ((top_inout).and.(myr2.eq.nr)) Then
       Call Top_Inflow_Outflow()
    Else
       Call Density_Terms_Radial()
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(10) = timings(10) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. nnodes-1) Then
          Print*, 'Density Radial'
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    If (magnetic) Then
       Call Induction_Terms()

       If (call_num .gt. 1) Then
          vars1(:,:,:,6:8) = maga
       EndIf

       If (Debug) Then
          If (myrank .eq. nnodes-1) Then
             Print*, 'Induction'
             Print*, 'minmax dvars1'
             Do iv=1,nv
                Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
             EndDo
          EndIf
       EndIf
       If (Do_Timings) Then
          tstop = MPI_WTIME()
          timings(11) = timings(11) + tstop-tstart
       EndIf
    EndIf

    If (xi_cs .lt. 1d0) dvars1(:,:,:,1) = (xi_cs*xi_cs)*dvars1(:,:,:,1)

    !Wave testing, avoid coupling modes!
    If (Unstratified) Then 
       If (magnetic) Then
          dvars1(:,:,:,1:2) = 0d0          
          dvars1(:,:,:,4:5) = 0d0
          dvars1(:,:,:,7:8) = 0d0
       Else
          dvars1(:,:,:,2) = 0d0
          dvars1(:,:,:,4) = 0d0
          dvars1(:,:,:,5) = 0d0
       EndIf
    EndIf
  End Subroutine RHS

  !----------------------------------------------------------------------
  Subroutine Source_Theta_Terms()
    Implicit None
    Real*8 :: worktb1(binfo(2,1)%bsize%bnt,mynphi)
    Real*8 :: worktb2(binfo(2,2)%bsize%bnt,mynphi)
    Real*8, Dimension(mynth,mynphi) :: work1,work2,work3,work4,work5,work6
    Real*8, Dimension(mynth,mynphi) :: rho,w,v,u,t,s,br,bt,bp,oorho
    Real*8, Dimension(mynth,mynphi) :: drho,dw,dv,du,ds
    Real*8 :: ri,rv,epsr,om2,ddiff,omsq
    Integer :: ir, it, ip

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: worktb1, worktb2, work1, work2, work3, work4, work5, work6
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho,w,v,u,t,s,br,bt,bp,oorho
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: drho,dw,dv,du,ds

    !  fill sine,cosine,cotangent work arrays
    Do ip=1,mynphi
       work2(:,ip)=sines(myth1:myth2)
       work3(:,ip)=cosines(myth1:myth2)
    EndDo
    !#DIR$ VECTOR ALIGNED
    work4=work3/work2 
    If (Laplacian) work5=1d0/work2**2
    om2 = 2d0*omega_0 !1

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do ir=1,mynr
       ri = r_inv(myr1+ir-1)
       rv = r_ind(myr1+ir-1)
       epsr = eps(myr1+ir-1)
       omsq = rv*omega_0**2

       !  useful aliases, overlap computation and memory ops
       rho  = vars1(:,:,ir,1)
       !#DIR$ VECTOR ALIGNED
       dw   = omsq*work2*work3
       !#DIR$ VECTOR ALIGNED
       du   = -gravity(myr1+ir-1) + omsq*work2*work2

       w    = vars1(:,:,ir,2)
       !#DIR$ VECTOR ALIGNED
       drho = -rho*w*work4*ri
       !#DIR$ VECTOR ALIGNED
       dv   = -om2*work3*w
       !#DIR$ VECTOR ALIGNED
       du   = du + ri*w*w

       !  dw/d(theta)
       worktb1 = binfo(2,1)%bdata(:,:,ir,2)
       worktb2 = binfo(2,2)%bdata(:,:,ir,2)
       Call dbydt(w,work1,worktb1,worktb2,2)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*ri
       !#DIR$ VECTOR ALIGNED
       drho  = drho - rho*work1
       !#DIR$ VECTOR ALIGNED
       dw    = dw - w*work1 

       If (Laplacian) Then
          !#DIR$ VECTOR ALIGNED
          work6 = mur(myr1+ir-1)*invrho(:,:,ir)

          ddiff = -fourthirds*ri*ri
          !#DIR$ VECTOR ALIGNED
          dw = dw + (ddiff*work5-ri*dlnmu(myr1+ir-1))*work6*w

          vars2(:,:,ir,1) = work1
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work6
          ddiff = -onethird*(7d0*ri+2d0*dlnmu(myr1+ir-1))
          !#DIR$ VECTOR ALIGNED
          du = du + ddiff*(ri*work6*work4*w+work1)
          ddiff = fourthirds*ri
          !#DIR$ VECTOR ALIGNED
          dw = dw + ddiff*work4*work1
          Call d2bydt2(w,work1,worktb1,worktb2,2) 
          !#DIR$ VECTOR ALIGNED
          dw = dw + ddiff*ri*work6*work1
       EndIf

       v  = vars1(:,:,ir,3)
       !#DIR$ VECTOR ALIGNED
       dw = dw + v*(ri*work4*v + om2*work3)
       !#DIR$ VECTOR ALIGNED
       dv = dv - ri*v*w*work4
       !#DIR$ VECTOR ALIGNED
       du = du + v*(ri*v + om2*work2)

       !  terms using dv/d(theta)
       worktb1 = binfo(2,1)%bdata(:,:,ir,3)
       worktb2 = binfo(2,2)%bdata(:,:,ir,3)
       Call dbydt(v,work1,worktb1,worktb2,3)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*ri
       !#DIR$ VECTOR ALIGNED
       dv    = dv - w*work1

       If (Laplacian) Then
          ddiff = ri*ri
          !#DIR$ VECTOR ALIGNED
          dv    = dv + ri*work6*(work4*work1 - (ri*work5+dlnmu(myr1+ir-1))*v)
          vars2(:,:,ir,2) = work1
          Call d2bydt2(v,work1,worktb1,worktb2,3)
          !#DIR$ VECTOR ALIGNED
          dv    = dv + ddiff*work6*work1
       EndIf

       u    = vars1(:,:,ir,4)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,ir,1) = drho - 2d0*ri*rho*u
       !#DIR$ VECTOR ALIGNED
       dv   = dv - u*(ri*v + om2*work2)
       !#DIR$ VECTOR ALIGNED
       dw   = dw - ri*u*w

       If (.not. do_lorentz) Then
          dvars1(:,:,ir,3) = dv
       EndIf

       !  now terms using du/d(theta)
       worktb1 = binfo(2,1)%bdata(:,:,ir,4)
       worktb2 = binfo(2,2)%bdata(:,:,ir,4)
       Call dbydt(u,work1,worktb1,worktb2,4)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*ri
       !#DIR$ VECTOR ALIGNED
       du    = du - w*work1

       If (Laplacian) Then
          vars2(:,:,ir,3) = work1
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work6
          ddiff = (2d0*fourthirds*ri+dlnmu(myr1+ir-1))
          !#DIR$ VECTOR ALIGNED
          dw = dw + ddiff*work1
          !#DIR$ VECTOR ALIGNED
          du = du + ri*work4*work1
          Call d2bydt2(u,work1,worktb1,worktb2,4)
          ddiff = ri*ri
          !#DIR$ VECTOR ALIGNED
          du = du + ddiff*work6*work1
       EndIf

       If (.not. do_lorentz) Then
          dvars1(:,:,ir,4) = du 
       EndIf

       t = temperature(:,:,ir)
       If (abs(epsr) .gt. 0d0) Then
          !#DIR$ VECTOR ALIGNED
          ds = epsr*invrho(:,:,ir)/t
       Else 
          ds = 0d0
       EndIf

       !  ds/d(theta)
       s  = vars1(:,:,ir,5) 
       worktb1 = binfo(2,1)%bdata(:,:,ir,5)
       worktb2 = binfo(2,2)%bdata(:,:,ir,5)
       Call dbydt(s,work1,worktb1,worktb2,5)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*ri

       If (do_lorentz) Then
          !#DIR$ VECTOR ALIGNED
          dw = dw - gamm1*t*work1
       Else
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,2) = dw - gamm1*t*work1
       EndIf
       !#DIR$ VECTOR ALIGNED
       ds    = ds - w*work1

       If ((Laplacian).or.(SLaplacian)) Then
          !#DIR$ VECTOR ALIGNED
          oorho = invrho(:,:,ir)/t

          ddiff = ri*ks(myr1+ir-1)

          !#DIR$ VECTOR ALIGNED
          work1 = ddiff*work1*oorho

          !#DIR$ VECTOR ALIGNED
          ds = ds + work4*work1

          Call d2bydt2(s,work1,worktb1,worktb2,5) 
          ddiff = ddiff*ri**2

          !#DIR$ VECTOR ALIGNED
          ds = ds + ddiff*work1*oorho
       EndIf

       If (Thermal_Drag) Then
          Do ip=1,mynphi
             dvars1(:,ip,ir,5) = ds(:,ip)+(sr1-s(:,ip))*inv_drag_time(ir)
          EndDo
       Else
          dvars1(:,:,ir,5) = ds
       EndIf

       If (do_lorentz) Then
          oorho = invrho(:,:,ir)
          !#DIR$ VECTOR ALIGNED
          oorho = oorho*one_over_4pi*ri

          bt = vars1(:,:,ir,6)
          !#DIR$ VECTOR ALIGNED
          du = du - bt*bt*oorho

          bp = vars1(:,:,ir,7)
          !#DIR$ VECTOR ALIGNED
          dw = dw - bp*bp*work4*oorho
          !#DIR$ VECTOR ALIGNED
          dv = dv + bp*bt*work4*oorho
          !#DIR$ VECTOR ALIGNED
          du = du - bp*bp*oorho

          !  dbp/d(theta)
          worktb1 = binfo(2,1)%bdata(:,:,ir,7)
          worktb2 = binfo(2,2)%bdata(:,:,ir,7)
          Call dbydt(bp,work1,worktb1,worktb2,8)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          dv = dv + bt*work1
          !#DIR$ VECTOR ALIGNED
          dw = dw - bp*work1

          br = vars1(:,:,ir,8)
          !  dbr/d(theta)
          worktb1 = binfo(2,1)%bdata(:,:,ir,8)
          worktb2 = binfo(2,2)%bdata(:,:,ir,8)
          Call dbydt(br,work1,worktb1,worktb2,9)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,2) = dw + br*(bt*oorho - work1)
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,3) = dv + br*bp*oorho
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,4) = du + bt*work1
       EndIf
    EndDo
  End Subroutine Source_Theta_Terms

  !----------------------------------------------------------------------
  Subroutine Phi_Terms()
    Implicit None
    Real*8 :: workpb1(mynth,binfo(3,1)%bsize%bnp)
    Real*8 :: workpb2(mynth,binfo(3,2)%bsize%bnp)
    Real*8, Dimension(mynth,mynphi) :: work1, work2
    Real*8, Dimension(mynth,mynphi) :: rho,w,v,u,t,s,br,bt,bp,oorho
    Real*8, Dimension(mynth,mynphi) :: drho,dw,dv,du,ds
    Real*8 :: oors(mynth,mynphi), cscth(mynth,mynphi), cotth(mynth,mynphi)
    Real*8 :: ddiff, ddiff2, ri
    Integer :: ir, it, ip

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: workpb1, workpb2, work1, work2, oors, cscth, cotth
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho,w,v,u,t,s,br,bt,bp,oorho
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: drho,dw,dv,du,ds

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do ip=1,mynphi
       !#DIR$ VECTOR ALIGNED
       cscth(:,ip) = 1d0/sines(myth1:myth2)
    EndDo

    If (Laplacian) Then
       !#DIR$ VECTOR ALIGNED
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          cotth(:,ip) = cosines(myth1:myth2)*cscth(:,ip)
       EndDo
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do ir=1,mynr
       ri = r_inv(myr1+ir-1)
       !  useful aliases, overlap computation and memory access
       rho  = vars1(:,:,ir,1)
       !#DIR$ VECTOR ALIGNED
       oors = ri*cscth(:,:)

       v    = vars1(:,:,ir,3)

       !drho/dphi
       workpb1 = binfo(3,1)%bdata(:,:,ir,1)
       workpb2 = binfo(3,2)%bdata(:,:,ir,1)
       Call dbydp(rho,work1,workpb1,workpb2,1)

       drho  = dvars1(:,:,ir,1)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*oors

       dv    = dvars1(:,:,ir,3)
       !#DIR$ VECTOR ALIGNED
       drho = drho - v*work1

       oorho = invrho(:,:,ir)
       t  = temperature(:,:,ir)
       !#DIR$ VECTOR ALIGNED
       dv = dv - pconst*work1*t*oorho 

       ds = dvars1(:,:,ir,5)
       If (do_radheat) Then
          w = work1*oorho
          workpb1 = tbcinfo(3,1)%bdata(:,:,ir,1)
          workpb2 = tbcinfo(3,2)%bdata(:,:,ir,1)
          Call dbydp(t,work1,workpb1,workpb2,6)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oors/t
          ds = ds + dvars2(:,:,ir,1)*(3d0*work1-w)*work1

          Call d2bydp2(t,work1,workpb1,workpb2,6)
          !#DIR$ VECTOR ALIGNED
          work1=work1*oors**2

          !#DIR$ VECTOR ALIGNED
          ds = ds + dvars2(:,:,ir,1)*work1/t
       EndIf

       !  dv/d(phi):
       workpb1 = binfo(3,1)%bdata(:,:,ir,3)
       workpb2 = binfo(3,2)%bdata(:,:,ir,3)
       Call dbydp(v,work1,workpb1,workpb2,3)

       w = vars1(:,:,ir,2)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*oors
       !#DIR$ VECTOR ALIGNED
       dv = dv - v*work1
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,ir,1) = drho - rho*work1

       dw = dvars1(:,:,ir,2)
       If (Laplacian) Then
          u = vars1(:,:,ir,4)
          du = dvars1(:,:,ir,4)
          !#DIR$ VECTOR ALIGNED
          work2 = oorho*mur(myr1+ir-1)
          vars2(:,:,ir,4) = work1
          !Viscous heating theta-theta, phi-phi terms
          !#DIR$ VECTOR ALIGNED
          ds = ds + fourthirds*work2*(vars2(:,:,ir,1)+u*ri)*(vars2(:,:,ir,1)+2d0*u*ri-work1-ri*cotth*w)/t
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work2
          ddiff = -seventhirds*ri
          ddiff2 = ddiff-twothirds*dlnmu(myr1+ir-1)
          !#DIR$ VECTOR ALIGNED
          dw = dw + ddiff*cotth*work1
          !#DIR$ VECTOR ALIGNED
          du = du + ddiff2*work1
          Call d2bydp2(v,work1,workpb1,workpb2,3)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oors**2
          !#DIR$ VECTOR ALIGNED
          dv = dv + fourthirds*work1*work2
       EndIf

       !  dw/d(phi):
       workpb1 = binfo(3,1)%bdata(:,:,ir,2)
       workpb2 = binfo(3,2)%bdata(:,:,ir,2)
       Call dbydp(w,work1,workpb1,workpb2,2)
       !#DIR$ VECTOR ALIGNED
       work1=work1*oors
       !#DIR$ VECTOR ALIGNED
       dw = dw - v*work1

       If (Laplacian) Then
          !Viscous heating theta-phi term
          !#DIR$ VECTOR ALIGNED
          ds = ds + work2*(vars2(:,:,ir,2)+work1-ri*cotth(:,:)*v)**2/t
          vars2(:,:,ir,2) = work1
          ddiff = seventhirds*ri
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work2
          !#DIR$ VECTOR ALIGNED
          dv = dv + ddiff*cotth(:,:)*work1
          Call d2bydp2(w,work1,workpb1,workpb2,2)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oors**2 
          !#DIR$ VECTOR ALIGNED
          dw = dw + work1*work2
       Else
          u = vars1(:,:,ir,4)
          du = dvars1(:,:,ir,4)
       EndIf

       If (.not. do_lorentz) Then
          dvars1(:,:,ir,2) = dw
       EndIf

       !  du/d(phi)
       workpb1 = binfo(3,1)%bdata(:,:,ir,4)
       workpb2 = binfo(3,2)%bdata(:,:,ir,4)
       Call dbydp(u,work1,workpb1,workpb2,4)
       !#DIR$ VECTOR ALIGNED
       work1=work1*oors
       !#DIR$ VECTOR ALIGNED
       du = du-v*work1 

       If (Laplacian) Then
          vars2(:,:,ir,5) = work1
          ddiff = (2d0*fourthirds*ri+dlnmu(myr1+ir-1))
          !#DIR$ VECTOR ALIGNED
          dv = dv + ddiff*work1*work2
          Call d2bydp2(u,work1,workpb1,workpb2,4)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oors**2 
          !#DIR$ VECTOR ALIGNED
          du = du + work1*work2
       EndIf

       If (.not. do_lorentz) Then
          dvars1(:,:,ir,4) = du
       EndIf

       !  ds/d(phi):
       s = vars1(:,:,ir,5)
       workpb1 = binfo(3,1)%bdata(:,:,ir,5)
       workpb2 = binfo(3,2)%bdata(:,:,ir,5)
       Call dbydp(s,work1,workpb1,workpb2,5)
       !#DIR$ VECTOR ALIGNED
       work1 = work1*oors
       !#DIR$ VECTOR ALIGNED
       ds    = ds - v*work1
       !#DIR$ VECTOR ALIGNED
       dv    = dv - gamm1*t*work1

       If (.not. do_lorentz) Then
          dvars1(:,:,ir,3) = dv
       EndIf

       If ((Laplacian).or.(SLaplacian)) Then
          Call d2bydp2(s,work1,workpb1,workpb2,5) !11*mynphi*mynphi+25*mynphi
          ddiff = ks(myr1+ir-1)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*(oors**2)
          !#DIR$ VECTOR ALIGNED
          ds = ds + ddiff*work1*(oorho/t)
       EndIf

       dvars1(:,:,ir,5) = ds

       If ((magnetic).and.(do_lorentz)) Then
          oorho = one_over_4pi*oorho*oors

          bt  = vars1(:,:,ir,6)
          !  dbt/d(phi)
          workpb1 = binfo(3,1)%bdata(:,:,ir,6)
          workpb2 = binfo(3,2)%bdata(:,:,ir,6)
          Call dbydp(bt,work1,workpb1,workpb2,7)
          bp = vars1(:,:,ir,7)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          dv = dv - bt*work1
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,2) = dw + bp*work1

          br  = vars1(:,:,ir,8)
          !  dbr/d(phi)
          workpb1 = binfo(3,1)%bdata(:,:,ir,8)
          workpb2 = binfo(3,2)%bdata(:,:,ir,8)
          Call dbydp(br,work1,workpb1,workpb2,9)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,3) = dv - br*work1
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,4) = du + bp*work1
       EndIf
    EndDo
  End Subroutine Phi_Terms

  !----------------------------------------------------------------------
  Subroutine Radial_Terms()
    Implicit None
    Real*8 :: workrb1(mynth,binfo(1,1)%bsize%bnr)
    Real*8 :: workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8, Dimension(mynth,mynr) :: work1,work2,work3,work4,work5,work6,work7,work8,work9,work10
    Real*8, Dimension(mynth,mynr) :: w,v,u,t,s,br,bt,bp,oorho
    Real*8, Dimension(mynth,mynr) :: dw,dv,du,ds
    Real*8 :: cotth(mynth)
    Integer :: ir, it, ip

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: workrb1, workrb2, work1, work2, work3, work4, work5, work6, work7, work8, work9, work10, cotth
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: w,v,u,t,s,br,bt,bp,oorho
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dw,dv,du,ds

    If (non_uniform) Then
       Do ip=1,mynth
          work3(ip,:)=dxdr(myr1:myr2)
          work9(ip,:)=d2xdr2(myr1:myr2)
       EndDo
       work9 = work9/work3
    EndIf

    If ((Laplacian).or.(SLaplacian)) Then
       !#DIR$ VECTOR ALIGNED
       cotth = cosines(myth1:myth2)/sines(myth1:myth2)
       !#DIR$ VECTOR ALIGNED
       Do ip=1,mynth
          work4(ip,:) = ks(myr1:myr2)
          !#DIR$ VECTOR ALIGNED
          work5(ip,:) = 2d0*r_inv(myr1:myr2)+dlnmu(myr1:myr2)
          !#DIR$ VECTOR ALIGNED
          work6(ip,:) = (2d0*r_inv(myr1:myr2)+dlnks(myr1:myr2))*ks(myr1:myr2)
          work10(ip,:) = r_inv(myr1:myr2)
       EndDo
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    !#DIR$ VECTOR ALIGNED
    Do ip=1,mynphi
       !  dw/dr:
       w   = vars1(:,ip,:,2)
       workrb1 = binfo(1,1)%bdata(:,ip,:,2)
       workrb2 = binfo(1,2)%bdata(:,ip,:,2)
       Call dbydr(w,work1,workrb1,workrb2,2)
       u = vars1(:,ip,:,4)
       If (non_uniform) Then
          work1 = work1*work3
       EndIf
       dw = dvars1(:,ip,:,2)
       !#DIR$ VECTOR ALIGNED
       dw = dw - u*work1

       t = temperature(:,ip,:)
       If (Laplacian) Then 
          du = dvars1(:,ip,:,4)
          ds = dvars1(:,ip,:,5)
          !#DIR$ VECTOR ALIGNED              
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             work2(:,ir) = mur(myr1+ir-1)*invrho(:,ip,ir)
             !#DIR$ VECTOR ALIGNED
             du(:,ir) = du(:,ir) + onethird*cotth*work10(:,ir)*work1(:,ir)
             !#DIR$ VECTOR ALIGNED
             dw(:,ir) = dw(:,ir) + work2(:,ir)*work5(:,ir)*work1(:,ir)
             !Viscous heating r-theta
             !#DIR$ VECTOR ALIGNED
             ds(:,ir) = ds(:,ir) + work2(:,ir)*(work1(:,ir)+vars2(:,ip,ir,3)-w(:,ir)*work10(:,ir))**2/t(:,ir)
          EndDo

          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             dw = dw + work9*work1*work2
          EndIf

          Call d2bydr2(w,work1,workrb1,workrb2,2)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work3**2
          EndIf
          !#DIR$ VECTOR ALIGNED
          dw = dw + work1*work2
       EndIf

       If (.not. do_lorentz) Then
          dvars1(:,ip,:,2) = dw
       EndIf

       !  dv/dr:
       v = vars1(:,ip,:,3)
       workrb1 = binfo(1,1)%bdata(:,ip,:,3)
       workrb2 = binfo(1,2)%bdata(:,ip,:,3)
       Call dbydr(v,work1,workrb1,workrb2,3)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work3
       EndIf
       dv = dvars1(:,ip,:,3)
       !#DIR$ VECTOR ALIGNED
       dv = dv - u*work1

       If (Laplacian) Then
          !Viscous heating r-theta
          !#DIR$ VECTOR ALIGNED
          ds = ds + work2*(work1+vars2(:,ip,:,5)-v*work10)**2/t
          !#DIR$ VECTOR ALIGNED
          dv = dv + work1*work5*work2
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             dv = dv + work9*work1*work2
          EndIf

          Call d2bydr2(v,work1,workrb1,workrb2,3)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work3**2
          EndIf
          !#DIR$ VECTOR ALIGNED
          dv = dv + work1*work2
       EndIf

       If (.not. do_lorentz) Then
          dvars1(:,ip,:,3) = dv
       EndIf

       If (.not. Laplacian) du = dvars1(:,ip,:,4)
       !  now terms involving du/dr:
       workrb1 = binfo(1,1)%bdata(:,ip,:,4)
       workrb2 = binfo(1,2)%bdata(:,ip,:,4)
       Call dbydr(u,work1,workrb1,workrb2,4)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work3
       EndIf
       !#DIR$ VECTOR ALIGNED
       du = du-u*work1
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ IVDEP
          !#DIR$ VECTOR ALIGNED
          Do it=1,mynth
             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1)-vars1(it,ip,ir,1)*work1(it,ir)
          EndDo
       EndDo

       If (Laplacian) Then
          !#DIR$ VECTOR ALIGNED
          du = du + fourthirds*work2*work5*(work1-u*work10)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             du = du + fourthirds*work9*work1*work2
          EndIf
          !Viscous heating r-r
          !#DIR$ VECTOR ALIGNED
          ds = ds + fourthirds*work2*work1*(work1-vars2(:,ip,:,1)-u*work10)/t
          !Viscous heating r-theta
          !#DIR$ VECTOR ALIGNED
          ds = ds + work2*(work1+vars2(:,ip,:,5)-v*work10)**2/t
          !Viscous heating r-phi
          !#DIR$ VECTOR ALIGNED
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             ds(:,ir) = ds(:,ir) + fourthirds*work2(:,ir)*(vars2(:,ip,ir,4)+(cotth*w(:,ir)+u(:,ir))*work10(:,ir))* &
                  & (vars2(:,ip,ir,4)+(cotth*w(:,ir)+u(:,ir))*work10(:,ir)-work1(:,ir))/t(:,ir)
          EndDo
          Call d2bydr2(u,work1,workrb1,workrb2,4)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work3**2
          EndIf
          !#DIR$ VECTOR ALIGNED
          du = du + fourthirds*work1*work2
       Else
          ds = dvars1(:,ip,:,5)
       EndIf

       !  ds/dr:
       s  = vars1(:,ip,:,5)
       workrb1 = binfo(1,1)%bdata(:,ip,:,5)
       workrb2 = binfo(1,2)%bdata(:,ip,:,5)
       Call dbydr(s,work1,workrb1,workrb2,5)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work1 = work1*work3
       EndIf
       !#DIR$ VECTOR ALIGNED
       du = du - gamm1*t*work1
       !#DIR$ VECTOR ALIGNED
       ds = ds - u*work1

       If (do_radheat) dvars2(:,ip,:,2) = work1

       If ((Laplacian).or.(SLaplacian)) Then
          !#DIR$ VECTOR ALIGNED
          oorho = invrho(:,ip,:)/t
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          ds = ds + work1*work6
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             ds = ds + work9*work1*work4
          EndIf
          Call d2bydr2(s,work1,workrb1,workrb2,5)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work3**2
          EndIf
          !#DIR$ VECTOR ALIGNED
          ds = ds + work4*work1*oorho
       EndIf

       If (.not. do_lorentz) Then
          dvars1(:,ip,:,4) = du
       EndIf
       dvars1(:,ip,:,5) = ds

       If ((magnetic).and.(do_lorentz)) Then
          !#DIR$ VECTOR ALIGNED
          oorho = invrho(:,ip,:)*one_over_4pi

          bt = vars1(:,ip,:,6)
          !  dbt/dr:
          workrb1 = binfo(1,1)%bdata(:,ip,:,6)
          workrb2 = binfo(1,2)%bdata(:,ip,:,6)
          Call dbydr(bt,work1,workrb1,workrb2,7)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work3
          EndIf
          br = vars1(:,ip,:,8)
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          du = du - bt*work1
          !#DIR$ VECTOR ALIGNED
          dvars1(:,ip,:,2) = dw + br*work1

          !  dbp/dr:
          bp = vars1(:,ip,:,7)
          workrb1 = binfo(1,1)%bdata(:,ip,:,7)
          workrb2 = binfo(1,2)%bdata(:,ip,:,7)
          Call dbydr(bp,work1,workrb1,workrb2,8)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work3
          EndIf
          !#DIR$ VECTOR ALIGNED
          work1 = work1*oorho
          !#DIR$ VECTOR ALIGNED
          dvars1(:,ip,:,3) = dv + br*work1
          !#DIR$ VECTOR ALIGNED
          dvars1(:,ip,:,4) = du - bp*work1
       EndIf
    EndDo
  End Subroutine Radial_Terms

  !----------------------------------------------------------------------
  Subroutine Density_Terms_Theta()
    Implicit None
    !  Now that everything else is defined we know the density
    !  gradient normal to the boundary. We can now evaluate the
    !  terms involving grad(rho)
    Real*8 :: worktb1(binfo(2,1)%bsize%bnt,mynphi)
    Real*8 :: worktb2(binfo(2,2)%bsize%bnt,mynphi)
    Real*8, Dimension(mynth,mynphi) :: work1, work2
    Real*8, Dimension(mynth,mynphi) :: rho, dw, t, ds
    Real*8 :: ri, rv
    Integer :: ir, it, ip

    If (do_radheat) Then
       Do ip=1,mynphi
          work2(:,ip) = cosines(myth1:myth2)/sines(myth1:myth2)
       EndDo
    EndIf

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: worktb1, worktb2, work1
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho,dw,t

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do ir=1,mynr
       rv  = r_ind(myr1+ir-1)*invpconst
       ri  = r_inv(myr1+ir-1)

       !  d(rho)/d(theta): set theta BC for d(rho)/d(theta) so that dw=0
       t = temperature(:,:,ir)
       rho = vars1(:,:,ir,1)
       dw  = dvars1(:,:,ir,2)
       If (myth1 .eq. 1) Then
          !#DIR$ VECTOR ALIGNED
          worktb1(1,:) = rv*dw(1,:)*rho(1,:)/t(1,:)
          worktb2 = binfo(2,2)%bdata(:,:,ir,1)
       ElseIf ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
          !#DIR$ VECTOR ALIGNED
          worktb2(1,:) = rv*dw(mynth,:)*rho(mynth,:)/t(mynth,:)
          worktb1 = binfo(2,1)%bdata(:,:,ir,1)
       Else
          worktb1 = binfo(2,1)%bdata(:,:,ir,1)
          worktb2 = binfo(2,2)%bdata(:,:,ir,1)
       EndIf

       Call dbydt(rho,work1,worktb1,worktb2,1)
       work1 = ri*work1
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1) - vars1(it,ip,ir,2)*work1(it,ip)
          EndDo
       EndDo


       work1 = work1*invrho(:,:,ir)
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,2) = dw(it,ip) - pconst*work1(it,ip)*t(it,ip)
          EndDo
       EndDo

       If (do_radheat) Then
          rho = work1
          ds = dvars1(:,:,ir,5)
          worktb1 = tbcinfo(2,1)%bdata(:,:,ir,1)
          worktb2 = tbcinfo(2,2)%bdata(:,:,ir,1)
          Call dbydt(t,work1,worktb1,worktb2,6) 
          work1 = ri*work1/t
          !#DIR$ VECTOR ALIGNED
          ds = ds + dvars2(:,:,ir,1)*work1*(3d0*work1 + ri*work2 - rho)

          Call d2bydt2(t,work1,worktb1,worktb2,6)
          !#DIR$ VECTOR ALIGNED
          dvars1(:,:,ir,5) = ds + (ri*ri)*dvars2(:,:,ir,1)*work1/t
       EndIf
    EndDo
  End Subroutine Density_Terms_Theta

!!$  !----------------------------------------------------------------------
!!$  Subroutine Theta_Inflow_Outflow()
!!$    Implicit None
!!$    Real*8 :: worktb1(binfo(2,1)%bsize%bnt,mynphi)
!!$    Real*8 :: worktb2(binfo(2,2)%bsize%bnt,mynphi)
!!$    Real*8, Dimension(mynth,mynphi) :: work1
!!$    Real*8, Dimension(mynth,mynphi) :: rho, dw, t
!!$    Real*8 :: ri, rv, meanu, meanw, meanv, means, umax, tau
!!$
!!$    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: worktb1, worktb2, work1
!!$    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho,dw,t
!!$
!!$    If (myth1 .eq. 1) Then
!!$       it = 1
!!$    ElseIf (myth2 .eq. nth) Then
!!$       it = nth
!!$    EndIf
!!$
!!$    ntot  = 0d0
!!$    meanv = 0d0
!!$    meanw = 0d0
!!$    Ltot  = 0d0
!!$    mrsint = 0d0
!!$    Do ir=1,mynr
!!$       rv = r_ind(myr1+ir-1)
!!$       ri = dr_ind(myr1+ir-1)
!!$       Do ip=1,mynphi
!!$          umax = vars1(it,ip,ir,2)
!!$          enter = min(0d0,umax)/min(-1d-6,umax)
!!$          ntot = ntot + enter
!!$          meanu = meanu + enter*vars1(it,ip,ir,2)
!!$          enter = enter*vars1(it,ip,ir,1)*ri*rv**3
!!$          Ltot = Ltot + enter*vars1(it,ip,ir,3)
!!$          mrsint = mrsint + enter
!!$       EndDo
!!$    EndDo
!!$    snd(1) = ntot
!!$    snd(2) = meanu
!!$    snd(3) = Ltot
!!$    snd(4) = mrsint
!!$    snd(5) = Sum(vars1(it,:,:,2))
!!$    snd(6) = Sum(vars1(it,:,:,5))
!!$    If (myth1.eq.1) Then
!!$       snd(7) = Abs(Minval(vars1(2,:,:,2)))
!!$       rcv = 0d0
!!$       Call Global_AllReduce(snd,rcv,theta_bot_comm,'csm')
!!$       ntot = Dble(nth*nphi)
!!$       meanu = -rcv(2)/(ntot-rcv(1))
!!$       meanv = -rcv(3)/rcv(4)
!!$       meanw = rcv(5)/ntot/deltat
!!$       means = sr(1)-rcv(6)/ntot
!!$       tau = rcv(7)/dr_ind(1)
!!$       meanv = tau*meanv
!!$       meanw = tau*meanw
!!$    ElseIf (myth2 .eq. nth) Then
!!$       snd(7) = Abs(Maxval(vars1(nth-1,:,:,2)))
!!$       rcv = 0d0
!!$       Call Global_AllReduce(snd,rcv,theta_top_comm,'csm')
!!$       ntot = Dble(nth*nphi)
!!$       meanw = -rcv(2)/(ntot-rcv(1))
!!$       meanv = -rcv(3)/rcv(4)
!!$       meanu = rcv(5)/ntot/deltat
!!$       means = sr(1)-rcv(6)/ntot
!!$       tau = rcv(7)/dr_ind(1)
!!$       meanv = tau*meanv
!!$       meanw = tau*meanw
!!$    EndIf
!!$
!!$    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
!!$    !#DIR$ VECTOR ALIGNED
!!$    Do ir=1,mynr
!!$       rv  = r_ind(myr1+ir-1)*invpconst
!!$       ri  = r_inv(myr1+ir-1)
!!$
!!$       !  d(rho)/d(theta): set theta BC for d(rho)/d(theta) so that dw=0
!!$       t = temperature(:,:,ir)
!!$       rho = vars1(:,:,ir,1)
!!$       dw  = dvars1(:,:,ir,2)
!!$       If (myth1 .eq. 1) Then
!!$          !#DIR$ VECTOR ALIGNED
!!$          worktb1(1,:) = rv*dw(1,:)*rho(1,:)/t(1,:)
!!$          worktb2 = binfo(2,2)%bdata(:,:,ir,1)
!!$       ElseIf ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
!!$          !#DIR$ VECTOR ALIGNED
!!$          worktb2(1,:) = rv*dw(mynth,:)*rho(mynth,:)/t(mynth,:)
!!$          worktb1 = binfo(2,1)%bdata(:,:,ir,1)
!!$       EndIf
!!$
!!$       Call dbydt(rho,work1,worktb1,worktb2,1)
!!$
!!$       !#DIR$ VECTOR ALIGNED
!!$       !#DIR$ IVDEP
!!$       Do ip=1,mynphi
!!$          !#DIR$ VECTOR ALIGNED
!!$          !#DIR$ IVDEP
!!$          Do it=1,mynth
!!$             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1) - ri*vars1(it,ip,ir,2)*work1(it,ip)
!!$             dvars1(it,ip,ir,2) = dw(it,ip) - ri*pconst*work1(it,ip)*temperature(it,ip,ir)*invrho(it,ip,ir)
!!$          EndDo
!!$       EndDo
!!$    EndDo
!!$  End Subroutine Theta_Inflow_Outflow

  Subroutine Density_Terms_Radial()
    Implicit None
    !  Now that everything else is defined we know the density
    !  gradient normal to the boundary. We can now evaluate the
    !  terms involving grad(rho)
    Real*8 :: tmpv
    Real*8 :: workrb1(mynth,binfo(1,1)%bsize%bnr)
    Real*8 :: workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8, Dimension(mynth,mynr) :: work1, work2, work3, work4
    Real*8, Dimension(mynth,mynr) :: rho,u,t,du,ds
    Integer :: ir, it, ip

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: workrb1, workrb2, work1, work2, work3
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho,u,t,du

    If (non_uniform) Then
       Do it=1,mynth
          work2(it,:)=dxdr(myr1:myr2)
          work3(it,:)=d2xdr2(myr1:myr2)
       EndDo
       work3 = work3/work2
    EndIf

    If (do_radheat) Then
       !#DIR$ VECTOR ALIGNED
       Do it=1,mynth
          !#DIR$ VECTOR ALIGNED
          work4(it,:) = 2d0*r_inv(myr1:myr2)-dlnkr(myr1:myr2)
       EndDo
    EndIf

    tmpv = (m0-mtot)/(vol*damping_time)
    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    Do ip=1,mynphi

       t   = temperature(:,ip,:)
       du  = dvars1(:,ip,:,4)

       rho = vars1(:,ip,:,1)
       u   = vars1(:,ip,:,4)
       If (myr1 .eq. 1) Then
          If ((lbcr_rho.Eq.2).Or.(lbcr_rho.Eq.4)) Then
             !  normal (radial) component
             !#DIR$ VECTOR ALIGNED
             workrb1(:,1) = du(:,1)*rho(:,1)*invpconst/t(:,1)
          Else
             workrb1(:,1) = rho(:,1)
          EndIf

          If (non_uniform) workrb1(:,1) = workrb1(:,1)/work2(:,1)
          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
       ElseIf (myr2 .eq. nr) Then
          If ((lbcr_rho.Eq.3).Or.(lbcr_rho.Eq.4)) Then
             !#DIR$ VECTOR ALIGNED
             workrb2(:,1) = du(:,mynr)*rho(:,mynr)*invpconst/t(:,mynr)
          Else
             workrb2(:,1) = rho(:,mynr)
          EndIf
          If (non_uniform) workrb2(:,1) = workrb2(:,1)/work2(:,mynr)
          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
       Else
          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
       EndIf

       Call dbydr(rho,work1,workrb1,workrb2,1)
       If (non_uniform) work1=work1*work2

       !If (myr1 .eq. 1) Then
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
       !   Do it=1,mynth
       !      dvars1(it,ip,1,1) = dvars1(it,ip,1,1) + tmpv
       !   EndDo
       !EndIf

       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth            
             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1) - u(it,ir)*work1(it,ir)
          EndDo
       EndDo

       work1 = work1*invrho(:,ip,:)
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,4) = du(it,ir) - pconst*t(it,ir)*work1(it,ir)
          EndDo
       EndDo

       !Drag upper boundary to desired entropy
       If ((myr2 .eq. nr).and.(.not. top_inout).and.(.not. Set_Stop).and.(.not. ell0_heating).and.(.not. do_radheat)) Then
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          ir = mynr
          Do it=1,mynth
             !dvars1(it,ip,mynr,5) = dvars1(it,ip,mynr,5) + (sr(1)+delta_s-vars1(it,ip,mynr,5))/damping_time
             !dvars1(it,ip,mynr,5) = dvars1(it,ip,mynr,5) + (sr(nr)+delta_s-vars1(it,ip,mynr,5))/damping_time
             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1) + rho(it,ir)*(T0-t(it,ir))/(damping_time*(gamma-2d0)*t(it,ir))
             dvars1(it,ip,ir,5) = dvars1(it,ip,ir,5) - Cv*(T0-t(it,ir))/(damping_time*(gamma-2d0)*t(it,ir))
          EndDo
       EndIf

       If (do_radheat) Then
          rho = work1
          ds = dvars1(:,ip,:,5)
          If (myr1 .eq. 1) Then
             workrb1(:,1) = t(:,1)*(gamm1*rho(:,1)+invCv*dvars2(:,ip,1,2))/dxdr(1)
          Else
             workrb1 = tbcinfo(1,1)%bdata(:,ip,:,1)
          EndIf
          If (myr2 .eq. nr) Then
             workrb2(:,1) = t(:,mynr)*(gamm1*rho(:,mynr)+invCv*dvars2(:,ip,mynr,2))/dxdr(nr)
          Else
             workrb2 = tbcinfo(1,2)%bdata(:,ip,:,1)
          EndIf
          Call dbydr(t,work1,workrb1,workrb2,6)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work2
          EndIf
          work1 = work1/t
          !#DIR$ VECTOR ALIGNED
          ds = ds + dvars2(:,ip,:,1)*work1*(work4+3d0*work1-rho)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             ds = ds + dvars2(:,ip,:,1)*work1*work3
          EndIf

          Call d2bydr2(t,work1,workrb1,workrb2,6)
          If (non_uniform) Then
             !#DIR$ VECTOR ALIGNED
             work1 = work1*work2**2
          EndIf
          !#DIR$ VECTOR ALIGNED
          dvars1(:,ip,:,5) = ds + dvars2(:,ip,:,1)*work1/t
       EndIf
    EndDo
  End Subroutine Density_Terms_Radial

  Subroutine Bottom_Inflow_Outflow()
    Implicit None
    Real*8 :: enter, leave, masstot, ntot, meanu, means, umax, tau, tconst, snd(4), rcv(4), bc1(mynth,1), bc2(mynth,4)
    Real*8, Dimension(mynth,mynr) :: rho, ur, temp, dwork, du
    Integer :: it, ip, ir

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: bc1, bc2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho,ur,temp,dwork,du

    If (.not. Bottom_Specified) Then
       ntot  = 0d0
       Do ip=1,mynphi
          Do it=1,mynth
             tau = vars1(it,ip,mynr,1)*sines(myth1+it-1)
             meanu = meanu + tau*vars1(it,ip,1,4)
             masstot = masstot + tau
          EndDo
       EndDo
       snd(1) = meanu
       snd(2) = masstot
       snd(3) = Sum(vars1(:,:,1,5))
       snd(4) = Abs(Minval(vars1(:,:,2,4)))
       rcv = 0d0
       Call Global_AllReduce(snd,rcv,bot_comm,'csm')
       ntot = Dble(nth*nphi)
       meanu = -rcv(1)/rcv(2)/deltat
       means = sr(1)-rcv(3)/ntot
       tau = rcv(4)/dr_ind(1)
       ntot = rhor(2)/rhor(1)
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    !#DIR$ VECTOR ALIGNED
    !#DIR$ IVDEP
    Do ip=1,mynphi
       temp = temperature(:,ip,:)
       du  = dvars1(:,ip,:,4)
       rho = vars1(:,ip,:,1)
       ur  = vars1(:,ip,:,4)
       !#DIR$ VECTOR ALIGNED
       bc1(:,1) = invpconst*rho(:,1)*du(:,1)/temp(:,1) !Sets rho, and du = 0 on r2
       bc2 = binfo(1,2)%bdata(:,ip,:,1)
       Call dbydr(rho,dwork,bc1,bc2,-4)

       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1) - ur(it,ir)*dwork(it,ir)
             dvars1(it,ip,ir,4) = du(it,ir) - pconst*temp(it,ir)*dwork(it,ir)*invrho(it,ip,ir)
          EndDo
       EndDo
    EndDo

    If (Bottom_Specified) Then
       tau = sqrt(gamm1*Cp*Tr(nr))/(20d0*dr_ind(1))
       tconst = Min(abs(sim_time)/1d5,1d0)
       meanu = tau*(m0-mtot)/vol
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             !enter = vars1(it,ip,2,4)
             !enter = min(0d0,enter)/min(-1d0,enter)
             !leave = tau*enter
             !enter = tau*(1d0-enter)
             dvars1(it,ip,1,1) = dvars1(it,ip,1,1) + meanu
             !dvars1(it,ip,1,2) = dvars1(it,ip,1,2) + (enter*tconst*ashvars(it,ip,2) + leave*vars1(it,ip,2,2) - tau*vars1(it,ip,1,2)) 
             !dvars1(it,ip,1,3) = dvars1(it,ip,1,3) + (enter*tconst*ashvars(it,ip,3) + leave*vars1(it,ip,2,3) - tau*vars1(it,ip,1,3)) 
             !dvars1(it,ip,1,4) = dvars1(it,ip,1,4) + (enter*tconst*ashvars(it,ip,1) + leave*vars1(it,ip,2,4) - tau*vars1(it,ip,1,4)) 
             !dvars1(it,ip,1,5) = dvars1(it,ip,1,5) + (enter*(sr1(it)+tconst*ash_fen_mod*ashvars(it,ip,4)) + leave*vars1(it,ip,2,5) - tau*vars1(it,ip,1,5))

             dvars1(it,ip,1,2) = dvars1(it,ip,1,2) + tau*(tconst*ashvars(it,ip,2) - vars1(it,ip,1,2)) 
             dvars1(it,ip,1,3) = dvars1(it,ip,1,3) + tau*(tconst*ashvars(it,ip,3) - vars1(it,ip,1,3)) 
             dvars1(it,ip,1,4) = dvars1(it,ip,1,4) + tau*(tconst*ashvars(it,ip,1) - vars1(it,ip,1,4)) 
             dvars1(it,ip,1,5) = dvars1(it,ip,1,5) + tau*((sr(1)+tconst*ash_fen_mod*ashvars(it,ip,4)) - vars1(it,ip,1,5))
          EndDo
       EndDo
    Else
       tconst = tau*(m0-mtot)/vol
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,1,1) = dvars1(it,ip,1,1) + tconst
             dvars1(it,ip,1,4) = dvars1(it,ip,1,4) + tau*(ntot*vars1(it,ip,2,4)-vars1(it,ip,1,4)) + meanu
             dvars1(it,ip,1,5) = dvars1(it,ip,1,5) + tau*(vars1(it,ip,2,5)-vars1(it,ip,1,5) + means)
          EndDo
       EndDo
    EndIf
  End Subroutine Bottom_Inflow_Outflow

  Subroutine Top_Inflow_Outflow()
    Implicit None
    Real*8 :: enter, leave, ntot, masstot, massin, meanw, meanv,  meanu, means, mdvel, umax, tau, taur(mynr)
    Real*8 :: snd(7), rcv(7), bc1(mynth,4), bc2(mynth,1)
    Real*8, Dimension(mynth,mynr) :: rho, ur, temp, dwork, du
    Integer :: ir, it, ip

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: bc1, bc2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rho, ur, temp, dwork, du

    meanw = 0d0
    meanv = 0d0
    meanu = 0d0
    massin = 0d0
    masstot = 0d0
    mdvel = 0d0
    If (do_entropy_rain) Then
       Do ip=1,mynphi
          Do it=1,mynth
             umax = vars1(it,ip,mynr,4)
             enter = 1d0-min(0d0,umax)/min(-1d-6,umax)
             leave = 1d0-enter
             tau = vars1(it,ip,mynr,1)*sines(myth1+it-1)
             meanw = meanw + enter*vars1(it,ip,mynr,2)*tau*sines(myth1+it-1)
             meanv = meanv + enter*vars1(it,ip,mynr,3)*tau*sines(myth1+it-1)
             meanu = meanu + umax*tau
             massin = massin + leave*tau*sines(myth1+it-1)
             masstot = masstot + tau
             mdvel = mdvel + tau*(dvel_rain(it,ip)+vars1(it,ip,mynr-1,4))
          EndDo
       EndDo
       taur = sqrt(gamm1*Cp*Tr(nr))*((r_ind(myr1:myr2)-r_ind(myr1))/(r2-r_ind(myr1)))**4/dr_ind(nr)/2d0
    Else
       Do ip=1,mynphi
          Do it=1,mynth
             umax = vars1(it,ip,mynr,4)
             enter = 1d0-min(0d0,umax)/min(-1d-6,umax)
             tau = vars1(it,ip,mynr,1)*sines(myth1+it-1)
             meanw = meanw + enter*vars1(it,ip,mynr,2)*tau
             meanv = meanv + enter*vars1(it,ip,mynr,3)*tau
             meanu = meanu + umax*tau
             massin = massin + (1d0-enter)*tau
             masstot = masstot + tau
          EndDo
       EndDo
    EndIf
    snd(1) = meanw
    snd(2) = meanv
    snd(3) = meanu
    snd(4) = massin
    snd(5) = masstot    
    snd(6) = Sum(vars1(:,:,mynr,5))
    snd(7) = mdvel
    rcv = 0d0
    Call Global_AllReduce(snd,rcv,top_comm,'sum')
    ntot = Dble(nth*nphi)
    meanw = -rcv(1)/max(rhor(nr),rcv(4))
    meanv = -rcv(2)/max(rhor(nr),rcv(4))
    meanu = -rcv(3)/rcv(5)/deltat
    means = sr(1)+delta_s-rcv(6)/ntot
    mdvel = -rcv(7)/rcv(5)
    tau = sqrt(gamm1*Cp*Tr(nr))/dr_ind(nr)/20d0
    meanw = tau*meanw
    meanv = tau*meanv
    means = tau*means
    ntot = rhor(nr-1)/rhor(nr)
    
    !#DIR$ VECTOR ALIGNED
    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    Do ip=1,mynphi
       rho = vars1(:,ip,:,1)
       bc1 = binfo(1,1)%bdata(:,ip,:,1)
       bc2(:,1) = 0d0
       Call dbydr(rho,dwork,bc1,bc2,-2)

       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,1) = dvars1(it,ip,ir,1) - vars1(it,ip,ir,4)*dwork(it,ir)
             dvars1(it,ip,ir,4) = dvars1(it,ip,ir,4) - pconst*temperature(it,ip,ir)*invrho(it,ip,ir)*dwork(it,ir)
          EndDo
       EndDo
    EndDo

    If (do_entropy_rain) Then
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             !dvars1(it,ip,1:mynr-1,4) = dvars1(it,ip,1:mynr-1,4) + (mdvel + dvel_rain(it,ip) - vars1(it,ip,mynr,4))*taur(1:mynr-1)
             dvars1(it,ip,mynr,4) = dvars1(it,ip,mynr,4) + (mdvel + dvel_rain(it,ip) + vars1(it,ip,mynr-1,4) - 2d0*vars1(it,ip,mynr,4))*taur(mynr) + meanu
             !dvars1(it,ip,1:mynr-1,5) = dvars1(it,ip,1:mynr-1,5) + (sr(1)+delta_s + dent_rain(it,ip) - vars1(it,ip,mynr,5))*taur(1:mynr-1)
             dvars1(it,ip,mynr,5) = dvars1(it,ip,mynr,5) + (sr(1)+delta_s + dent_rain(it,ip) + vars1(it,ip,mynr-1,5) - 2d0*vars1(it,ip,mynr,5))*taur(mynr)
          EndDo
       EndDo
    Else
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,mynr,4) = dvars1(it,ip,mynr,4) + tau*(ntot*vars1(it,ip,mynr-1,4)-vars1(it,ip,mynr,4)) + meanu
             dvars1(it,ip,mynr,5) = dvars1(it,ip,mynr,5) + tau*(vars1(it,ip,mynr-1,5)-vars1(it,ip,mynr,5)) + means
          EndDo
       EndDo
    EndIf
  End Subroutine Top_Inflow_Outflow

  !Computes d2urdrdt, d2utdrdt, d2urdrdp, d2updrdp
  Subroutine Laplacian_Cross_Terms_Radial()
    Implicit None
    Real*8 :: workrb1(mynth,binfo(1,1)%bsize%bnr), workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: radrb1(mynth,binfo(1,1)%bsize%bnr), radrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8, Dimension(mynth,mynr) :: work1, work2, work3, work4, work5, work6, du
    Integer :: ir, it, ip

    Do ip=1,mynth
       work3(ip,:) = r_ind(myr1:myr2)
       work4(ip,:) = r_inv(myr1:myr2)
    EndDo

    If (myr1 .eq. 1) Then
       radrb1(:,1) = r1
    Else
       Do ip=1,mynth
          radrb1(ip,:) = r_ind(myr1-4:myr1-1)
       EndDo
    EndIf

    If (myr2 .eq. nr) Then
       radrb2(:,1) = r2
    Else
       Do ip=1,mynth
          radrb2(ip,:) = r_ind(myr2+1:myr2+4)
       EndDo
    EndIf

    If (non_uniform) Then
       Do ip=1,mynth
          work5(ip,:)=dxdr(myr1:myr2)
       EndDo
    EndIf

    invrho = vars2(:,:,:,2)
    vars2(:,:,:,2) = vars2(:,:,:,5)
    vars2(:,:,:,5) = invrho
    !#DIR$ VECTOR ALIGNED
    invrho =1d0/vars1(:,:,:,1)

    Call Boundary_Exchange_Lapl(1)

    !Set Boundaries
    If ((myr1 .eq. 1).and.(Bottom_Stress_Free)) Then
       !#DIR$ VECTOR ALIGNED
       dbinfo(1,1)%bdata(:,:,1,1) = vars2(:,:,1,1)/r1/dxdr(1)
       !#DIR$ VECTOR ALIGNED
       dbinfo(1,1)%bdata(:,:,1,1) = vars2(:,:,1,1)/r1/dxdr(1)
    EndIf

    If ((myr2 .eq. nr).and.(Bottom_Stress_Free)) Then
       !#DIR$ VECTOR ALIGNED
       dbinfo(1,2)%bdata(:,:,1,4) = vars2(:,:,mynr,4)/r2/dxdr(nr)
       !#DIR$ VECTOR ALIGNED
       dbinfo(1,2)%bdata(:,:,1,4) = vars2(:,:,mynr,4)/r2/dxdr(nr)
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    !#DIR$ VECTOR ALIGNED
    Do ip=1,mynphi
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynth
          !#DIR$ VECTOR ALIGNED
          work6(ir,:) = onethird*mur(myr1:myr2)*r_inv(myr1:myr2)*invrho(it+ir-1,ip,:)
       EndDo
       !d2utdrdt
       workrb1 = radrb1*dbinfo(1,1)%bdata(:,ip,:,1)
       workrb2 = radrb2*dbinfo(1,2)%bdata(:,ip,:,1)
       !#DIR$ VECTOR ALIGNED
       work1 = work3*vars2(:,ip,:,1)
       Call dbydr(work1,work2,workrb1,workrb2,2)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work2=work2*work5
       EndIf
       du = dvars1(:,ip,:,4)
       !#DIR$ VECTOR ALIGNED
       du = du + work6*work2

       !d2urdrdp
       workrb1 = radrb1*dbinfo(1,1)%bdata(:,ip,:,2)
       workrb2 = radrb2*dbinfo(1,2)%bdata(:,ip,:,2)
       !#DIR$ VECTOR ALIGNED
       work1 = work3*vars2(:,ip,:,2)
       Call dbydr(work1,work2,workrb1,workrb2,4)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work2=work2*work5
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,3) = dvars1(it,ip,ir,3) + work6(it,ir)*work2(it,ir)
          EndDo
       EndDo

       !d2urdrdt
       workrb1 = radrb1*dbinfo(1,1)%bdata(:,ip,:,3)
       workrb2 = radrb2*dbinfo(1,2)%bdata(:,ip,:,3)
       !#DIR$ VECTOR ALIGNED
       work1 = work3*vars2(:,ip,:,3)
       Call dbydr(work1,work2,workrb1,workrb2,4)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work2=work2*work5
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,2) = dvars1(it,ip,ir,2) + work6(it,ir)*work2(it,ir)
          EndDo
       EndDo

       !d2updrdp
       !#DIR$ VECTOR ALIGNED
       workrb1 = radrb1*dbinfo(1,1)%bdata(:,ip,:,4)
       !#DIR$ VECTOR ALIGNED
       workrb2 = radrb2*dbinfo(1,2)%bdata(:,ip,:,4)
       !#DIR$ VECTOR ALIGNED
       work1 = work3*vars2(:,ip,:,4)
       Call dbydr(work1,work2,workrb1,workrb2,3)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          work2=work2*work5
       EndIf
       !#DIR$ VECTOR ALIGNED
       dvars1(:,ip,:,4) = du + work6*work2
    EndDo

  End Subroutine Laplacian_Cross_Terms_Radial

  Subroutine Laplacian_Cross_Terms_Theta()
    Implicit None
    Real*8 :: sintb1(binfo(2,1)%bsize%bnt,mynphi), worktb1(binfo(2,1)%bsize%bnt,mynphi)
    Real*8 :: sintb2(binfo(2,2)%bsize%bnt,mynphi), worktb2(binfo(2,2)%bsize%bnt,mynphi)
    Real*8, Dimension(mynth,mynphi) :: work1, work2, work3, work4
    Real*8 :: ri
    Integer :: ir, it, ip

    Do ip=1,mynphi
       work3(:,ip) = sines(myth1:myth2)
    EndDo

    If (myth1.eq.1) Then
       sintb1(1,:) = sines(1)
    Else
       Do ip=1,mynphi
          sintb1(:,ip) = sines(myth1-4:myth1-1)
       EndDo
    EndIf

    If (myth2.eq.nth) Then
       If (Theta_Symmetric) Then
          Do ip=1,mynphi
             sintb2(:,ip) = sin(th2+dth*(/(Dble(it),it=1,4)/))
          EndDo
       Else
          sintb2(1,:) = sines(nth)
       EndIf
    Else
       Do ip=1,mynphi
          sintb2(:,ip) = sines(myth2+1:myth2+4)
       EndDo
    EndIf

    Call Boundary_Exchange_Lapl(2)

    !Set Boundaries
    If (myth1 .eq. 1) Then
       !#DIR$ VECTOR ALIGNED
       dbinfo(2,1)%bdata(1,:,:,4) = vars2(1,:,:,4)*cosines(1)/sines(1)
    EndIf

    If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
       !#DIR$ VECTOR ALIGNED
       dbinfo(2,2)%bdata(1,:,:,4) = vars2(mynth,:,:,4)*cosines(nth)/sines(nth)
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do ir=1,mynr
       ri = r_inv(myr1+ir-1)*mur(myr1+ir-1)
       !d2updtdp
       !#DIR$ VECTOR ALIGNED
       worktb1 = sintb1*dbinfo(2,1)%bdata(:,:,ir,4)
       !#DIR$ VECTOR ALIGNED
       worktb2 = sintb2*dbinfo(2,2)%bdata(:,:,ir,4)
       !#DIR$ VECTOR ALIGNED
       work1 = work3*vars2(:,:,ir,4)
       Call dbydt(work1,work2,worktb1,worktb2,3)
       !#DIR$ VECTOR ALIGNED
       work4 = onethird*ri*invrho(:,:,ir)/work3
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,2) = dvars1(it,ip,ir,2) + work4(it,ip)*work2(it,ip)
          EndDo
       EndDo

       !d2utdtdp
       !#DIR$ VECTOR ALIGNED
       worktb1 = sintb1*dbinfo(2,1)%bdata(:,:,ir,5)
       !#DIR$ VECTOR ALIGNED
       worktb2 = sintb2*dbinfo(2,2)%bdata(:,:,ir,5)
       !#DIR$ VECTOR ALIGNED
       work1 = work3*vars2(:,:,ir,5)
       Call dbydt(work1,work2,worktb1,worktb2,2)
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ip=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do it=1,mynth
             dvars1(it,ip,ir,3) = dvars1(it,ip,ir,3) + work4(it,ip)*work2(it,ip)
          EndDo
       EndDo
    EndDo
  End Subroutine Laplacian_Cross_Terms_Theta

  Subroutine Induction_Terms()
    Implicit None
    Real*8 :: jconst
    Integer :: ir, it, ip, iv

    If (Laplacian) Then
       jconst = 4d0*pi/c_light**2
       !CurlB
       Call VectorCurl()
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          !#DIR$ VECTOR ALIGNED
          vars2(:,:,ir,1) = -etar(myr1+ir-1)*dvars1(:,:,ir,6)
          !#DIR$ VECTOR ALIGNED
          vars2(:,:,ir,2) = -etar(myr1+ir-1)*dvars1(:,:,ir,7)
          !#DIR$ VECTOR ALIGNED
          vars2(:,:,ir,3) = -etar(myr1+ir-1)*dvars1(:,:,ir,8)
          !Ohmic Heating = 4 pi eta j^2/c^2
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do ip=1,mynphi
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do it=1,mynth
                dvars1(it,ip,ir,5) = dvars1(it,ip,ir,5) + (jconst*etar(myr1+ir-1))*(dvars1(it,ip,ir,6)**2+dvars1(it,ip,ir,7)**2+dvars1(it,ip,ir,8)**2)*invrho(it,ip,ir)/temperature(it,ip,ir)
             EndDo
          EndDo
       EndDo
    EndIf

    If (mod(istep,skip) .eq. 0) Then 
       !Compute div.B
       Call VectorDivergence()
       divBtot = 0d0
       Do ir=1,mynr
          Do ip=1,mynphi
             divBtot = divBtot + r_ind(myr1+ir-1)**2*dr_ind(myr1+ir-1)*Sum(sines(myth1:myth2)*abs(invrho(:,ip,ir)))
          EndDo
       EndDo
       divBmax = Maxval(abs(invrho))
    EndIf

    !Compute vxB
    !#DIR$ VECTOR ALIGNED
    dvars1(:,:,:,6) = vars1(:,:,:,3)*vars1(:,:,:,8)-vars1(:,:,:,4)*vars1(:,:,:,7)
    !#DIR$ VECTOR ALIGNED
    dvars1(:,:,:,7) = vars1(:,:,:,4)*vars1(:,:,:,6)-vars1(:,:,:,2)*vars1(:,:,:,8)
    !#DIR$ VECTOR ALIGNED
    dvars1(:,:,:,8) = vars1(:,:,:,2)*vars1(:,:,:,7)-vars1(:,:,:,3)*vars1(:,:,:,6)
    If (Laplacian) Then
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do iv=1,3
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do it=1,mynth
                   dvars1(it,ip,ir,iv+5) = dvars1(it,ip,ir,iv+5)+vars2(it,ip,ir,iv)
                EndDo
             EndDo
          EndDo
       EndDo
    EndIf
  End Subroutine Induction_Terms

  Subroutine Drag_Toward_Solar()
    Implicit None
    Real*8, Dimension(:), Allocatable :: mvp,snd,rcv
    Integer :: it
    Allocate(mvp(nth),snd(nth),rcv(nth))

    mvp = 0d0
    snd = 0d0
    rcv = 0d0
    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=64
    Do it=1,mynth
       mvp(myth1+it-1) = Sum(vars1(it,:,1,3))
    EndDo
    snd = mvp
    Call Global_AllReduce(snd,rcv,bot_comm,'sum')
    mvp = rcv/Dble(nphi)

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=64
    Do it=1,mynth
       dvars1(it,:,1,3) = dvars1(it,:,1,3) + (vphi_prof(myth1+it-1)-mvp(myth1+it-1))/damping_time
    EndDo

    Deallocate(mvp,snd,rcv)
  End Subroutine Drag_Toward_Solar

  !========================================================================

  Subroutine Bound()
    Implicit None
    Integer :: ip
    !  set velocities on latitudinal boundaries
    If (myth1 .eq. 1) Then
       vars1(1,:,:,2) = 0d0
       If (.not. Theta_Stress_Free) vars1(1,:,:,4) = 0d0
    EndIf

    If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
       vars1(mynth,:,:,2)=0d0
       If (.not. Theta_Stress_Free) vars1(mynth,:,:,4) = 0d0
    EndIf

    !  set velocities on lower radial boundary
    If (myr1.Eq.1) Then
       !  meridional (latitudinal) component
       If (.not. Bottom_Stress_Free) Then
          vars1(:,:,1,2)=0d0
          vars1(:,:,1,3)=0d0
       EndIf

       If (.not. bot_inout) Then
          vars1(:,:,1,4)=0d0
          If (Bottom_Constant_Temp) Then
             Do ip=1,mynphi
                vars1(:,ip,1,5)=sr1
             EndDo
          EndIf
       EndIf
    EndIf

    !  set velocities on upper radial boundary
    If (myr2.Eq.nr) Then
       If (.not. top_inout) Then
          !  meridional (latitudinal) component
          If (.not. Top_Stress_Free) Then
             vars1(:,:,mynr,2)=0d0
             vars1(:,:,mynr,3)=0d0
          EndIf

          vars1(:,:,mynr,4) = 0d0
          If ((Top_Constant_Temp).and.(Set_Stop)) Then
             Do ip=1,mynphi
                vars1(:,ip,mynr,5)=sr2
             EndDo
          EndIf
       EndIf
    EndIf

    If (magnetic) Then
       !  zero normal mag field on latitudinal boundaries       
       If (Theta_Perf_Conductor) Then
          If (myth1 .eq. 1) Then
             If (.not.Unstratified) vars1(1,:,:,7)=0d0
             vars1(1,:,:,8)=0d0
          EndIf
          If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
             If (.not.Unstratified) vars1(mynth,:,:,7)=0d0
             vars1(mynth,:,:,8)=0d0
          EndIf
       Else
          !  zero tangential field on theta boundaries
          If (myth1 .eq. 1) Then
             vars1(1,:,:,6)=0d0
          EndIf
          If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
             vars1(mynth,:,:,6)=0d0
          EndIf
       EndIf

       !  zero tangential field on radial boundaries
       If (myr1.Eq.1) Then
          If (.not. bot_inout) Then
             If (Bottom_Radial_Field) Then
                If (Unstratified) vars1(:,:,1,6) = 0d0
                vars1(:,:,1,8)=0d0
             ElseIf (Bottom_Perfect_Cond) Then
                vars1(:,:,1,6) = 0d0
                vars1(:,:,1,7) = 0d0
             EndIf
          EndIf
       EndIf
       If (myr2 .eq. nr) Then
          If (.not. top_inout) then
             If (Top_Radial_Field) Then
                If (Unstratified) vars1(:,:,mynr,6) = 0d0
                vars1(:,:,mynr,8)=0d0
             ElseIf (Top_Perfect_Cond) Then
                vars1(:,:,mynr,6) = 0d0
                vars1(:,:,mynr,7) = 0d0
             EndIf
          EndIf
       EndIf

    EndIf
  End Subroutine Bound

End Module Physics
