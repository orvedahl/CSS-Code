!-------------------------------------------------------------------------
!
!        SPHERICAL COMPRESSIBLE MHD CODE 
!
!    K. Augustson, N. Hurlburt, M. DeRosa  2009
!
!-------------------------------------------------------------------------
Program CSS2

  Use Time_Step_RK4
  Use Report
  Use Validation
  Use BesselFunction

  Implicit None

  Real*8 :: main_tstart, main_tend,ttime,stime,ftime, eps_exp !Timing variables
  Real*8, Allocatable, Dimension(:,:) :: full_timings, timings_buff
  Integer :: namelist_unit=20 !Input namelist
  Logical :: update_timestep=.False.

  Call Initialize_Parallel_Part1()

  Call Setup()

  Call cpu_time(main_tstart)
  Call Compute_Step_Size(deltat)
  Do istep=begiter+1,begiter+niters

     !  step forward
     deltat_old=deltat
     If (update_timestep.and.(mod(istep,skip) .eq. 0)) Then 
        Call Compute_Step_Size(deltat)
     EndIf

     If (RK2TVD) Then
        Call step_rk2_tvd()
     ElseIf (RK3TVD) Then
        Call step_rk3_tvd()
     Else
        Call step_rk4()
     EndIf
     sim_time=sim_time+deltat
     Call Output_Check(.False.)

     !Call Barrier(node_comm)
     If ((myrank .eq. 0).and.(mod(istep,skip) .eq. 0)) Print*, 'Iteration', istep, 'at time', sim_time, 'with time step', deltat
  Enddo
  Call cpu_time(main_tend)
  ttime=main_tend-main_tstart

  If (Do_Timings) Then 
     Allocate(full_timings(17,nnodes+1),timings_buff(17,nnodes+1))
     full_timings = 0d0
     full_timings(1:16,myrank+2) = timings
     full_timings(17,myrank+2) = ttime
     Call MPI_REDUCE(full_timings(1,1),timings_buff(1,1),(nnodes+1)*17,MPI_REAL8,MPI_SUM,0,node_comm,ierr)
     full_timings = timings_buff
     If (myrank .eq. 0) Then !Write out timings to timings.dat
        full_timings(1,1) = nnodes
        full_timings(2:17,1) = 0d0
        Inquire(iolength=irec_azavg) full_timings
        Open(namelist_unit, file=Trim(Trim(perm_dir)//'timings.dat'), form='unformatted', access='direct', status='unknown',recl=irec_azavg)
        Write(namelist_unit, rec=1) full_timings
        Close(namelist_unit)
     EndIf
     Deallocate(full_timings,timings_buff,timings)
  EndIf

  !  print closing messages
  If (myrank.eq.0) Then
    Print*
    Print*,'Calculation completed successfully!'
    Print*,'  Number of processors      = ',nnodes
    Print*,'  Time elapsed (sec)         = ',ttime
    Print*,'  ',niters,' steps at ',niters/ttime,' step/sec'
    Print*
  Endif

  Call Finalize_Diffusion()
  Call Finalize_Step_RK4()
  Call Finalize_Physics()
  Call Finalize_IO()
  Call Finalize_Derivatives()
  Call Finalize_Parallel()
  Call Graceful_Exit()

Contains

  !------------------------------------------------------------------------

  Subroutine Setup()
    Implicit None
    Real*8 :: input_time, para_time, phys_time, der_time, cond_time, io_time(3), sbr8(8), rbr8(8)

    !  read in namelists and set boundary conditions
    Call cpu_time(stime)
    Call Read_Input()
    Call cpu_time(ftime)

    input_time = ftime-stime
    
    !Print*,'Input Read.'
    !  initialize variables within various modules
    Call cpu_time(stime)
    Call Initialize_Parallel_Part2()
    Call cpu_time(ftime)
    para_time = ftime-stime

    !Print*,'Parallel Module Initialized.'
    Call cpu_time(stime)
    Call Initialize_Derivatives()
    Call cpu_time(ftime)
    der_time = ftime-stime

    !Print*,'Derivatives Module Initialized.'
    Call cpu_time(stime)
    Call Initialize_Physics()
    Call cpu_time(ftime)
    phys_time = ftime-stime

    Call Initialize_Step_RK4()

    !Print*,'Phyics Module Initialized.'
    Call cpu_time(stime)
    Call Initialize_IO()
    Call cpu_time(ftime)
    io_time(1) = ftime-stime


    !Print*,'I/O Initialized.'
       
    !  print some stuff
    If (myrank.eq.0) Then
       Print*
       Print*,'Problem Size and Array Bounds of Current Run'
       Print*,'  nr,nth,nphi,nv = ',nr,nth,nphi,nv
       Print*,'  r1,r2 =',r1,r2
       Print*,'  th1,th2 =',th1,th2
       Print*,'  phi1,phi2 =',phi1,phi2
       Print*
       Print*,'Boundary Conditions of Current Run'
       Select Case (lbcr_v)
       Case (1)  ;  Print*,'    bottom = no-slip,       top = no-slip'
       Case (2)  ;  Print*,'    bottom = stress-free,   top = no-slip'
       Case (3)  ;  Print*,'    bottom = no-slip,       top = stress-free'
       Case (4)  ;  Print*,'    bottom = stress-free,   top = stress-free'
       Case Default  !  shouldn't be able to get here
          Print*,'  ERROR in Setup: boundary conditions set incorrectly'
          Call Graceful_Exit()
       EndSelect
       Print*,'  Radial Temperature BC'
       Select Case (lbcr_t)
       Case (1)  ;  Print*,'    bottom = constant temp, top = constant temp'
       Case (2)  ;  Print*,'    bottom = constant flux, top = constant temp'
       Case (3)  ;  Print*,'    bottom = constant temp, top = constant flux'
       Case (4)  ;  Print*,'    bottom = constant flux, top = constant flux'
       Case Default  !  shouldn't be able to get here
          Print*,'  ERROR in Setup: boundary conditions set incorrectly'
          Call Graceful_Exit()
       EndSelect
       Print*,'  Latitudinal Velocity BC'
       Print*,'    bottom = impenetrable,  top = impenetrable'
       Print*,'    bottom = stress-free,   top = stress-free'
       Print*,'  Latitudinal Temperature BC'
       Print*,'    bottom = conducting,    top = conducting'
       Print*
       Print*,'Control Parameters of Current Run'
       Print*,'  nsteps = ',niters
       Print*,'  cfl_safety = ',cfl_safety
       Print*
       If (magnetic) Then
          Print*,'  Do Lorentz forces =',do_lorentz
       Else
          Print*,'Solution is nonmagnetic'
       Endif
    Endif

    If (Validate) Then
       Call TestDerivatives()
    EndIf

    !  either set input conditions or read checkpoint file to initialize
    If (Start_From_Checkpoint) Then 
       If (myrank.Eq.0) Then
          Print*,'Restarting from ',Trim(input_checkpoint_filename)
       Endif
       
       !  in case input file is nonmagnetic, set initial magnetic conditions
       Call cpu_time(stime)
       Call Checkpoint_Input(input_checkpoint_filename)

       !  read in from previous checkpoint
       Call Read_Stellar_Model()
       Call cpu_time(ftime)
       io_time(2) = ftime-stime
             
       Call cpu_time(stime)       
       !Initialize the thermal rain bc 
       Call Initialize_Flow(Start_From_Checkpoint)

       If (magnetic) Then
          Call Initialize_Mag(Start_From_Checkpoint)
       EndIf

       If ((do_entropy_rain).and.(myr2.eq.nr)) Then
          Call Initialize_Rain_Squalls()
       End If
          
       Call cpu_time(ftime)
       cond_time = ftime-stime
       
       Call cpu_time(stime)
       Call Initialize_Diffusion()
       istep = begiter
       Call Output_Check(Start_From_Checkpoint)
       Call cpu_time(ftime)
       io_time(3) = ftime-stime
    Else
       If (myrank.Eq.0) Then
          Print*,'        '
          Print*,'Applying initial conditions'
       Endif
       
       Call cpu_time(stime)
       !  initialize rho,t,u,v,w, and br,bt,bp if magnetic
       Call Read_Stellar_Model()
       If (myrank .eq. 0) Then
          Print*,'        '
          Print*,'Stellar Model Read.'
       EndIf
       Call cpu_time(ftime)
       io_time(2) = ftime-stime
       
       Call cpu_time(stime)
       !Initialize the thermal rain bc 
       Call Initialize_Flow(Start_From_Checkpoint)

       If (magnetic) Call Initialize_Mag(Start_From_Checkpoint)

       If ((do_entropy_rain).and.(myr2.eq.nr)) Then
          Call Initialize_Rain_Squalls()
       End If
       
       If (myrank .eq. 0) Then
          Print*,'        '
          Print*,'Flow initialized.'
       EndIf
       Call cpu_time(ftime)
       cond_time = ftime-stime

       begiter=0
       istep=0
       Call cpu_time(stime)
       Call Initialize_Diffusion()
       !  do an initial checkpoint
       Call Output_Check(Start_From_Checkpoint)
       Call cpu_time(ftime)
       io_time(3) = ftime-stime
       If (myrank.eq.0) Print*, 'Output Checked'
    EndIf

    sbr8 = (/input_time,para_time,phys_time,der_time,cond_time,io_time(1),io_time(2),io_time(3)/)
    Call Global_Reduce(sbr8,rbr8,0,node_comm,'max')
    input_time = Sum(rbr8)

    If (Do_Timings) Then
       Allocate(timings(16))
       timings = 0d0
       timings(1) = Sum(sbr8)
    EndIf

    If (myrank .eq. 0) Then
       Print*, 'Start up times (read input, init MPI cart, physics, derivatives, init conds., I/O):', rbr8
       Print*, 'Total startup time:', input_time
    EndIf
  End Subroutine Setup

  !------------------------------------------------------------------------

  Subroutine Read_Stellar_Model()
    Implicit None
    Real*8 :: bc1, bc2, Gconst
    Real*8 :: beta, temp, drin, hrho, alp, temparr(nr)
    Integer :: unit = 17, model_nr
    Real*8, Allocatable, Dimension(:) :: din, dout
    Real*8, Allocatable, Dimension(:,:) :: model
    Integer :: ir

    gamm1 = gamma-1d0
    Cv = Cp/gamma
    If (Isothermal) Then
       gravity = g0
       hrho = g0/(gamm1*Cv*T0)
       rhor = rho0*dexp(-(r_ind-r1)*hrho)
       rhob1(4) = rho0*dexp(dr_ind(1)*hrho)
       rhob1(3) = rho0*dexp(2d0*dr_ind(1)*hrho)
       rhob1(2) = rho0*dexp(3d0*dr_ind(1)*hrho)
       rhob1(1) = rho0*dexp(4d0*dr_ind(1)*hrho)

       rhob2(1) = rho0*dexp(-(r2-r1+dr_ind(nr))*hrho)
       rhob2(2) = rho0*dexp(-(r2-r1+2d0*dr_ind(nr))*hrho)
       rhob2(3) = rho0*dexp(-(r2-r1+3d0*dr_ind(nr))*hrho)
       rhob2(4) = rho0*dexp(-(r2-r1+4d0*dr_ind(nr))*hrho)
       tr = T0
       pr = gamm1*Cv*T0*rhor
       sr = Cv*(dlog(gamm1*Cv*T0)-gamm1*dlog(rhor))
       dsdrm = g0/T0
       kr = 0d0
       dlnkr = 0d0
       eps = 0d0
       ks = kappa_value
       dlnks = 0d0
       mur = kappa_value*rhor
       dlnmu = -hrho
       etar = ks       
    ElseIf (Isentropic) Then
       Gconst = 6.67d-8
       rhor2 = rho0
       Tr2 = T0
       Pr2 = gamm1*Cv*rhor2*Tr2
       Sr2 = Cv*(dlog(Pr2)-gamma*dlog(rhor2))
       alp = Gconst*mstar/(Cp*T0*r2)

       rhor(1) = rho0*(1d0+alp*(r2/r1-1d0))**(1d0/gamm1)

       If (non_uniform) Then
          beta = 0.23d0
          r_ind = alp*r2/(alp - 1d0 + ((1d0-(rhor(1)/rho0)**beta)*x_ind + (rhor(1)/rho0)**beta)**(gamm1/beta))          
          r_inv = 1d0/r_ind

          dr_ind(1) = r_ind(1)-alp*r2/(alp - 1d0 + (-(1d0-(rhor(1)/rho0)**beta)/Dble(nr-1) + (rhor(1)/rho0)**beta)**(gamm1/beta))
          Do ir=2,nr
             dr_ind(ir) = r_ind(ir)-r_ind(ir-1)
          EndDo
       EndIf

       rhor = rho0*(1d0+alp*(r2*r_inv-1d0))**(1d0/gamm1)

       If (myrank .eq. 0) Print*, 'Density Contrast:', rhor(1)/rho0

       gravity = Gconst*mstar*r_inv**2
       Tr = T0*(1d0+alp*(r2*r_inv-1d0))
       Pr = gamm1*Cv*T0*rho0*(1d0+alp*(r2*r_inv-1d0))**(gamma/gamm1)
       Sr = Sr2(1) !Adiabat
       inv_drag_time = thermal_drag_time*(1d0+alp*(r2*r_inv(myr1:myr2)-1d0))*(r_ind(myr1:myr2)/r1)**2/(1d0+alp*(r2/r1-1d0))
       inv_drag_time = 1d0/inv_drag_time

       ks = kappa_value*rhor*Tr*(rhor/rhor(nr))**g0

       etar = ks
       mur = sld_pr*kappa_value*rhor
       dlnks = -(gamma+g0)*alp*r2*r_inv**2/(1d0+alp*(r2*r_inv-1d0))/gamm1
       dlnmu = -alp*r2*r_inv**2/(1d0+alp*(r2*r_inv-1d0))/gamm1

       !initialize model_deriv
       dsdrm = 0d0
       drhodrm = -alp*r2*r_inv**2/(1d0+alp*(r2*r_inv-1d0))/gamm1

       If (non_uniform) Then
          dxdr = beta*(rhor**beta)*drhodrm/(rho0**beta-rhor(1)**beta)
          d2xdr2  = beta*(rhor**beta)*((beta-gamma+1d0)*drhodrm +2d0*r_inv)*drhodrm/(rho0**beta-rhor(1)**beta)
       EndIf

       If (do_radheat) Then
          Allocate(din(nr),dout(nr))
          eps = 0d0
          drin = Dble(nr-1)/(x_ind(nr)-x_ind(1))
          Lrad = lumr1 + (lumr2-lumr1)*((r_ind-r1)/(r2-r1))**eps_exp
          dlnkr = tr*(dsdrm/Cp-gamm1*rhor*gravity/(gamma*pr))
          kr = -16d0*pi*a_const*c_light*(r_ind**2)*(tr**3)*dlnkr/(3d0*rhor*Lrad)
          din = dlog(kr)
          bc1 = din(1)/drin
          bc2 = din(nr)/drin
          dout = 0d0
          Call dbyd1_dt0(din,dout,bc1,bc2,nr,1)
          dlnkr = drin*dout
          If (non_uniform) dlnkr = dlnkr*dxdr

          kr = 1d0/kr !Inverse rosseland opacity.

          k0 = Lrad/(4d0*pi*r_ind**2)
          If (myrank .eq. 0) Then
             Print*, 'kr=[',kr,']'
             Print*, 'dlnkr=[',dlnkr,']'
          EndIf
       EndIf

       rhob1(4) = rho0*(1d0+alp*(r2/(r1-dr_ind(1))-1d0))**(1d0/gamm1)
       rhob1(3) = rho0*(1d0+alp*(r2/(r1-2d0*dr_ind(1))-1d0))**(1d0/gamm1)
       rhob1(2) = rho0*(1d0+alp*(r2/(r1-3d0*dr_ind(1))-1d0))**(1d0/gamm1)
       rhob1(1) = rho0*(1d0+alp*(r2/(r1-4d0*dr_ind(1))-1d0))**(1d0/gamm1)

       rhob2(1) = rho0*(1d0+alp*(r2/(r2+dr_ind(nr))-1d0))**(1d0/gamm1)
       rhob2(2) = rho0*(1d0+alp*(r2/(r2+2d0*dr_ind(nr))-1d0))**(1d0/gamm1)
       rhob2(3) = rho0*(1d0+alp*(r2/(r2+3d0*dr_ind(nr))-1d0))**(1d0/gamm1)
       rhob2(4) = rho0*(1d0+alp*(r2/(r2+4d0*dr_ind(nr))-1d0))**(1d0/gamm1)

    ElseIf (Unstratified) Then
       rhor = rho0
       Tr = T0
       Pr = g0*rho0
       Sr = 2d0*Cv*dlog(g0)
       gravity = 0d0
       dsdrm = 0d0
       rhob1 = rho0
       rhob2 = rho0

       If (non_uniform) Then
          r_ind = r1 + 5d0*(r2-r1)*x_ind*(1d0-9d0*x_ind*(1d0-2d0*x_ind/3d0)/5d0)/2d0
          r_inv = 1d0/r_ind

          dr_ind(1) = r_ind(2)-r_ind(1)
          Do ir=2,nr
             dr_ind(ir) = r_ind(ir)-r_ind(ir-1)
          EndDo
          dxdr = 2d0/((1d0-18d0*x_ind*(1d0-x_ind)/5d0)*(r2-r1)*5d0)
          d2xdr2 = 1d0/(9d0*(r2-r1)*(2d0*x_ind-1d0))

          If (myrank .eq. 0) Then
             Print*, 'r=[',r_ind,']'
             Print*, 'dxdr=[',dxdr,']'
          EndIf
       EndIf
    Else
       Allocate(model(nr,6))
       
       If (myrank .eq. 0) Then
          ! read in the model
          open(unit,file=Trim(stellar_model_filename),status='old',form='formatted')
          
          read(unit,101) model_nr
          
          If (model_nr .ne. nr) Then
             Print*, 'Model background size differs from input size'
             close(unit)
             Call Graceful_Exit()
          End If
          
          Do ir=1,nr
             read(unit, 102) temp
             gravity(ir) = temp
          End Do
          Do ir=1,nr
             read(unit, 102) temp
             rhor(ir) = temp
          End Do
          Do ir=1,nr
             read(unit, 102) temp
             sr(ir) = temp
          End Do
          Do ir=1,nr
             read(unit, 102) temp
             dsdrm(ir) = temp
          End Do
          Do ir=1,nr
             read(unit,102) temp
             Lrad(ir) = temp
          EndDo
       
          Do ir=1,4
             read(unit,102) temp
             rhob1(ir) = temp
          EndDo

          Do ir=1,4
             read(unit,102) temp
             rhob2(ir) = temp
          EndDo

          close(unit)
          
101       format(i5)
102       format(e16.9)
       EndIf
       
       !Distribute background
       model = 0d0
       If (myrank .eq. 0) Then
          model(:,1) = gravity(:)
          model(:,2) = rhor(:)
          model(:,3) = sr(:)
          model(:,4) = dsdrm(:)
          model(:,5) = Lrad(:)
          model(1:4,6) = rhob1
          model(5:8,6) = rhob2
          
          Do ir=1,nnodes-1
             Call Send(model,ir,ir,node_comm)
          EndDo
       Else
          Call Receive(model,0,myrank,node_comm)
          gravity = model(:,1)
          rhor    = model(:,2)
          sr      = model(:,3)
          dsdrm   = model(:,4)
          Lrad    = model(:,5)
          rhob1 = model(1:4,6)
          rhob2 = model(5:8,6)
       EndIf
       
       Deallocate(model)
       
       gamm1 = gamma-1d0
       Cv = Cp/gamma
       tr = gravity
       gravity = gravity - r_ind*omega_0**2
       rhor = Compute_Density(sr)
       gravity = tr

       tr = rhor**gamm1*dexp(sr/Cv)/gamm1/Cv
       pr = gamm1*Cv*rhor*tr

       drin = Dble(nr-1)/(x_ind(nr)-x_ind(1))
       Allocate(din(nr),dout(nr))

       Lrad = Lrad*lumr1/Lrad(1)

       din = Lrad
       bc1 = Lrad(1)/drin
       bc2 = Lrad(nr)/drin
       dout = 0d0
       Call dbyd1_dt0(din,dout,bc1,bc2,nr,1)
       dlnkr = drin*dout
       dlnkr = -dlnkr/(4d0*pi*r_ind**2)
       If (non_uniform) dlnkr = dlnkr*dxdr
       eps = dlnkr
       eps(nr-8:nr) = 0d0 !zero out spurious wiggles near upper boundary
       k0 = Lrad/(4d0*pi*r_ind**2)

       If (lrad_update) Lrad_new = Lrad

       If (do_radheat) Then
          eps = 0d0
          Lrad = Lrad + lumr2*((r_ind-r1)/(r2-r1))**eps_exp
          dlnkr = tr*(dsdrm/Cp-gamm1*rhor*gravity/(gamma*pr))
          kr = -16d0*pi*a_const*c_light*(r_ind**2)*(tr**3)*dlnkr/(3d0*rhor*Lrad)
          din = dlog(kr)
          bc1 = din(1)/drin
          bc2 = din(nr)/drin
          dout = 0d0
          Call dbyd1_dt0(din,dout,bc1,bc2,nr,1)
          dlnkr = drin*dout
          If (non_uniform) dlnkr = dlnkr*dxdr

          kr = 1d0/kr !Inverse rosseland opacity.

          k0 = Lrad/(4d0*pi*r_ind**2)
          If (myrank .eq. 0) Then
             Print*, 'kr=[',kr,']'
             Print*, 'dlnkr=[',dlnkr,']'
          EndIf
       EndIf

       Deallocate(din,dout)
    EndIf

    If (ell0_heating) Then
       eps = 0d0
       k0 = 0d0
       alp = eps_exp
       k0(1) = -lumr2/(4d0*pi*(r2-r1)*((r2-r1)**2/(alp+3d0)+r1*(r2-r1)/(alp+2d0)+r1*r1/(alp+1d0)))
       !Cooling carrying flux out
       eps = dlnkr + k0(1)*((r_ind-r1)/(r2-r1))**alp
       k0 = -eps*((r_ind-r1)**3/(alp+3d0)+r1*(r_ind-r1)**2/(alp+2d0)+r1*r1*(r_ind-r1)/(alp+1d0))/r_ind**2
       k0 = k0 + Lrad/(4d0*pi*r_ind**2)
       !Solar-like radiative heating for deep shell
       !eps = eps + 3.55d0*lumr1*((r2-r_ind)/(r2-r1))**2.55d0/(4d0*pi*r_ind**2*(r2-r1))
       !k0 = (lumr2*((r_ind-r1)/(r2-r1))**alp + lumr1*((r2-r_ind)/(r2-r1))**3.55d0)/(4d0*pi*r_ind**2)
    EndIf

    !initialize model_deriv
    dsdr1 = dsdrm(1)
    dsdr2 = dsdrm(nr)

    If (SLaplacian) Then
       If (Bottom_Constant_Flux) Then
          dsdr1 = -lumr1*one_over_4pi/(ks(1)*r1**2)
       EndIf

       If (Top_Constant_Flux) Then
          dsdr2 = -lumr2*one_over_4pi/(ks(nr)*r2**2)
       EndIf
    EndIf


!!$    If (SLaplacian) Then
!!$       temparr = (k0*4d0*pi*r_ind**2-lstar)/(4d0*pi*ks*r_ind**2)
!!$       dsdr1 = temparr(1)
!!$       dsdr2 = temparr(nr)
!!$
!!$       If (Top_Constant_Temp) Then
!!$          If (Bottom_Constant_Flux) Then
!!$             Do ir=nr-1,1,-1
!!$                sr(ir) = Sum(temparr(ir:nr)*dr_ind(ir:nr))
!!$             EndDo
!!$             sr(1:nr-1) = sr(nr)-sr(1:nr-1)
!!$          Else
!!$             Do ir=nr-1,1,-1
!!$                sr(ir) = Sum(temparr(ir:nr)*dr_ind(ir:nr))
!!$             EndDo
!!$
!!$             Do ir=2,nr
!!$                sr(ir) = sr(ir) + Sum(temparr(1:ir)*dr_ind(1:ir))
!!$             EndDo
!!$             sr = 0.5d0*sr
!!$             sr1 = Cv*(dlog(pr(1))-gamma*dlog(rhor(1)))
!!$             sr2 = delta_s+sr1
!!$             sr = (sr1(1)-sr(1)*sr2(1)/sr(nr))/(1d0-sr(1)/sr(nr)) + (sr1(1)-sr2(1))*sr/(sr(1)-sr(nr))
!!$          EndIf
!!$       Else
!!$         If (Bottom_Constant_Flux) Then
!!$             Do ir=nr-1,1,-1
!!$                sr(ir) = Sum(temparr(ir:nr)*dr_ind(ir:nr))
!!$             EndDo
!!$
!!$             Do ir=2,nr
!!$                sr(ir) = sr(ir) + Sum(temparr(1:ir)*dr_ind(1:ir))
!!$             EndDo
!!$             sr = 0.5d0*sr
!!$
!!$             sr(1:nr-1) = sr2(1)+sr(1:nr-1)
!!$          Else
!!$             Do ir=2,nr
!!$                sr(ir) = Sum(temparr(1:ir)*dr_ind(1:ir))
!!$             EndDo
!!$             sr(2:nr) = sr(1) + sr(2:nr)
!!$          EndIf
!!$       EndIf
!!$       MeanV(:,5) = sr
!!$       rhor = Compute_Density()
!!$    EndIf

    rhor1 = rhor(1)
    rhor2 = rhor(nr)
    Tr1 = tr(1)
    Tr2 = tr(nr)
    Sr1 = sr(1)
    Sr2 = sr(nr)
    Pr1 = pr(1)
    Pr2 = pr(nr)
   
  End Subroutine Read_Stellar_Model

  Subroutine Initialize_Flow(checkpoint_start)
    Implicit None
    Integer :: seed(1)
    Real*8 :: xvar, radius, drad, sinth, Etot, av, bv, cvar, Gconst
    Real*8 :: sendb(5), recvb(5), omtmp(mynth), iomtmp(mynth), diomtmp(mynth), rand(mynth,mynphi)
    Real*8 :: rhot(mynphi), wt(mynphi), vt(mynphi), ut(mynphi), tt(mynphi), yvar(mynphi), temparr(nr), density(mynth,nr)
    Logical, Intent(In) :: checkpoint_start
    Integer :: ir, it, ip, iv

    !Generate small scale perturbations in entropy
    vars1(:,:,:,1:5) = 0d0

    Call system_clock(seed(1))
    Call random_seed(put=seed)

    If (Initialize_MRI) Then    
       Do ir=1,mynr
          xvar = 2d0*sqrt(pi*rhor(myr1+ir-1))*(omega_0*r1**2)*r_inv(myr1+ir-1)/bamp
          Do ip=1,mynphi
             vars1(:,ip,ir,5) = cos(xvar*sines(myth1:myth2)*(theta(myth1:myth2)-0.5d0*pi))*(1d0-cos((r_ind(myr1+ir-1)-r1)/(r2-r1)*pi)**128)
          EndDo
       EndDo

       Do ir=1,mynr
          Do ip=1,mynphi
             xvar = r_ind(myr1+ir-1)*((omega_0*(r1*r_inv(myr1+ir-1))**2)**2)
             Etot = rhor(myr1+ir-1)*r_ind(myr1+ir-1)*xvar/pr(myr1+ir-1)/gamma
             xvar = xvar*(4d0+r_ind(myr1+ir-1)*gravity(myr1+ir-1)*dexp(-sr(1)/Cv)/rhor(myr1+ir-1)**gamm1/gamma)/gravity(myr1+ir-1)
             vars1(:,ip,ir,1) = rhor(myr1+ir-1)*(1d0 + 0.5d0*xvar*sines(myth1:myth2)**2)
             vars1(:,ip,ir,3) = omega_0*r_ind(myr1+ir-1)*sines(myth1:myth2)*((r1*r_inv(myr1+ir-1))**2-1d0)
             vars1(:,ip,ir,5) = sr(myr1+ir-1)*(1d0 + samp*vars1(:,ip,ir,5))+0.5d0*Cp*(Etot-xvar)*sines(myth1:myth2)**2
          EndDo          
       EndDo
    ElseIf (.not. solar_prof) Then !Rigid rotation
       If (top_inout) Then
          If (bot_inout) Then
             Do ir=1,mynr
                Call random_number(harvest=rand)
                Do it=1,mynth
                   xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                   vars1(it,:,ir,5) = (2d0*rand(it,:)-1d0)*(1d0-cos(xvar)**16)
                EndDo
             EndDo
          Else
             Do ir=1,mynr
                Call random_number(harvest=rand)
                yvar = pi*(r_ind(myr1+ir-1)-r1)/(r2-r1)/2d0
                Do it=1,mynth
                   xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                   vars1(it,:,ir,5) = (2d0*rand(it,:)-1d0)*(1d0-cos(xvar)**16)*(1d0-cos(yvar)**16)
                EndDo
             EndDo
          EndIf
       Else
          If (bot_inout) Then
             Do ir=1,mynr
                Call random_number(harvest=rand)
                yvar = pi*(r_ind(myr1+ir-1)-r1)/(r2-r1)/2d0
                Do it=1,mynth
                   xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                   vars1(it,:,ir,5) = (2d0*rand(it,:)-1d0)*(1d0-cos(xvar)**16)*(1d0-sin(yvar)**16)
                EndDo
             EndDo
          Else
             Do ir=1,mynr
                Call random_number(harvest=rand)
                yvar = pi*(r_ind(myr1+ir-1)-r1)/(r2-r1)
                Do it=1,mynth
                   xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                   vars1(it,:,ir,5) = (2d0*rand(it,:)-1d0)*(1d0-cos(xvar)**16)*(1d0-cos(yvar)**16)
                EndDo
             EndDo
          EndIf
       EndIf

       If (SLaplacian) Then
          Gconst = 6.67d-8
          av = (Gconst*mstar/(Cp*T0*r2))
          bv = 0.5d0*(omega_0*r2)**2/(Cp*T0)

          Do ir=1,mynr
             cvar = (rhor(myr1+ir-1)/rho0)**(gamm1) !1d0 + av*(r2*r_inv(myr1+ir-1)-1d0)
             Do ip=1,mynphi
                !rotation compenstated density.
                vars1(:,ip,ir,1) = rho0*(cvar + bv*((r_ind(myr1+ir-1)*sines(myth1:myth2)/r2)**2-1d0))**(1d0/gamm1)
                vars1(:,ip,ir,4) = vamp*vars1(:,ip,ir,5)
                vars1(:,ip,ir,5) = (1d0+samp*vars1(:,ip,ir,5))*sr(myr1+ir-1)
             EndDo
          EndDo

       Else
          If ((Unstratified).and.(.not.magnetic)) Then 
             !Initialize acoustic wave
             cvar = sqrt(g0)
             rhot = cos(pi*(theta(myth1:myth2)-0.5d0*(th1+th2))/(th2-th1))**42
             Do ir=1,mynr
                bv = samp*cos(pi*(r_ind(myr1+ir-1)-0.5d0*(r1+r2))/(r2-r1))**42
                Do ip=1,mynphi
                   vars1(:,ip,ir,1) = rho0*(1d0 + bv*rhot*cos(nmodes*pi*phi(myphi1+ip-1)/phi2))
                   vars1(:,ip,ir,3) = cvar*bv*rhot*cos(nmodes*pi*phi(myphi1+ip-1)/phi2)
                   vars1(:,ip,ir,5) = sr(1)
                EndDo
             EndDo
          ElseIf (Isentropic) Then
             Gconst = 6.67d-8
             
             av = (Gconst*mstar/(Cp*T0*r2))
             bv = 0.5d0*(omega_0*r2)**2/(Cp*T0)
             
             !Initialize prho and ps
             Do ir=1,mynr
                cvar = 1d0 + av*(r2*r_inv(myr1+ir-1)-1d0)
                Do ip=1,mynphi
                   !rotation compenstated density.
                   vars1(:,ip,ir,1) = rho0*(cvar + bv*((r_ind(myr1+ir-1)*sines(myth1:myth2)/r2)**2-1d0))**(1d0/gamm1)
                   vars1(:,ip,ir,2) = 0d0
                   vars1(:,ip,ir,3) = 0d0
                   vars1(:,ip,ir,4) = vamp*vars1(:,ip,ir,5)
                   vars1(:,ip,ir,5) = (1d0+samp*vars1(:,ip,ir,5))*sr(myr1+ir-1) + delta_s*((r_ind(myr1+ir-1)-r1)/(r2-r1))**10
                EndDo
             EndDo
             If (do_entropy_rain) delta_s = 0d0
          Else

             temparr = gravity             
             Do it=1,mynth
                gravity = temparr - r_ind*(sines(myth1+it-1)*omega_0)**2
                density(it,:) = Compute_Density(sr)
             EndDo
             gravity = temparr

             !Initialize prho and ps
             Do ir=1,mynr
                Do ip=1,mynphi
                   !rotation compenstated density.
                   vars1(:,ip,ir,1) = density(:,myr1+ir-1)
                   vars1(:,ip,ir,2) = 0d0
                   vars1(:,ip,ir,3) = 0d0
                   vars1(:,ip,ir,4) = vamp*vars1(:,ip,ir,5)
                   vars1(:,ip,ir,5) = (1d0+samp*vars1(:,ip,ir,5))*sr(myr1+ir-1) + delta_s*((r_ind(myr1+ir-1)-r1)/(r2-r1))**10
                EndDo
             EndDo
          EndIf
       EndIf
    EndIf

    Allocate(vphi_prof(nth))
    vphi_prof = 0d0
    If (solar_prof) Then
       !scaled to match equatorial value at 0.89 Rsun in inversion   
       av = 467.425d0
       bv = -51.25d0*cvar/455.425d0
       cvar = -81.525d0*cvar/455.425d0
       av = 2d0*pi*av/1d9
       bv = 2d0*pi*bv/1d9
       cvar = 2d0*pi*cvar/1d9

       vphi_prof = av+bv*cosines**2+cvar*cosines**4 !at 0.89 

       Do ir=1,mynr
          Do it=1,mynth             
             xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
             yvar = pi*nmodes*(phi(myphi1:myphi2)-phi1)/(phi2-phi1)
             vars1(it,:,ir,5) = sin(xvar)**2*sin(pi*(r_ind(myr1+ir-1)-r1)/(r2-r1))**2*sin(yvar)*sin(nmodes*xvar)
          EndDo
       EndDo
       
       !Initialize prho and ps
       Do ir=1,mynr
          Do ip=1,mynphi
             xvar = r_ind(myr1+ir-1)/r1
             omtmp = (vphi_prof(myth1:myth2)*xvar**2+omega_0*(1d0-xvar**2))**2 !Omega^2
             iomtmp = -0.5d0*(av**2*xvar**4+omega_0**2*(1d0+xvar**4))*cosines(myth1:myth2)**2 - &
                  & 0.5d0*bv*xvar**2*(av*xvar**2+omega_0*(1d0-xvar**2))*cosines(myth1:myth2)**4 - &
                  & xvar**2*((bv*xvar)**2+2d0*bv*((av-omega_0)*xvar**2+omega_0))*cosines(myth1:myth2)**6/6d0 - & 
                  & 0.25d0*bv*cvar*xvar**4*cosines(myth1:myth2)**8 - 0.1d0*cvar**2*xvar**4*cosines(myth1:myth2)**10 + &
                  & omega_0*xvar**2*(omega_0-av*(1d0-xvar**2))*cos(2d0*theta(myth1:myth2))  !Int(sint*cost*Omega^2 dt)
             diomtmp = -2d0*(av**2+omega_0**2)*xvar**3*cosines(myth1:myth2)**2-bv*xvar*(omega_0+2d0*(av-omega_0)*xvar**2)*cosines(myth1:myth2)**4 - &
                  & 2d0*xvar*(cvar*omega_0+xvar**2*(bv**2+(2d0*av-omega_0)*cvar))*cosines(myth1:myth2)**6/3d0 - &
                  & -bv*cvar*xvar**3*cosines(myth1:myth2)**8 - 0.4d0*cvar**2*xvar**3*cosines(myth1:myth2)**10 - &
                  & +omega_0*xvar*(av*(2d0*xvar**2-1d0)+omega_0)*cos(2d0*theta(myth1:myth2)) !d/dr Int(sint*cost*Omega^2 dt)
             diomtmp = diomtmp/r1
             
             vars1(:,ip,ir,1) = rhor(myr1+ir-1)*r_ind(myr1+ir-1)*(omtmp*sines(myth1:myth2)**2-r_ind(myr1+ir-1)*diomtmp-&
                  & (2d0-gravity(myr1+ir-1)*r_ind(myr1+ir-1)*dexp(-sr(1)/Cv)/rhor(myr1+ir-1)**gamm1/gamma)*iomtmp)/gravity(myr1+ir-1)
             vars1(:,ip,ir,5) = (1d0+samp*vars1(:,ip,ir,5))*sr(myr1+ir-1) + Cp*rhor(myr1+ir-1)*r_ind(myr1+ir-1)**2*iomtmp/pr(myr1+ir-1)/gamma-Cp*vars1(:,ip,ir,1)/rhor(myr1+ir-1)
             vars1(:,ip,ir,1) = vars1(:,ip,ir,1) + rhor(myr1+ir-1)
          EndDo
       EndDo

       vphi_prof = r1*(vphi_prof-omega_0)
    EndIf

    If (Unstratified) Then
      If ((myth1 .eq. 1).or.(myth2 .eq. nth)) Then
          rhot1 = rho0
          st1 = sr(1)
          pt1 = Pr(1)
          tt1 = T0
       EndIf
       
       If (myr1 .eq. 1) Then
          rhor1 = rho0
          sr1 = sr(1)
          pr1 = pr(1)
          Tr1 = T0
       EndIf
       
       If (myr2 .eq. nr) Then
          rhor2 = rho0
          sr2 = sr(1)
          pr2 = pr(1)
          Tr2 = T0
       EndIf

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
             rhot = vars1(it,:,ir,1)
             mtot = mtot + radius**2*sinth*Sum(rhot)*drad*dth*dphi !nops = 7*mynphi
             mrsint = mrsint + radius**3*sinth**2*Sum(rhot)*drad*dth*dphi
             wt = vars1(it,:,ir,2)**2 !2*mynphi
             vt = vars1(it,:,ir,3)
             ut = vars1(it,:,ir,4)**2 !2*mynphi
             vol = vol + sinth*radius**2*drad*dth*dphi
             Ltot = Ltot + radius**3*sinth*Sum(rhot*Sqrt(wt+sinth**2*(vt+radius*sinth*omega_0)**2))*drad*dth*dphi !8+29*mynphi
             vt = vt*vt !2*mynphi
             Etot = Etot + radius**2*sinth*Sum(rhot*(Cv*T0+0.5d0*(ut+wt+vt)))*drad*dth*dphi !6+6*mynphi
          EndDo !71*mynphi+21 (hydro)
       EndDo !Total ops = mynr*mynth*(71*mynphi+21)
       L0 = Ltot
       E0 = Etot
       m0 = mtot
    Else
       If ((myth1 .eq. 1).or.(myth2 .eq. nth)) Then
          rhot1 = vars1(1,1,:,1)
          st1 = vars1(1,1,:,5)
          pt1 = rhot1**gamma*dexp(st1/Cv)
          tt1 = pt1/rhot1/gamm1/Cv
       EndIf
       
       If (myr1 .eq. 1) Then
          rhor1 = vars1(:,1,1,1)
          sr1 = vars1(:,1,1,5)
          pr1 = rhor1**gamma*dexp(sr1/Cv)
          Tr1 =pr1/rhor1/gamm1/Cv
       EndIf
       
       If (myr2 .eq. nr) Then
          rhor2 = vars1(:,1,mynr,1)
          sr2 = vars1(:,1,mynr,5)
          pr2 = rhor2**gamma*dexp(sr2/Cv)
          Tr2 = pr2/rhor2/gamm1/Cv
       EndIf

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
             rhot = vars1(it,:,ir,1)
             mtot = mtot + radius**2*sinth*Sum(rhot)*drad*dth*dphi !nops = 7*mynphi
             mrsint = mrsint + radius**3*sinth**2*Sum(rhot)*drad*dth*dphi
             wt = vars1(it,:,ir,2)**2 !2*mynphi
             vt = vars1(it,:,ir,3)
             ut = vars1(it,:,ir,4)**2 !2*mynphi
             tt = (rhot**gamm1)*exp(vars1(it,:,ir,5)/Cv)/gamm1/Cv
             vol = vol + sinth*radius**2*drad*dth*dphi
             Ltot = Ltot + radius**3*sinth*Sum(rhot*Sqrt(wt+sinth**2*(vt+radius*sinth*omega_0)**2))*drad*dth*dphi !8+29*mynphi
             vt = vt*vt !2*mynphi
             Etot = Etot + radius**2*sinth*Sum(rhot*(Cv*tt+0.5d0*(ut+wt+vt)))*drad*dth*dphi !6+6*mynphi
          EndDo !71*mynphi+21 (hydro)
       EndDo !Total ops = mynr*mynth*(71*mynphi+21)
       L0 = Ltot
       E0 = Etot
       m0 = mtot
    EndIf
    sendb(1) = vol
    sendb(2) = m0
    sendb(3) = E0
    sendb(4) = L0
    sendb(5) = mrsint
    
    recvb = 0d0
    Call Global_AllReduce(sendb,recvb,node_comm,'sum')
    vol = recvb(1)
    m0 = recvb(2)
    E0 = recvb(3)
    L0 = recvb(4)
    mrsint = recvb(5)

    If (.Not. checkpoint_start) Then
       vars0 = vars1
    EndIf
    
    If (startbio .and. checkpoint_start .and. (myr1.eq.1)) Then
       vars0(:,:,1,4) = vamp*vars0(:,:,2,4)
       ir = 1
       sendb(1) = Sum(vars0(:,:,ir,4))
       recvb = 0d0
       Call Global_AllReduce(sendb,recvb,bot_comm,'sum')
       recvb = recvb/Dble(nth*nphi)
       vars0(:,:,ir,4) = vars0(:,:,ir,4) - recvb(1)
    EndIf
    
    If (starttio .and. checkpoint_start .and. (myr2.eq.nr)) Then
       ir = mynr
       vars0(:,:,ir,4) = vars0(:,:,ir-1,4)*vamp
       sendb(1) = Sum(vars0(:,:,ir,4))
       recvb = 0d0
       Call Global_AllReduce(sendb,recvb,top_comm,'sum')
       recvb = recvb/Dble(nth*nphi)
       vars0(:,:,ir,4) = vars0(:,:,ir,4) - recvb(1)
       vars0(:,:,ir,5) = vars0(:,:,ir-1,5)
    EndIf

    vars1 = vars0

    If (bot_inout) Then
       mtot = 0d0
       Ltot = 0d0
       mrsint = 0d0
       !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
          !#DIR$ VECTOR ALIGNED
          Do ip=1,mynphi
             mtot = mtot + Sum(vars1(:,ip,ir,1)*sines(myth1:myth2))*dr_ind(myr1+ir-1)*r_ind(myr1+ir-1)**2
             mrsint = mrsint + Sum(vars1(:,ip,ir,1)*sines(myth1:myth2)**3)*dr_ind(myr1+ir-1)*r_ind(myr1+ir-1)**3
             Ltot = Ltot + Sum(vars1(:,ip,ir,1)*vars1(:,ip,ir,3)*sines(myth1:myth2)**2)*dr_ind(myr1+ir-1)*r_ind(myr1+ir-1)**3
          EndDo
       EndDo
       sendb(1) = mtot*dth*dphi
       sendb(2) = Ltot*dth*dphi
       sendb(3) = mrsint*dth*dphi
       recvb = 0d0
       Call Global_AllReduce(sendb,recvb,node_comm,'sum')
       mtot = recvb(1)
       Ltot = recvb(2)
       mrsint = recvb(3)
    EndIf

    If (Bottom_Specified) Then
       ash_fen_mod = 0d0
       If (myr1 .eq. 1) Then
          Call Initialize_ASH_Inflow()
          Call ASH_Enthalpy_Flux()
       EndIf
       ash_fen_mod = 2d0*ash_fen_mod
       Call MPI_Bcast(ash_fen_mod,1,MPI_REAL8,0,node_comm,ierr)
    EndIf
  End Subroutine Initialize_Flow

  Subroutine Initialize_Mag(checkpoint_start)
    Implicit None
    Logical, Intent(In) :: checkpoint_start
    Real*8 :: sendb, recvb, Etot, xvar, yvar, arg, work(mynth,mynr), workrb1(mynth,1), workrb2(mynth,4), dwork(mynth,mynr), din(nr), dout(nr), drin
    Real*8 :: bjpoo6(nr), bjtmp(0:10), dbjtmp(0:10), bytmp(0:10), dbytmp(0:10), rand(mynth,mynphi,mynr), bc1, bc2
    Integer :: ir, it, ip, iv, seed(1)

    If ((checkpoint_start).and.(Magnetic_In)) Then
       vars1(:,:,:,6:8) = vars0(:,:,:,6:8)
    Else
       If (Initialize_MRI) Then
          vars1(:,:,:,6) = 0d0
          Do ir=1,mynr
             Do ip=1,mynphi
                vars1(:,ip,ir,7) = -0.5d0*bamp*r_ind(myr1+ir-1)/sines(myth1:myth2)
             EndDo
          EndDo
          vars1(:,:,:,8) = 0d0
       ElseIf (Unstratified) Then

          !Assumes r1=1d0 & r2 = R ~ 4.973(fifth zero of besselJ, for roughly seven peaks)
          sendb = 0d0
          bjtmp = 0d0
          dbjtmp = 0d0
          bytmp = 0d0
          dbytmp = 0d0
          recvb = 1d0/6d0
          arg = (r1**3)/3d0
          call jyv(recvb,arg,sendb,bjtmp,dbjtmp,bytmp,dbytmp)
          xvar = dcos(pi/6d0)*bjtmp(0)-dsin(pi/6d0)*bytmp(0)
          yvar = bjtmp(0)
          Etot = xvar/yvar
          xvar = dcos(pi/6d0)
          yvar = dsin(pi/6d0)

          Do ir=1,nr
             sendb = 0d0
             bjtmp = 0d0
             dbjtmp = 0d0
             bytmp = 0d0
             dbytmp = 0d0
             arg = (r_ind(ir)**3)/3d0
             call jyv(recvb,arg,sendb,bjtmp,dbjtmp,bytmp,dbytmp)
             bjpoo6(ir) = (xvar-Etot)*bjtmp(0)-yvar*bytmp(0)
          EndDo

          xvar = bamp*samp/sqrt(4d0*pi*rho0)
          !Velocity wave
          bjpoo6 = xvar*(r_ind**(1.5d0))*bjpoo6

          vars1(:,:,:,1) = rho0
          vars1(:,:,:,2:4) = 0d0
          vars1(:,:,:,5) = sr(1)
          vars1(:,:,:,6:8) = 0d0
          xvar = sqrt(4d0*pi*rho0)
          Do ir=1,mynr
             Do it=1,mynth
                !yvar = -xvar*(r_inv(myr1+ir-1)**2)/sines(myth1+it-1)
                Etot = -bamp*cosines(myth1+it-1)*r_inv(myr1+ir-1)*r2**2/sines(myth1+it-1)
                vars1(it,:,ir,3) = bjpoo6(myr1+ir-1)/sines(myth1+it-1)
                !vars1(it,:,ir,6) = yvar*bjpoo6(myr1+ir-1)
                vars1(it,:,ir,7) = Etot
             EndDo
          EndDo

          If (myrank .eq. 0) Print*,'bjpoo6=[',bjpoo6,']'

       Else

          Call system_clock(seed(1))
          Call random_seed(put=seed)
          Call random_number(harvest=rand)
          rand = bamp*r2*(2d0*rand-1d0)
          Do ir=1,mynr
             yvar = pi*(r_ind(myr1+ir-1)-r1)/(r2-r1)
             Do it=1,mynth
                xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                vars1(it,:,ir,6) = rand(it,:,ir)*(1d0-cos(xvar)**2)*(1d0-cos(yvar)**2)
             EndDo
          EndDo


          seed = seed + 1010101
          Call random_seed(put=seed)
          Call random_number(harvest=rand)
          rand = bamp*r2*(2d0*rand-1d0)
          Do ir=1,mynr
             yvar = pi*(r_ind(myr1+ir-1)-r1)/(r2-r1)
             Do it=1,mynth
                xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                vars1(it,:,ir,7) = rand(it,:,ir)*(1d0-cos(xvar)**2)*(1d0-cos(yvar)**2)
             EndDo
          EndDo

          seed = seed + 1010101
          Call random_seed(put=seed)
          Call random_number(harvest=rand)
          rand = bamp*r2*(2d0*rand-1d0)
          Do ir=1,mynr
             yvar = pi*(r_ind(myr1+ir-1)-r1)/(r2-r1)
             Do it=1,mynth
                xvar = pi*(theta(myth1+it-1)-th1)/(th2-th1)
                vars1(it,:,ir,8) = rand(it,:,ir)*(1d0-cos(xvar)**2)*(1d0-cos(yvar)**2)
             EndDo
          EndDo

          Do ir=1,mynr
             Do ip=1,mynphi
                vars1(:,ip,ir,8) = vars1(:,ip,ir,8) - bamp*r2*((4d0*(r_ind(myr1+ir-1)-r1)*(r2-r_ind(myr1+ir-1))/(r2-r1)**2)**6)* &
                     & ((4d0*(theta(myth1:myth2)-th1)*(th2-theta(myth1:myth2))/(th2-th1)**2)**11)*(th2-th1)/11d0
             EndDo
          EndDo
       EndIf
       
       vars0 = vars1
    EndIf

  End Subroutine Initialize_Mag

  !------------------------------------------------------------------------

  Subroutine Read_Input()
    Implicit None
    Integer :: input_az_quantities(max_qs), input_shell_avg_quantities(max_qs)
    Integer :: input_shell_slice_quantities(max_qs), input_shell_slice_radii(max_qs)
    Integer :: iq

    Namelist /ProblemSpec_Namelist/ niters, cfl_safety, nr, nth, nphi, nnr, nnt, nnp, r1, r2, &
         & th1, th2, phi1, phi2, non_uniform, x_1, x_2, perm_dir, input_checkpoint_filename, &
         & Start_From_Checkpoint, Magnetic_In, stellar_model_filename, samp, bamp, vamp, nmodes, &
         & skip, Validate, Debug, update_timestep, Do_Timings, RK3TVD, RK2TVD

    !Physics Flags
    Namelist /Physics_Namelist/ magnetic, do_radheat, do_entropy_rain, & 
         & solar_prof, Initialize_MRI, Isothermal, Unstratified, SLDiffusion, Laplacian, Isentropic, Thermal_Drag, SLaplacian, &
         !Physics Parameters
         & gamma, Cp, omega_0, kappa_value, g0, rho0, T0, thermal_drag_time, eps_exp, ell0_heating, mstar, lrad_update, lrad_radius, lrad_tanh_coef, &
         !Diffusion Parameters
         & slope_limiter, sld_coef, sld_speedr, sld_speedh, sld_exp, sld_pr, sld_prm, Mass_Diffusion, update_background, sld_tanh_coef, &
         & Magnetic_Diffusion, dstype, do_ohmic_heating, do_viscous_heating, flux_limiter, flux_cutoff, mag_flux_limiter, ent_flux_limiter, &
         !Background State Parameters
         & ufreq, k0expon, lstar, lumr1, lumr2, damping_time, delta_s, &
         !Entropy Rain Parameters
         & start_rain, nsqualls, meandur, meansize, rainstats, rhotopin, utopin, &
         & stopin, tstruct, sstruct, placetype, ncells, build_lanes, &
         !Modify Sound Speed
         & xi_cs

    !Top Boundary Conditions
    Namelist /Boundary_Namelist/ Top_Impenetrable, Top_Radial_Outflow, Top_Stress_Free, &
         & Top_No_Slip, Top_Constant_Temp, Top_Constant_Flux, Set_Stop, &
         !Bottom Boundary Conditions
         & Bottom_Impenetrable, Bottom_Stress_Free, Bottom_No_Slip, Bottom_Constant_Temp, Bottom_Constant_Flux, Bottom_Specified, &
         !Theta Boundary Conditions
         & Theta_Symmetric, Theta_Open, Theta_Stress_Free, &
         !Magnetic Boundary Conditions
         & Top_Radial_Field, Top_Perfect_Cond, Bottom_Radial_Field, Bottom_Perfect_Cond, Theta_Perf_Conductor, Theta_Normal_Field, &
         !Inflow Outflow
         & bot_inout, startbio, top_inout, starttio, ashnth, ashnphi

    Namelist /IO_Namelist/ scalar_freq, shell_avg_freq, shell_avg_recnum, shell_slice_freq, &
         & shell_slice_recnum, az_avg_freq, az_avg_recnum, checkpoint_freq, input_az_quantities, &
         & input_shell_avg_quantities, input_shell_slice_quantities, input_shell_slice_radii

    input_az_quantities = 0
    input_shell_avg_quantities = 0
    input_shell_slice_quantities = 0
    input_shell_slice_radii = 0

    !  read in namelists
    Open(namelist_unit,file='input',status='old')
    Read(namelist_unit,nml=ProblemSpec_Namelist)
    Read(namelist_unit,nml=Physics_Namelist)
    Read(namelist_unit,nml=Boundary_Namelist)
    Read(namelist_unit,nml=IO_Namelist)
    Close(namelist_unit)

    az_qs = 0
    shell_avg_qs = 0
    shell_slice_qs = 0
    shell_slice_rads = 0

    Do iq = 1, max_qs
       If (input_shell_avg_quantities(iq) .gt. 0) Then
          shell_avg_qs = shell_avg_qs + 1
       End If  
       If (input_az_quantities(iq) .gt. 0) Then
          az_qs = az_qs + 1
       End If
       If (input_shell_slice_quantities(iq) .gt. 0) Then
          shell_slice_qs = shell_slice_qs + 1
       End If
       If (input_shell_slice_radii(iq) .gt. 0) Then
          shell_slice_rads = shell_slice_rads + 1
       End If
    End Do

    If (az_qs .eq. 0) Then 
       az_qs=max_qs
    End If

    If (shell_avg_qs .eq. 0) Then 
       shell_avg_qs=max_qs
    End If
    
    If (shell_slice_qs .eq. 0) Then 
       shell_slice_qs=max_qs
    End If

    If (shell_slice_rads .eq. 0) Then
       shell_slice_rads=max_rads
    End If
    
    Allocate(az_quantities(az_qs),shell_avg_quantities(shell_avg_qs),shell_slice_quantities(shell_slice_qs),shell_slice_radii(shell_slice_rads))

    If (az_qs .lt. max_qs) Then
       Do iq=1,az_qs
          az_quantities(iq) = input_az_quantities(iq)
       End Do
    Else
       Do iq=1,az_qs
          az_quantities(iq) = iq
       End Do
    End If

    If (shell_avg_qs .lt. max_qs) Then
       Do iq=1,shell_avg_qs
          shell_avg_quantities(iq) = input_shell_avg_quantities(iq)
       End Do
    Else
       Do iq=1,shell_avg_qs
          shell_avg_quantities(iq) = iq
       End Do
    End If

    If (shell_slice_qs .lt. max_qs) Then
       Do iq=1,shell_slice_qs
          shell_slice_quantities(iq) = input_shell_slice_quantities(iq)
       End Do
    Else
       Do iq=1,shell_slice_qs
          shell_slice_quantities(iq) = iq
       End Do
    End If

    If (shell_slice_rads .lt. max_rads) Then
       Do iq=1,shell_slice_rads
          shell_slice_radii(iq) = input_shell_slice_radii(iq)
       End Do
    Else
       Do iq=1,shell_slice_rads
          shell_slice_radii(iq) = iq*nr/(shell_slice_rads+1)
       End Do
    End If

    Call Check_Quantities()

    If ((input_checkpoint_filename .eq. '0').and.(Start_From_Checkpoint)) Then
       !Read checkpoint.log
       Open(namelist_unit, file=Trim(Trim(perm_dir)//'checkpoint.log'), status='unknown', access='sequential', position='append')
       Backspace(namelist_unit)
       Read(namelist_unit, fmt='(i7)') iq
       Close(namelist_unit)
       istep = iq
       input_checkpoint_filename=''
       Call Make_Filename(input_checkpoint_filename,input_checkpoint_filename)
    EndIf

    If (magnetic) do_lorentz=.True.

    ! Constants
    gamm1 = gamma-1d0
    Cv = Cp/gamma

    !  convert input angles from degrees to radians    
    th1=th1*deg2rad
    th2=th2*deg2rad
    phi1=phi1*deg2rad
    phi2=phi2*deg2rad

    !  set number of variables
    If (magnetic) Then 
      nv=8
    Else 
      nv=5
    Endif

    !  check for conflicting boundary specifications
    If ((Top_No_Slip).Eqv.(Top_Stress_Free)) Then
      Print*,'  ERROR:  Choose one: Top_No_Slip/Top_Stress_Free'
      Stop
    Elseif ((Top_Constant_Temp).Eqv.(Top_Constant_Flux)) Then
      Print*,'  ERROR:  Choose one: Bottom_Constant_Temp/Bottom_Constant_Flux'
      Stop
    Elseif ((Bottom_No_Slip).Eqv.(Bottom_Stress_Free)) Then
      Print*,'  ERROR:  Choose one  Bottom_No_Slip/Bottom_Stress_Free'
      Stop
    Elseif ((Bottom_Constant_Temp).Eqv.(Bottom_Constant_Flux)) Then
      Print*,'  ERROR:  Choose one: Bottom_Constant_Temp/Bottom_Constant_Flux'
      Stop
    Endif

    !  set tangential velocity BC on radial boundaries
    If ((Top_No_Slip).And.(Bottom_No_Slip)) Then
       lbcr_v=1
       lbcr_w=1
    Elseif ((Top_No_Slip).And.(Bottom_Stress_Free)) Then
       lbcr_v=2
       lbcr_w=2
    Elseif ((Top_Stress_Free).And.(Bottom_No_Slip)) Then
       lbcr_v=3
       lbcr_w=3
    Elseif ((Top_Stress_Free).And.(Bottom_Stress_Free)) Then
       lbcr_v=4
       lbcr_w=4
    Else
      Print*,'  ERROR:  problem with tangential velocity BC on radial surfaces'
      Stop
    Endif

    If ((Top_Constant_Temp).And.(Bottom_Constant_Temp)) Then
      lbcr_s=1
    ElseIf ((Top_Constant_Temp).And.(Bottom_Constant_Flux)) Then
      lbcr_s=2
    ElseIf ((Top_Constant_Flux).And.(Bottom_Constant_Temp)) Then
      lbcr_s=3
    ElseIf ((Top_Constant_Flux).And.(Bottom_Constant_Flux)) Then
      lbcr_s=4
    Else
      Print*,'  ERROR:  problem with temperature BC on radial surfaces'
      Stop
    Endif

!    If (bot_inout) Then
!       If (top_inout) Then
!          lbcr_u = 0
!       Else
!          lbcr_u = 2
!       EndIf
!    Else
       If (top_inout) Then
          lbcr_u = 3
       Else
          lbcr_u = 1
       EndIf
!    EndIf

    lbcr_rho = 4
    lbcr_t   = 4 !lbcr_rho

    !  set theta temperature & density BC
    lbct_rho = 4
    lbct_t = 0 !lbct_s

    !Theta boundaries
    !Set normal velocity BC on latitudinal boundaries
    lbct_w = 1

    !Set tangential velocity BC on latitudinal boundaries
    If (Theta_Stress_Free) Then
       lbct_u = 0
    Else
       lbct_u = 1
    EndIf
    lbct_v = 4
    lbct_s = 0

    !  magnetic stuff: boundary conditions and initial field configuration
    If (magnetic) Then

       !  set radial magnetic field BC: assumes bth=bph=0 on radial boundary,
       !  which implies that dbr/dr=-2br/r on boundary (from div.B=0)     
       
       If ((Top_Radial_Field).And.(Bottom_Radial_Field)) Then
          If ((bot_inout) .and. (top_inout)) Then
             lbcr_br = 1
             lbcr_bt = 1
             lbcr_bp = 1

             lbcr_ar = 1
             lbcr_at = 1
             lbcr_ap = 1
          ElseIf ((bot_inout) .and. (.not. top_inout)) Then
             lbcr_br = 3
             lbcr_bt = 1
             lbcr_bp = 1

             lbcr_ar = 1
             lbcr_at = 3
             lbcr_ap = 3
          ElseIf ((.not. bot_inout) .and. (top_inout)) Then
             lbcr_br = 2
             lbcr_bt = 1
             lbcr_bp = 1

             lbcr_ar = 1
             lbcr_at = 2
             lbcr_ap = 2
          Else
             lbcr_br = 4
             lbcr_bt = 1
             If (Unstratified) Then
                lbcr_bp = 4
             Else
                lbcr_bp = 1
             EndIf

             lbcr_ar = 1
             If (Unstratified) Then
                lbcr_at = 1
             Else
                lbcr_at = 0
             EndIf
             lbcr_ap = 0
          EndIf
       Endif

       If ((Top_Perfect_Cond).and.(Bottom_Perfect_Cond)) Then
          lbcr_br = 1
          lbcr_bt = 1
          lbcr_bp = 1

          lbcr_ar = 1
          lbcr_at = 1
          lbcr_ap = 1
       EndIf

       If (Theta_Normal_Field) Then
          lbct_br = 1
          lbct_bt = 4
          lbct_bp = 1

          lbct_ar = 0
          lbct_at = 1
          lbct_ap = 0
       EndIf
       
       !  set latitudinal magnetic field BC
       If (Theta_Perf_Conductor) Then
          lbct_br = 1
          lbct_bt = 1
          lbct_bp = 1

          lbct_ar = 1
          lbct_at = 1
          lbct_ap = 1
       Endif
    Endif

  End Subroutine Read_Input

  Subroutine Check_Quantities()
     Integer :: iq, nqs

     nqs = Size(shell_avg_quantities)

     Do iq=1,nqs
        If (shell_avg_quantities(iq) .eq. 8) Then
           do_ke = .True.
        End If

        If (shell_avg_quantities(iq) .eq. 13) Then
           do_me = .True.
        End If

        If (shell_avg_quantities(iq) .eq. 14) Then
           do_bpol = .True.
        End If
     End Do

     nqs = Size(az_quantities)

     Do iq=1,nqs
        If (az_quantities(iq) .eq. 8) Then
           do_ke = .True.
        End If

        If (az_quantities(iq) .eq. 12) Then
           do_me = .True.
        End If

        If (az_quantities(iq) .eq. 13) Then
           do_bpol = .True.
        End If

        If ((az_quantities(iq) .eq. 18).or.(az_quantities(iq).eq.19).or.(az_quantities(iq) .eq. 24).or.(az_quantities(iq).eq.25)) Then
           Do_SLD_Output = .True.
        EndIf
     End Do
  End Subroutine Check_Quantities

End Program CSS2


  !------------------------------------------------------------------------

  Integer Function Sum_Max_Op(inv, inoutv, len, type)
    Implicit None
    Integer :: len, type, ii
    Real*8 :: inv(len), inoutv(len)
    Do ii=1,len-1
       inoutv(ii) = inoutv(ii) + inv(ii)
    EndDo
    inoutv(len) = max(inoutv(len),inv(len))
    Sum_Max_Op = 0
    Return
  End Function Sum_Max_Op
