Module Constants

  Implicit None

  !Basic Constants
  Real*8, Parameter :: Pi=3.1415926535897932384626433832795028841972d0, one_over_4pi=1d0/pi/4d0, deg2rad=pi/180d0
  Real*8, Parameter :: onethird=1d0/3d0, onesixth=1d0/6d0, a_const = 7.56577d-15, c_light = 2.99792458d10
  Real*8, Parameter :: twothirds=2d0/3d0, fourthirds=4d0/3d0, seventhirds=7d0/3d0

  !Problem size
  Integer :: nr, nth, nphi, nv

  !Number of subdomains in each direction
  Integer :: nnr, nnt, nnp 

  !Physical Parameters
  Real*8 :: gamma, gamm1, Cp, Cv, invCp, invCv, invpconst, pconst, invtconst, omega_0, sim_time, deltat, tiny, dtcur, mstar

  !Background state and its boundary values
  Real*8, Allocatable, Dimension(:) :: rhor1, rhor2, Tr1, Tr2, Pr1, Pr2, Sr1, Sr2, dsdr1, dsdr2, rhot1, pt1, tt1, st1, inv_drag_time
  Real*8 :: sm1, sp1, rhom1, rhop1, rho0, T0, g0, lumr1, lumr2, lstar, kappa_value
  Real*8 :: rhob1(4), rhob2(4)
  Real*8, Allocatable, Dimension(:) :: sr, rhor, tr, pr, drhodrm
  Real*8, Allocatable, Dimension(:) :: gravity, dsdrm, vphi_prof, Lrad, opac0, opac, dopdlnd, dopdlnT !Initial background state variables
  Real*8, Allocatable, Dimension(:) :: kr, dlnkr, ks, dlnks, mur, dlnmu, etar, eps, k0 !Initial background state variables

  !Initialization filenames
  Character(256) :: stellar_model_filename, input_checkpoint_filename

  !Intialization Flags
  Logical :: Magnetic_In=.False., Start_From_Checkpoint=.False., non_uniform=.False.

  !Initialization amplitudes for entropy (samp), magnetic field (bamp) and velocity (vamp)
  Real*8 :: samp, bamp, vamp, nmodes, ash_fen_mod

  !Integer Boundary Condition Flags for theta (lbct) and radius (lbcr) for each variable
  Integer :: lbct_u, lbcr_u, lbct_v, lbcr_v, lbct_w, lbcr_w, lbct_t, lbcr_t, lbct_s, lbcr_s, lbcr_rho, lbct_rho
  Integer :: lbct_br, lbcr_br, lbct_bt, lbcr_bt, lbct_bp, lbcr_bp, lbct_ar, lbcr_ar, lbct_at, lbcr_at, lbct_ap, lbcr_ap

  !  default values for boundary namelist variables
  Logical :: Top_Impenetrable       = .False.  ! / pick \
  Logical :: Top_Radial_Outflow     = .False.  ! \  one /
  Logical :: Top_Stress_Free        = .False.  ! / pick \
  Logical :: Top_No_Slip            = .False.  ! \  one /
  Logical :: Top_Constant_Temp      = .True.   ! / pick \
  Logical :: Top_Constant_Flux      = .True.   ! \  one /
  Logical :: Top_Radial_Field       = .True.
  Logical :: Top_Perfect_Cond       = .False.
  Logical :: Set_Stop               = .False.
  Logical :: Bottom_Impenetrable    = .False.  ! / pick \
  Logical :: Bottom_Stress_Free     = .False.  ! / pick \
  Logical :: Bottom_No_Slip         = .False.  ! \  one /
  Logical :: Bottom_Constant_Temp   = .False.  ! / pick \
  Logical :: Bottom_Constant_Flux   = .False.  ! \  one /
  Logical :: Bottom_Specified       = .False.
  Logical :: Bottom_Radial_Field    = .True.
  Logical :: Bottom_Perfect_Cond    = .False.
  Logical :: Theta_Perf_Conductor   = .False.
  Logical :: Theta_Normal_Field     = .True.
  Logical :: Theta_Symmetric        = .False.
  Logical :: Theta_Open             = .False.
  Logical :: Theta_Stress_Free      = .True.
  
  !Physics Controls
  Logical :: magnetic=.False., do_lorentz=.False., Initialize_MRI=.False., Thermal_Drag=.False.
  Logical :: solar_prof=.False., SLDiffusion=.True.,  do_radheat=.False., Laplacian=.False., SLaplacian=.False.
  Logical :: Isothermal=.False., Isentropic=.False., Unstratified=.False.

  !Inflow/Outflow and Entropy Rain Variables
  Real*8 :: rhotopin, utopin, stopin, T2avg, thermal_drag_time
  Logical :: top_inout=.False., bot_inout=.False., startbio=.False., starttio=.False.
  Logical :: do_entropy_rain=.False., start_rain=.False.


  !Input quantities specifying how many iterations between records and how many records per file
  Integer :: scalar_freq, az_avg_freq, shell_avg_freq, shell_slice_freq, checkpoint_freq, az_avg_recnum, shell_avg_recnum, shell_slice_recnum
  Logical :: Do_SLD_Output
  !Loop variables
  Integer :: istep, niters, begiter, skip=1, iashtime, ufreq

  !Timing & Debug Variables
  Logical :: Do_Timings=.False., Validate=.False., Debug=.False.
  Real*8 :: tstart, tstop
  Real*8, Allocatable, Dimension(:) :: timings
  
End Module Constants
