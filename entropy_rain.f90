!|==============================================================================|
!| Kyle Augustson @ JILA version 2010.1                                         |
!| This module contains all the routines to build the sub-grid scale model      |
!| that mimics the flow structures seen in surface convection simulations.      |
!| =============================================================================|
!| A rain squall has a lifetime of "Duration" during which it cools the         | 
!| surface layer of gas.  Alternatively, volumetric cooling can be activiated   |
!| that adds a radial extent to the cooling region. They can have either a      |
!| fixed horizontal extent or chosen from a random distribution. The center     |
!| point of the rain squall is advected by the horizontal velocities at the     |
!| surface. The mean duration determines how much energy from the flow is       |
!| extracted by the cooling region on average. This amount of energy is         |
!| matched to the radiative energy lost per unit surface area on a time average.|
!| =============================================================================|

Module Entropy_Rain

  Use Constants
  Use Derivatives
  Use Voronoi

  Implicit None

  Type Rain_Squall !Properties of the squall
     Integer :: loct, locp
     Real*8 :: Theta_loc, Phi_loc
     Real*8 :: Duration, SigmaT, Start_Time, Size, Dels, Scaling, Vel, Delr
  End Type Rain_Squall

  Type Voronoi_Cell
     Integer :: loct, locp, nvert !Location of centroid on grid
     Real*8 :: Theta_loc, Phi_loc !Location on sphere
     Real*8 :: Duration, Start_Time, Scaling
     Real*8, Allocatable, Dimension(:,:) :: vertices
  End Type Voronoi_Cell

  Logical, Public :: build_lanes
  Integer, Public :: ncells, nsqualls, tstruct, sstruct, placetype, rainstats, isz_max
  Real*8,  Public :: meandur, meansize, meands, meandrho, meanvel, sigt, sigds, sigv
  Real*8,  Public, Allocatable, Dimension(:,:) :: dent_rain, dvel_rain, squall_density
  Type(Rain_Squall), Public, Allocatable, Dimension(:) :: rain_squalls
  Type(Voronoi_Cell), Public, Allocatable, Dimension(:) :: conv_cells

Contains

  !Top Level Routines
  Subroutine Initialize_Rain_Squalls()
    Implicit None
    Integer :: is, ic
    Integer*4 :: v_num, g_deg(3*ncells), g_start(3*ncells), g_face(18*ncells), seed(1)
    Real*8 :: v_xy(2,6*ncells)
    Real*8, Allocatable, Dimension(:,:)  :: params, ptmp

    meandrho = rhotopin*rhor(nr)
    meands = stopin 
    meanvel = utopin
    sigt = 3d0*dr_ind(nr)/abs(utopin)
    isz_max = Int(Ceiling(meansize*sqrt(-log(1d-5))/(r2*Min(dth,dphi*Minval(sines)))))
    
    Call system_clock(seed(1))
    Call random_seed(put=seed)

    If (start_rain) Then
       Allocate(dent_rain(mynth,mynphi),dvel_rain(mynth,mynphi))

       dvel_rain = 0d0
       dent_rain = 0d0

       If (build_lanes) Then
          Allocate(params(3*ncells,3),ptmp(ncells,3),conv_cells(3*ncells))
          
          If (mytop_rank .eq. 0) Then
             Call Generate_New_Cell_Parameters(ptmp,ncells)
          EndIf
          
          Call MPI_Bcast(ptmp,ncells*3,MPI_REAL8,0,top_comm,ierr)
          
          params(1:ncells,:) = ptmp
          Deallocate(ptmp)
          
          !Ensure periodicity
          !Left
          params(ncells+1:2*ncells,1) = params(1:ncells,1)-(phi2-phi1)
          params(ncells+1:2*ncells,2) = params(1:ncells,2)
          params(ncells+1:2*ncells,3) = params(1:ncells,3)
          !Right
          params(2*ncells+1:3*ncells,1) = (phi2-phi1)+params(1:ncells,1)
          params(2*ncells+1:3*ncells,2) = params(1:ncells,2)
          params(2*ncells+1:3*ncells,3) = params(1:ncells,3)

          Do is = 1, 3*ncells
             conv_cells(is)%Phi_loc = params(is,1)
             conv_cells(is)%Theta_loc = params(is,2)
             conv_cells(is)%Duration = params(is,3)
             conv_cells(is)%Start_Time = sim_time-1d0
             conv_cells(is)%Scaling = 0d0
             Call Find_Cell_Position(is)
          EndDo
          
          Call Compute_Voronoi(params(:,1),params(:,2),g_deg,g_start,g_face,v_num,v_xy,ierr)
          
          Deallocate(params)
          Do is = 1, 3*ncells
             conv_cells(is)%nvert = g_deg(is)
             Allocate(conv_cells(is)%vertices(2,g_deg(is)))
             Do ic = 1, g_deg(is)
                conv_cells(is)%vertices(:,ic) = v_xy(:,g_face(g_start(is)+ic-1))
             EndDo
          EndDo
          Call Build_Downflow_Lanes()
       EndIf

       If (nsqualls .gt. 0) Then
          Allocate(params(nsqualls,7),rain_squalls(nsqualls),squall_density(mynth,mynphi))
          squall_density = 0d0
          If (mytop_rank .eq. 0) Then
             Call Generate_New_Squall_Parameters(params,nsqualls)
             
             Print*, 'Mean properties of squalls'
             Print*, 'nsqualls=', nsqualls, 'meansize=', Sum(params(:,2))/nsqualls, 'meandur=', Sum(params(:,3))/nsqualls
             Print*, 'meands=', Sum(params(:,1))/nsqualls, 'meanvel=',  Sum(params(:,6))/nsqualls, 'meandrho=', Sum(params(:,7))/nsqualls
             Print*, 'minmax ds', Minval(params(:,1)), Maxval(params(:,1))
             Print*, 'minmax size', Minval(params(:,2)), Maxval(params(:,2))
             Print*, 'minmax dur', Minval(params(:,3)), Maxval(params(:,3))
             Print*, 'minmax vel', Minval(params(:,6)), Maxval(params(:,6))
             Print*, 'minmax drho', Minval(params(:,7)), Maxval(params(:,7))
             Print*, 'isz_max=', isz_max
             Print*, 'minmax theta_loc', Minval(params(:,4)), Maxval(params(:,4))
             Print*, 'minmax phi_loc', Minval(params(:,5)), Maxval(params(:,5))
          EndIf

          Call MPI_Bcast(params,nsqualls*7,MPI_REAL8,0,top_comm,ierr)

          Do is=1,nsqualls
             !Generate initial storm
             rain_squalls(is)%Dels       = params(is,1)
             rain_squalls(is)%Size       = params(is,2)
             rain_squalls(is)%Duration   = params(is,3)
             rain_squalls(is)%Theta_loc  = params(is,4)
             rain_squalls(is)%Phi_loc    = params(is,5)
             rain_squalls(is)%Vel        = params(is,6)
             rain_squalls(is)%Delr       = params(is,7)
             rain_squalls(is)%SigmaT     = sigt
             rain_squalls(is)%Start_Time = sim_time-1d0
             rain_squalls(is)%Scaling    = 0d0
          
             !Find nearest grid point to center.
             Call Find_Position(is)
          EndDo

          Deallocate(params)
       EndIf
    EndIf
  End Subroutine Initialize_Rain_Squalls

  Subroutine Finalize_Rain_Squalls()
    Implicit None
    Integer :: ic
    If (myr2 .eq. nr) Then
       Deallocate(dent_rain,dvel_rain)
       If (nsqualls .gt. 0) Deallocate(rain_squalls,squall_density)
       If (build_lanes) Then
          Do ic=1,3*ncells
             Deallocate(conv_cells(ic)%vertices)
          EndDo
          Deallocate(conv_cells)
       EndIf
    EndIf
  End Subroutine Finalize_Rain_Squalls

  Subroutine Update_Rain_Squalls()
    Implicit None
    Integer :: is, ic, it, nup, lct, lcp
    Integer*4 :: v_num, g_deg(3*ncells), g_start(3*ncells), g_face(18*ncells)
    Real*8 :: delphi, v_xy(2,6*ncells)
    Real*8, Allocatable, Dimension(:) :: aloctmp, alocrcv, iup, vs, vt
    Real*8, Allocatable, Dimension(:,:) :: params

    dvel_rain = 0d0
    dent_rain = 0d0

    If (build_lanes) Then
       !Search for dead convective cells, elsewise update
       Allocate(aloctmp(5*ncells), alocrcv(5*ncells),iup(ncells),vs(2*nth),vt(2*nth))

       vs = 0d0
       vt = 0d0
       Do it=1,mynth
          vs(myth1+it-1) = Sum(vars0(it,:,:,2))
          vs(myth1+it-1+nth) = Sum(vars0(it,:,:,3))
       EndDo
       
       Call Global_AllReduce(vs,vt,top_comm,'sum')
       vs = deltat*vt/Dble(mynr*nphi)/r2
       Deallocate(vt)
       
       iup = 0
       nup = 0
       aloctmp = -1d0
       Do ic=1,ncells
          If ((conv_cells(ic)%Start_Time + conv_cells(ic)%Duration) .le. sim_time) Then !Generate new cells if necessary
             iup(ic) = 1
             nup = nup+1
          Else !Determine scaling and update position
             lct = conv_cells(ic)%loct
             lcp = conv_cells(ic)%locp
             Call Temporal_Structure_Cell(ic)
             If ((lct.ge.myth1).and.(lct.le.myth2).and.(lcp.ge.myphi1).and.(lcp.le.myphi2)) Then
                !Update Position, theta based on average subdomain velocity
                conv_cells(ic)%Theta_loc = conv_cells(ic)%Theta_loc + vs(lct)
                If (conv_cells(ic)%Theta_loc .le. th1) Then
                   conv_cells(ic)%Start_Time = 0d0 !Kill on next time step
                   nup = nup + 1
                   iup(ic) = 1
                EndIf
                If (conv_cells(ic)%Theta_loc .ge. th2) Then
                   conv_cells(ic)%Start_Time = 0d0 !Kill on next time step
                   nup = nup + 1
                   iup(ic) = 1
                EndIf
                
                !Update Position, phi based on average subdomain velocity
                conv_cells(ic)%Phi_loc = conv_cells(ic)%Phi_loc + vs(lct+nth)
                If (conv_cells(ic)%Phi_loc .le. phi1) Then
                   delphi = conv_cells(ic)%Phi_loc-phi1
                   conv_cells(ic)%Phi_loc = phi2+delphi  !Wrap around
                EndIf
                If (conv_cells(ic)%Phi_loc .ge. phi2) Then
                   delphi = conv_cells(ic)%Phi_loc-phi2
                   conv_cells(ic)%Phi_loc = phi1+delphi !Wrap around
                EndIf
                Call Find_Cell_Position(ic)
             EndIf
          EndIf
          aloctmp(ic) = Dble(conv_cells(ic)%loct)
          aloctmp(ic+ncells) = Dble(conv_cells(ic)%locp)
          aloctmp(ic+2*ncells) = conv_cells(ic)%Theta_loc
          aloctmp(ic+3*ncells) = conv_cells(ic)%Phi_loc
          aloctmp(ic+4*ncells) = Dble(iup(ic))
       EndDo
       
       !Reduce updated positions
       alocrcv=-1d2
       Call Global_AllReduce(aloctmp,alocrcv,top_comm,'max')
       
       Do ic=1,ncells
          conv_cells(ic)%loct = Int(alocrcv(ic))
          conv_cells(ic)%locp = Int(alocrcv(ic+ncells))
          conv_cells(ic)%Theta_loc = alocrcv(ic+2*ncells)
          conv_cells(ic)%Phi_loc = alocrcv(ic+3*ncells)
          iup(ic) = Int(alocrcv(ic+4*ncells))
       EndDo
       Deallocate(aloctmp,alocrcv,vs)
       
       nup = Sum(iup)
       
       !If necessary create new cells
       If (nup .gt. 0) Then
          Allocate(params(nup,3))
          params = -1d0
          If (mytop_rank .eq. 0) Then
             Print*, 'Generating ', nup, 'new cells.'
             Call Generate_New_Cell_Parameters(params,nup)
          EndIf
          
          Call MPI_Bcast(params,nup*3,MPI_REAL8,0,top_comm,ierr)
          
          ic = 1
          Do is=1,ncells
             If ((iup(is) .eq. 1).and.(ic.le.nup)) Then
                conv_cells(is)%Phi_loc    = params(ic,1)
                conv_cells(is)%Theta_loc  = params(ic,2)
                conv_cells(is)%Duration   = params(ic,3)
                conv_cells(is)%Start_Time = sim_time
                conv_cells(is)%Scaling = 0d0
                Call Find_Cell_Position(is)
                ic = ic+1
             EndIf
          EndDo
          Deallocate(params)
       EndIf

       Allocate(params(3*ncells,2))
       Do is = 1, ncells
          params(is,1) = conv_cells(is)%Phi_loc
          params(is,2) = conv_cells(is)%Theta_loc
          conv_cells(is+ncells)%Scaling = conv_cells(is)%Scaling
          conv_cells(is+2*ncells)%Scaling = conv_cells(is)%Scaling
       EndDo
       
       !Ensure periodicity
       !Left
       params(ncells+1:2*ncells,1) = params(1:ncells,1)-(phi2-phi1)
       params(ncells+1:2*ncells,2) = params(1:ncells,2)
       
       !Right
       params(2*ncells+1:3*ncells,1) = (phi2-phi1)+params(1:ncells,1)
       params(2*ncells+1:3*ncells,2) = params(1:ncells,2)

       ierr = 0
       Call Compute_Voronoi(params(:,1),params(:,2),g_deg,g_start,g_face,v_num,v_xy,ierr)
       If (ierr .ne. 0) Print*, 'ierr=', ierr
       
       Deallocate(params,iup)
       
       Do is = 1, 3*ncells
          Deallocate(conv_cells(is)%vertices)
          conv_cells(is)%nvert = g_deg(is)
          Allocate(conv_cells(is)%vertices(2,g_deg(is)))
          Do ic = 1, g_deg(is)
             conv_cells(is)%vertices(:,ic) = v_xy(:,g_face(g_start(is)+ic-1))
          EndDo
       EndDo

       Call Build_Downflow_Lanes()
    EndIf

    If (nsqualls .gt. 0) Then
       !Search for dead rain squalls, elsewise update
       Allocate(aloctmp(5*nsqualls), alocrcv(5*nsqualls),iup(nsqualls))
       iup = 0
       nup = 0

       aloctmp = -1d0
       Do is=1,nsqualls
          If ((rain_squalls(is)%Start_Time + rain_squalls(is)%Duration) .le. sim_time) Then        !Generate new squalls if necessary
             iup(is) = 1
             nup = nup+1
          Else !Determine scaling and update position
             lct = rain_squalls(is)%loct
             lcp = rain_squalls(is)%locp

             Call Temporal_Structure(is)
             If ((lct.ge.myth1).and.(lct.le.myth2).and.(lcp.ge.myphi1).and.(lcp.le.myphi2)) Then
                !Update Position, theta
                rain_squalls(is)%Theta_loc = rain_squalls(is)%Theta_loc + deltat*vars0(lct-myth1+1,lcp-myphi1+1,mynr,2)/r2
                If (rain_squalls(is)%Theta_loc .le. th1) Then
                   rain_squalls(is)%Start_Time = 0d0 !Kill on next time step
                   nup = nup + 1
                   iup(is) = 1
                EndIf
                If (rain_squalls(is)%Theta_loc .ge. th2) Then
                   rain_squalls(is)%Start_Time = 0d0 !Kill on next time step
                   nup = nup + 1
                   iup(is) = 1
                EndIf

                !Update Position, phi
                rain_squalls(is)%Phi_loc = rain_squalls(is)%Phi_loc + deltat*vars0(lct-myth1+1,lcp-myphi1+1,mynr,3)/r2
                If (rain_squalls(is)%Phi_loc .le. phi1) Then
                   delphi = rain_squalls(is)%Phi_loc-phi1
                   rain_squalls(is)%Phi_loc = phi2+delphi  !Wrap around
                EndIf
                If (rain_squalls(is)%Phi_loc .ge. phi2) Then
                   delphi = rain_squalls(is)%Phi_loc-phi2
                   rain_squalls(is)%Phi_loc = phi1+delphi !Wrap around
                EndIf
                Call Find_Position(is)
             EndIf
          EndIf
          aloctmp(is) = Dble(rain_squalls(is)%loct)
          aloctmp(is+nsqualls) = Dble(rain_squalls(is)%locp)
          aloctmp(is+2*nsqualls) = rain_squalls(is)%Theta_loc
          aloctmp(is+3*nsqualls) = rain_squalls(is)%Phi_loc
          aloctmp(is+4*nsqualls) = Dble(iup(is))
       EndDo

       !Reduce updated positions
       alocrcv=-1d2
       Call Global_AllReduce(aloctmp,alocrcv,top_comm,'max')

       Do is=1,nsqualls
          rain_squalls(is)%loct = Int(alocrcv(is))
          rain_squalls(is)%locp = Int(alocrcv(is+nsqualls))
          rain_squalls(is)%Theta_loc = alocrcv(is+2*nsqualls)
          rain_squalls(is)%Phi_loc = alocrcv(is+3*nsqualls)
          iup(is) = Int(alocrcv(is+4*nsqualls))
       EndDo
       Deallocate(aloctmp,alocrcv)

       nup = Sum(iup)
       squall_density = 0d0
       If (sstruct .eq. 1) Then
          Call Build_Gaussian_Squalls()
       ElseIf (sstruct .eq. 0) Then
          Call Build_Polynomial_Squalls()
       EndIf
       
       !If necessary create new squalls
       If (nup .gt. 0) Then          
          Allocate(params(nup,7))
          If (mytop_rank .eq. 0) Then
             params = -1d0
             Print*, 'Generating ', nup, 'new squalls.'
             Call Generate_New_Squall_Parameters(params,nup)
          EndIf

          Call MPI_Bcast(params,nup*7,MPI_REAL8,0,top_comm,ierr)

          ic = 1
          Do is=1,nsqualls
             If (iup(is) .eq. 1) Then
                rain_squalls(is)%Dels       = params(ic,1)
                rain_squalls(is)%Size       = params(ic,2)
                rain_squalls(is)%Duration   = params(ic,3)
                rain_squalls(is)%Theta_loc  = params(ic,4)
                rain_squalls(is)%Phi_loc    = params(ic,5)
                rain_squalls(is)%Vel        = params(ic,6)
                rain_squalls(is)%Delr       = params(ic,7)
                rain_squalls(is)%SigmaT     = sigt
                rain_squalls(is)%Start_Time = sim_time
                rain_squalls(is)%Scaling    = 0d0
                Call Find_Position(is)
                ic = ic+1
             EndIf
          EndDo

          Deallocate(params,iup)
       EndIf
    EndIf
  End Subroutine Update_Rain_Squalls

  Subroutine Build_Gaussian_Squalls()
    Implicit None
    Integer :: is, it, ip, lct, lcp
    Real*8 :: tmpval, thtmp(mynth), exptmp(mynth)
    Do is=1,nsqualls
       !Evaluate at new position
       lct = rain_squalls(is)%loct
       lcp = rain_squalls(is)%locp
       If (((lct+mynth).ge.myth1).and.((lct-mynth).le.myth2)) Then
          If (((lcp+mynphi).ge.myphi1).and.((lcp-mynphi).le.myphi2)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = -(thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).gt.(-11.5d0)) Then !log(1d-5)
                      tmpval = rain_squalls(is)%Scaling*exp(exptmp(it))
                      squall_density(it,ip) = squall_density(it,ip) + tmpval
                   EndIf
                EndDo
             EndDo
          EndIf

          !Periodicity
          If ((lcp + mynphi .gt. nphi).and.(myphi1 .eq. 1)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = -(thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc+phi2-phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).gt.(-11.5d0)) Then !log(1d-5)
                      tmpval = rain_squalls(is)%Scaling*exp(exptmp(it))
                      squall_density(it,ip) = squall_density(it,ip) + tmpval
                   EndIf
                EndDo
             EndDo
          EndIf

          If ((lcp - mynphi .lt. 1).and.(myphi2 .eq. nphi)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = -(thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc-phi2+phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).gt.(-11.5d0)) Then !log(1d-5)
                      tmpval = rain_squalls(is)%Scaling*exp(exptmp(it))
                      squall_density(it,ip) = squall_density(it,ip) + tmpval
                   EndIf
                EndDo
             EndDo
          EndIf
       EndIf
    EndDo
    squall_density = 1d0/(squall_density + 1d0)

    Do is=1,nsqualls
       !Evaluate at new position
       lct = rain_squalls(is)%loct
       lcp = rain_squalls(is)%locp
       If (((lct+mynth).ge.myth1).and.((lct-mynth).le.myth2)) Then
          If (((lcp+mynphi).ge.myphi1).and.((lcp-mynphi).le.myphi2)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = -(thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).gt.(-11.5d0)) Then !log(1d-5)
                      tmpval = rain_squalls(is)%Scaling*exp(exptmp(it))*squall_density(it,ip)
                      dvel_rain(it,ip) = dvel_rain(it,ip) + rain_squalls(is)%Vel*tmpval
                      dent_rain(it,ip) = dent_rain(it,ip) + rain_squalls(is)%Dels*tmpval
                   EndIf
                EndDo
             EndDo
          EndIf
          !Periodicity
          If ((lcp + mynphi .gt. nphi).and.(myphi1 .eq. 1)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = -(thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc+phi2-phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).gt.(-11.5d0)) Then !log(1d-5)
                      tmpval = rain_squalls(is)%Scaling*exp(exptmp(it))*squall_density(it,ip)
                      dvel_rain(it,ip) = dvel_rain(it,ip) + rain_squalls(is)%Vel*tmpval
                      dent_rain(it,ip) = dent_rain(it,ip) + rain_squalls(is)%Dels*tmpval
                   EndIf
                EndDo
             EndDo
          EndIf

          If ((lcp - mynphi .lt. 1).and.(myphi2 .eq. nphi)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = -(thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc-phi2+phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).gt.(-11.5d0)) Then !log(1d-5)
                      tmpval = rain_squalls(is)%Scaling*exp(exptmp(it))*squall_density(it,ip)
                      dvel_rain(it,ip) = dvel_rain(it,ip) + rain_squalls(is)%Vel*tmpval
                      dent_rain(it,ip) = dent_rain(it,ip) + rain_squalls(is)%Dels*tmpval
                   EndIf
                EndDo
             EndDo
          EndIf
       EndIf
    EndDo
  End Subroutine Build_Gaussian_Squalls

  Subroutine Build_Polynomial_Squalls()
    Implicit None
    Integer :: is, it, ip, lct, lcp
    Real*8 :: tmpval, xpow(3), thtmp(mynth), exptmp(mynth)
    Do is=1,nsqualls
       !Evaluate at new position
       lct = rain_squalls(is)%loct
       lcp = rain_squalls(is)%locp
       If (((lct+mynth).ge.myth1).and.((lct-mynth).le.myth2)) Then
          If (((lcp+mynphi).ge.myphi1).and.((lcp-mynphi).le.myphi2)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = (thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it) .le. 1d0) Then
                      tmpval = sqrt(exptmp(it))
                      xpow(1) = tmpval**2
                      xpow(2) = xpow(1)*tmpval
                      xpow(3) = xpow(2)*tmpval
                      tmpval = 1d0-18d0*xpow(1)+32d0*xpow(2)-15d0*xpow(3)
                      tmpval = rain_squalls(is)%Scaling*tmpval
                      squall_density(it,ip) = squall_density(it,ip) + abs(tmpval)
                   EndIf
                EndDo
             EndDo
          EndIf

          !Periodicity
          If ((lcp + mynphi .gt. nphi).and.(myphi1 .eq. 1)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = (thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc+phi2-phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it) .le. 1d0) Then
                      tmpval = sqrt(exptmp(it))
                      xpow(1) = tmpval**2
                      xpow(2) = xpow(1)*tmpval
                      xpow(3) = xpow(2)*tmpval
                      tmpval = 1d0-18d0*xpow(1)+32d0*xpow(2)-15d0*xpow(3)
                      tmpval = rain_squalls(is)%Scaling*tmpval
                      squall_density(it,ip) = squall_density(it,ip) + abs(tmpval)
                   EndIf
                EndDo
             EndDo
          EndIf

          If ((lcp - mynphi .lt. 1).and.(myphi2 .eq. nphi)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = (thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc-phi2+phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it) .le. 1d0) Then
                      tmpval = sqrt(exptmp(it))
                      xpow(1) = tmpval**2
                      xpow(2) = xpow(1)*tmpval
                      xpow(3) = xpow(2)*tmpval
                      tmpval = 1d0-18d0*xpow(1)+32d0*xpow(2)-15d0*xpow(3)
                      tmpval = rain_squalls(is)%Scaling*tmpval
                      squall_density(it,ip) = squall_density(it,ip) + abs(tmpval)
                   EndIf
                EndDo
             EndDo
          EndIf
       EndIf
    EndDo
    squall_density = 1d0/(squall_density + 1d0)

    Do is=1,nsqualls
       !Evaluate at new position
       lct = rain_squalls(is)%loct
       lcp = rain_squalls(is)%locp
       If (((lct+mynth).ge.myth1).and.((lct-mynth).le.myth2)) Then
          If (((lcp+mynphi).ge.myphi1).and.((lcp-mynphi).le.myphi2)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = (thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).le. 1d0) Then
                      tmpval = sqrt(exptmp(it))
                      xpow(1) = tmpval**2
                      xpow(2) = xpow(1)*tmpval
                      xpow(3) = xpow(2)*tmpval
                      tmpval = 1d0-18d0*xpow(1)+32d0*xpow(2)-15d0*xpow(3)
                      tmpval = rain_squalls(is)%Scaling*squall_density(it,ip)*tmpval
                      dvel_rain(it,ip) = dvel_rain(it,ip) + rain_squalls(is)%Vel*tmpval
                      dent_rain(it,ip) = dent_rain(it,ip) + rain_squalls(is)%Dels*tmpval
                   EndIf
                EndDo
             EndDo
          EndIf
          !Periodicity
          If ((lcp + mynphi .gt. nphi).and.(myphi1 .eq. 1)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = (thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc+phi2-phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).le. 1d0) Then
                      tmpval = sqrt(exptmp(it))
                      xpow(1) = tmpval**2
                      xpow(2) = xpow(1)*tmpval
                      xpow(3) = xpow(2)*tmpval
                      tmpval = 1d0-18d0*xpow(1)+32d0*xpow(2)-15d0*xpow(3)
                      tmpval = rain_squalls(is)%Scaling*squall_density(it,ip)*tmpval
                      dvel_rain(it,ip) = dvel_rain(it,ip) + rain_squalls(is)%Vel*tmpval
                      dent_rain(it,ip) = dent_rain(it,ip) + rain_squalls(is)%Dels*tmpval
                   EndIf
                EndDo
             EndDo
          EndIf

          If ((lcp - mynphi .lt. 1).and.(myphi2 .eq. nphi)) Then
             thtmp = (theta(myth1:myth2)-rain_squalls(is)%Theta_loc)**2
             Do ip=1,mynphi
                exptmp = (thtmp+(phi(myphi1+ip-1)-rain_squalls(is)%Phi_loc-phi2+phi1)**2)*(r2/rain_squalls(is)%Size)**2
                Do it=1,mynth
                   If (exptmp(it).le. 1d0) Then
                      tmpval = sqrt(exptmp(it))
                      xpow(1) = tmpval**2
                      xpow(2) = xpow(1)*tmpval
                      xpow(3) = xpow(2)*tmpval
                      tmpval = 1d0-18d0*xpow(1)+32d0*xpow(2)-15d0*xpow(3)
                      tmpval = rain_squalls(is)%Scaling*squall_density(it,ip)*tmpval
                      dvel_rain(it,ip) = dvel_rain(it,ip) + rain_squalls(is)%Vel*tmpval
                      dent_rain(it,ip) = dent_rain(it,ip) + rain_squalls(is)%Dels*tmpval
                   EndIf
                EndDo
             EndDo
          EndIf
       EndIf
    EndDo
  End Subroutine Build_Polynomial_Squalls

  Subroutine Build_Downflow_Lanes()
    Implicit None
    Integer :: ic, ip, it, ie, my, nmax
    Real*8 :: dvx, dvy, x1, x2, y1, y2, tmp, xv, yv, xsz, ysz, rtmp, thtmp(mynth)
    !Make lanes along Voronoi edges
    nmax = Int(Sqrt(Dble(nth**2+nphi**2)))
    !xsz = lane_size*dphi
    !ysz = lane_size*dth
    xsz = isz_max*dphi
    ysz = isz_max*dth
    thtmp = theta(myth1:myth2)
    rtmp = (r2/meansize)**2
    Do ic=1,3*ncells
       Do ie=1,conv_cells(ic)%nvert
          If (ie .ne. conv_cells(ic)%nvert) Then
             x1 = conv_cells(ic)%vertices(1,ie)
             x2 = conv_cells(ic)%vertices(1,ie+1)

             y1 = conv_cells(ic)%vertices(2,ie)
             y2 = conv_cells(ic)%vertices(2,ie+1)
          Else
             x1 = conv_cells(ic)%vertices(1,ie)
             x2 = conv_cells(ic)%vertices(1,1)

             y1 = conv_cells(ic)%vertices(2,ie)
             y2 = conv_cells(ic)%vertices(2,1)
          EndIf

          dvx = x2-x1
          dvy = y2-y1

          x2 = dvx/dphi
          y2 = dvy/dth
          
          x2 = x2**2+y2**2
          my = Min(Int(Sqrt(x2)),nmax)

          Do it=0,my-1
             tmp = Dble(it)/Dble(my)
             xv = x1 + dvx*tmp
             yv = y1 + dvy*tmp
             If (((xv.ge.phi(myphi1)-xsz).and.(xv.le.phi(myphi2)+xsz)) .and. &
                  & ((yv.ge.theta(myth1)-ysz).and.(yv.le.theta(myth2)+ysz))) Then
                Do ip=1,mynphi
                   dvel_rain(:,ip) = dvel_rain(:,ip) + conv_cells(ic)%Scaling*exp(-rtmp*((thtmp-yv)**2+(phi(myphi1+ip-1)-xv)**2)) !
                EndDo
             EndIf
          EndDo
       EndDo
    EndDo

    x1 = Maxval(dvel_rain)
    Call Global_AllReduce(x1,x2,top_comm,'max')
    dvel_rain = dvel_rain/Max(x2,1d0)
    dent_rain = meands*dvel_rain
    dvel_rain = meanvel*dvel_rain
  End Subroutine Build_Downflow_Lanes

  Subroutine Generate_New_Cell_Parameters(params,nsqs)
    Implicit None
    Integer, Intent(In) :: nsqs
    Real*8, Intent(InOut), Dimension(:,:) :: params
    Real*8, Allocatable, Dimension(:) :: rand
    Real*8 :: dur, ph0, th0
    Integer :: is

    Allocate(rand(nsqs))
    Call Random_Exponential(nsqs,rand,0.3d0,4d0)
    ! Loop over squalls
    Do is=1, nsqs
       dur = 1d12 !3600d0*6d0*rand(is)
       Call Place_Cell(th0,ph0)
       
       params(is,1) = ph0
       params(is,2) = th0
       params(is,3) = dur
    EndDo
    Deallocate(rand)
  End Subroutine Generate_New_Cell_Parameters

  Subroutine Generate_New_Squall_Parameters(params,nsqs)
    Implicit None
    Integer, Intent(In) :: nsqs
    Real*8, Intent(InOut), Dimension(:,:) :: params
    Real*8 :: dels, dur, sz, ph0, th0, vel, delr
    Integer :: is
    Real*8, Allocatable, Dimension(:) :: rand

    Select Case (rainstats)
    Case(1)
       Allocate(rand(nsqs))
       Call Random_Exponential(nsqs,rand,0.3d0,4d0)
       ! Loop over squalls
       Do is=1, nsqs
          !Duration and size are exponentially distributed
          dur = meandur*rand(is)
          sz  = meansize*rand(is)
          vel = meanvel*rand(is)
          dels = meands*rand(is)
          delr = meandrho*rand(is)
          
          Call Place_Squall(th0,ph0)
          
          params(is,1) = dels
          params(is,2) = sz
          params(is,3) = dur
          params(is,4) = th0
          params(is,5) = ph0
          params(is,6) = vel
          params(is,7) = delr
       EndDo
       Deallocate(rand)
    Case (0)
       Allocate(rand(nsqs))
       Call Random_Exponential(nsqs,rand,0.3d0,4d0)
       Do is=1, nsqs
          !Duration and size are exponentially distributed
          dur = meandur*rand(is)
          sz  = meansize
          vel = meanvel
          dels = meands
          delr = meandrho
          
          Call Place_Squall(th0,ph0)
          
          params(is,1) = dels
          params(is,2) = sz
          params(is,3) = dur
          params(is,4) = th0
          params(is,5) = ph0
          params(is,6) = vel
          params(is,7) = delr
       EndDo
       Deallocate(rand)
    Case Default
       Allocate(rand(nsqs))
       Call Random_Exponential(nsqs,rand,0.3d0,4d0)
       Do is=1, nsqs
          !Duration and size are exponentially distributed
          dur = meandur*rand(is)
          sz  = meansize
          vel = meanvel
          dels = meands
          delr = meandrho
          
          Call Place_Squall(th0,ph0)
          
          params(is,1) = dels
          params(is,2) = sz
          params(is,3) = dur
          params(is,4) = th0
          params(is,5) = ph0
          params(is,6) = vel
          params(is,7) = delr
       EndDo
       Deallocate(rand)
    EndSelect
  End Subroutine Generate_New_Squall_Parameters

  !Mid Level Routines
  Subroutine Find_Position(isq)
    Implicit None
    Integer,Intent(In) :: isq
    
    rain_squalls(isq)%loct = Max(Min(Floor((rain_squalls(isq)%Theta_loc-th1)/dth)+1,nth),1)
    rain_squalls(isq)%locp = Max(Min(Floor((rain_squalls(isq)%Phi_loc-phi1)/dphi)+1,nphi),1)

  End Subroutine Find_Position

  !Mid Level Routines
  Subroutine Find_Cell_Position(isq)
    Implicit None
    Integer,Intent(In) :: isq
    
    conv_cells(isq)%loct = Max(Min(Floor((conv_cells(isq)%Theta_loc-th1)/dth)+1,nth),1)
    conv_cells(isq)%locp = Max(Min(Floor((conv_cells(isq)%Phi_loc-phi1)/dphi)+1,nphi),1)

  End Subroutine Find_Cell_Position

  Subroutine Temporal_Structure_Cell(isq)
    Implicit None
    Integer, Intent(In) :: isq
    Real*8 :: t, t0, ft

    t = sim_time
    t0 = conv_cells(isq)%Start_Time
    If ((t-t0).le.1d4) Then
       ft = sin(pi*(t-t0)/2d4)**2
    Else
       ft = 1d0
    EndIf
    conv_cells(isq)%Scaling = ft
  End Subroutine Temporal_Structure_Cell

  Subroutine Temporal_Structure(isq)
    Implicit None
    Integer, Intent(In) :: isq
    Real*8 :: t, t0, Dur, SigmaT, ft, t1

    t = sim_time
    t0 = rain_squalls(isq)%Start_Time
    Dur = rain_squalls(isq)%Duration
    SigmaT = rain_squalls(isq)%SigmaT

    t1 = t0+Dur

    Select Case(tstruct)
    Case (1) !Ramp up and cut off abruptly
       If ((t .gt. t0) .and. (t .lt. t1)) Then
          ft = 1d0
       Else
          ft = 0d0
       EndIf
    Case Default !Gradually taper down
       If ((t-t0).le.SigmaT) Then
          ft = sin(0.5d0*pi*(t-t0)/SigmaT)**2
       ElseIf ((t-t0).ge.(Dur-SigmaT)) Then
          ft = sin(0.5d0*pi*(Dur-t-t0)/SigmaT)**2
       Else
          ft = 1d0
       EndIf
    End Select
    rain_squalls(isq)%Scaling = ft
  End Subroutine Temporal_Structure

  !Bottom Level Routines
  Subroutine Place_Squall(th0,ph0)
    Implicit None
    Real*8, Intent(InOut) :: th0, ph0
    Call Place_Cell(th0,ph0)
!!$    Real*8 :: z, z2
!!$    Integer :: ic, ie
!!$    !Randomly choose a cell, uniform
!!$    Call random_number(harvest=z)
!!$    z = 2d0*z-1d0
!!$    If (z .gt. 0d0) Then
!!$       Call random_number(harvest=z2)
!!$       ic = Ceiling(z2*ncells)
!!$    Else
!!$       Call random_number(harvest=z2)
!!$       ic = Max(Floor(z2*ncells),1)
!!$    EndIf
!!$    
!!$    !Randomly choose an edge of the cell, uniform
!!$    Call random_number(harvest=z)
!!$    z = 2d0*z-1d0
!!$    If (z .gt. 0d0) Then
!!$       Call random_number(harvest=z2)
!!$       ie = Ceiling(z2*conv_cells(ic)%nvert)
!!$    Else
!!$       Call random_number(harvest=z2)
!!$       ie = Max(Floor(z2*conv_cells(ic)%nvert),1)
!!$    EndIf
!!$
!!$    !Randomly choose a position on the edge, uniform
!!$    Call random_number(harvest=z)
!!$    If (ie .ne. conv_cells(ic)%nvert) Then
!!$       ph0 = conv_cells(ic)%vertices(1,ie)+(conv_cells(ic)%vertices(1,ie+1)-conv_cells(ic)%vertices(1,ie))*z
!!$       th0 = conv_cells(ic)%vertices(2,ie)+(conv_cells(ic)%vertices(2,ie+1)-conv_cells(ic)%vertices(2,ie))*z
!!$    Else
!!$       ph0 = conv_cells(ic)%vertices(1,ie)+(conv_cells(ic)%vertices(1,1)-conv_cells(ic)%vertices(1,ie))*z
!!$       th0 = conv_cells(ic)%vertices(2,ie)+(conv_cells(ic)%vertices(2,1)-conv_cells(ic)%vertices(2,ie))*z
!!$    EndIf
!!$
!!$    !Randomize around edge
!!$    Call random_number(harvest=z)
!!$    z = 2d0*z**2-1d0
!!$    ph0 = ph0 + isz_max*z*dphi
!!$
!!$    Call random_number(harvest=z)
!!$    z = 2d0*z**2-1d0
!!$    th0 = th0 + isz_max*z*dth
  End Subroutine Place_Squall

  Subroutine Place_Cell(th0,ph0)
    Implicit None
    Real*8, Intent(Out) :: th0, ph0
    Real*8 :: z
    !uniformly random
    th0 = 2d0*th2
    Do While ((th0 .gt. th2) .or. (th0 .lt. th1))
       Call random_number(harvest=z)
       z = 1.1d0*z-0.1d0
       th0 = th1+(th2-th1)*z
    EndDo

    ph0 = 2d0*phi2
    Do While ((ph0 .gt. phi2) .or. (ph0 .lt. phi1))
       Call random_number(harvest=z)
       z = 1.1d0*z-0.1d0
       ph0 = phi1+(phi2-phi1)*z
    EndDo
  End Subroutine Place_Cell

  Subroutine Random_Exponential(nrand,g_rand,g0,g1)
    Implicit None
    !This routine provides exponentially distributed random numbers
    Integer, Intent(In) :: nrand
    Real*8, Intent(In) :: g0,g1
    Real*8, Intent(Out), Dimension(:) :: g_rand
    Real*8 :: rand_tmp
    Integer :: ii
    ii = 1
    Do While (ii .le. nrand) 
       Call random_number(harvest=rand_tmp)
       rand_tmp = -dlog(rand_tmp)
       If ((rand_tmp .lt. g1).and.(rand_tmp .gt. g0)) Then
          g_rand(ii) = rand_tmp
          ii=ii+1
       EndIf
    EndDo
  End Subroutine Random_Exponential

  Subroutine Random_Gaussian(nrand,g_rand,g0,g1)
    Implicit None
    !This routine uses the Marsaglia polar method to build normally distributed random numbers
    !The probability density function is P(g) = exp(-g^2/2)   for g >= g0,
    !                                         = 0             for g <  g0.
    Integer, Intent(In) :: nrand
    Integer :: i 
    Real*8, Intent(In) :: g0, g1
    Real*8, Intent(Out), Dimension(:) :: g_rand
    Real*8 :: x, y, s, tmp
    i=0
    Do While (i .lt. nrand)
       Call random_number(harvest=x)
       Call random_number(harvest=y)
! Uniform unity square
       x = 2d0*(x-0.5d0)
       y = 2d0*(y-0.5d0)
       s = x*x+y*y
! Cut out uniform unity circle
       If (s .lt. 1d0) Then
          tmp = x*sqrt(-2d0*log(s)/s)
          If ((tmp .ge. g0).and.(tmp.le.g1)) Then             
             i=i+1
             g_rand(i) = tmp
          EndIf
       EndIf
    EndDo
  End Subroutine Random_Gaussian

End Module Entropy_Rain
