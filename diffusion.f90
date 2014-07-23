Module Diffusion
   
  Use Constants
  Use Derivatives

  Implicit None

  Integer :: slope_limiter, dstype, flux_limiter=-1, nsld, mag_flux_limiter=-1, ent_flux_limiter=-1
  Integer, Allocatable, Dimension(:) :: varlist
  Real*8 :: sld_coef, sld_speedr, sld_speedh, sld_pr, sld_prm, sld_exp, flux_cutoff, sld_tanh_coef
  Logical :: Mass_Diffusion=.True., Magnetic_Diffusion=.True., do_ohmic_heating=.True., do_viscous_heating=.True., Do_SLD_Report=.False., update_background=.False.
  Real*8, Allocatable, Dimension(:) :: diffusion_speed, tinyvar
  Real*8, Allocatable, Dimension(:,:) :: meanV, smoothing_function
  Real*8, Allocatable, Dimension(:,:,:) :: vcc
  !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: meanV, meanA, diffusion_speed, smoothing_function, vcc

Contains

  Subroutine Initialize_Diffusion()
    Implicit None
    Real*8 :: norm
    Integer :: ir
    Allocate(tinyvar(nv), meanV(nr,2), diffusion_speed(0:mynr+1), smoothing_function(nr,nr), vcc(0:mynth+1,0:mynphi+1,0:mynr+1))

    MeanV = 0d0

    smoothing_function = 0d0
    !smoothing_function(1:3,1) = (/0.25d0,0.5d0,0.25d0/)
    smoothing_function(1:3,1) = (/0.5d0,1d0/3d0,1d0/6d0/)   
    !smoothing_function(1:7,1) = (/-29d0/20d0,6d0,-15d0/2d0,20d0/3d0,-15d0/4d0,6d0/5d0,-1d0/6d0/)
    norm = Sum(smoothing_function(:,1)*dr_ind)
    !#DIR$ VECTOR ALIGNED
    smoothing_function(:,1) = smoothing_function(:,1)/norm

    Do ir = 2, nr-1
       !smoothing_function(:,ir) = exp(-onethird*((r_ind(:)-r_ind(ir))/dr_ind(ir))**2)
       smoothing_function(ir-1:ir+1,ir) = (/0.25d0,0.5d0,0.25d0/)

       norm = Sum(smoothing_function(:,ir)*dr_ind)
       !#DIR$ VECTOR ALIGNED
       smoothing_function(:,ir) = smoothing_function(:,ir)/norm
    EndDo

    smoothing_function(nr-2:nr,nr) = (/1d0/6d0,1d0/3d0,0.5d0/)
    !smoothing_function(nr-2:nr,nr) = (/0.25d0,0.5d0,0.5d0/)
    
    norm = Sum(smoothing_function(:,nr)*dr_ind)
    !#DIR$ VECTOR ALIGNED
    smoothing_function(:,nr) = smoothing_function(:,nr)/norm

    If (myr1 .eq. 1) Then
       diffusion_speed(0) = sld_speedh*(rhor(1)/rhor(nr))**(sld_exp) + sld_speedr
       diffusion_speed(1:mynr+1) = sld_speedh*(rhor(myr1:myr2+1)/rhor(nr))**(sld_exp) + sld_speedr
    ElseIf (myr2 .eq. nr) Then
       diffusion_speed(0:mynr) = sld_speedh*(rhor(myr1-1:myr2)/rhor(nr))**(sld_exp) + sld_speedr
       diffusion_speed(mynr+1) = sld_speedh + sld_speedr
    Else
       diffusion_speed = sld_speedh*(rhor(myr1-1:myr2+1)/rhor(nr))**(sld_exp) + sld_speedr
    EndIf
    diffusion_speed = 0.5d0*diffusion_speed*(1d0+tanh(sld_tanh_coef*(r_ind-flux_cutoff)/(r2-r1)))
    diffusion_speed = diffusion_speed**2

  End Subroutine Initialize_Diffusion

  Subroutine Finalize_Diffusion()
    Implicit None
    Deallocate(meanV, diffusion_speed, smoothing_function, vcc)
  End Subroutine Finalize_Diffusion

  Subroutine Compute_Diffusion(call_num)
    Implicit None
    Real*8 :: Vtmp(nr,2), Atmp(nr,2), norm
    Integer, Intent(In) :: call_num
    Integer :: ir, it, ip, iv

    If (magnetic) Then
       tinyvar = 1d-8*(/rhor(nr),1d0,1d0,1d0,1d5,1d8,1d8,1d8/)
    Else
       tinyvar = 1d-8*(/rhor(nr),1d0,1d0,1d0,1d5/)
    EndIf

    If (Do_Timings) Then
       tstart = MPI_WTIME()
    EndIf

    If (Allocated(varlist)) Deallocate(varlist)
    If (Do_SLD_Report) Then
       nsld = 2
       Allocate(varlist(nsld))
       varlist(1) = 3
       varlist(2) = 5
    Else
       If (Magnetic_Diffusion) Then
          If (Mass_Diffusion) Then
             nsld = 8
             Allocate(varlist(nsld))
             varlist = (/1,2,3,4,5,6,7,8/)
          Else
             nsld = 7
             Allocate(varlist(nsld))
             varlist = (/2,3,4,5,6,7,8/)
          EndIf
       Else
          If (Mass_Diffusion) Then
             nsld = 5
             Allocate(varlist(nsld))
             varlist = (/1,2,3,4,5/)
          Else
             nsld = 4
             Allocate(varlist(nsld))
             varlist = (/2,3,4,5/)
          EndIf
       EndIf
    EndIf

    !#DIR$ VECTOR ALIGNED
    invrho = 1d0/vars1(:,:,:,1)
    If (Unstratified) Then
       temperature = T0
    Else
       !#DIR$ VECTOR ALIGNED
       temperature = (vars1(:,:,:,1)**(-gamm1))*dexp(-vars1(:,:,:,5)*invCv)/invtconst !inverse temperature
    EndIf

    If (((mod(istep,ufreq).eq.0).or.(istep .eq. begiter+1)).and.(call_num .eq. 1)) Then

       If (update_background) Then
          Atmp = 0d0
          Do ir=1,mynr
             Do ip=1,mynphi
                Do it=1,mynth
                   Atmp(myr1+ir-1,1) = Atmp(myr1+ir-1,1) + vars1(it,ip,ir,1)
                   Atmp(myr1+ir-1,2) = Atmp(myr1+ir-1,2) + vars1(it,ip,ir,5)
                EndDo
             EndDo
          EndDo
          
          Vtmp = 0d0
          Call Global_AllReduce(Atmp,Vtmp,node_comm,'sum')
          !#DIR$ VECTOR ALIGNED
          Vtmp = Vtmp/Dble(nphi*nth)
          
          If (istep .eq. begiter+1) meanV = Vtmp

          iv=2
          ir=1
          Atmp(ir,iv) = Sum((/sr(ir)-dr_ind(ir)*dsdrm(ir),Vtmp(ir,iv),Vtmp(ir+1,iv)/)*smoothing_function(1:3,ir)*dr_ind(1:3))
          Do iv=1,2
             Do ir=2,nr
                Atmp(ir,iv) = Sum(Vtmp(:,iv)*smoothing_function(:,ir)*dr_ind)
             EndDo
          EndDo

          meanV = 0.99d0*meanV + 0.01d0*Atmp

          !meanV(:,2) = Vtmp(:,2)
          !meanV = Vtmp

          meanV(:,1) = Compute_Density(meanV(:,2))
          !meanV(:,2) = sr(1)
       Else
          meanV(:,1) = rhor
          meanV(:,2) = Sr
          !meanV = 0d0
       EndIf
    EndIf

    Do ir=1,mynr
       Do ip=1,mynphi
          Do it=1,mynth
             vars1(it,ip,ir,1) = vars1(it,ip,ir,1) - meanV(myr1+ir-1,1)
             vars1(it,ip,ir,5) = vars1(it,ip,ir,5) - meanV(myr1+ir-1,2)
          EndDo
       EndDo
    EndDo

    Do ir=0,mynr+1
       vcc(:,:,ir) = diffusion_speed(ir)
    EndDo

    If (dstype .eq. 1) Then
       !#DIR$ VECTOR ALIGNED
       Do iv=2,4
          If (myr1 .ne. 1) Then
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                vcc(1:mynth,ip,0) = vcc(1:mynth,ip,0) + (binfo(1,1)%bdata(:,ip,4,iv))**2
             EndDo
          EndIf

          !#DIR$ VECTOR ALIGNED
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             vcc(1:mynth,0,ir) = vcc(1:mynth,0,ir) + (binfo(3,1)%bdata(:,4,ir,iv))**2

             !#DIR$ VECTOR ALIGNED
             Do ip=1,mynphi
                If (myth1 .ne. 1) Then
                   vcc(0,ip,ir) = vcc(0,ip,ir) + (binfo(2,1)%bdata(4,ip,ir,iv))**2
                EndIf

                !#DIR$ VECTOR ALIGNED
                Do it=1,mynth
                   vcc(it,ip,ir) = vcc(it,ip,ir) + (vars1(it,ip,ir,iv))**2
                EndDo

                If (myth2 .ne. nth) Then
                   vcc(mynth+1,ip,ir) = vcc(mynth+1,ip,ir) + (binfo(2,2)%bdata(1,ip,ir,iv))**2
                EndIf
             EndDo

             !#DIR$ VECTOR ALIGNED
             vcc(1:mynth,mynphi+1,ir) = vcc(1:mynth,mynphi+1,ir) + (binfo(3,2)%bdata(:,1,ir,iv))**2
          EndDo

          If (myr2 .ne. nr) Then
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                vcc(1:mynth,ip,mynr+1) = vcc(1:mynth,ip,mynr+1) + (binfo(1,2)%bdata(:,ip,1,iv))**2
             EndDo
          EndIf
       EndDo
    EndIf

    If (magnetic) Then
       If (dstype .eq. 1) Then
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do iv=6,8
             If (myr1 .ne. 1) Then
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      vcc(it,ip,0) = vcc(it,ip,0) + one_over_4pi*(binfo(1,1)%bdata(it,ip,4,iv))**2/binfo(1,1)%bdata(it,ip,4,1)
                   EndDo
                EndDo
             EndIf

             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ir=1,mynr
                !#DIR$ VECTOR ALIGNED
                vcc(1:mynth,0,ir) = vcc(1:mynth,0,ir) + one_over_4pi*(binfo(3,1)%bdata(:,4,ir,iv))**2/binfo(3,1)%bdata(:,4,ir,1)

                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   If (myth1 .ne. 1) Then
                      !#DIR$ VECTOR ALIGNED
                      vcc(0,ip,ir) = vcc(0,ip,ir) + one_over_4pi*(binfo(2,1)%bdata(4,ip,ir,iv))**2/binfo(2,1)%bdata(4,ip,ir,1)
                   EndIf

                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      vcc(it,ip,ir) = vcc(it,ip,ir) + one_over_4pi*invrho(it,ip,ir)*(vars1(it,ip,ir,iv))**2
                   EndDo

                   If (myth2 .ne. nth) Then
                      !#DIR$ VECTOR ALIGNED
                      vcc(mynth+1,ip,ir) = vcc(mynth+1,ip,ir) + one_over_4pi*(binfo(2,2)%bdata(1,ip,ir,iv))**2/binfo(2,2)%bdata(1,ip,ir,1)
                   Endif
                EndDo

                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                vcc(1:mynth,mynphi+1,ir) = vcc(1:mynth,mynphi+1,ir) + one_over_4pi*(binfo(3,2)%bdata(:,1,ir,iv))**2/binfo(3,2)%bdata(:,1,ir,1)
             EndDo

             If (myr2 .ne. nr) Then
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      vcc(it,ip,mynr+1) = vcc(it,ip,mynr+1) + one_over_4pi*(binfo(1,2)%bdata(it,ip,1,iv))**2/binfo(1,2)%bdata(it,ip,1,1)
                   EndDo
                EndDo
             EndIf
          EndDo
       EndIf
    EndIf

    If (myr1 .eq. 1) Then
       vcc(:,:,0) = vcc(:,:,1)
    EndIf

    If (myr2 .eq. nr) Then
       vcc(:,:,mynr+1) = vcc(:,:,mynr)
    EndIf

    If (myth1 .eq. 1) Then
       vcc(0,:,:) = vcc(1,:,:)
    EndIf

    If (myth2 .eq. nth) Then
       vcc(mynth+1,:,:) = vcc(mynth,:,:)
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(12) = timings(12) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    !Do r direction
    Call Diffusion_R()

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(13) = timings(13) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. 7) Then
          Print*, 'Diffusion R'
          Print*, 'minmax vars1'
          Do iv=1,nv
             Print*, Minval(vars1(:,:,:,iv)), Maxval(vars1(:,:,:,iv))
          EndDo
          !Do iv=1,2
          !   Print*, 'MeanV(:,iv)=[', MeanV(:,iv),']'
          !EndDo
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    !Do theta direction
    Call Diffusion_Th()

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(14) = timings(14) + tstop-tstart
       tstart = MPI_WTIME()
    EndIf

    If (Debug) Then
       If (myrank .eq. 7) Then
          Print*, 'Diffusion T'
          Print*, 'minmax dvars1'
          Do iv=1,nv
             Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
          EndDo
       EndIf
    EndIf

    If (.not. Do_SLD_Output) Then
       !Do phi direction
       Call Diffusion_Phi()

       If (Do_Timings) Then
          tstop = MPI_WTIME()
          timings(15) = timings(15) + tstop-tstart
          tstart = MPI_WTIME()
       EndIf

       If (Debug) Then
          If (myrank .eq. 7) Then
             Print*, 'Diffusion P'
             Print*, 'minmax dvars1'
             Do iv=1,nv
                Print*, Minval(dvars1(:,:,:,iv)), Maxval(dvars1(:,:,:,iv))
             EndDo
          EndIf
       EndIf

       If (Do_Viscous_Heating) Then
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do iv=2,4
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ir=1,mynr
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      dvars1(it,ip,ir,5) = dvars1(it,ip,ir,5) - vars1(it,ip,ir,iv)*dvars1(it,ip,ir,iv)*temperature(it,ip,ir) !Viscous heating
                   EndDo
                EndDo
             EndDo
          EndDo
       EndIf

       If (Mass_Diffusion) Then

          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do iv=2,4
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ir=1,mynr
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      dvars1(it,ip,ir,iv) = dvars1(it,ip,ir,iv) - invrho(it,ip,ir)*vars1(it,ip,ir,iv)*dvars1(it,ip,ir,1) !Diffuse density, but do not change momentum
                   EndDo
                EndDo
             EndDo
          EndDo

          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do ir=1,mynr
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do it=1,mynth
                   dvars1(it,ip,ir,5) = dvars1(it,ip,ir,5) - (gamm1*Cv)*invrho(it,ip,ir)*dvars1(it,ip,ir,1) !Diffuse density, but do not change energy
                EndDo
             EndDo
          EndDo
       EndIf

       If (Do_Ohmic_Heating) Then
          !Curl dA/dt
          temperature = (0.5d0*one_over_4pi)*invrho*temperature
          maga(:,:,:,1) = dvars1(:,:,:,6)
          maga(:,:,:,2) = dvars1(:,:,:,7)
          maga(:,:,:,3) = dvars1(:,:,:,8)
          
          Call VectorCurlA()

          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do iv=6,8
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ir=1,mynr
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      dvars1(it,ip,ir,5) = dvars1(it,ip,ir,5) - vars1(it,ip,ir,iv)*dvars1(it,ip,ir,iv)*temperature(it,ip,ir) !Viscous heating
                   EndDo
                EndDo
             EndDo
          EndDo

          !Restore dA/dt terms
          dvars1(:,:,:,6) = maga(:,:,:,1)
          dvars1(:,:,:,7) = maga(:,:,:,2)
          dvars1(:,:,:,8) = maga(:,:,:,3)
       EndIf
    EndIf

    If (Do_Timings) Then
       tstop = MPI_WTIME()
       timings(12) = timings(12) + tstop-tstart
    EndIf

    If (Unstratified) Then
       dvars1(:,:,:,1:5) = 0d0
       If (magnetic) dvars1(:,:,:,7:8) = 0d0
    EndIf

  End Subroutine Compute_Diffusion

  Subroutine Diffusion_R()
    Implicit None
    Real*8, Dimension(mynth,mynr) :: tmp, flux
    Real*8, Dimension(mynth,mynr+1) :: vtmp, delue, rph, qmh, qph, qp3h, phiL, phiR
    Real*8 :: b1(mynth,4), b2(mynth,4), rtmp
    Integer :: ir, it, ip, iv, isl

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: tmp, flux, vtmp, delue, rph, b1, b2, qmh, qph, qp3h, phiL, phiR

    !#DIR$ VECTOR ALIGNED
    Do isl=1,nsld
       iv = varlist(isl)
       !#DIR$ VECTOR ALIGNED
       Do ip=1,mynphi

          If (Do_SLD_Output) Then
             If (iv .eq. 3) Then
                !#DIR$ VECTOR ALIGNED
                tmp = vars1(:,ip,:,1)*vars1(:,ip,:,3)
                If (myr1 .ne. 1) Then
                   Do ir=1,4
                      !#DIR$ VECTOR ALIGNED
                      b1(:,ir) = (binfo(1,1)%bdata(:,ip,ir,1)-meanV(myr1-5+ir,1))*(binfo(1,1)%bdata(:,ip,ir,3))
                   EndDo
                EndIf
                If (myr2 .ne. nr) Then
                   Do ir=1,4
                      !#DIR$ VECTOR ALIGNED
                      b2(:,ir) = (binfo(1,2)%bdata(:,ip,ir,1)-meanV(myr2+ir,1))*(binfo(1,2)%bdata(:,ip,ir,3))
                   EndDo
                EndIf
             ElseIf (iv .eq. 5) Then
                tmp = vars1(:,ip,:,5)
                If (myr1 .ne. 1) Then
                   Do ir=1,4
                      !#DIR$ VECTOR ALIGNED
                      b1(:,ir) = binfo(1,1)%bdata(:,ip,ir,5)-meanV(myr1-5+ir,2)
                   EndDo
                EndIf
                If (myr2 .ne. nr) Then
                   Do ir=1,4
                      !#DIR$ VECTOR ALIGNED
                      b2(:,ir) = binfo(1,2)%bdata(:,ip,ir,5)-meanV(myr2+ir,2)
                   EndDo
                EndIf
             EndIf
          Else
             If (iv .gt. 5) Then
                tmp = maga(:,ip,:,iv-5)
                If (myr1 .ne. 1) Then
                   Do ir=1,4
                      !#DIR$ VECTOR ALIGNED
                      b1(:,ir) = dbinfo(1,1)%bdata(:,ip,ir,iv-5)
                   EndDo
                EndIf
                If (myr2 .ne. nr) Then
                   Do ir=1,4
                      !#DIR$ VECTOR ALIGNED
                      b2(:,ir) = dbinfo(1,2)%bdata(:,ip,ir,iv-5)
                   EndDo
                EndIf
             Else
                tmp = vars1(:,ip,:,iv)
                If (iv .eq. 1) Then
                   If (myr1 .ne. 1) Then
                      Do ir=1,4
                         !#DIR$ VECTOR ALIGNED
                         b1(:,ir) = binfo(1,1)%bdata(:,ip,ir,iv)-meanV(myr1-5+ir,iv)
                      EndDo
                   EndIf
                   If (myr2 .ne. nr) Then
                      Do ir=1,4
                         !#DIR$ VECTOR ALIGNED
                         b2(:,ir) = binfo(1,2)%bdata(:,ip,ir,iv)-meanV(myr2+ir,iv)
                      EndDo
                   EndIf
                ElseIf (iv .eq. 5) Then
                   If (myr1 .ne. 1) Then
                      Do ir=1,4
                         !#DIR$ VECTOR ALIGNED
                         b1(:,ir) = binfo(1,1)%bdata(:,ip,ir,iv)-meanV(myr1-5+ir,2)
                      EndDo
                   EndIf
                   If (myr2 .ne. nr) Then
                      Do ir=1,4
                         !#DIR$ VECTOR ALIGNED
                         b2(:,ir) = binfo(1,2)%bdata(:,ip,ir,iv)-meanV(myr2+ir,2)
                      EndDo
                   EndIf
                Else
                   If (myr1 .ne. 1) Then
                      Do ir=1,4
                         !#DIR$ VECTOR ALIGNED
                         b1(:,ir) = binfo(1,1)%bdata(:,ip,ir,iv)
                      EndDo
                   EndIf
                   If (myr2 .ne. nr) Then
                      Do ir=1,4
                         !#DIR$ VECTOR ALIGNED
                         b2(:,ir) = binfo(1,2)%bdata(:,ip,ir,iv)
                      EndDo
                   EndIf
                EndIf
             EndIf
          EndIf
 
          If (myr1 .eq. 1) Then
             Select Case (iv)
             Case (1)
                b1 = 0d0
                !b1(:,1) = rhob1(1)*tmp(:,1)/rhor(1)
                !b1(:,2) = rhob1(2)*tmp(:,1)/rhor(1)
                !b1(:,3) = rhob1(3)*tmp(:,1)/rhor(1)
                !b1(:,4) = rhob1(4)*tmp(:,1)/rhor(1)
                !b1(:,1) = tmp(:,1)
                !b1(:,2) = tmp(:,1)
                !b1(:,3) = tmp(:,1)
                !b1(:,4) = tmp(:,1)
             Case (2)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                !b1(:,3) = tmp(:,1)
                !b1(:,4) = tmp(:,1)
                b1(:,4) = -10d0*onethird*tmp(:,1)/r_ind(1)+6d0*tmp(:,2)/r_ind(2)-2d0*tmp(:,3)/r_ind(3)+onethird*tmp(:,4)/r_ind(4)
                b1(:,3) = 8d0*b1(:,4)/(r_ind(1)-dr_ind(1))-8d0*tmp(:,2)/r_ind(2)+tmp(:,3)/r_ind(3)
             Case (3)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                !b1(:,3) = tmp(:,1)
                !b1(:,4) = tmp(:,1)
                b1(:,4) = -10d0*onethird*tmp(:,1)/r_ind(1)+6d0*tmp(:,2)/r_ind(2)-2d0*tmp(:,3)/r_ind(3)+onethird*tmp(:,4)/r_ind(4)
                b1(:,3) = 8d0*b1(:,4)/(r_ind(1)-dr_ind(1))-8d0*tmp(:,2)/r_ind(2)+tmp(:,3)/r_ind(3)
             Case (4)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                b1(:,3) = tmp(:,1)
                b1(:,4) = tmp(:,1)
             Case (5)
                If ((Bottom_Constant_Flux).and.(Laplacian.or.SLaplacian)) Then
                   b1(:,1) = tmp(:,1) + 4d0*dr_ind(1)*lumr1*one_over_4pi/(r1**2*ks(1)*dxdr(1))
                   b1(:,2) = tmp(:,1) + 3d0*dr_ind(1)*lumr1*one_over_4pi/(r1**2*ks(1)*dxdr(1))
                   b1(:,3) = tmp(:,1) + 2d0*dr_ind(1)*lumr1*one_over_4pi/(r1**2*ks(1)*dxdr(1))
                   b1(:,4) = tmp(:,1) + dr_ind(1)*lumr1*one_over_4pi/(r1**2*ks(1)*dxdr(1))
                Else
                   b1 = 0d0
                   !b1(:,1) = sr(1)-4d0*dr_ind(1)*dsdrm(1)
                   !b1(:,2) = sr(1)-3d0*dr_ind(1)*dsdrm(1)
                   !b1(:,3) = sr(1)-2d0*dr_ind(1)*dsdrm(1)
                   !b1(:,4) = sr(1)-dr_ind(1)*dsdrm(1)
                   !b1(:,1) = tmp(:,1)
                   !b1(:,2) = tmp(:,1)
                   !b1(:,3) = tmp(:,1)
                   !b1(:,4) = tmp(:,1)
                EndIf
             Case (6)
                If (Bottom_Perfect_Cond) Then
                   b1 = 0d0
                Else
                   b1(:,1) = tmp(:,1)
                   b1(:,2) = tmp(:,1)
                   b1(:,3) = tmp(:,1)
                   b1(:,4) = tmp(:,1)
                EndIf
             Case (7)
                If (Bottom_Perfect_Cond) Then
                   b1 = 0d0
                Else
                   b1(:,1) = tmp(:,1)
                   b1(:,2) = tmp(:,1)
                   b1(:,3) = tmp(:,1)
                   b1(:,4) = tmp(:,1)
                EndIf
             Case (8)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                b1(:,3) = tmp(:,1)
                b1(:,4) = tmp(:,1)
             End Select
          EndIf

          If (myr2 .eq. nr) Then
             Select Case (iv)
             Case (1)
                !b2(:,1) = -10d0*onethird*tmp(:,mynr) + 6d0*tmp(:,mynr-1) - 2d0*tmp(:,mynr-2) + onethird*tmp(:,mynr-3) - (4d0*dr*invpconst*gravity(nr)/tr(nr))*tmp(:,mynr)
                !b2(:,2) = tmp(:,mynr-2) - 8d0*tmp(:,mynr-1) + 8d0*b2(:,1) + (12d0*dr*invpconst*gravity(nr)/tr(nr))*tmp(:,mynr)
                !b2(:,1) = rhob2(1)*tmp(:,mynr)/rhor(nr)
                !b2(:,2) = rhob2(2)*tmp(:,mynr)/rhor(nr)
                !b2(:,3) = rhob2(3)*tmp(:,mynr)/rhor(nr)
                !b2(:,4) = rhob2(4)*tmp(:,mynr)/rhor(nr)
                b2 = 0d0
             Case (2)
                b2(:,1) = -10d0*onethird*tmp(:,mynr)/r_ind(nr)+6d0*tmp(:,mynr-1)/r_ind(nr-1)-2d0*tmp(:,mynr-2)/r_ind(nr-2)+onethird*tmp(:,mynr-3)/r_ind(nr-3)
                b2(:,2) = tmp(:,mynr-2)/r_ind(nr-2)-8d0*tmp(:,mynr-1)/r_ind(nr-1)+8d0*b2(:,1)/(r_ind(nr)+dr_ind(nr))
                !b2(:,1) = tmp(:,mynr)
                !b2(:,2) = tmp(:,mynr)
                b2(:,3) = tmp(:,mynr)
                b2(:,4) = tmp(:,mynr)
                !b2 = 0d0
             Case (3)
                b2(:,1) = -10d0*onethird*tmp(:,mynr)/r_ind(nr)+6d0*tmp(:,mynr-1)/r_ind(nr-1)-2d0*tmp(:,mynr-2)/r_ind(nr-2)+onethird*tmp(:,mynr-3)/r_ind(nr-3)
                b2(:,2) = tmp(:,mynr-2)/r_ind(nr-2)-8d0*tmp(:,mynr-1)/r_ind(nr-1)+8d0*b2(:,1)/(r_ind(nr)+dr_ind(nr))
                !b2(:,1) = tmp(:,mynr)
                !b2(:,2) = tmp(:,mynr)
                b2(:,3) = tmp(:,mynr)
                b2(:,4) = tmp(:,mynr)
                !b2 = 0d0
             Case (4)
                b2(:,1) = tmp(:,mynr)
                b2(:,2) = tmp(:,mynr)
                b2(:,3) = tmp(:,mynr)
                b2(:,4) = tmp(:,mynr)
             Case (5)
                If ((Top_Constant_Flux).and.(Laplacian.or.SLaplacian)) Then
                   b2(:,1) = tmp(:,mynr) - dr_ind(1)*lumr2*one_over_4pi/(r2**2*ks(nr)*dxdr(nr))
                   b2(:,2) = tmp(:,mynr) - 2d0*dr_ind(1)*lumr2*one_over_4pi/(r2**2*ks(nr)*dxdr(nr))
                   b2(:,3) = tmp(:,mynr) - 3d0*dr_ind(1)*lumr2*one_over_4pi/(r2**2*ks(nr)*dxdr(nr))
                   b2(:,4) = tmp(:,mynr) - 4d0*dr_ind(1)*lumr2*one_over_4pi/(r2**2*ks(nr)*dxdr(nr))
                Else !Extrapolate
                   !b2(:,1) = sr(nr)+dr_ind(nr)*dsdrm(nr)
                   !b2(:,2) = sr(nr)+2d0*dr_ind(nr)*dsdrm(nr)
                   !b2(:,3) = sr(nr)+3d0*dr_ind(nr)*dsdrm(nr)
                   !b2(:,4) = sr(nr)+4d0*dr_ind(nr)*dsdrm(nr)
                   !b2(:,1) = tmp(:,mynr)
                   !b2(:,2) = tmp(:,mynr)
                   !b2(:,3) = tmp(:,mynr)
                   !b2(:,4) = tmp(:,mynr)
                   b2 = 0d0
                EndIf
             Case (6)
                If (Top_Perfect_Cond) Then
                   b2 = 0d0
                Else
                   b2(:,1) = tmp(:,mynr)
                   b2(:,2) = tmp(:,mynr)
                   b2(:,3) = tmp(:,mynr)
                   b2(:,4) = tmp(:,mynr)
                EndIf
             Case (7)
                If (Top_Perfect_Cond) Then
                   b2 = 0d0
                Else
                   b2(:,1) = tmp(:,mynr)
                   b2(:,2) = tmp(:,mynr)
                   b2(:,3) = tmp(:,mynr)
                   b2(:,4) = tmp(:,mynr)
                EndIf
             Case (8)
                b2(:,1) = tmp(:,mynr)
                b2(:,2) = tmp(:,mynr)
                b2(:,3) = tmp(:,mynr)
                b2(:,4) = tmp(:,mynr)
             End Select
          EndIf

          !Gradients
          !#DIR$ VECTOR ALIGNED
          qmh(:,1) = b1(:,4) - b1(:,3)
          qmh(:,2) = tmp(:,1) - b1(:,4)
          !#DIR$ VECTOR ALIGNED
          qmh(:,3:mynr+1) = tmp(:,2:mynr) - tmp(:,1:mynr-1)

          !#DIR$ VECTOR ALIGNED
          qph(:,1:mynr) = qmh(:,2:mynr+1)
          qph(:,mynr+1) = b2(:,1) - tmp(:,mynr)

          !#DIR$ VECTOR ALIGNED
          qp3h(:,1:mynr) = qph(:,2:mynr+1)
          qp3h(:,mynr+1) = b2(:,2) - b2(:,1)

          !Van Albada
          !#DIR$ VECTOR ALIGNED
          Do ir=1,mynr+1
             Do it=1,mynth
                rtmp = abs(qmh(it,ir))
                rtmp = 1d0/max(rtmp,tinyvar(iv))
                rtmp = qmh(it,ir)*rtmp**2
                rtmp = rtmp*qph(it,ir)
                rtmp = min(1d0,rtmp) !rtmp*(rtmp+1d0)/(rtmp*rtmp+1d0)
                phiL(it,ir) = qmh(it,ir)*max(0d0,rtmp)

                rtmp = abs(qp3h(it,ir))
                rtmp = 1d0/max(rtmp,tinyvar(iv))
                rtmp = qp3h(it,ir)*rtmp**2
                rtmp = rtmp*qph(it,ir)
                rtmp = min(1d0,rtmp) !rtmp*(rtmp+1d0)/(rtmp*rtmp+1d0)
                phiR(it,ir) = qp3h(it,ir)*max(0d0,rtmp)
             EndDo
          EndDo

          !#DIR$ VECTOR ALIGNED
          delue = qph - 0.5d0*(phiR+phiL)

          If ((iv .lt. 5).and.(flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          ElseIf ((iv .eq. 5).and.(ent_flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**ent_flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          ElseIf ((iv .gt. 5).and.(mag_flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**mag_flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          EndIf

          If (myr1 .eq. 1) Then
             !#DIR$ VECTOR ALIGNED
             rph(:,1) = r_ind(1)-0.5d0*dr_ind(1)
             !#DIR$ VECTOR ALIGNED
             delue(:,1) = 0d0 !delue(:,1)*rph(:,1)**2
          Else
             !#DIR$ VECTOR ALIGNED
             rph(:,1) = 0.5d0*(r_ind(myr1-1)+r_ind(myr1))
             !#DIR$ VECTOR ALIGNED
             delue(:,1) = delue(:,1)*rph(:,1)**2
          EndIf
          
          !#DIR$ VECTOR ALIGNED
          Do ir=2,mynr
             !#DIR$ VECTOR ALIGNED
             rph(:,ir) = 0.5d0*(r_ind(myr1+ir-1)+r_ind(myr1+ir-2))
          EndDo

          !#DIR$ VECTOR ALIGNED
          delue(:,2:mynr) = delue(:,2:mynr)*rph(:,2:mynr)**2
          
          If (myr2 .eq. nr) Then
             !#DIR$ VECTOR ALIGNED
             rph(:,mynr+1) = r_ind(myr2)+0.5d0*dr_ind(myr2)
             !#DIR$ VECTOR ALIGNED
             delue(:,mynr+1) = 0d0 !delue(:,mynr+1)*rph(:,mynr+1)**2
          Else
             !#DIR$ VECTOR ALIGNED
             rph(:,mynr+1) = 0.5d0*(r_ind(myr2+1)+r_ind(myr2))
             !#DIR$ VECTOR ALIGNED
             delue(:,mynr+1) = delue(:,mynr+1)*rph(:,mynr+1)**2
          EndIf

          !#DIR$ VECTOR ALIGNED
          vtmp = sld_coef*sqrt(0.5d0*(vcc(1:mynth,ip,1:mynr+1) + vcc(1:mynth,ip,0:mynr)))

          !#DIR$ VECTOR ALIGNED
          delue = vtmp*delue

          !#DIR$ VECTOR ALIGNED
          tmp = rph(:,2:mynr+1)**3-rph(:,1:mynr)**3

          !#DIR$ VECTOR ALIGNED
          tmp = 3d0/tmp

          !#DIR$ VECTOR ALIGNED
          flux = tmp*(delue(:,2:mynr+1)-delue(:,1:mynr))
          
          If ((sld_pr.ne.1d0).or.(sld_prm.ne.1d0)) Then
             If (iv .le. 3) Then
                !#DIR$ VECTOR ALIGNED
                flux = sld_pr*flux
             ElseIf (iv .ge. 5) Then
                !#DIR$ VECTOR ALIGNED
                flux = sld_prm*flux
             EndIf
          EndIf

          If (Do_SLD_Output) Then
             If (iv .eq. 3) Then
                !#DIR$ VECTOR ALIGNED
                Do ir=1,mynr
                   !#DIR$ VECTOR ALIGNED
                   dvars1(:,ip,ir,1) = flux(:,ir)*r_ind(myr1+ir-1)*sines(myth1:myth2) !Flux at cell center
                EndDo
             ElseIf (iv .eq. 5) Then
                !#DIR$ VECTOR ALIGNED
                Do ir=1,mynr
                   !#DIR$ VECTOR ALIGNED
                   dvars1(:,ip,ir,3) = flux(:,ir) !Flux at cell center
                EndDo
             EndIf
          Else
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ir=1,mynr
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do it=1,mynth
                   dvars1(it,ip,ir,iv) = dvars1(it,ip,ir,iv) + flux(it,ir)
                EndDo
             EndDo
          EndIf
       EndDo
    EndDo
  End Subroutine Diffusion_R

  Subroutine Diffusion_Th()
    Implicit None
    Real*8, Dimension(mynphi,mynth) :: tmp, flux
    Real*8, Dimension(mynphi,mynth+1) :: vtmp, delue, sph, cph, qmh, qph, qp3h, phiL, phiR
    Real*8 :: b1(mynphi,4), b2(mynphi,4)
    Real*8 :: stmp, rtmp
    Integer :: ir, it, ip, iv, isl

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: tmp, flux, vtmp, delue, sph, cph, b1, b2, qmh, qph, qp3h, phiL, phiR

    !#DIR$ VECTOR ALIGNED
    Do isl=1,nsld
       iv = varlist(isl)
       !#DIR$ VECTOR ALIGNED
       Do ir=1,mynr
          If (Do_SLD_Output) Then
             If (iv .eq. 3) Then
                !#DIR$ VECTOR ALIGNED
                tmp = transpose(vars1(:,:,ir,1)*vars1(:,:,ir,3))
                If (myth1 .ne. 1) Then
                   Do ip=1,mynphi
                      !#DIR$ VECTOR ALIGNED
                      b1(ip,:) = (binfo(2,1)%bdata(:,ip,ir,1)-meanV(myr1+ir-1,1))*(binfo(2,1)%bdata(:,ip,ir,3))
                   EndDo
                EndIf
                If (myth2 .ne. nth) Then
                   Do ip=1,mynphi
                      !#DIR$ VECTOR ALIGNED
                      b2(ip,:) = (binfo(2,2)%bdata(:,ip,ir,1)-meanV(myr1+ir-1,1))*(binfo(2,2)%bdata(:,ip,ir,3))
                   EndDo
                EndIf
             ElseIf (iv .eq. 5) Then
                tmp = transpose(vars1(:,:,ir,5))
                If (myth1 .ne. 1) Then
                   Do ip=1,mynphi
                      !#DIR$ VECTOR ALIGNED
                      b1(ip,:) = binfo(2,1)%bdata(:,ip,ir,5)-meanV(myr1+ir-1,2)
                   EndDo
                EndIf
                If (myth2 .ne. nth) Then
                   Do ip=1,mynphi
                      !#DIR$ VECTOR ALIGNED
                      b2(ip,:) = binfo(2,2)%bdata(:,ip,ir,5)-meanV(myr1+ir-1,2)
                   EndDo
                EndIf
             EndIf
          Else
             If (iv .gt. 5) Then
                !#DIR$ VECTOR ALIGNED
                tmp = transpose(maga(:,:,ir,iv-5))
                If (myth1 .ne. 1) Then
                   !#DIR$ VECTOR ALIGNED
                   Do ip=1,mynphi
                      !#DIR$ VECTOR ALIGNED
                      b1(ip,:) = dbinfo(2,1)%bdata(:,ip,ir,iv-5)
                   EndDo
                EndIf
                If (myth2 .ne. nth) Then
                   !#DIR$ VECTOR ALIGNED
                   Do ip=1,mynphi
                      !#DIR$ VECTOR ALIGNED
                      b2(ip,:) = dbinfo(2,2)%bdata(:,ip,ir,iv-5)
                   EndDo
                EndIf
             Else
                tmp = transpose(vars1(:,:,ir,iv))
                If (iv .eq. 1) Then
                   If (myth1 .ne. 1) Then
                      !#DIR$ VECTOR ALIGNED
                      Do ip=1,mynphi
                         !#DIR$ VECTOR ALIGNED
                         b1(ip,:) = binfo(2,1)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,iv)
                      EndDo
                   EndIf
                   If (myth2 .ne. nth) Then
                      !#DIR$ VECTOR ALIGNED
                      Do ip=1,mynphi
                         !#DIR$ VECTOR ALIGNED
                         b2(ip,:) = binfo(2,2)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,iv)
                      EndDo
                   EndIf
                ElseIf (iv .eq. 5) Then
                  If (myth1 .ne. 1) Then
                      !#DIR$ VECTOR ALIGNED
                      Do ip=1,mynphi
                         !#DIR$ VECTOR ALIGNED
                         b1(ip,:) = binfo(2,1)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,2)
                      EndDo
                   EndIf
                   If (myth2 .ne. nth) Then
                      !#DIR$ VECTOR ALIGNED
                      Do ip=1,mynphi
                         !#DIR$ VECTOR ALIGNED
                         b2(ip,:) = binfo(2,2)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,2)
                      EndDo
                   EndIf
                Else
                  If (myth1 .ne. 1) Then
                      !#DIR$ VECTOR ALIGNED
                      Do ip=1,mynphi
                         !#DIR$ VECTOR ALIGNED
                         b1(ip,:) = binfo(2,1)%bdata(:,ip,ir,iv)
                      EndDo
                   EndIf
                   If (myth2 .ne. nth) Then
                      !#DIR$ VECTOR ALIGNED
                      Do ip=1,mynphi
                         !#DIR$ VECTOR ALIGNED
                         b2(ip,:) = binfo(2,2)%bdata(:,ip,ir,iv)
                      EndDo
                   EndIf
                EndIf
             EndIf
          EndIf

          If (myth1 .eq. 1) Then
             Select Case (iv)
             Case (1)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                b1(:,3) = tmp(:,1)
                b1(:,4) = tmp(:,1)
             Case (2)
                b1 = 0d0
             Case (3)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                !b1(:,3) = tmp(:,1)
                !b1(:,4) = tmp(:,1)
                b1(:,4) = -10d0*onethird*tmp(:,1)/sines(1)+6d0*tmp(:,2)/sines(2)-2d0*tmp(:,3)/sines(3)+onethird*tmp(:,4)/sines(4)
                b1(:,3) = 8d0*b1(:,4)/sin(theta(1)-dth)-8d0*tmp(:,2)/sines(2)+tmp(:,3)/sines(3)
             Case (4)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                b1(:,3) = tmp(:,1)
                b1(:,4) = tmp(:,1)
             Case (5)
                b1(:,1) = tmp(:,1)
                b1(:,2) = tmp(:,1)
                b1(:,3) = tmp(:,1)
                b1(:,4) = tmp(:,1)
             Case (6)
                If (Theta_Normal_Field) Then
                   b1 = 0d0
                Else
                   b1(:,1) = tmp(:,1)
                   b1(:,2) = tmp(:,1)
                   b1(:,3) = tmp(:,1)
                   b1(:,4) = tmp(:,1)
                EndIf
             Case (7)
                If (Theta_Perf_Conductor) Then
                   b1 = 0d0
                Else
                   b1(:,1) = tmp(:,1)
                   b1(:,2) = tmp(:,1)
                   b1(:,3) = tmp(:,1)
                   b1(:,4) = tmp(:,1)
                EndIf
             Case (8)
                If (Theta_Perf_Conductor) Then
                   b1 = 0d0
                Else
                   b1(:,1) = tmp(:,1)
                   b1(:,2) = tmp(:,1)
                   b1(:,3) = tmp(:,1)
                   b1(:,4) = tmp(:,1)
                EndIf
             End Select
          EndIf

          If (myth2 .eq. nth) Then
             Select Case (iv)
             Case (1)
                b2(:,1) = tmp(:,mynth)
                b2(:,2) = tmp(:,mynth)
                b2(:,3) = tmp(:,mynth)
                b2(:,4) = tmp(:,mynth)
             Case (2)
                b2 = 0d0
             Case (3)
                b2(:,1) = -10d0*onethird*tmp(:,mynth)/sines(nth)+6d0*tmp(:,mynth-1)/sines(nth-1)-2d0*tmp(:,mynth-2)/sines(nth-2)+onethird*tmp(:,mynth-3)/sines(nth-3)
                b2(:,2) = tmp(:,mynth-2)/sines(nth-2)-8d0*tmp(:,mynth-1)/sines(nth-1)+8d0*b2(:,1)/sin(theta(nth)+dth)
                !b2(:,1) = tmp(:,mynth)
                !b2(:,2) = tmp(:,mynth)
                b2(:,3) = tmp(:,mynth)
                b2(:,4) = tmp(:,mynth)
             Case (4)
                b2(:,1) = tmp(:,mynth)
                b2(:,2) = tmp(:,mynth)
                b2(:,3) = tmp(:,mynth)
                b2(:,4) = tmp(:,mynth)
             Case (5)
                b2(:,1) = tmp(:,mynth)
                b2(:,2) = tmp(:,mynth)
                b2(:,3) = tmp(:,mynth)
                b2(:,4) = tmp(:,mynth)
             Case (6)
                If (Theta_Normal_Field) Then
                   b2 = 0d0
                Else
                   b2(:,1) = tmp(:,mynth)
                   b2(:,2) = tmp(:,mynth)
                   b2(:,3) = tmp(:,mynth)
                   b2(:,4) = tmp(:,mynth)
                EndIf
             Case (7)
                If (Theta_Perf_Conductor) Then
                   b2 = 0d0
                Else
                   b2(:,1) = tmp(:,mynth)
                   b2(:,2) = tmp(:,mynth)
                   b2(:,3) = tmp(:,mynth)
                   b2(:,4) = tmp(:,mynth)
                EndIf
             Case (8)
                If (Theta_Perf_Conductor) Then
                   b2 = 0d0
                Else
                   b2(:,1) = tmp(:,mynth)
                   b2(:,2) = tmp(:,mynth)
                   b2(:,3) = tmp(:,mynth)
                   b2(:,4) = tmp(:,mynth)
                EndIf
             End Select
          EndIf

          !Gradients
          !#DIR$ VECTOR ALIGNED
          qmh(:,1) = b1(:,4) - b1(:,3)
          qmh(:,2) = tmp(:,1) - b1(:,4)
          !#DIR$ VECTOR ALIGNED
          qmh(:,3:mynth+1) = tmp(:,2:mynth) - tmp(:,1:mynth-1)

          !#DIR$ VECTOR ALIGNED
          qph(:,1:mynth) = qmh(:,2:mynth+1)
          qph(:,mynth+1) = b2(:,1) - tmp(:,mynth)

          !#DIR$ VECTOR ALIGNED
          qp3h(:,1:mynth) = qph(:,2:mynth+1)
          qp3h(:,mynth+1) = b2(:,2) - b2(:,1)

          !Van Albada
          !#DIR$ VECTOR ALIGNED
          Do it=1,mynth+1
             Do ip=1,mynphi
                rtmp = abs(qmh(ip,it))
                rtmp = 1d0/max(rtmp,tinyvar(iv))
                rtmp = qmh(ip,it)*rtmp**2
                rtmp = rtmp*qph(ip,it)
                rtmp = min(1d0,rtmp) !rtmp*(rtmp+1d0)/(rtmp*rtmp+1d0)
                phiL(ip,it) = qmh(ip,it)*max(0d0,rtmp)

                rtmp = abs(qp3h(ip,it))
                rtmp = 1d0/max(rtmp,tinyvar(iv))
                rtmp = qp3h(ip,it)*rtmp**2
                rtmp = rtmp*qph(ip,it)
                rtmp = min(1d0,rtmp) !rtmp*(rtmp+1d0)/(rtmp*rtmp+1d0)
                phiR(ip,it) = qp3h(ip,it)*max(0d0,rtmp)
             EndDo
          EndDo

          !#DIR$ VECTOR ALIGNED
          delue = qph - 0.5d0*(phiR+phiL)

          If ((iv .lt. 5).and.(flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          ElseIf ((iv .eq. 5).and.(ent_flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**ent_flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          ElseIf ((iv .gt. 5).and.(mag_flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**mag_flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          EndIf

          If (myth1 .eq. 1) Then
             delue(:,1) = 0d0
             !#DIR$ VECTOR ALIGNED
             sph(:,1) = abs(sin(theta(1)-0.5d0*dth))
             !#DIR$ VECTOR ALIGNED
             cph(:,1) = cos(theta(1)-0.5d0*dth)
          Else
             !#DIR$ VECTOR ALIGNED
             sph(:,1) = sin(0.5d0*(theta(myth1-1)+theta(myth1)))
             !#DIR$ VECTOR ALIGNED
             cph(:,1) = cos(0.5d0*(theta(myth1-1)+theta(myth1)))
          EndIf
          
          !#DIR$ VECTOR ALIGNED
          Do it=2,mynth
             !#DIR$ VECTOR ALIGNED
             sph(:,it) = sin(0.5d0*(theta(myth1+it-1)+theta(myth1+it-2)))
             !#DIR$ VECTOR ALIGNED
             cph(:,it) = cos(0.5d0*(theta(myth1+it-1)+theta(myth1+it-2)))
          EndDo
          
          If (myth2 .eq. nth) Then
             delue(:,mynth+1) = 0d0
             !#DIR$ VECTOR ALIGNED
             sph(:,mynth+1) = sin(theta(myth2)+0.5d0*dth)
             !#DIR$ VECTOR ALIGNED
             cph(:,mynth+1) = cos(theta(myth2)+0.5d0*dth)
          Else
             !#DIR$ VECTOR ALIGNED
             sph(:,mynth+1) = sin(0.5d0*(theta(myth2)+theta(myth2+1)))
             !#DIR$ VECTOR ALIGNED
             cph(:,mynth+1) = cos(0.5d0*(theta(myth2)+theta(myth2+1)))
          EndIf

          !#DIR$ VECTOR ALIGNED          
          vtmp = sld_coef*transpose(sqrt(0.5d0*(vcc(1:mynth+1,1:mynphi,ir)+vcc(0:mynth,1:mynphi,ir))))

          !#DIR$ VECTOR ALIGNED
          delue = vtmp*sph*delue

          If ((myr1.eq.1).and.(ir.eq.1)) Then
             stmp = r_ind(1)-0.5d0*dr_ind(1)
          Else
             stmp = 0.5d0*(r_ind(myr1+ir-1)+r_ind(myr1+ir-2))
          EndIf

          If ((myr2.eq.nr).and.(ir.eq.mynr)) Then
             rtmp = r_ind(nr)+0.5d0*dr_ind(nr)
          Else
             rtmp = 0.5d0*(r_ind(myr1+ir)+r_ind(myr1+ir-1))
          EndIf
          rtmp = 1.5d0*(rtmp**2-stmp**2)/(rtmp**3-stmp**3)

          !#DIR$ VECTOR ALIGNED
          tmp = rtmp/(cph(:,1:mynth)-cph(:,2:mynth+1))
          
          !#DIR$ VECTOR ALIGNED
          flux = tmp*(delue(:,2:mynth+1)-delue(:,1:mynth))

          If ((sld_pr.ne.1d0).or.(sld_prm.ne.1d0)) Then
             If (iv.lt.5) Then
                !#DIR$ VECTOR ALIGNED
                flux = sld_pr*flux
             ElseIf (iv .gt. 5) Then
                !#DIR$ VECTOR ALIGNED
                flux = sld_prm*flux
             EndIf
          EndIf

          If (Do_SLD_Output) Then
             If (iv .eq. 3) Then
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      dvars1(it,ip,ir,2) = r_ind(myr1+ir-1)*sines(myth1+it-1)*flux(ip,it)
                   EndDo
                EndDo
             ElseIf (iv .eq. 5) Then
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do ip=1,mynphi
                   !#DIR$ VECTOR ALIGNED
                   !#DIR$ IVDEP
                   Do it=1,mynth
                      dvars1(it,ip,ir,4) = flux(ip,it)
                   EndDo
                EndDo
             EndIf
          Else
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do ip=1,mynphi
                !#DIR$ VECTOR ALIGNED
                !#DIR$ IVDEP
                Do it=1,mynth
                   dvars1(it,ip,ir,iv) = dvars1(it,ip,ir,iv) + flux(ip,it)
                EndDo
             EndDo
          EndIf
       EndDo
    EndDo

  End Subroutine Diffusion_Th

  Subroutine Diffusion_Phi()
    Implicit None
    Real*8, Dimension(mynth,mynphi) :: tmp, flux
    Real*8, Dimension(mynth,mynphi+1) :: vtmp, delue, qmh, qph, qp3h, phiR, phiL
    Real*8 :: b1(mynth,4), b2(mynth,4), geom(mynth)
    Real*8 :: stmp, rtmp, ctmp
    Integer :: ir, it, ip, iv, isl

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: tmp, flux, vtmp, delue, b1, b2, qmh, qph, qp3h, phiR, phiL

    Do isl=1,nsld
       iv = varlist(isl)
       Do ir=1,mynr
          Do it=myth1,myth2
             If (it.eq.1) Then
                rtmp = theta(1)-0.5d0*dth
                stmp = cos(rtmp)
             Else
                rtmp = 0.5d0*(theta(it)+theta(it-1))
                stmp = cos(rtmp)
             EndIf

             If (it .eq. nth) Then
                ctmp = theta(it)+0.5d0*dth
                rtmp = ctmp - rtmp
                ctmp = cos(ctmp)
             Else
                ctmp = 0.5d0*(theta(it+1)+theta(it))
                rtmp = ctmp - rtmp
                ctmp = cos(ctmp)
             EndIf
             ctmp = rtmp/(stmp-ctmp)
             geom(it-myth1+1) = ctmp
          EndDo

          If ((myr1.eq.1).and.(ir.eq.1)) Then
             stmp = r_ind(1)-0.5d0*dr_ind(1)
          Else
             stmp = 0.5d0*(r_ind(myr1+ir-1)+r_ind(myr1+ir-2))
          EndIf
          If ((myr2.eq.nr).and.(ir.eq.mynr)) Then
             rtmp = r_ind(nr)+0.5d0*dr_ind(nr)
          Else
             rtmp = 0.5d0*(r_ind(myr1+ir)+r_ind(myr1+ir-1))
          EndIf
          rtmp = 1.5d0*(rtmp**2-stmp**2)/(rtmp**3-stmp**3)/dphi

          !#DIR$ VECTOR ALIGNED
          geom = geom*rtmp

          If (iv .gt. 5) Then
             !#DIR$ VECTOR ALIGNED
             tmp = maga(:,:,ir,iv-5)
             Do ip=1,4
                !#DIR$ VECTOR ALIGNED
                b1(:,ip) = dbinfo(3,1)%bdata(:,ip,ir,iv-5)
                !#DIR$ VECTOR ALIGNED
                b2(:,ip) = dbinfo(3,2)%bdata(:,ip,ir,iv-5)
             EndDo
          Else
             tmp = vars1(:,:,ir,iv)
             If (iv .eq. 1) Then
                Do ip=1,4
                   !#DIR$ VECTOR ALIGNED
                   b1(:,ip) = binfo(3,1)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,iv)
                   !#DIR$ VECTOR ALIGNED
                   b2(:,ip) = binfo(3,2)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,iv)
                EndDo
             ElseIf (iv .eq. 5) Then
                Do ip=1,4
                   !#DIR$ VECTOR ALIGNED
                   b1(:,ip) = binfo(3,1)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,2)
                   !#DIR$ VECTOR ALIGNED
                   b2(:,ip) = binfo(3,2)%bdata(:,ip,ir,iv)-meanV(myr1+ir-1,2)
                EndDo
             Else
                Do ip=1,4
                   !#DIR$ VECTOR ALIGNED
                   b1(:,ip) = binfo(3,1)%bdata(:,ip,ir,iv)
                   !#DIR$ VECTOR ALIGNED
                   b2(:,ip) = binfo(3,2)%bdata(:,ip,ir,iv)
                EndDo
             EndIf
          EndIf

          !Gradients
          !#DIR$ VECTOR ALIGNED
          qmh(:,1) = b1(:,4) - b1(:,3)
          qmh(:,2) = tmp(:,1) - b1(:,4)
          !#DIR$ VECTOR ALIGNED
          qmh(:,3:mynphi+1) = tmp(:,2:mynphi) - tmp(:,1:mynphi-1)

          !#DIR$ VECTOR ALIGNED
          qph(:,1:mynphi) = qmh(:,2:mynphi+1)
          qph(:,mynphi+1) = b2(:,1) - tmp(:,mynphi)

          !#DIR$ VECTOR ALIGNED
          qp3h(:,1:mynphi) = qph(:,2:mynphi+1)
          qp3h(:,mynphi+1) = b2(:,2) - b2(:,1)

          !Van Albada
          !#DIR$ VECTOR ALIGNED
          Do ip=1,mynphi+1
             Do it=1,mynth
                rtmp = abs(qmh(it,ip))
                rtmp = 1d0/max(rtmp,tinyvar(iv))
                rtmp = qmh(it,ip)*rtmp**2
                rtmp = rtmp*qph(it,ip)
                rtmp = min(1d0,rtmp) !rtmp*(rtmp+1d0)/(rtmp*rtmp+1d0)
                phiL(it,ip) = qmh(it,ip)*max(0d0,rtmp)

                rtmp = abs(qp3h(it,ip))
                rtmp = 1d0/max(rtmp,tinyvar(iv))
                rtmp = qp3h(it,ip)*rtmp**2
                rtmp = rtmp*qph(it,ip)
                rtmp = min(1d0,rtmp) !rtmp*(rtmp+1d0)/(rtmp*rtmp+1d0)
                phiR(it,ip) = qp3h(it,ip)*max(0d0,rtmp)
             EndDo
          EndDo

          !#DIR$ VECTOR ALIGNED
          delue = qph - 0.5d0*(phiR+phiL)

          If ((iv .lt. 5).and.(flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          ElseIf ((iv .eq. 5).and.(ent_flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**ent_flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          ElseIf ((iv .gt. 5).and.(mag_flux_limiter .gt. 0)) Then
             where (qph*delue .le. 0d0) 
                phiR = 0d0
             elsewhere
                phiR = (delue/qph)**mag_flux_limiter
             endwhere

             qph = phiR
             where(qph .gt. 1d0)
                phiR = 1d0
             endwhere
             delue = phiR*delue
          EndIf

          !#DIR$ VECTOR ALIGNED
          vtmp = sld_coef*sqrt(0.5d0*(vcc(1:mynth,0:mynphi,ir) + vcc(1:mynth,1:mynphi+1,ir)))

          !#DIR$ VECTOR ALIGNED
          delue = vtmp*delue

          !#DIR$ VECTOR ALIGNED
          flux = delue(:,2:mynphi+1)-delue(:,1:mynphi)

          If ((sld_pr.ne.1d0).or.(sld_prm.ne.1d0)) Then
             If (iv.lt.5) Then
                !#DIR$ VECTOR ALIGNED
                flux = sld_pr*flux
             ElseIf (iv .gt. 5) Then
                !#DIR$ VECTOR ALIGNED
                flux = sld_prm*flux
             EndIf
          EndIf

          !Compute Diffusive Flux:
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do ip=1,mynphi
             !#DIR$ VECTOR ALIGNED
             !#DIR$ IVDEP
             Do it=1,mynth
                dvars1(it,ip,ir,iv) = dvars1(it,ip,ir,iv) + geom(it)*flux(it,ip)
             EndDo
          EndDo
       EndDo
    EndDo
  End Subroutine Diffusion_Phi

  Function Compute_Density(ent) Result(rho)  
    Real*8, Intent(In) :: ent(nr)
    Real*8 :: matrix(nr,nr), rhs(nr), rho(nr), egrad(nr)
    Real*8 :: drho(nr), c1(nr), c2(nr), rho_temp(nr)
    Real*8 :: errvars, errfuns, tol, gam1, gam, bc1, bc2, drin
    Integer :: pivot(nr), iteration, rr, max_iter

    drin = Dble(nr-1)/(x_ind(nr)-x_ind(1))
    rhs = ent
    bc1 = rhs(1)/drin
    bc2 = rhs(nr)/drin
    egrad = 0d0
    Call dbyd1_dt0(rhs,egrad,bc1,bc2,nr,1)
    If (non_uniform) Then
       egrad = dxdr*egrad
    EndIf
    egrad = drin*egrad

    ! use the solar model interpolated onto our grid as an initial guess
    gam1 = gamma
    rho = rhor

    ! initialize a few constants
    gam = 1d0-gam1
    c1 = egrad/Cp
    c2 = (gravity / gam1) * exp( - gam1 * rhs / Cp)
    If (non_uniform) Then
       c1 = c1/dxdr
       c2 = c2/dxdr
    EndIf

    ! now iterate until convergence
    iteration = 1
    max_iter = 100
    tol = 1d-13
    errvars = 1d0
    errfuns = 1d0
    Do While ((iteration .le. max_iter) .and. ((errvars .gt. tol) .and. (errfuns .gt. tol))) 

       ! first set up the coefficient matrix
       matrix = 0d0
       matrix(1,1:7) = (/-147d0,360d0,-450d0,400d0,-225d0,72d0,-10d0/)
       matrix(2,1:7) = (/-10d0,-77d0,150d0,-100d0,50d0,-15d0,2d0/)
       matrix(3,1:7) = (/2d0,-24d0,-35d0,80d0,-30d0,8d0,-1d0/)
       Do rr=4,nr-3
          matrix(rr,rr-3:rr+3) = (/-1d0,9d0,-45d0,0d0,45d0,-9d0,1d0/)
       EndDo
       matrix(nr-2,nr-6:nr) = (/1d0,-8d0,30d0,-80d0,35d0,24d0,-2d0/)
       matrix(nr-1,nr-6:nr) = (/-2d0,15d0,-50d0,100d0,-150d0,77d0,10d0/)
       matrix(nr,nr-6:nr)   = (/10d0,-72d0,225d0,-400d0,450d0,-360d0,147d0/)
       matrix = matrix*drin/60d0

       ! now set up the right-hand side
       rho_temp = log(rho)
       drho = matmul(matrix,rho_temp)
       rhs = - ( drho + c1 + c2 * rho**gam )

       ! Add the diagonal terms
       Do rr = 1, nr
          matrix(rr,rr)= matrix(rr,rr) + gam*c2(rr)*rho(rr)**gam ! the diagonal
       Enddo

       ! This is the Boundary Condition
       matrix(nr,:) = 0d0
       matrix(nr,nr) = 1d0
       rhs(nr) = 0d0

       !matrix(1,:) = 0d0
       !matrix(1,1) = 1d0
       !rhs(1) = 0d0

       errfuns = sum(abs(rhs)) 
       errfuns = errfuns / sum(abs(drho))

       ! now solve for the variations,...           
       call lu_solve_1D(matrix, rhs, pivot)
       drho = rhs

       ! add in the variations
       rho = exp(log(rho) + drho)

       ! check for convergence
       errvars = sum(abs(drho))/(sum(abs(rho)))

       If (myrank .eq. 0) Print*, 'Compute_Density (vars,funcs) : ', iteration, errvars, errfuns

       iteration = iteration + 1

    End Do

  End Function Compute_Density

  Subroutine LU_Solve_1D(mat, rhs, pvt, na, nb)
    Real*8,intent(in) :: mat(:,:)
    Real*8,Intent(inout) :: rhs(:)
    Integer,intent(inout) :: pvt(:)
    Integer, Optional :: na, nb
    Integer :: n, info

    If (Present(na)) Then
       n = na
    Else
       n = Size(mat,1)
    End If
    
    Call dgesv(n, 1, mat, Size(mat,1), pvt, rhs, Size(rhs,1), info)
    
  End Subroutine LU_Solve_1D

End Module Diffusion
