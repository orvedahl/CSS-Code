Module Report
  
  Use Checkpoint
  Use Diffusion
  Implicit None
  
  Character(256) :: az_avg_name, shell_avg_name, shell_slice_name, scalar_filename, az_avg_filename, shell_avg_filename, shell_slice_filename
  !Number of quantities, radii and current record counter
  Integer :: az_qs, shell_avg_qs, shell_slice_qs, shell_slice_rads, az_record=0, sh_avg_record=0, sh_slice_record=0 
  !Integer Record lengths
  Integer :: irec_azavg, irec_shlavg, irec_shlslice
  !Integer quantities selecting output quantities
  Integer, Allocatable, Dimension(:) :: az_quantities, shell_avg_quantities, shell_slice_quantities, shell_slice_radii 
  !Node look up table for shell slices
  Integer, Allocatable, Dimension(:,:) :: ilupr 
  !Flags to report on particular quantities
  Logical :: do_ke, do_me, do_bpol, do_mach 
  !Derivative work array for angular momentum output
  Real*8 :: lrad_radius, lrad_tanh_coef
  Real*8, Target, Allocatable, Dimension(:,:,:,:) :: derivs
  Real*8, Allocatable, Dimension(:) :: avgP, avgT, avgS, avgrho

Contains

  Subroutine Initialize_IO()
    checkpoint_name = Trim(perm_dir)//'Checkpoints/' 
    az_avg_name = Trim(perm_dir)//'AZ_Avgs/az_avg_'
    shell_avg_name = Trim(perm_dir)//'Shell_Avgs/shell_avg_'
    shell_slice_name = Trim(perm_dir)//'Shell_Slices/shell_slice_'
    scalar_filename = Trim(perm_dir)//'scalar.dat'
    Call Make_Inverse_Lookup(shell_slice_radii)
    If (do_radheat) Then
       If (Do_SLD_Output) Then
          Allocate(derivs(mynth,mynphi,mynr,8))
       Else
          Allocate(derivs(mynth,mynphi,mynr,4))
       EndIf
    Else
       If (Do_SLD_Output) Then
          Allocate(derivs(mynth,mynphi,mynr,7))
       Else
          Allocate(derivs(mynth,mynphi,mynr,3))
       EndIf
    EndIf
    Allocate(avgT(nr),avgrho(nr),avgP(nr),avgS(nr))
    irec_azavg = -1
    irec_shlavg = -1
    irec_shlslice = -1
  End Subroutine Initialize_IO
  
  Subroutine Finalize_IO()
    Deallocate(az_quantities,shell_avg_quantities,shell_slice_quantities,shell_slice_radii,ilupr,derivs,avgT,avgrho,avgP,avgS)
  End Subroutine Finalize_IO
  
  Subroutine Make_Filename(filename,fout)
    Character(256), Intent(In) :: filename
    Character(256), Intent(Out) :: fout
    Character(10) :: temp, pretemp
    
    pretemp = ' '
    Write(temp,fmt='(i7.7)') istep
    fout = Trim(filename)//Trim(pretemp)//Trim(temp)
    fout = Trim(fout)

  End Subroutine Make_Filename
  
  !Check output frequencies
  Subroutine Output_Check(checkpoint_start)
    Character(256) :: temp,temph
    Logical :: checkpoint_start
    
    If (checkpoint_start) Then
       !Shell_Avgs
       If ((istep .eq. 0) .or. (Mod(istep,shell_avg_freq*shell_avg_recnum) .eq. 0)) Then
          Call Make_Filename(shell_avg_name,temp)
          shell_avg_filename = Trim(temp)
          sh_avg_record=1 
       End If
       
       If (Mod(istep,shell_avg_freq) .eq. 0) Then
          If (magnetic) Call Compute_B(0)
          Call Compute_Derivatives()
          Call Shell_Average_Output(shell_avg_quantities)
          sh_avg_record=sh_avg_record+1 
       End If
       
       !AZ_Avgs
       If ((istep .eq. 0) .or. (Mod(istep,az_avg_freq*az_avg_recnum) .eq. 0)) Then             
          Call Make_Filename(az_avg_name,temp)
          az_avg_filename = Trim(temp)
          az_record=1
       End If
       
       If (Mod(istep,az_avg_freq) .eq. 0) Then
          Call Azimuthal_Average_Output(az_quantities)
          az_record=az_record+1
       End If
       
       !Shell_Slices
       If ((istep .eq. 0) .or. (Mod(istep,shell_slice_freq*shell_slice_recnum) .eq. 0)) Then
          Call Make_Filename(shell_slice_name,temp)
          shell_slice_filename = Trim(temp)
          sh_slice_record=1
       End If
       
       If (Mod(istep,shell_slice_freq) .eq. 0) Then
          Call Shell_Slice_Output(shell_slice_quantities,shell_slice_radii,0)
          sh_slice_record=sh_slice_record+1
       End If

       If ((Mod(istep,shell_avg_freq) .eq. 0).and.magnetic) Then
          vars0(:,:,:,6:8) = maga
       EndIf
    Else          
       !Scalar
       If ((istep .eq. 0) .or. (Mod(istep,scalar_freq) .eq. 0)) Then
          If (magnetic) Call Compute_B(0)
          Call Scalar_Output()
       End If
       
       !Shell_Avgs
       If ((istep .eq. 0) .or. (Mod(istep,shell_avg_freq*shell_avg_recnum) .eq. 0)) Then
          Call Make_Filename(shell_avg_name,temp)
          shell_avg_filename = Trim(temp)
          sh_avg_record=1 
       End If
       
       If (Mod(istep,shell_avg_freq) .eq. 0) Then
          Call Compute_Derivatives()
          Call Shell_Average_Output(shell_avg_quantities)
          sh_avg_record=sh_avg_record+1 
       End If
       
       !AZ_Avgs
       If ((istep .eq. 0) .or. (Mod(istep,az_avg_freq*az_avg_recnum) .eq. 0)) Then
          Call Make_Filename(az_avg_name,temp)
          az_avg_filename = Trim(temp)
          az_record=1
       End If
       
       If (Mod(istep,az_avg_freq) .eq. 0) Then
          Call Azimuthal_Average_Output(az_quantities)
          az_record=az_record+1
       End If
       
       !Shell_Slices
       If ((istep .eq. 0) .or. (Mod(istep,shell_slice_freq*shell_slice_recnum) .eq. 0)) Then
          Call Make_Filename(shell_slice_name,temp)
          shell_slice_filename = Trim(temp)
          sh_slice_record=1
       End If
       
       If (Mod(istep,shell_slice_freq) .eq. 0) Then
          Call Shell_Slice_Output(shell_slice_quantities,shell_slice_radii,0)
          sh_slice_record=sh_slice_record+1
       End If

       If (((istep .eq. 0) .or. (Mod(istep,scalar_freq) .eq. 0)) .and. magnetic) Then
          vars0(:,:,:,6:8) = maga
       EndIf

       !Checkpoint
       If ((istep .eq. 0) .or. (Mod(istep,checkpoint_freq) .eq. 0)) Then
          Call Make_Checkpoint_Filename(checkpoint_name,temp,temph)
          checkpoint_filename = Trim(temp)
          checkpoint_header_filename = Trim(temph)
          Call Checkpoint_Output()
       End If
    EndIf
    
  End Subroutine Output_Check
  
  Subroutine Scalar_Output()
    Real*8, Pointer, Dimension(:,:,:) :: prho,pw,pv,pu,ps,pbr,pbt,pbp
    Integer :: nsize
    Real*8 :: umach,kinener,cke,mcke,drke,magener,bmax,snd(5),rcv(5)
    Real*8, Allocatable, Dimension(:,:,:) :: usqr,bsqr
    Integer :: scalar_unit=10   

    !  useful aliases
    Nullify(prho,pw,pv,pu,ps,pbr,pbt,pbp)
    prho => vars0(:,:,:,1)
    pw   => vars0(:,:,:,2)
    pv   => vars0(:,:,:,3)
    pu   => vars0(:,:,:,4)
    ps   => vars0(:,:,:,5)
    If (magnetic) Then
       pbt => vars0(:,:,:,6)
       pbp => vars0(:,:,:,7)
       pbr => vars0(:,:,:,8)
    Endif
    
    nsize=nr*nth*nphi
    
    !  compute total ke
    Allocate(usqr(mynth,mynphi,mynr))
    If (magnetic) Then
       !  compute magnetic energy here
       Allocate(bsqr(mynth,mynphi,mynr))  !  bsqr holds B squared
       bsqr=pbr**2+pbt**2+pbp**2
       magener=Sum(bsqr)/(8*Pi*nsize)
    Else
       magener = 0d0
    Endif

    If (Unstratified) Then
       temperature = g0
    Else
       temperature = (prho**gamm1)*dexp(ps*invCv)*invtconst*gamm1*Cp
    EndIf
    usqr = pu**2
    cke = Sum(prho*usqr)/Dble(2*nsize)
    usqr = pw**2
    mcke = cke+Sum(prho*usqr)/Dble(2*nsize)
    drke = Sum(prho*pv**2)/Dble(2*nsize)
    usqr = usqr+pv**2+pu**2
    If (magnetic) Then
       umach = Sqrt(Maxval(usqr/(2d0*bsqr/prho+usqr+temperature)))
    Else
       umach = Sqrt(Maxval(usqr/(usqr+temperature)))
    EndIf
    kinener = Sum(prho*usqr)/Dble(2*nsize)
    Deallocate(usqr)

    !  reduce to get global values of above
    snd = 0d0
    snd(1) = umach
    If (magnetic) Then
       snd(2) = Maxval(sqrt(bsqr))*2d0*sqrt(2d0*pi)
    Else
       snd(2) = 0d0
    EndIf
    rcv    = 0d0
    Call Global_Reduce(snd,rcv,0,node_comm,'max')
    umach = rcv(1)
    bmax  = rcv(2)
    
    rcv = 0d0
    snd(1) = kinener
    snd(2) = cke
    snd(3) = mcke
    snd(4) = drke
    snd(5) = magener
    Call Global_Reduce(snd,rcv,0,node_comm,'sum')
    kinener=rcv(1)
    cke=rcv(2)
    mcke=rcv(3)
    drke=rcv(4)
    magener = rcv(5)
    
    If (myrank.eq.0) Then
       Open(scalar_unit,file=trim(scalar_filename),status='unknown',access='sequential',position='append')
       If (magnetic) Then
          Write(scalar_unit,fmt='(i10,9d15.7)') istep,sim_time,deltat,kinener,cke,mcke,drke,umach,magener,bmax
       Else
          Write(scalar_unit,fmt='(i10,7d15.7)') istep,sim_time,deltat,kinener,cke,mcke,drke,umach
       End If
       Close(scalar_unit)
    Endif
    If (Allocated(bsqr)) Deallocate(bsqr)
    Nullify(prho,pw,pv,pu,ps,pbr,pbt,pbp)
  End Subroutine Scalar_Output
  
  Subroutine Compute_Derivatives()
    Implicit None
    Real*8, Allocatable, Dimension(:,:) :: rho,s,v,dsdr,drvdr,work1,work2
    Real*8 :: workrb1(mynth,binfo(1,1)%bsize%bnr)
    Real*8 :: workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: worktb1(binfo(2,1)%bsize%bnt,mynphi)
    Real*8 :: worktb2(binfo(2,2)%bsize%bnt,mynphi)
    Real*8 :: tempv
    Integer :: iq, nqs, ir, it, ip, iv

    derivs = 0d0
    vars1 = vars0

    !Exchange boundary info
    Call Boundary_Exchange_Init(1)

    !Do radial derivatives
    Allocate(rho(mynth,mynr),v(mynth,mynr),s(mynth,mynr),dsdr(mynth,mynr),drvdr(mynth,mynr),work2(mynth,mynr))
    If (non_uniform) Then
       Allocate(work1(mynth,mynr))
       Do ip=1,mynth
          work1(ip,:) = dxdr(myr1:myr2)
       EndDo
    EndIf

    Call Boundary_Exchange_Finish(1)
    Call Boundary_Exchange_Init(2)

    Do ip=1,mynphi
       rho = vars0(:,ip,:,1)
       v   = vars0(:,ip,:,3)
       s   = vars0(:,ip,:,5)

       If (myr1 .eq. 1) Then
          If ((lbcr_s .eq. 2).or.(lbcr_s.eq.4)) Then
             workrb1(:,1) = dsdr1(:)/dxdr(1)
          Else
             workrb1(:,1) = s(:,1)
          EndIf
       Else
          workrb1 = binfo(1,1)%bdata(:,ip,:,5)
       EndIf

       If (myr2 .eq. nr) Then
          If ((lbcr_s .eq. 3).or.(lbcr_s.eq.4)) Then
             workrb2(:,1) = dsdr2(:)/dxdr(nr)
          Else
             workrb2(:,1) = s(:,mynr)
          EndIf
       Else
          workrb2 = binfo(1,2)%bdata(:,ip,:,5)
       EndIf

       Call dbydr(s,dsdr,workrb1,workrb2,5)
       If (non_uniform) Then
          dsdr = dsdr*work1
       EndIf

       If (myr1 .eq. 1) Then
          workrb1(:,1) = rho(:,1)*v(:,1)/r1/dxdr(1)
       Else
          workrb1 = binfo(1,1)%bdata(:,ip,:,1)*binfo(1,1)%bdata(:,ip,:,3)
       EndIf

       If (myr2 .eq. nr) Then
          workrb2(:,1) = rho(:,mynr)*v(:,mynr)/r2/dxdr(nr)
       Else
          workrb2 = binfo(1,2)%bdata(:,ip,:,1)*binfo(1,2)%bdata(:,ip,:,3)
       EndIf
       work2 = rho*v
       Call dbydr(work2,drvdr,workrb1,workrb2,3)
       If (non_uniform) Then
          drvdr = drvdr*work1
       EndIf
       derivs(:,ip,:,1) = dsdr
       derivs(:,ip,:,2) = drvdr
    EndDo
    Deallocate(rho,v,s,dsdr,drvdr,work2)

    If (non_uniform) Deallocate(work1)

    Call Boundary_Exchange_Finish(2)

    Allocate(rho(mynth,mynphi),v(mynth,mynphi),drvdr(mynth,mynphi),work2(mynth,mynphi))
    Do ir=1,mynr
       rho = vars0(:,:,ir,1)
       v   = vars0(:,:,ir,3)
       If (myth1 .eq. 1) Then
          worktb1(1,:) = cosines(1)*rho(1,:)*v(1,:)/sines(1)
       Else
          worktb1 = binfo(2,1)%bdata(:,:,ir,1)*binfo(2,1)%bdata(:,:,ir,3)
       EndIf
       If (myth2 .eq. nth) Then
          worktb2(1,:) = cosines(nth)*rho(mynth,:)*v(mynth,:)/sines(nth)
       Else
          worktb2 = binfo(2,2)%bdata(:,:,ir,1)*binfo(2,2)%bdata(:,:,ir,3)
       EndIf
       work2 = rho*v
       Call dbydt(work2,drvdr,worktb1,worktb2,3)
       derivs(:,:,ir,3) = drvdr
    EndDo
    Deallocate(rho,v,drvdr,work2)

    nqs = size(az_quantities)
    Do iq=1,nqs
       If ((az_quantities(iq) .eq. 18).or.(az_quantities(iq).eq.19).or.(az_quantities(iq).eq.24).or.(az_quantities(iq).eq.25)) Then
          Do_SLD_Output = .True.
       EndIf
    EndDo

    nqs = size(shell_avg_quantities)
    Do iq=1,nqs
       If (shell_avg_quantities(iq).eq.19) Then
          Do_SLD_Output = .True.
       EndIf
    EndDo

    If (Do_SLD_Output) Then
       dvars1 = 0d0
       Call Compute_Diffusion(2)
       If (do_radheat) Then
          derivs(:,:,:,5) = dvars1(:,:,:,1)
          derivs(:,:,:,6) = dvars1(:,:,:,2)
          derivs(:,:,:,7) = dvars1(:,:,:,3)
          derivs(:,:,:,8) = dvars1(:,:,:,4)
       Else
          derivs(:,:,:,4) = dvars1(:,:,:,1)
          derivs(:,:,:,5) = dvars1(:,:,:,2)
          derivs(:,:,:,6) = dvars1(:,:,:,3)
          derivs(:,:,:,7) = dvars1(:,:,:,4)
       EndIf
       Do_SLD_Output = .False.
    EndIf

    temperature = (vars0(:,:,:,1)**gamm1)*dexp(vars0(:,:,:,5)*invCv)*invtconst

    If (do_radheat) Then
       Allocate(dsdr(mynth,mynr),drvdr(mynth,mynr),work2(mynth,mynr))
       If (non_uniform) Then
          Allocate(work1(mynth,mynr))
          Do ip=1,mynth
             work1(ip,:) = dxdr(myr1:myr2)
          EndDo
       EndIf
       Do ip=1,mynphi
          !work2 = vars0(:,ip,:,1)
          !workrb1 = binfo(1,1)%bdata(:,ip,:,1)
          !workrb2 = binfo(1,2)%bdata(:,ip,:,1)
          !Call dbydr(work2,dsdr,workrb1,workrb2,-1)
          !If (non_uniform) Then
          !   dsdr = dsdr*work1
          !EndIf

          If (myr1 .eq. 1) Then
             workrb1(:,1) = -3d0*vars0(:,ip,1,1)*lumr1/(16d0*pi*a_const*c_light*(r1**2)*kr(1)*temperature(:,ip,1)**3)/dxdr(1) !*(gamm1*dsdr(:,1)+invCv*derivs(:,ip,1,1))/dxdr(1)
          Else
             workrb1 = invtconst*(binfo(1,1)%bdata(:,ip,:,1)**gamm1)*dexp(invCv*binfo(1,1)%bdata(:,ip,:,5))
          EndIf
          If (myr2 .eq. nr) Then
             workrb2(:,1) =  -3d0*vars0(:,ip,mynr,1)*lumr2/(16d0*pi*a_const*c_light*(r2**2)*kr(nr)*temperature(:,ip,mynr)**3)/dxdr(nr) !*(gamm1*dsdr(:,mynr)+invCv*derivs(:,ip,mynr,1))/dxdr(nr)
          Else
             workrb2 = invtconst*(binfo(1,2)%bdata(:,ip,:,1)**gamm1)*dexp(invCv*binfo(1,2)%bdata(:,ip,:,5))
          EndIf

          work2 = temperature(:,ip,:)
          Call dbydr(work2,drvdr,workrb1,workrb2,6)
          If (non_uniform) Then
             drvdr = drvdr*work1
          EndIf
          derivs(:,ip,:,4) = drvdr
       EndDo

       Deallocate(dsdr,drvdr,work2)       
       If (non_uniform) Deallocate(work1)
    EndIf
  End Subroutine Compute_Derivatives

  Subroutine Azimuthal_Average_Output(quantities)
    Implicit None
    Real*8, Pointer, Dimension(:,:,:) :: prho,pw,pv,pu,ps,pbr,pbt,pbp
    Real*8, Allocatable, Dimension(:,:,:) :: averages,usqr,bsqr,bpol,buff,rbuff, avgvs
    Real*8 :: rtmp
    Integer, Intent(In) :: quantities(:)
    Integer :: iq, ir, it, nq, az_avg_unit=20    

    !  useful aliases
    Nullify(prho,pw,pv,pu,ps,pbr,pbt,pbp)
    prho => vars0(:,:,:,1)
    pw   => vars0(:,:,:,2)
    pv   => vars0(:,:,:,3)
    pu   => vars0(:,:,:,4)
    ps   => vars0(:,:,:,5)
    If (magnetic) Then
       Nullify(pbr,pbt,pbp)
       pbt => vars0(:,:,:,6)
       pbp => vars0(:,:,:,7)
       pbr => vars0(:,:,:,8)
    End If

    nq = size(quantities)

    Allocate(buff(mynth,mynr,nq))

    Allocate(averages(nth,nr+2,nq))
    averages = 0d0
    If (irec_azavg .eq. -1) Then
       Inquire(iolength=irec_azavg) averages
    EndIf

    If (magnetic) Then
       Allocate(avgvs(nth,nr,7),rbuff(nth,nr,7))
    Else
       Allocate(avgvs(nth,nr,4),rbuff(nth,nr,4))
    EndIf

    avgvs = 0d0
    Do ir=1,mynr
       Do it=1,mynth
          avgvs(myth1+it-1,myr1+ir-1,1) = Sum(prho(it,:,ir))/Dble(nphi)
          avgvs(myth1+it-1,myr1+ir-1,2) = Sum(pw(it,:,ir))/Dble(nphi)
          avgvs(myth1+it-1,myr1+ir-1,3) = Sum(pv(it,:,ir))/Dble(nphi)
          avgvs(myth1+it-1,myr1+ir-1,4) = Sum(pu(it,:,ir))/Dble(nphi)
       EndDo
    EndDo

    If (magnetic) Then
       Do ir=1,mynr
          Do it=1,mynth
             avgvs(myth1+it-1,myr1+ir-1,5) = Sum(pbt(it,:,ir))/Dble(nphi)
             avgvs(myth1+it-1,myr1+ir-1,6) = Sum(pbp(it,:,ir))/Dble(nphi)
             avgvs(myth1+it-1,myr1+ir-1,7) = Sum(pbr(it,:,ir))/Dble(nphi)
          EndDo
       EndDo
    EndIf

    Call Global_Allreduce(avgvs,rbuff,node_comm,'sum')
    avgvs = rbuff
    Deallocate(rbuff)

    !Kinetic Energy
    If (do_ke) Then
       Allocate(usqr(mynth,mynphi,mynr))
       usqr=pu**2+pv**2+pw**2
       usqr=prho*usqr/2
    End If

    If (magnetic) Then
       !Magnetic Energy
       If (do_me) Then
          Allocate(bsqr(mynth,mynphi,mynr))
          bsqr=pbr**2+pbt**2+pbp**2
          bsqr=bsqr/pi/8
       End If

       !Poloidal Field
       If (do_bpol) Then
          Allocate(bpol(mynth,mynphi,mynr))
          bpol = pbr**2+pbt**2
          bpol = Sqrt(bpol)
       End If
    End If

    Do iq=1,nq
       Select Case (quantities(iq))
       Case (1)
          !  Rho
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(prho(it,:,ir)-avgrho(myr1+ir-1))/Dble(nphi)
             EndDo
          EndDo
       Case (2)
          !  Temperature
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(temperature(it,:,ir)-avgT(myr1+ir-1))/Dble(nphi)
             EndDo
          EndDo
       Case (3)
          !  Entropy
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(ps(it,:,ir)-avgS(myr1+ir-1))/Dble(nphi)
             EndDo
          EndDo
       Case (4) 
          ! Pressure
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(gamm1*Cv*prho(it,:,ir)*temperature(it,:,ir)-avgP(myr1+ir-1))/Dble(nphi)
             EndDo
          EndDo
       Case (5)
          !  Vr
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = avgvs(myth1+it-1,myr1+ir-1,4)
             EndDo
          EndDo
       Case (6)
          !  Vtheta
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = avgvs(myth1+it-1,myr1+ir-1,2)
             EndDo
          EndDo
       Case (7)
          !  Vphi
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = avgvs(myth1+it-1,myr1+ir-1,3)
             EndDo
          EndDo
       Case (8)
          !  Kinetic Energy
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(usqr(it,:,ir))/Dble(nphi)
             EndDo
          EndDo
       Case (9)
          !  Br
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = avgvs(myth1+it-1,myr1+ir-1,7)
             EndDo
          EndDo
       Case (10)
          !  Btheta
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = avgvs(myth1+it-1,myr1+ir-1,5)
             EndDo
          EndDo
       Case (11)
          !  Bphi
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = avgvs(myth1+it-1,myr1+ir-1,6)
             End Do
          End Do
       Case (12)
          !  Magnetic Energy
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(bsqr(it,:,ir))/Dble(nphi)
             End Do
          End Do
       Case (13)
          !  Poloidal Field Magnitude
          Do ir=1,mynr
             Do it=1,mynth
                buff(it,ir,iq) = Sum(bpol(it,:,ir))/Dble(nphi)
             End Do
          End Do
       Case (14)
          ! Radial MC Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)
                buff(it,ir,iq) = rtmp*Sum(prho(it,:,ir)*pu(it,:,ir))*(avgvs(myth1+it-1,myr1+ir-1,3)+rtmp*omega_0)/Dble(nphi)
             EndDo
          EndDo
       Case (15)
          ! Theta MC Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)
                buff(it,ir,iq) = rtmp*Sum(prho(it,:,ir)*pw(it,:,ir))*(avgvs(myth1+it-1,myr1+ir-1,3)+rtmp*omega_0)/Dble(nphi)
             EndDo
          EndDo
       Case (16)
          ! Radial RS Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)
                buff(it,ir,iq) = rtmp*avgvs(myth1+it-1,myr1+ir-1,1)*Sum((pu(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,4))*(pv(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,3)))/Dble(nphi)
             EndDo
          EndDo
       Case (17)
          ! Theta RS Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)
                buff(it,ir,iq) = rtmp*avgvs(myth1+it-1,myr1+ir-1,1)*Sum((pw(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,2))*(pv(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,3)))/Dble(nphi)
             EndDo
          EndDo
       Case (18)
          ! Radial Diff Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                If (do_radheat) Then
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,5))/Dble(nphi)
                Else
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,4))/Dble(nphi)
                EndIf
             EndDo
          EndDo
       Case (19)
          ! Theta Diff Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                If (do_radheat) Then
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,6))/Dble(nphi)
                Else
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,5))/Dble(nphi)
                EndIf
             EndDo
          EndDo
       Case (20)
          ! Radial Mean B Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)*one_over_4pi
                buff(it,ir,iq) = rtmp*avgvs(myth1+it-1,myr1+ir-1,7)*avgvs(myth1+it-1,myr1+ir-1,6)
             EndDo
          EndDo
       Case (21)
          ! Theta Mean B Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)*one_over_4pi
                buff(it,ir,iq) = rtmp*avgvs(myth1+it-1,myr1+ir-1,5)*avgvs(myth1+it-1,myr1+ir-1,6)
             EndDo
          EndDo
       Case (22)
          ! Radial MS Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)*one_over_4pi
                buff(it,ir,iq) = rtmp*Sum((pbr(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,7))*(pbp(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,6)))
             EndDo
          EndDo
       Case (23)
          ! Theta MS Amom Flux
          Do ir=1,mynr
             Do it=1,mynth
                rtmp = r_ind(myr1+ir-1)*sines(myth1+it-1)*one_over_4pi
                buff(it,ir,iq) = rtmp*Sum((pbt(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,5))*(pbp(it,:,ir)-avgvs(myth1+it-1,myr1+ir-1,6)))
             EndDo
          EndDo
       Case (24)
          ! Radial Entropy Diff
          Do ir=1,mynr
             Do it=1,mynth
                If (do_radheat) Then
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,7))/Dble(nphi)
                Else
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,6))/Dble(nphi)
                EndIf
             EndDo
          EndDo
       Case (25)
          ! Theta Entropy Diff
          Do ir=1,mynr
             Do it=1,mynth
                If (do_radheat) Then
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,8))/Dble(nphi)
                Else
                   buff(it,ir,iq) = Sum(derivs(it,:,ir,7))/Dble(nphi)
                EndIf
             EndDo
          EndDo
       Case Default
          ! Error
          Write(*,*) 'Unrecognized Quantity In Azimuthal_Average_Output'
       EndSelect
    End Do

    averages(myth1:myth2,myr1+2:myr2+2,:) = buff
    Allocate(rbuff(nth,nr+2,nq))
    rbuff = 0d0

    Call MPI_Reduce(averages(1,1,1),rbuff(1,1,1),nth*(nr+2)*nq,MPI_REAL8,MPI_SUM,nnodes/2,node_comm,ierr)

    If (myrank .eq. nnodes/2) Then
       averages = rbuff

       averages(1,1,1) = nth
       averages(2,1,1) = nr
       averages(3,1,1) = nq
       averages(4,1,1) = r1
       averages(5,1,1) = r2
       averages(6,1,1) = th1
       averages(7,1,1) = th2
       averages(8:nq+7,1,1) = quantities
       averages(1:nr,2,1) = r_ind

       !If (nth .gt. nr) Then
       !Else !assumes nr<2*nth-nq-8
       !   averages(1:nth,2,1) = r_ind(1:nth)
       !   averages(nq+8:nq+7+nr-nth,1,1) = r_ind(nth+1:nr)
       !EndIf

       !Write out averages
       Open(az_avg_unit, file=Trim(az_avg_filename), form='unformatted', access='direct', status='unknown',recl=irec_azavg)
       Write(az_avg_unit, rec=az_record) averages
       Close(az_avg_unit)
    EndIf

    If (Allocated(averages)) Deallocate(averages)
    If (Allocated(rbuff)) Deallocate(rbuff)
    If (Allocated(buff)) Deallocate(buff)
    If (Allocated(avgvs)) Deallocate(avgvs)

    If (do_ke) Then
       Deallocate(usqr)
    EndIf
    If (do_me) Then
       Deallocate(bsqr)
    EndIf
    If (do_bpol) Then
       Deallocate(bpol)
    EndIf

    Nullify(prho,pw,pv,pu,ps)
    If (magnetic) Nullify(pbr,pbt,pbp)
  End Subroutine Azimuthal_Average_Output
  
  Subroutine Shell_Average_Output(quantities)
    Implicit None
    Real*8, Pointer, Dimension(:,:,:) :: prho,pw,pv,pu,ps,pbr,pbt,pbp
    Real*8, Allocatable, Dimension(:,:,:) :: usqr,bsqr,bpol
    Real*8, Allocatable, Dimension(:,:) :: averages,buff,rbuff
    Real*8, Allocatable, Dimension(:) :: dtavg,avgvr,rcv,snd
    Real*8 :: bc1, bc2, drt
    Integer, Intent(In) :: quantities(:)
    Integer :: iq,ir,it,nsize,nq,shell_avg_unit=30
    
    !  useful aliases
    Nullify(prho,pw,pv,pu,ps)
    prho => vars0(:,:,:,1)
    pw   => vars0(:,:,:,2)
    pv   => vars0(:,:,:,3)
    pu   => vars0(:,:,:,4)
    ps   => vars0(:,:,:,5)
    If (magnetic) Then
       Nullify(pbr,pbt,pbp)
       pbt => vars0(:,:,:,6)
       pbp => vars0(:,:,:,7)
       pbr => vars0(:,:,:,8)
    End If

    nq = size(quantities)
    nsize=nphi*nth
    
    Allocate(buff(mynr,nq))
    
    Allocate(dtavg(nr),avgvr(nr),snd(5*nr),rcv(5*nr))
    avgT = 0d0
    avgP = 0d0
    avgrho = 0d0
    avgS = 0d0
    avgvr = 0d0
    Do ir=1,mynr
       avgS(myr1+ir-1) = Sum(ps(:,:,ir))
       avgP(myr1+ir-1) = gamm1*Cv*Sum(prho(:,:,ir)*temperature(:,:,ir))
       avgT(myr1+ir-1) = Sum(temperature(:,:,ir))
       avgrho(myr1+ir-1) = Sum(prho(:,:,ir))
       avgvr(myr1+ir-1) = Sum(pu(:,:,ir))
    EndDo
    rcv = 0d0
    snd(1:nr) = avgT
    snd(nr+1:2*nr) = avgrho
    snd(2*nr+1:3*nr) = avgP
    snd(3*nr+1:4*nr) = avgvr
    snd(4*nr+1:5*nr) = avgS
    Call Global_AllReduce(snd,rcv,node_comm,'sum')
    rcv = rcv/Dble(nth*nphi)
    avgT = rcv(1:nr)
    avgrho = rcv(nr+1:2*nr)
    avgP = rcv(2*nr+1:3*nr)
    avgvr = rcv(3*nr+1:4*nr)
    avgS = rcv(4*nr+1:5*nr)

    drt = Dble(nr-1)/(x_ind(nr)-x_ind(1))
    bc1 = avgT(1)/drt
    bc2 = avgT(nr)/drt
    Call dbyd1_dt0(avgT,dtavg,bc1,bc2,nr,1)
    dtavg = drt*dtavg
    If (non_uniform) dtavg = dtavg*dxdr

    Allocate(averages(nr,nq+2))
    averages = 0d0
    If (irec_shlavg .eq. -1) Then
       Inquire(iolength=irec_shlavg) averages
    EndIf
    
    !Kinetic Energy

    If (do_ke) Then
       Allocate(usqr(mynth,mynphi,mynr))
       usqr=pu**2+pv**2+pw**2
       usqr=0.5d0*prho*usqr
    End If
    
    If (magnetic) Then
       !Magnetic Energy
       If (do_me) Then
          Allocate(bsqr(mynth,mynphi,mynr))
          bsqr=pbr**2+pbt**2+pbp**2
          bsqr=0.125d0*bsqr/pi
       End If
       
       !Poloidal Field
       If (do_bpol) Then
          Allocate(bpol(mynth,mynphi,mynr))
          bpol = pbr**2+pbt**2
          bpol = Sqrt(bpol)
       End If
    End If
    
    Do iq=1,nq
       Select Case (quantities(iq))
       Case (1) !  Rho
          Do ir=1,mynr
             buff(ir,iq) = Sum(prho(:,:,ir))/Dble(nsize)
          EndDo
       Case (2) !  Temperature
          Do ir=1,mynr
             buff(ir,iq) = Sum(temperature(:,:,ir))/Dble(nsize)
          EndDo
       Case (3) !  Entropy
          Do ir=1,mynr
             buff(ir,iq) = Sum(ps(:,:,ir))/Dble(nsize)
          EndDo
       Case (4) ! Pressure
          Do ir=1,mynr
             buff(ir,iq) = gamm1*Cv*Sum(prho(:,:,ir)*temperature(:,:,ir))/Dble(nsize)
          EndDo
       Case (5) !  Vr
          Do ir=1,mynr
             buff(ir,iq) = avgvr(myr1+ir-1)
          EndDo
       Case (6) !  Vtheta
          Do ir=1,mynr
             buff(ir,iq) = Sum(pw(:,:,ir))/Dble(nsize)
          EndDo
       Case (7) !  Vphi
          Do ir=1,mynr
             buff(ir,iq) = Sum(pv(:,:,ir))/Dble(nsize)
          EndDo
       Case (8) !  Kinetic Energy
          Do ir=1,mynr
             buff(ir,iq) = Sum(usqr(:,:,ir))/Dble(nsize)
          EndDo
       Case (9) !  ds/dr
          Do ir=1,mynr
             buff(ir,iq) = Sum(derivs(:,:,ir,1))/Dble(nsize)
          EndDo
       Case (10) !  Br
          Do ir=1,mynr
             buff(ir,iq) = Sum(pbr(:,:,ir))/Dble(nsize)
          EndDo
       Case (11) !  Btheta
          Do ir=1,mynr
             buff(ir,iq) = Sum(pbt(:,:,ir))/Dble(nsize)
          EndDo
       Case (12) !  Bphi
          Do ir=1,mynr
             buff(ir,iq) = Sum(pbp(:,:,ir))/Dble(nsize)
          EndDo
       Case (13) !  Magnetic Energy
          Do ir=1,mynr
             buff(ir,iq) = Sum(bsqr(:,:,ir))/Dble(nsize)
          EndDo
       Case (14) !  Poloidal Field Magnitude
          Do ir=1,mynr
             buff(ir,iq) = Sum(bpol(:,:,ir))/Dble(nsize)
          EndDo
       Case (15) ! F_ks (kappa0 and kappar terms)
          Do ir=1,mynr
             If (do_radheat) Then
                buff(ir,iq) = -(4d0*a_const*c_light*kr(myr1+ir-1)/3d0)*Sum((temperature(:,:,ir)**3)*derivs(:,:,ir,4)/prho(:,:,ir))/Dble(nsize)
                If ((Laplacian).or.(SLaplacian)) Then
                   buff(ir,iq) = buff(ir,iq)-ks(myr1+ir-1)*Sum(derivs(:,:,ir,1))/Dble(nsize)
                EndIf
             Else
                If ((Laplacian).or.(SLaplacian)) Then
                   buff(ir,iq) = -ks(myr1+ir-1)*Sum(derivs(:,:,ir,1))/Dble(nsize)+k0(myr1+ir-1)/Dble(nnp*nnt)
                Else
                   buff(ir,iq) = k0(myr1+ir-1)/Dble(nnp*nnt)
                EndIf
             EndIf
          EndDo
       Case (16) !F_en (enthalpy flux)
          Do ir=1,mynr
             buff(ir,iq) = Cp*Sum((pu(:,:,ir)-avgvr(myr1+ir-1))*prho(:,:,ir)*(temperature(:,:,ir)-avgT(myr1+ir-1)))/Dble(nsize)
          EndDo          
       Case (17) !F_ke (kinetic energy flux)
          Do ir=1,mynr
             buff(ir,iq) = Sum(usqr(:,:,ir)*(pu(:,:,ir)-avgvr(myr1+ir-1)))/Dble(nsize)
          EndDo
       Case (18) !F_ac (Acoustic flux)
          Do ir=1,mynr
             buff(ir,iq) = Sum((gamm1*Cv*prho(:,:,ir)*temperature(:,:,ir)-avgP(myr1+ir-1))*(pu(:,:,ir)-avgvr(myr1+ir-1)))/Dble(nsize)
          EndDo
       Case (19)
          ! Radial Entropy Diff Flux
          Do ir=1,mynr
             If (do_radheat) Then
                buff(ir,iq) = Sum(derivs(:,:,ir,7))/Dble(nsize)
             Else
                buff(ir,iq) = Sum(derivs(:,:,ir,6))/Dble(nsize)
             EndIf
          EndDo
       Case Default
          ! Error
          Write(*,*) 'Unrecognized Quantity In Shell_Average_Output'
       EndSelect
    End Do

    averages(myr1:myr2,3:nq+2) = buff
    
    Allocate(rbuff(nr,nq+2))
    rbuff = 0d0

    Call MPI_Reduce(averages(1,1),rbuff(1,1),nr*(nq+2),MPI_REAL8,MPI_SUM,nnodes/3,node_comm,ierr)    

    If (myrank .eq. nnodes/3) Then
       averages = rbuff
       averages(1,1) = nr
       averages(2,1) = nq
       averages(3,1) = r1
       averages(4,1) = r2
       averages(5:nq+4,1) = quantities
       averages(:,2) = r_ind
       !Write out averages
       Open(shell_avg_unit, file=Trim(shell_avg_filename), form='unformatted', access='direct', status='unknown',recl=irec_shlavg)
       
       Write(shell_avg_unit, rec=sh_avg_record) averages
       Close(shell_avg_unit)

       If ((lrad_update).and.(.not.do_radheat)) Then
          !dlnkr = (Lrad - 4d0*pi*(averages(:,nq-1)+averages(:,nq))*exp(-k0expon*((r_ind-flux_cutoff)/(r2-r1))**8)*r_ind**2)
          !dlnkr = 4d0*pi*(averages(:,nq-1)+averages(:,nq))*((Lrad/Maxval(Lrad))**k0expon)*(r_ind**2)
          drt = Dble(nr-1)/(x_ind(nr)-x_ind(1))
          bc1 = tr(1)/drt
          bc2 = tr(nr)/drt
          Call dbyd1_dt0(tr,avgvr,bc1,bc2,nr,1)
          avgvr = drt*avgvr
          If (non_uniform) avgvr = avgvr*dxdr

          !kappa_rosseland
          avgvr = 16d0*pi*a_const*c_light*(r_ind**2)*(tr**3)*avgvr/(3d0*rhor*Lrad)
          kr = avgvr

          dlnkr = 16d0*pi*a_const*c_light*(r_ind**2)*(avgT**3)*dtavg/(3d0*rhor*avgvr)
          Lrad_new = 0.99d0*Lrad_new + 0.01d0*dlnkr
          Lrad_new(1) = Lrad(1)

          avgvr = Lrad_new
          bc1 = Lrad_new(1)/drt
          bc2 = Lrad_new(nr)/drt
          dtavg=0d0
          Call dbyd1_dt0(avgvr,dtavg,bc1,bc2,nr,1)
          dlnkr = drt*dtavg
          If (non_uniform) dlnkr = dlnkr*dxdr

          dlnkr = -dlnkr/(4d0*pi*r_ind**2)
          k0 = Lrad_new/(4d0*pi*r_ind**2)
       EndIf

       If ((lrad_update).and.(do_radheat)) Then
          dlnkr = tr*(dsdrm/Cp-gamm1*rhor*gravity/(gamma*pr))

          !oringinal kappa rosseland
          avgvr = -16d0*pi*a_const*c_light*(r_ind**2)*(tr**3)*dlnkr/(3d0*rhor*Lrad)

          !updated kappa rosseland
          kr = -16d0*pi*a_const*c_light*(r_ind**2)*(avgT**3)*dtavg/(3d0*avgrho*Lrad)

          !Allow surface to change opacity, but depth to adjust luminosity.
          kr = kr*(1d0+tanh(lrad_tanh_coef*(r_ind-lrad_radius)/(r2-r1))) + avgvr*(1d0-tanh(lrad_tanh_coef*(r_ind-lrad_radius)/(r2-r1)))
          kr = 0.5d0*kr
          
          dlnkr = 16d0*pi*a_const*c_light*(r_ind**2)*(avgT**3)*dtavg/(3d0*avgrho*kr)

          Lrad_new = 0.99d0*Lrad_new + 0.01d0*dlnkr
          Lrad_new(1) = Lrad(1)

          avgvr = dlog(kr)
          bc1 = avgvr(1)/drt
          bc2 = avgvr(nr)/drt
          Call dbyd1_dt0(avgvr,dlnkr,bc1,bc2,nr,1)
          dlnkr = drt*dlnkr
          If (non_uniform) dlnkr = dlnkr*dxdr

          kr = 1d0/kr !Inverse rosseland opacity.
          k0 = Lrad_new/(4d0*pi*r_ind**2)
       EndIf
    EndIf
    
    If (lrad_update) Then
       Deallocate(avgvr)
       Allocate(avgvr(3*nr))
       avgvr(1:nr) = k0
       avgvr(nr+1:2*nr) = kr
       avgvr(2*nr+1:) = dlnkr
       Call MPI_BCAST(avgvr,3*nr,MPI_REAL8,nnodes/3,node_comm,ierr)
       k0 = avgvr(1:nr)
       kr = avgvr(nr+1:2*nr)
       dlnkr = avgvr(2*nr+1:)
       If (.not. do_radheat) Then
          eps = dlnkr
          eps(nr-8:nr) = 0d0
       EndIf
    EndIf

    If (Allocated(dtavg)) Deallocate(dtavg)
    If (Allocated(avgvr)) Deallocate(avgvr)
    If (Allocated(snd)) Deallocate(snd)
    If (Allocated(rcv)) Deallocate(rcv)
    If (Allocated(rbuff)) Deallocate(rbuff)
    If (Allocated(buff)) Deallocate(buff)
    If (Allocated(averages)) Deallocate(averages)
    
    If (do_ke) Then
       Deallocate(usqr)
    EndIf
    If (do_me) Then
       Deallocate(bsqr)
    EndIf
    If (do_bpol) Then
       Deallocate(bpol)
    EndIf
    
    Nullify(prho,pw,pv,pu,ps)
    If (magnetic) Nullify(pbr,pbt,pbp)

  End Subroutine Shell_Average_Output
  
  Subroutine Shell_Slice_Output(quantities,radii,io_proc)
    Implicit None
    Real*8, Pointer, Dimension(:,:,:) :: prho,pw,pv,pu,ps,pbr,pbt,pbp
    Real*8, Allocatable, Dimension(:,:,:,:) :: averages,buff
    Real*8, Allocatable, Dimension(:,:,:) :: usqr,bsqr,bpol,rbuff,sbuff
    Integer, Intent(In) :: quantities(:),radii(:),io_proc
    Integer :: iq,ir,ip,it,nip,nq,myir,nradii,ip_nt1,ip_nt2,ip_np1,ip_np2,ibuff(4),shell_unit=40
    
    !  useful aliases
    Nullify(prho,pw,pv,pu,ps)
    prho => vars0(:,:,:,1)
    pw   => vars0(:,:,:,2)
    pv   => vars0(:,:,:,3)
    pu   => vars0(:,:,:,4)
    ps   => vars0(:,:,:,5)
    If (magnetic) Then
       Nullify(pbr,pbt,pbp)
       pbt => vars0(:,:,:,6)
       pbp => vars0(:,:,:,7)
       pbr => vars0(:,:,:,8)
    End If
    
    nq = size(quantities)
    nradii = size(radii)
    
    nip=0
    Do ir=1,nradii
       If ((radii(ir) .ge. myr1) .and. (radii(ir) .le. myr2)) Then 
          nip=nip+1
       End If
    End Do
    
    Allocate(buff(mynth,mynphi,nip,nq),sbuff(mynth,mynphi,nq))
    
    !temperature = prho**gamm1*dexp(ps/Cv)/gamm1/Cv
    
    If (myrank.Eq.io_proc) Then
       Allocate(averages(nth,nphi+1,nradii,nq))
       averages = 0d0
       averages(1,1,1,1) = nradii
       averages(2,1,1,1) = nr
       averages(3,1,1,1) = nth
       averages(4,1,1,1) = nphi
       averages(5,1,1,1) = nq
       averages(6,1,1,1) = r1
       averages(7,1,1,1) = r2
       averages(8,1,1,1) = th1
       averages(9,1,1,1) = th2
       averages(10,1,1,1) = phi1
       averages(11,1,1,1) = phi2
       averages(12:nq+11,1,1,1) = quantities
       averages(nq+12:nradii+nq+11,1,1,1) = radii
       If (irec_shlslice .eq. -1) Then
          Inquire(iolength=irec_shlslice) averages
       EndIf
    End If
    
    !Kinetic Energy
    If (do_ke) Then
       Allocate(usqr(mynth,mynphi,mynr))
       usqr=pu**2+pv**2+pw**2
       usqr=prho*usqr/2
    End If
    
    If (magnetic) Then
       !Magnetic Energy
       If (do_me) Then
          Allocate(bsqr(mynth,mynphi,mynr))
          bsqr=pbr**2+pbt**2+pbp**2
          bsqr=bsqr/(8d0*pi)
       End If
       
       !Poloidal Field
       If (do_bpol) Then
          Allocate(bpol(mynth,mynphi,mynr))
          bpol = pbr**2+pbt**2
          bpol = Sqrt(bpol)
       End If
    End If
    
    Do iq=1,nq
       it=0
       Do ir=1,nradii
          If ((radii(ir) .ge. myr1) .and. (radii(ir) .le. myr2)) Then 
             it=it+1
             myir = radii(ir)-myr1+1
             If (quantities(iq).Eq.1) Then
                !  Rho
                buff(:,:,it,iq) = prho(:,:,myir)
             Else If (quantities(iq).Eq.2) Then
                !  Temperature
                buff(:,:,it,iq) = temperature(:,:,myir)
             Else If (quantities(iq).Eq.3) Then
                !  Entropy
                buff(:,:,it,iq) = ps(:,:,myir)
             Else If (quantities(iq).Eq.4) Then
                !  Pressure
                buff(:,:,it,iq) = gamm1*Cv*prho(:,:,myir)*temperature(:,:,myir)
             Else If (quantities(iq).Eq.5) Then
                !  Vr
                buff(:,:,it,iq) = pu(:,:,myir)
             Else If (quantities(iq).Eq.6) Then
                !  Vtheta
                buff(:,:,it,iq) = pw(:,:,myir)
             Else If (quantities(iq).Eq.7) Then
                !  Vphi
                buff(:,:,it,iq) = pv(:,:,myir)
             Else If (quantities(iq).Eq.8) Then
                !  Kinetic Energy
                buff(:,:,it,iq) = usqr(:,:,myir)
             Else If (quantities(iq).Eq.10) Then
                !  Br
                buff(:,:,it,iq) = pbr(:,:,myir)
             Else If (quantities(iq).Eq.11) Then
                !  Btheta
                buff(:,:,it,iq) = pbt(:,:,myir)
             Else If (quantities(iq).Eq.12) Then
                !  Bphi
                buff(:,:,it,iq) = pbp(:,:,myir)
             Else If (quantities(iq).Eq.13) Then
                !  Magnetic Energy
                buff(:,:,it,iq) = bsqr(:,:,myir)
             Else If (quantities(iq).Eq.14) Then
                !  Poloidal Field Magnitude
                buff(:,:,it,iq) = bpol(:,:,myir)
             Else
                ! Error
                Write(*,*) 'Unrecognized Quantity In Shell_Slice_Output'
             End If
          End If
       End Do
    End Do
    
    If (myrank .eq. io_proc) Then
       ip=0
       Do ir=1,nradii
          Do it=1,nnodes
             If ((ilupr(ir,it) .ne. -1).and.(ilupr(ir,it) .ne. io_proc)) Then
                ip_nt1 = 0
                ip_nt2 = 0
                ip_np1 = 0
                ip_np2 = 0
                ibuff = 0
                Call Receive(ibuff,ilupr(ir,it),ilupr(ir,it)+ir,node_comm)
                ip_nt1=ibuff(1)
                ip_nt2=ibuff(2)
                ip_np1=ibuff(3)
                ip_np2=ibuff(4)
                Allocate(rbuff(ip_nt2-ip_nt1+1,ip_np2-ip_np1+1,nq))
                rbuff = 0d0
                Call Receive(rbuff,ilupr(ir,it),ilupr(ir,it)+2*ir,node_comm)
                averages(ip_nt1:ip_nt2,ip_np1+1:ip_np2+1,ir,:) = rbuff
                Deallocate(rbuff)
             ElseIf ((ilupr(ir,it) .ne. -1) .and. (ilupr(ir,it) .eq. io_proc)) Then
                ip=ip+1
                averages(myth1:myth2,myphi1+1:myphi2+1,ir,:) = buff(:,:,ip,:)
             EndIf
          End Do
       End Do
       !Write out averages
       Open(shell_unit, file=Trim(shell_slice_filename), form='unformatted', access='direct', status='unknown',recl=irec_shlslice)
       
       Write(shell_unit, rec=sh_slice_record) averages
       Close(shell_unit)
    End If
    
    ! send buffer to IO node
    If (myrank .ne. io_proc) Then
       it=0
       Do ir = 1, nradii
          If ((radii(ir) .ge. myr1) .and. (radii(ir) .le. myr2)) Then
             it=it+1
             Call Send((/myth1,myth2,myphi1,myphi2/),io_proc,myrank+ir,node_comm)
             sbuff = buff(:,:,it,:)
             Call Send(sbuff,io_proc,myrank+2*ir,node_comm)
          End If
       End Do
    EndIf

    If (Allocated(averages)) Deallocate(averages)
    If (Allocated(buff)) Deallocate(buff)
    If (Allocated(sbuff)) Deallocate(sbuff)
    
    If (do_ke) Then
       Deallocate(usqr)
    EndIf
    If (do_me) Then
       Deallocate(bsqr)
    EndIf
    If (do_bpol) Then
       Deallocate(bpol)
    EndIf
    
    Nullify(prho,pw,pv,pu,ps)
    If (magnetic) Nullify(pbr,pbt,pbp)
  End Subroutine Shell_Slice_Output
  
  Subroutine Make_Inverse_Lookup(radii)
    Implicit None
    Integer, Intent(In) :: radii(:)
    Integer :: ii,ir,nrad,irem,rrank,ynr,myd1(nnr),myd2(nnr),mynd(nnr),temp(nnodes,2),ccoords(nnodes),isnd(nnodes),ircv(nnodes)
    
    nrad = Size(radii)
    
    ynr = nr/nnr
    
    irem = Mod(nr,nnr) ! guaranteed to be less than nnr, so add one to each node less than nnr
    !Distrbute any remaining points over iremr nodes for load balancing
    
    If (irem .gt. 0) Then
       myd1(1) = 1
       myd2(1) = ynr+1
       mynd(1) = ynr+1
       Do ii=2,nnr
          If (ii .le. irem) Then
             mynd(ii) = ynr+1
          Else
             mynd(ii) = ynr
          EndIf
          myd1(ii) = myd2(ii-1)+1
          myd2(ii) = myd1(ii)+mynd(ii)-1
       EndDo
    Else
       mynd = ynr
       myd1(1) = 1
       myd2(1) = ynr
       Do ii=2,nnr
          myd1(ii) = myd2(ii-1)+1
          myd2(ii) = myd1(ii)+mynd(ii)-1
       EndDo
    EndIf
    
    ccoords = 0
    ccoords(myrank+1) = cartcoords(3)
    isnd = ccoords
    ircv = 0
    Call Global_AllReduce(isnd,ircv,node_comm,'sum')
    ccoords = ircv
    
    Do ir=0,nnodes-1
       !Determine node data bounds
       rrank = ccoords(ir+1)
       temp(ir+1,1)=myd1(rrank+1)
       temp(ir+1,2)=myd2(rrank+1)
    EndDo
    
    Allocate(ilupr(nrad,nnodes))
    
    Do ii=1,nrad
       Do ir=0,nnodes-1
          If ((radii(ii) .ge. temp(ir+1,1)) .and. (radii(ii) .le. temp(ir+1,2))) Then
             ilupr(ii,ir+1)=ir
          Else
             ilupr(ii,ir+1)=-1
          EndIf
       EndDo
    EndDo
  End Subroutine Make_Inverse_Lookup

End Module Report

