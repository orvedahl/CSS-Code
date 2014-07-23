Module Validation

  Use Physics

Contains

  Subroutine TestDerivatives()
    Real*8 :: work1(mynth,mynr), work2(mynth,mynr)
    Real*8 :: workt(mynth,mynphi), workp(mynth,mynphi), workr(mynth,mynr)
    Real*8 :: dworkt(mynth,mynphi), dworkp(mynth,mynphi), dworkr(mynth,mynr)
    Real*8 :: worktb1(binfo(2,1)%bsize%bnt,mynphi), workpb1(mynth,binfo(3,1)%bsize%bnp), workrb1(mynth,binfo(1,1)%bsize%bnr)
    Real*8 :: worktb2(binfo(2,2)%bsize%bnt,mynphi), workpb2(mynth,binfo(3,2)%bsize%bnp), workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: invr, testd1(2), testd2(2), testd3(2), snd(3), rcv(3), nflop, stime, ftime, flops
    Integer :: itst

    !Setup function and exact derivatives
    Do ir=1,mynr
       Do ip=1,mynphi
          vars1(:,ip,ir,1) = cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))
          vars1(:,ip,ir,2) = -pi*nmodes*cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*sin(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))/(th2-th1)
          vars1(:,ip,ir,3) = pi*nmodes*cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &cos(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))/(phi2-phi1)
          vars1(:,ip,ir,4) = -pi*nmodes*sin(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))/(r2-r1)
       EndDo
    EndDo

    Do iv=0,4
       Call cpu_time(stime)
       !Take Gradient
       !Theta
       Call Boundary_Exchange_Init(2)
       If (myth1 .eq. 1) Then
          Select Case(iv)
          Case(0)
             binfo(2,1)%bdata(1,:,:,1) = 0d0
          Case(1)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,1)
          Case(2)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,2)
          Case(3)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,1)
          Case(4)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,2)
          EndSelect
       EndIf
       If (myth2 .eq. nth) Then
          Select Case(iv)
          Case(0)
             binfo(2,2)%bdata(1,:,:,1) = 0d0
          Case(1)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,1)
          Case(2)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,1)
          Case(3)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,2)
          Case(4)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,2)
          EndSelect
       EndIf

       Call Boundary_Exchange_Finish(2)
       Call Boundary_Exchange_Init(3)
       Do ir=1,mynr
          workt = vars1(:,:,ir,1)
          worktb1 = binfo(2,1)%bdata(:,:,ir,1)
          worktb2 = binfo(2,2)%bdata(:,:,ir,1)
          Call dbydt(workt,dworkt,worktb1,worktb2,-iv) !Nflops = (11*mynth*+25)*mynr*mynphi
          dvars1(:,:,ir,1) = dworkt
       EndDo

       testd2(1) = Sum(abs(vars1(:,:,:,2)-dvars1(:,:,:,1)))/Sum(abs(vars1(:,:,:,2))) !Nflops=5*mynth*mynr*mynphi
       testd2(2) = Maxval(abs(vars1(:,:,:,2)-dvars1(:,:,:,1)))/Maxval(abs(vars1(:,:,:,2))) !Nflops=3*mynth*mynr*mynphi

       !Phi
       Call Boundary_Exchange_Finish(3)
       Call Boundary_Exchange_Init(1)
       Do ir=1,mynr
          workp = vars1(:,:,ir,1)
          workpb1 = binfo(3,1)%bdata(:,:,ir,1)
          workpb2 = binfo(3,2)%bdata(:,:,ir,1)
          Call dbydp(workp,dworkp,workpb1,workpb2,-iv) !Nflops = (11*mynphi*+25)*mynr*mynth
          dvars1(:,:,ir,2) = dworkp
       EndDo

       testd3(1) = Sum(abs(vars1(:,:,:,3)-dvars1(:,:,:,2)))/Sum(abs(vars1(:,:,:,3))) !Nflops=5*mynth*mynr*mynphi
       testd3(2) = Maxval(abs(vars1(:,:,:,3)-dvars1(:,:,:,2)))/Maxval(abs(vars1(:,:,:,3))) !Nflops=3*mynth*mynr*mynphi

       If (non_uniform) Then
          Do ip=1,mynth
             work1(ip,:) = dxdr(myr1:myr2)
          EndDo
       EndIf

       !Radius
       Call Boundary_Exchange_Finish(1)
       If (myr1 .eq. 1) Then
          Select Case(iv)
          Case(0)
             binfo(1,1)%bdata(:,:,1,1) = 0d0
          Case(1)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,1)
          Case(2)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,4)
          Case(3)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,1)
          Case(4)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,4)
          EndSelect
       EndIf

       If (myr2 .eq. nr) Then
          Select Case(iv)
          Case(0)
             binfo(1,2)%bdata(:,:,1,1) = 0d0
          Case(1)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,1)
          Case(2)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,1)
          Case(3)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,4)
          Case(4)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,4)
          EndSelect
       EndIf

       Do ip=1,mynphi
          workr = vars1(:,ip,:,1)
          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
          Call dbydr(workr,dworkr,workrb1,workrb2,-iv) !Nflops = (11*mynr*+25)*mynth*mynphi
          If (non_uniform) dworkr = dworkr*work1             
          dvars1(:,ip,:,3) = dworkr
       EndDo

       testd1(1) = Sum(abs(vars1(:,:,:,4)-dvars1(:,:,:,3)))/Sum(abs(vars1(:,:,:,4))) !Nflops=5*mynth*mynr*mynphi
       testd1(2) = Maxval(abs(vars1(:,:,:,4)-dvars1(:,:,:,3)))/Maxval(abs(vars1(:,:,:,4))) !Nflops=3*mynth*mynr*mynphi

       Call cpu_time(ftime)

       nflop = (11d0*Dble(mynr)+25d0)*Dble(mynth)*Dble(mynphi)+(11d0*Dble(mynphi)+25d0)*Dble(mynth)*Dble(mynr)+(11d0*Dble(mynth)+25d0)*Dble(mynr)*Dble(mynth)+24d0*Dble(mynr*mynth*mynphi)
       ftime = ftime-stime
       flops = nflop/ftime

       snd(1) = testd1(1)
       snd(2) = testd2(1)
       snd(3) = testd3(1)
       rcv = 0d0
       Call Global_AllReduce(snd,rcv,node_comm,'sum')
       testd1(1) = rcv(1)
       testd2(1) = rcv(2)
       testd3(1) = rcv(3)

       snd(1) = testd1(2)
       snd(2) = testd2(2)
       snd(3) = testd3(2)
       rcv = 0d0
       Call Global_AllReduce(snd,rcv,node_comm,'max')
       testd1(2) = rcv(1)
       testd2(2) = rcv(2)
       testd3(2) = rcv(3)

       snd(1) = flops
       Call Global_AllReduce(snd,rcv,node_comm,'sum')
       flops = rcv(1)

       If (myrank .eq. 0) Then
          Print*, 'Validation of derivatives:'
          Print*, '                      '

          Print*, 'Relative error in radial derivatives for ibc=', iv
          Print*, 'minmax df_ex', Minval(vars1(:,:,:,4)), Maxval(vars1(:,:,:,4))
          Print*, 'minmax df_cdf', Minval(dvars1(:,:,:,3)), Maxval(dvars1(:,:,:,3))
          Print*, 'minmax binfo(1,1)', Minval(binfo(1,1)%bdata(:,:,:,1)), Maxval(binfo(1,1)%bdata(:,:,:,1))
          Print*, 'minmax binfo(1,2)', Minval(binfo(1,2)%bdata(:,:,:,1)), Maxval(binfo(1,2)%bdata(:,:,:,1))
          Print*, '||df_ex-df_cfd||_L1 =', testd1(1)
          Print*, '||df_ex-df_cfd||_Linf =', testd1(2)

          Print*, '                      '
          Print*, 'Relative error in theta derivatives for ibc=', iv
          Print*, 'minmax df_ex', Minval(vars1(:,:,:,2)), Maxval(vars1(:,:,:,2))
          Print*, 'minmax df_cdf', Minval(dvars1(:,:,:,1)), Maxval(dvars1(:,:,:,1))
          Print*, 'minmax binfo(2,1)', Minval(binfo(2,1)%bdata(:,:,:,1)), Maxval(binfo(2,1)%bdata(:,:,:,1))
          Print*, 'minmax binfo(2,2)', Minval(binfo(2,2)%bdata(:,:,:,1)), Maxval(binfo(2,2)%bdata(:,:,:,1))
          Print*, '||df_ex-df_cfd||_L1 =', testd2(1)
          Print*, '||df_ex-df_cfd||_Linf =', testd2(2)

          Print*, '                      '
          Print*, 'Relative error in phi derivatives for ibc=', iv
          Print*, 'minmax df_ex', Minval(vars1(:,:,:,3)), Maxval(vars1(:,:,:,3))
          Print*, 'minmax df_cdf', Minval(dvars1(:,:,:,2)), Maxval(dvars1(:,:,:,2))
          Print*, 'minmax binfo(3,1)', Minval(binfo(3,1)%bdata(:,:,:,1)), Maxval(binfo(3,1)%bdata(:,:,:,1))
          Print*, 'minmax binfo(3,2)', Minval(binfo(3,2)%bdata(:,:,:,1)), Maxval(binfo(3,2)%bdata(:,:,:,1))
          Print*, '||df_ex-df_cfd||_L1 =', testd3(1)
          Print*, '||df_ex-df_cfd||_Linf =', testd3(2)
          Print*, '                      '
          Print*, 'ex-fd in r = [', vars1(mynth/2,mynphi/2,:,4)-dvars1(mynth/2,mynphi/2,:,3), ']'
          Print*, '                      '
          Print*, 'ex-fd in th = [', vars1(:,mynphi/2,mynr/2,2)-dvars1(:,mynphi/2,mynr/2,1), ']'
          Print*, '                      '
          Print*, 'ex-fd in phi = [', vars1(mynth/2,:,mynr/2,3)-dvars1(mynth/2,:,mynr/2,2), ']'
          Print*, '                      '
          !Print*, '                      '       
          !Print*, 'Total MFLOPS      =', flops*1d-6
          !Print*, 'Total MFLOPS/core =', flops*1d-6/Dble(nnodes)
       EndIf

       Call Barrier(node_comm)       
       invr = Max(testd1(2),testd2(2),testd3(2))
       If (invr .gt. 1d-1) Then
          Print*, 'Derivative Validation Failed with max error=', invr, ' with boundary type', iv
          Call Graceful_Exit()
       EndIf
    EndDo

    !Speed test
    iv = 4
    ir = 1
    workt = vars1(:,:,ir,1)
    worktb1 = binfo(2,1)%bdata(:,:,ir,1)
    worktb2 = binfo(2,2)%bdata(:,:,ir,1)
    Call cpu_time(stime)
    Do itst=1,nphi*nphi
       !Do ir=1,mynr
       !workt = vars1(:,:,ir,1)
       !worktb1 = binfo(2,1)%bdata(:,:,ir,1)
       !worktb2 = binfo(2,2)%bdata(:,:,ir,1)
       Call dbydt(workt,dworkt,worktb1,worktb2,-iv) !Nflops = (11*mynth*+25)*mynr*mynphi
       !dvars1(:,:,ir,1) = dworkt

!!$       workt = vars1(:,:,ir,1)
!!$       worktb1 = binfo(2,1)%bdata(:,:,ir,1)
!!$       worktb2 = binfo(2,2)%bdata(:,:,ir,1)
!!$       Call dbydt(workt,dworkt,worktb1,worktb2,-iv) !Nflops = (11*mynth*+25)*mynr*mynphi
!!$       dvars1(:,:,ir,2) = dworkt
!!$
!!$       workt = vars1(:,:,ir,1)
!!$       worktb1 = binfo(2,1)%bdata(:,:,ir,1)
!!$       worktb2 = binfo(2,2)%bdata(:,:,ir,1)
!!$       Call dbydt(workt,dworkt,worktb1,worktb2,-iv) !Nflops = (11*mynth*+25)*mynr*mynphi
!!$       dvars1(:,:,ir,3) = dworkt
!!$
!!$       workt = vars1(:,:,ir,1)
!!$       worktb1 = binfo(2,1)%bdata(:,:,ir,1)
!!$       worktb2 = binfo(2,2)%bdata(:,:,ir,1)
!!$       Call dbydt(workt,dworkt,worktb1,worktb2,-iv) !Nflops = (11*mynth*+25)*mynr*mynphi
!!$       dvars1(:,:,ir,4) = dworkt
       !EndDo
    EndDo

    Call cpu_time(ftime)
    !nflop = 4d0*(11d0*Dble(mynth)+25d0)*Dble(mynr)*Dble(mynth)
    nflop = (11d0*Dble(mynth)+20d0)*Dble(mynphi)
    nflop = nflop*Dble(nphi*nphi)
    ftime = ftime-stime
    flops = nflop/ftime

    snd(1) = flops
    Call Global_AllReduce(snd,rcv,node_comm,'sum')
    flops = rcv(1)

    If (myrank .eq. 0) Then
       Print*, '                      '
       Print*, 'Theta Derivatives'
       Print*, 'Total MFLOPS      =', flops*1d-6
       Print*, 'Total MFLOPS/core =', flops*1d-6/Dble(nnodes)
    EndIf

    !Phi
    ir = 1
    workp = vars1(:,:,ir,1)
    workpb1 = binfo(3,1)%bdata(:,:,ir,1)
    workpb2 = binfo(3,2)%bdata(:,:,ir,1)
    Call cpu_time(stime)
    Do itst=1,nphi*nphi
       !Do ir=1,mynr

       !workp = vars1(:,:,ir,1)
       !workpb1 = binfo(3,1)%bdata(:,:,ir,1)
       !workpb2 = binfo(3,2)%bdata(:,:,ir,1)
       Call dbydp(workp,dworkp,workpb1,workpb2,-iv) !Nflops = (11*mynphi*+25)*mynr*mynth
       !dvars1(:,:,ir,1) = dworkp

!!$       workp = vars1(:,:,ir,1)
!!$       workpb1 = binfo(3,1)%bdata(:,:,ir,1)
!!$       workpb2 = binfo(3,2)%bdata(:,:,ir,1)
!!$       Call dbydp(workp,dworkp,workpb1,workpb2,-iv) !Nflops = (11*mynphi*+25)*mynr*mynth
!!$       dvars1(:,:,ir,2) = dworkp
!!$
!!$       workp = vars1(:,:,ir,1)
!!$       workpb1 = binfo(3,1)%bdata(:,:,ir,1)
!!$       workpb2 = binfo(3,2)%bdata(:,:,ir,1)
!!$       Call dbydp(workp,dworkp,workpb1,workpb2,-iv) !Nflops = (11*mynphi*+25)*mynr*mynth
!!$       dvars1(:,:,ir,3) = dworkp
!!$
!!$       workp = vars1(:,:,ir,1)
!!$       workpb1 = binfo(3,1)%bdata(:,:,ir,1)
!!$       workpb2 = binfo(3,2)%bdata(:,:,ir,1)
!!$       Call dbydp(workp,dworkp,workpb1,workpb2,-iv) !Nflops = (11*mynphi*+25)*mynr*mynth
!!$       dvars1(:,:,ir,4) = dworkp
!!$
!!$       workp = vars1(:,:,ir,1)
!!$       workpb1 = binfo(3,1)%bdata(:,:,ir,1)
!!$       workpb2 = binfo(3,2)%bdata(:,:,ir,1)
!!$       Call dbydp(workp,dworkp,workpb1,workpb2,-iv) !Nflops = (11*mynphi*+25)*mynr*mynth
!!$       dvars1(:,:,ir,5) = dworkp
       !EndDo
    EndDo

    Call cpu_time(ftime)
    !nflop = 5d0*(11d0*Dble(mynphi)+25d0)*Dble(mynth)*Dble(mynr)
    nflop = (11d0*Dble(mynphi)+20d0)*Dble(mynth)
    nflop = nflop*Dble(nphi*nphi)
    ftime = ftime-stime
    flops = nflop/ftime

    snd(1) = flops
    Call Global_AllReduce(snd,rcv,node_comm,'sum')
    flops = rcv(1)

    If (myrank .eq. 0) Then
       Print*, '                      '
       Print*, 'Phi Derivatives'
       Print*, 'Total MFLOPS      =', flops*1d-6
       Print*, 'Total MFLOPS/core =', flops*1d-6/Dble(nnodes)
    EndIf

    !Radius
    ip = 1
    workr = vars1(:,ip,:,1)
    workrb1 = binfo(1,1)%bdata(:,ip,:,1)
    workrb2 = binfo(1,2)%bdata(:,ip,:,1)
    Call cpu_time(stime)
    Do itst=1,nphi*nphi
       !Do ip=1,mynphi
          !workr = vars1(:,ip,:,1)
          !workrb1 = binfo(1,1)%bdata(:,ip,:,1)
          !workrb2 = binfo(1,2)%bdata(:,ip,:,1)
          Call dbydr(workr,dworkr,workrb1,workrb2,-iv) !Nflops = (11*mynr*+25)*mynth*mynphi
          !If (non_uniform) dworkr = dworkr*work1             
          !dvars1(:,ip,:,1) = dworkr

!!$          workr = vars1(:,ip,:,1)
!!$          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
!!$          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
!!$          Call dbydr(workr,dworkr,workrb1,workrb2,-iv) !Nflops = (11*mynr*+25)*mynth*mynphi
!!$          If (non_uniform) dworkr = dworkr*work1             
!!$          dvars1(:,ip,:,2) = dworkr
!!$
!!$          workr = vars1(:,ip,:,1)
!!$          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
!!$          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
!!$          Call dbydr(workr,dworkr,workrb1,workrb2,-iv) !Nflops = (11*mynr*+25)*mynth*mynphi
!!$          If (non_uniform) dworkr = dworkr*work1             
!!$          dvars1(:,ip,:,3) = dworkr
!!$
!!$          workr = vars1(:,ip,:,1)
!!$          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
!!$          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
!!$          Call dbydr(workr,dworkr,workrb1,workrb2,-iv) !Nflops = (11*mynr*+25)*mynth*mynphi
!!$          If (non_uniform) dworkr = dworkr*work1             
!!$          dvars1(:,ip,:,4) = dworkr
!!$       EndDo
    EndDo
    Call cpu_time(ftime)
    !nflop = 4d0*(11d0*Dble(mynr)+25d0)*Dble(mynth)*Dble(mynphi)
    nflop = (11d0*Dble(mynr)+20d0)*Dble(mynth)
    nflop = nflop*Dble(nphi*nphi)
    ftime = ftime-stime
    flops = nflop/ftime

    snd(1) = flops
    Call Global_AllReduce(snd,rcv,node_comm,'sum')
    flops = rcv(1)

    If (myrank .eq. 0) Then
       Print*, '                      '
       Print*, 'Radial Derivatives'
       Print*, 'Total MFLOPS      =', flops*1d-6
       Print*, 'Total MFLOPS/core =', flops*1d-6/Dble(nnodes)
    EndIf

    !Setup function and exact derivatives
    Do ir=1,mynr
       Do ip=1,mynphi
          vars1(:,ip,ir,1) = cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))
          vars2(:,ip,ir,2) = -(pi*nmodes)**2*cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))/(th2-th1)**2
          vars2(:,ip,ir,3) = -(pi*nmodes)**2*cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))/(phi2-phi1)**2
          vars2(:,ip,ir,4) = -(pi*nmodes)**2*cos(pi*nmodes*(r_ind(myr1+ir-1)-r1)/(r2-r1))*cos(pi*nmodes*(theta(myth1:myth2)-th1)/(th2-th1))*&
               &sin(pi*nmodes*(phi(myphi1+ip-1)-phi1)/(phi2-phi1))/(r2-r1)**2
       EndDo
    EndDo

    Do iv=0,4
       Call cpu_time(stime)
       !Take Gradient
       !Theta
       Call Boundary_Exchange_Init(2)
       If (myth1 .eq. 1) Then
          Select Case(iv)
          Case(0)
             binfo(2,1)%bdata(1,:,:,1) = 0d0
          Case(1)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,1)
          Case(2)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,2)
          Case(3)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,1)
          Case(4)
             binfo(2,1)%bdata(1,:,:,1) = vars1(1,:,:,2)
          EndSelect
       EndIf
       If (myth2 .eq. nth) Then
          Select Case(iv)
          Case(0)
             binfo(2,2)%bdata(1,:,:,1) = 0d0
          Case(1)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,1)
          Case(2)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,1)
          Case(3)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,2)
          Case(4)
             binfo(2,2)%bdata(1,:,:,1) = vars1(mynth,:,:,2)
          EndSelect
       EndIf
       Call Boundary_Exchange_Finish(2)
       Call Boundary_Exchange_Init(3)
       Do ir=1,mynr
          workt = vars1(:,:,ir,1)
          worktb1 = binfo(2,1)%bdata(:,:,ir,1)
          worktb2 = binfo(2,2)%bdata(:,:,ir,1)
          Call d2bydt2(workt,dworkt,worktb1,worktb2,-iv)
          dvars1(:,:,ir,1) = dworkt
       EndDo

       testd2(1) = Sum(abs(vars2(:,:,:,2)-dvars1(:,:,:,1)))/Sum(abs(vars2(:,:,:,2)))
       testd2(2) = Maxval(abs(vars2(:,:,:,2)-dvars1(:,:,:,1)))/Maxval(abs(vars2(:,:,:,2)))

       !Phi
       Call Boundary_Exchange_Finish(3)
       Do ir=1,mynr
          workp = vars1(:,:,ir,1)
          workpb1 = binfo(3,1)%bdata(:,:,ir,1)
          workpb2 = binfo(3,2)%bdata(:,:,ir,1)
          Call d2bydp2(workp,dworkp,workpb1,workpb2,-iv)
          dvars1(:,:,ir,2) = dworkp
       EndDo


       Call Boundary_Exchange_Init(1)
       testd3(1) = Sum(abs(vars2(:,:,:,3)-dvars1(:,:,:,2)))/Sum(abs(vars2(:,:,:,3)))
       testd3(2) = Maxval(abs(vars2(:,:,:,3)-dvars1(:,:,:,2)))/Maxval(abs(vars2(:,:,:,3)))

       If (non_uniform) Then
          Do ip=1,mynth
             work1(ip,:) = dxdr(myr1:myr2)
             work2(ip,:) = d2xdr2(myr1:myr2)
          EndDo
       EndIf
       !Radius

       If (myr1 .eq. 1) Then
          Select Case(iv)
          Case(0)
             binfo(1,1)%bdata(:,:,1,1) = 0d0
          Case(1)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,1)
          Case(2)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,4)
          Case(3)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,1)
          Case(4)
             binfo(1,1)%bdata(:,:,1,1) = vars1(:,:,1,4)
          EndSelect
       EndIf

       If (myr2 .eq. nr) Then
          Select Case(iv)
          Case(0)
             binfo(1,2)%bdata(:,:,1,1) = 0d0
          Case(1)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,1)
          Case(2)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,1)
          Case(3)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,4)
          Case(4)
             binfo(1,2)%bdata(:,:,1,1) = vars1(:,:,mynr,4)
          EndSelect
       EndIf
       dvars1(:,:,:,3) = 0d0
       Call Boundary_Exchange_Finish(1)
       Do ip=1,mynphi
          workr = vars1(:,ip,:,1)
          workrb1 = binfo(1,1)%bdata(:,ip,:,1)
          workrb2 = binfo(1,2)%bdata(:,ip,:,1)
          If (non_uniform) Then
             Call dbydr(workr,dworkr,workrb1,workrb2,-iv)
             dworkr = dworkr*work2
             dvars1(:,ip,:,3) = dworkr
          EndIf
          Call d2bydr2(workr,dworkr,workrb1,workrb2,-iv)
          If (non_uniform) Then
             dworkr = dworkr*work1**2
          EndIf
          dvars1(:,ip,:,3) = dvars1(:,ip,:,3) + dworkr
       EndDo

       testd1(1) = Sum(abs(vars2(:,:,:,4)-dvars1(:,:,:,3)))/Sum(abs(vars2(:,:,:,4)))
       testd1(2) = Maxval(abs(vars2(:,:,:,4)-dvars1(:,:,:,3)))/Maxval(abs(vars2(:,:,:,4)))

       Call cpu_time(ftime)

       nflop = (11d0*Dble(mynr)+25d0)*Dble(mynth)*Dble(mynphi)+(11d0*Dble(mynphi)+25d0)*Dble(mynth)*Dble(mynr)+(11d0*Dble(mynth)+25d0)*Dble(mynr)*Dble(mynth)+24d0*Dble(mynr*mynth*mynphi)
       ftime = ftime-stime
       flops = nflop/ftime

       snd(1) = testd1(1)
       snd(2) = testd2(1)
       snd(3) = testd3(1)
       rcv = 0d0
       Call Global_AllReduce(snd,rcv,node_comm,'sum')
       testd1(1) = rcv(1)
       testd2(1) = rcv(2)
       testd3(1) = rcv(3)

       snd(1) = testd1(2)
       snd(2) = testd2(2)
       snd(3) = testd3(2)
       rcv = 0d0
       Call Global_AllReduce(snd,rcv,node_comm,'max')
       testd1(2) = rcv(1)
       testd2(2) = rcv(2)
       testd3(2) = rcv(3)

       snd(1) = flops
       Call Global_AllReduce(snd,rcv,node_comm,'sum')
       flops = rcv(1)

       If (myrank .eq. 0) Then
          Print*, 'Validation of 2nd derivatives:'
          Print*, '                      '
          Print*, 'Relative error in radial 2nd derivatives for ibc=', iv
          Print*, 'minmax df_ex', Minval(vars2(:,:,:,4)), Maxval(vars2(:,:,:,4))
          Print*, 'minmax df_cdf', Minval(dvars1(:,:,:,3)), Maxval(dvars1(:,:,:,3))
          Print*, 'minmax binfo(1,1)', Minval(binfo(1,1)%bdata(:,:,:,1)), Maxval(binfo(1,1)%bdata(:,:,:,1))
          Print*, 'minmax binfo(1,2)', Minval(binfo(1,2)%bdata(:,:,:,1)), Maxval(binfo(1,2)%bdata(:,:,:,1))
          Print*, '||df_ex-df_cfd||_L1 =', testd1(1)
          Print*, '||df_ex-df_cfd||_Linf =', testd1(2)

          Print*, '                      '
          Print*, 'Relative error in theta 2nd derivatives for ibc=', iv
          Print*, 'minmax df_ex', Minval(vars2(:,:,:,2)), Maxval(vars2(:,:,:,2))
          Print*, 'minmax df_cdf', Minval(dvars1(:,:,:,1)), Maxval(dvars1(:,:,:,1))
          Print*, 'minmax binfo(2,1)', Minval(binfo(2,1)%bdata(:,:,:,1)), Maxval(binfo(2,1)%bdata(:,:,:,1))
          Print*, 'minmax binfo(2,2)', Minval(binfo(2,2)%bdata(:,:,:,1)), Maxval(binfo(2,2)%bdata(:,:,:,1))
          Print*, '||df_ex-df_cfd||_L1 =', testd2(1)
          Print*, '||df_ex-df_cfd||_Linf =', testd2(2)

          Print*, '                      '
          Print*, 'Relative error in phi 2nd derivatives for ibc=', iv
          Print*, 'minmax df_ex', Minval(vars2(:,:,:,3)), Maxval(vars2(:,:,:,3))
          Print*, 'minmax df_cdf', Minval(dvars1(:,:,:,2)), Maxval(dvars1(:,:,:,2))
          Print*, 'minmax binfo(3,1)', Minval(binfo(3,1)%bdata(:,:,:,1)), Maxval(binfo(3,1)%bdata(:,:,:,1))
          Print*, 'minmax binfo(3,2)', Minval(binfo(3,2)%bdata(:,:,:,1)), Maxval(binfo(3,2)%bdata(:,:,:,1))
          Print*, '||df_ex-df_cfd||_L1 =', testd3(1)
          Print*, '||df_ex-df_cfd||_Linf =', testd3(2)

          !Print*, '                      '
          !Print*, 'Total MFLOPS      =', flops*1d-6
          !Print*, 'Total MFLOPS/core =', flops*1d-6/Dble(nnodes)
       EndIf
       Call Barrier(node_comm)       
       invr = Max(testd1(2),testd2(2),testd3(2))
       If (invr .gt. 1d-1) Then
          Print*, 'Second Derivative Validation Failed with max error=', invr, ' with boundary type', iv
          Call Graceful_Exit()
       EndIf
    EndDo
  End Subroutine TestDerivatives
End Module Validation
