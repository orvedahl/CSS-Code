Module Checkpoint

  Use IFPORT
  Use Physics
  Use Interpolation
  
  Implicit None
  
  Integer, Parameter :: max_qs=30, max_rads=20
  Integer :: irec_rain=-1

Contains

  Subroutine Make_Directory(dir,good)
    Implicit None
    Character(256), Intent(In) :: dir
    Logical(4), Intent(Out) :: good

    Print*, 'Trim(dir)=',Trim(dir)
    Inquire(Directory=Trim(dir),Exist=good)
    If (.not. good) Then
       good = makedirqq(Trim(dir))
       If (.not. good) Then
          Print*, 'Directory not created, attempting again:'
          good = makedirqq(Trim(dir))
          If (.not. good) Then
             Print*, 'Directory not created, attempting again:'
             good = makedirqq(Trim(dir))
             If (.not. good) Then
                Print*, 'Directory not created, filesystem error. Exiting...'
             EndIf
          EndIf
       EndIf
    EndIf

  End Subroutine Make_Directory
  
  Subroutine Make_Checkpoint_Filename(filename,fout,fhead,interp)
    Implicit None
    Character(256), Intent(In) :: filename
    Character(256), Intent(Out) :: fout,fhead
    Logical, Intent(In), Optional :: interp
    Character(10) :: temp, pretemp
    Integer :: k
    Logical(4) :: good
    
    pretemp = ' '
    Write(temp,fmt='(i7.7)') istep
    fout = Trim(filename)//Trim(temp)
    fout = Trim(fout)

    If (Present(interp)) Then
       fout = Trim(fout)//'_interp/'
    Else
       fout = Trim(fout)//'/'
    EndIf
    
    !Make directory for output
    k = 0
    If (myrank .eq. 0) Then
       Call Make_Directory(fout,good)
       If (.not. good) Then
          k = 1
       EndIf
    EndIf

    Call MPI_Bcast(k,1,MPI_INTEGER,0,node_comm,ierr)

    If (k .eq. 1) Then
       Call Graceful_Exit()
    EndIf

    fhead = Trim(fout)//'header'
    fout = Trim(fout)//'checkpoint.'    
  End Subroutine Make_Checkpoint_Filename

  Subroutine Read_Rain_Squalls(fname)
    Implicit None
    Character(256), Intent(In) :: fname
    Real*8, Allocatable :: read_buff(:)
    Integer :: checkpoint_unit, lct, lcp, is, ic, ir

    If (myr2 .eq. nr) Then
       If (rainstats .eq. 1) Then
          isz_max = Int(Ceiling(meansize*sqrt(-log(1d-5))/(r2*Min(dth,dphi*Minval(sines)))))
       EndIf

       Allocate(dent_rain(mynth,mynphi),dvel_rain(mynth,mynphi))

       dvel_rain = 0d0
       dent_rain = 0d0

       checkpoint_unit=12
       If (build_lanes) Then
          Allocate(read_buff(7),conv_cells(3*ncells))
          
          Inquire(iolength=irec_rain) read_buff
          Open(checkpoint_unit,file=fname,status='old',form='unformatted', access='direct',recl=irec_rain)
          Do ir=1,ncells
             Read(checkpoint_unit,rec=ir) read_buff
             conv_cells(ir)%loct = Int(read_buff(1))
             conv_cells(ir)%locp = Int(read_buff(2))
             conv_cells(ir)%Theta_loc = read_buff(3)
             conv_cells(ir)%Phi_loc = read_buff(4)
             conv_cells(ir)%Duration = read_buff(5)
             conv_cells(ir)%Start_Time = read_buff(6)
             conv_cells(ir)%Scaling = read_buff(7)
          EndDo
          Close(checkpoint_unit)
          Deallocate(read_buff)
       EndIf

       If (nsqualls .gt. 0) Then
          Allocate(rain_squalls(nsqualls),read_buff(12),squall_density(mynth,mynphi)) 
          
          Inquire(iolength=irec_rain) read_buff
          Open(checkpoint_unit,file=fname,status='old',form='unformatted', access='direct',recl=irec_rain)
          Do ir=1,nsqualls
             Read(checkpoint_unit,rec=ir) read_buff
             
             rain_squalls(ir)%Duration = read_buff(1)
             rain_squalls(ir)%SigmaT = read_buff(2)
             rain_squalls(ir)%Start_Time = read_buff(3)
             rain_squalls(ir)%Size = read_buff(4)
             rain_squalls(ir)%Dels = read_buff(5)
             rain_squalls(ir)%Theta_loc = read_buff(6)
             rain_squalls(ir)%Phi_loc = read_buff(7)
             rain_squalls(ir)%Scaling = read_buff(8)
             rain_squalls(ir)%loct = Int(read_buff(9))
             rain_squalls(ir)%locp = Int(read_buff(10))
             rain_squalls(ir)%Vel = read_buff(11)
             rain_squalls(ir)%Delr = read_buff(12)
          EndDo
       EndIf
       Close(checkpoint_unit)
       Deallocate(read_buff)
       If (mytop_rank .eq. 0) Print*, 'rain_squalls read in.'

       Call Update_Rain_Squalls()

    EndIf
    Call Barrier(node_comm)
    
  End Subroutine Read_Rain_Squalls
     
  Subroutine Checkpoint_Input(filename)
    Implicit None
    Character(256), Intent(In) :: filename
    Integer :: checkpoint_unit, nr_in,nth_in,nphi_in,nv_in,istep_in
    Real*8 :: r1_in, r2_in, th1_in, th2_in, phi1_in, phi2_in
    Integer :: cinfo(6), sizes(3), subsizes(3), starts(3), icheck, ir, iv
    Integer :: fhandle, finfo, fmode, ftype, ierr, status(MPI_STATUS_SIZE)
    Integer(kind=MPI_OFFSET_KIND) :: disp, offset
    Real*8, Allocatable, Dimension(:,:,:) :: data, buf
    Character(32), Dimension(8) :: names
    Character(256) :: fname, fhead, temp, temph
    names(1) = 'rho'
    names(2) = 'w'
    names(3) = 'v'
    names(4) = 'u'
    names(5) = 's'
    names(6) = 'bt'
    names(7) = 'bp'
    names(8) = 'br'

    fhead = Trim(checkpoint_name)//Trim(filename)//'/header'
    checkpoint_unit = 11
             
    ! Get the header
    Open(checkpoint_unit, file=Trim(fhead),status='old')

    If ((do_entropy_rain).and.(.not.start_rain)) Then
       If (Bottom_Specified) Then
          Read(checkpoint_unit,fmt='(7i9)') nr_in,nth_in,nphi_in,nv_in,istep_in,nsqualls,iashtime
          Read(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1_in,r2_in,th1_in,th2_in,phi1_in,phi2_in
       Else
          Read(checkpoint_unit,fmt='(6i9)') nr_in,nth_in,nphi_in,nv_in,istep_in,nsqualls
          Read(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1_in,r2_in,th1_in,th2_in,phi1_in,phi2_in
       EndIf
    Else
       If (Bottom_Specified) Then
          Read(checkpoint_unit,fmt='(6i9)') nr_in,nth_in,nphi_in,nv_in,istep_in,iashtime
          Read(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1_in,r2_in,th1_in,th2_in,phi1_in,phi2_in
       Else
          Read(checkpoint_unit,fmt='(5i9)') nr_in,nth_in,nphi_in,nv_in,istep_in
          Read(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1_in,r2_in,th1_in,th2_in,phi1_in,phi2_in
       EndIf
    EndIf
    Close(checkpoint_unit)

    begiter = istep_in
    
    ! Check sizes
    icheck = Int(abs(th1-th1_in)*180d0/pi)
    icheck = Max(icheck,Int(abs(phi1-phi1_in)*180d0/pi))
    icheck = Max(icheck,Int(abs(th2-th2_in)*180d0/pi))
    icheck = Max(icheck,Int(abs(phi2-phi2_in)*180d0/pi))

    !If ((nr .ne. nr_in).or.(icheck .gt. 0)) Then
    If (nr .ne. nr_in) Then
       Print*, 'nr=', nr, 'nr_in=', nr_in
       Print*, 'Use IDL routines to resize radial domain.'
       Call Graceful_Exit()
    ElseIf ((nth .ne. nth_in).or.(nphi .ne. nphi_in).or.(icheck .ne. 0)) Then
       Print*, 'Doing horizontal interpolation.'
       Call Checkpoint_Interpolate(filename,nphi_in,nth_in)
       deltat = Min(dth/odth,dphi/odphi)*deltat
       Call Make_Checkpoint_Filename(checkpoint_name,temp,temph,.True.)
       checkpoint_filename = Trim(temp)
       checkpoint_header_filename = Trim(temph)
       Call Checkpoint_Output()
    Else
       If (mycheck_rank .eq. 0) Then
          Allocate(data(mynth,nphi,nr))
       Else
          Allocate(data(mynth,mynphi,mynr))
       EndIf
    
       Do iv=1,min(nv,nv_in)
          Call Barrier(check_comm)
          If (mycw_rank .ne. MPI_PROC_NULL) Then
             fname = Trim(checkpoint_name)//Trim(filename)//'/checkpoint.'//Trim(names(iv))
             fmode = MPI_MODE_RDONLY
             finfo = MPI_INFO_NULL
             sizes = (/nth,nphi,nr/)
             !subsizes = (/nth/nnt,nphi,nr/)
             !starts = (/mycw_rank*nth/nnt,0,0/)
             subsizes = (/mynth,nphi,nr/)
             starts = (/myth1-1,0,0/)
             Call MPI_FILE_OPEN(cw_comm,Trim(fname),fmode,finfo,fhandle,ierr)
             Call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,ftype,ierr)
             Call MPI_TYPE_COMMIT(ftype,ierr)
             disp = 0
             Call MPI_FILE_SET_VIEW(fhandle,disp,MPI_REAL8,ftype,'native',finfo,ierr)
             offset = disp
             Call MPI_FILE_READ_AT_ALL(fhandle,offset,data,nr*mynth*nphi,MPI_REAL8,status,ierr)
             CALL MPI_FILE_CLOSE(fhandle,ierr)
             Call Barrier(cw_comm)
          EndIf

          If (mycheck_rank .ne. 0) Then
             !Send chunk size info to rank 1
             cinfo = (/myr1,myr2,mynr,myphi1,myphi2,mynphi/)
             Call Send(cinfo,0,mycheck_rank,check_comm)
          EndIf
          
          If (mycheck_rank .eq. 0) Then
             vars0(:,:,:,iv) = data(:,myphi1:myphi2,myr1:myr2)
             Do ir=1,nnr*nnp-1
                cinfo = 0
                Call Receive(cinfo,ir,ir,check_comm)
                Allocate(buf(mynth,cinfo(6),cinfo(3)))
                buf = data(:,cinfo(4):cinfo(5),cinfo(1):cinfo(2))
                Call Send(buf,ir,ir,check_comm)
                Deallocate(buf)
             EndDo
          EndIf
          
          If (mycheck_rank .ne. 0) Then
             data = 0d0
             Call Receive(data,0,mycheck_rank,check_comm)
             vars0(:,:,:,iv) = data
          EndIf
       EndDo

       Deallocate(data)
    EndIf
    If ((do_entropy_rain) .and. (.not. start_rain)) Then
       fname = Trim(checkpoint_name)//Trim(filename)//'/checkpoint.rain_squalls'
       Call Read_Rain_Squalls(fname)
    EndIf

  End Subroutine Checkpoint_Input

  Subroutine Checkpoint_Output()
    Implicit None
    Integer :: checkpoint_unit, cinfo(6), sizes(3), subsizes(3), starts(3), iv, ir
    Integer :: fhandle, finfo, fmode, ftype, ierr, status(MPI_STATUS_SIZE)
    Integer(kind=MPI_OFFSET_KIND) :: disp, offset
    Real*8, Allocatable, Dimension(:) :: writebuf
    Real*8, Allocatable, Dimension(:,:,:) :: data, buf
    Character(32), Dimension(8) :: names
    Character(256) :: fname
    names(1) = 'rho'
    names(2) = 'w'
    names(3) = 'v'
    names(4) = 'u'
    names(5) = 's'
    names(6) = 'bt'
    names(7) = 'bp'
    names(8) = 'br'

    checkpoint_unit = 11
    If (Bottom_Specified) Then
       Call MPI_Bcast(iashtime,1,MPI_INTEGER,0,node_comm,ierr)
    EndIf
    If (mytop_rank .eq. 0) Then
       ! Output a header
       Open(checkpoint_unit, file=Trim(checkpoint_header_filename),status='unknown')
       If (do_entropy_rain) Then
          If (Bottom_Specified) Then
             Write(checkpoint_unit,fmt='(7i9)') nr,nth,nphi,nv,istep,nsqualls,iashtime
             Write(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1,r2,th1,th2,phi1,phi2
          Else
             Write(checkpoint_unit,fmt='(6i9)') nr,nth,nphi,nv,istep,nsqualls
             Write(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1,r2,th1,th2,phi1,phi2
          EndIf
       Else
          If (Bottom_Specified) Then
             Write(checkpoint_unit,fmt='(6i9)') nr,nth,nphi,nv,istep,iashtime
             Write(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1,r2,th1,th2,phi1,phi2
          Else
             Write(checkpoint_unit,fmt='(5i9)') nr,nth,nphi,nv,istep
             Write(checkpoint_unit,fmt='(8e14.7)') deltat,sim_time,r1,r2,th1,th2,phi1,phi2
          EndIf
       EndIf
       Close(checkpoint_unit)
    EndIf

    If (mycheck_rank .eq. 0) Then
       Allocate(data(mynth,nphi,nr))
    Else
       Allocate(data(mynth,mynphi,mynr))
    EndIf

    Do iv=1,nv
       If (mycheck_rank .ne. 0) Then
          !Send chunk size info to rank 1
          cinfo = (/myr1,myr2,mynr,myphi1,myphi2,mynphi/)
          Call Send(cinfo,0,mycheck_rank,check_comm)
          data = vars0(:,:,:,iv)
          Call Send(data,0,mycheck_rank+nnodes,check_comm)          
       EndIf

       If (mycheck_rank .eq. 0) Then
          data = 0d0
          data(:,myphi1:myphi2,myr1:myr2) = vars0(:,:,:,iv)
          Do ir=1,nnr*nnp-1
             cinfo = 0
             Call Receive(cinfo,ir,ir,check_comm)
             Allocate(buf(mynth,cinfo(6),cinfo(3)))
             buf = 0d0
             Call Receive(buf,ir,ir+nnodes,check_comm)
             data(:,cinfo(4):cinfo(5),cinfo(1):cinfo(2)) = buf
             Deallocate(buf)
          EndDo
       EndIf

       If (mycw_rank .ne. MPI_PROC_NULL) Then
          fname = Trim(checkpoint_filename)//Trim(names(iv))
          fmode = ior(MPI_MODE_WRONLY,MPI_MODE_CREATE)
          finfo = MPI_INFO_NULL
          Call MPI_FILE_OPEN(cw_comm,Trim(fname),fmode,finfo,fhandle,ierr)
          sizes = (/nth,nphi,nr/)
          subsizes = (/mynth,nphi,nr/)
          starts = (/myth1-1,0,0/)
          Call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,ftype,ierr)
          Call MPI_TYPE_COMMIT(ftype,ierr)
          disp = 0
          Call MPI_FILE_SET_VIEW(fhandle,disp,MPI_REAL8,ftype,'native',finfo,ierr)
          offset = disp
          Call MPI_FILE_WRITE_AT_ALL(fhandle,offset,data,mynth*nphi*nr,MPI_REAL8,status,ierr)
          CALL MPI_FILE_CLOSE(fhandle,ierr)
          Call Barrier(cw_comm)
       EndIf
       Call Barrier(node_comm)
    EndDo

    Deallocate(data)
    
    If ((mytop_rank .eq. 0) .and. (do_entropy_rain)) Then
       If (nsqualls .gt. 0) Then
          Allocate(writebuf(12))
          writebuf = 0d0
          Inquire(iolength=irec_rain) writebuf
          Open(checkpoint_unit,file=Trim(Trim(checkpoint_filename)//'rain_squalls'),status='unknown',form='unformatted', access='direct',recl=irec_rain)
          Do ir=1,nsqualls
             writebuf(1)=rain_squalls(ir)%Duration
             writebuf(2)=rain_squalls(ir)%SigmaT
             writebuf(3)=rain_squalls(ir)%Start_Time
             writebuf(4)=rain_squalls(ir)%Size
             writebuf(5)=rain_squalls(ir)%Dels
             writebuf(6)=rain_squalls(ir)%Theta_loc
             writebuf(7)=rain_squalls(ir)%Phi_loc
             writebuf(8)=rain_squalls(ir)%Scaling
             writebuf(9)=Dble(rain_squalls(ir)%loct)
             writebuf(10)=Dble(rain_squalls(ir)%locp)
             writebuf(11)=rain_squalls(ir)%Vel
             writebuf(12)=rain_squalls(ir)%Delr
             Write(checkpoint_unit,rec=ir) writebuf
          EndDo
          Close(checkpoint_unit)
          Deallocate(writebuf)
       EndIf

       If (build_lanes) Then
          Allocate(writebuf(7))
          writebuf = 0d0
          Inquire(iolength=irec_rain) writebuf
          Open(checkpoint_unit,file=Trim(Trim(checkpoint_filename)//'conv_cells'),status='unknown',form='unformatted', access='direct',recl=irec_rain)
          Do ir=1,ncells
             writebuf(1)=Dble(conv_cells(ir)%loct)
             writebuf(2)=Dble(conv_cells(ir)%locp)
             writebuf(3)=conv_cells(ir)%Theta_loc
             writebuf(4)=conv_cells(ir)%Phi_loc
             writebuf(5)=conv_cells(ir)%Duration
             writebuf(6)=conv_cells(ir)%Start_Time
             writebuf(7)=conv_cells(ir)%Scaling
             Write(checkpoint_unit,rec=ir) writebuf
          EndDo
          Close(checkpoint_unit)
          Deallocate(writebuf)
       EndIf
    EndIf

    If (mytop_rank .eq. 0) Then
       ! Update the log
       Open(checkpoint_unit, file=Trim(Trim(perm_dir)//'checkpoint.log'), status='unknown', access='sequential', position='append')
       Write(checkpoint_unit, fmt='(i7)') istep
       Close(checkpoint_unit)
    EndIf

  End Subroutine Checkpoint_Output

End Module Checkpoint
