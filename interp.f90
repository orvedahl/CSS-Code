Module Interpolation
  Use Derivatives
  Implicit None 
  Integer :: onth, onphi !Old size
  Real*8 :: odphi, odth
  Real*8, Allocatable, Dimension(:) :: oldtheta, oldphi !Old grid
  Real*8, Allocatable, Dimension(:,:,:) :: olddata, data
  Character(256) :: perm_dir, checkpoint_name, checkpoint_filename, checkpoint_header_filename

Contains

  Subroutine Checkpoint_Interpolate(filename,xin,yin)
    Character(256), Intent(In) :: filename
    Integer, Intent(In) :: xin, yin
    Integer :: cinfo(6), sizes(3), subsizes(3), starts(3), iv, ir, it, ip
    Integer :: fhandle, finfo, fmode, ftype, ierr, status(MPI_STATUS_SIZE)
    Integer(kind=MPI_OFFSET_KIND) :: disp, offset
    Real*8, Allocatable, Dimension(:,:,:) :: buf
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

    onphi = xin
    onth = yin

    If (mycheck_rank .eq. 0) Then
       Allocate(olddata(onth,onphi,mynr),data(nth,nphi,mynr))
    Else
       Allocate(data(mynth,mynphi,mynr))
    EndIf
    
    Allocate(oldtheta(onth),oldphi(onphi))

    odth = (th2-th1)/Dble(onth-1)
    odphi = (phi2-phi1)/Dble(onphi)
    oldtheta = odth*(/(Dble(it),it=0,onth-1)/)+th1
    oldphi = odphi*(/(Dble(ip),ip=0,onphi-1)/)+phi1

    Do iv=1,nv
       olddata = 0d0
       data = 0d0
       Call Barrier(interp_comm)
       
       !Read in checkpoint
       If (myintr_rank .ne. MPI_PROC_NULL) Then
          fname = Trim(checkpoint_name)//Trim(filename)//'/checkpoint.'//Trim(names(iv))
          fmode = MPI_MODE_RDONLY
          finfo = MPI_INFO_NULL
          sizes = (/onth,onphi,nr/)
          subsizes = (/onth,onphi,mynr/)
          starts = (/0,0,myr1-1/)
          Call MPI_FILE_OPEN(intread_comm,Trim(fname),fmode,finfo,fhandle,ierr)
          If (mycw_rank .lt. nnr) Then
             Call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,ftype,ierr)
             Call MPI_TYPE_COMMIT(ftype,ierr)
             disp = 0
             Call MPI_FILE_SET_VIEW(fhandle,disp,MPI_REAL8,ftype,'native',finfo,ierr)
             offset = disp
             Call MPI_FILE_READ_AT_ALL(fhandle,offset,olddata,mynr*onth*onphi,MPI_REAL8,status,ierr)
          EndIf
          Call MPI_FILE_CLOSE(fhandle,ierr)
          
          Do ir=1,mynr
             Call Interp2d(data(:,:,ir),olddata(:,:,ir),oldphi,oldtheta,phi,theta)
          EndDo
          Call Barrier(intread_comm)
       EndIf

       If (myint_rank .ne. 0) Then
          !Send chunk size info to rank 1
          cinfo = (/myth1,myth2,mynth,myphi1,myphi2,mynphi/)
          Call Send(cinfo,0,myint_rank,interp_comm)
       EndIf
          
       If (myint_rank .eq. 0) Then
          vars0(:,:,:,iv) = data(myth1:myth2,myphi1:myphi2,:)
          Do ir=1,nnt*nnp-1
             cinfo = 0
             Call Receive(cinfo,ir,ir,interp_comm)
             Allocate(buf(mynth,cinfo(6),cinfo(3)))
             buf = data(:,cinfo(4):cinfo(5),cinfo(1):cinfo(2))
             Call Send(buf,ir,ir,interp_comm)
             Deallocate(buf)
          EndDo
       EndIf
          
       If (myint_rank .ne. 0) Then
          data = 0d0
          Call Receive(data,0,myint_rank,interp_comm)
          vars0(:,:,:,iv) = data
       EndIf
    EndDo
    Deallocate(olddata,data,oldtheta,oldphi)
  End Subroutine Checkpoint_Interpolate

  Subroutine Interp2d(newarr, arr, ox, oy, nx, ny)
    Real*8, Dimension(:,:), Intent(InOut) :: newarr !Interpolated data
    Real*8, Dimension(:,:), Intent(In) :: arr !Old data
    Real*8, Dimension(:), Intent(In) :: ox, oy !Old grid
    Real*8, Dimension(:), Intent(In) :: nx, ny !New grid
    Integer :: onx, ony !Size of ox and oy
    Integer :: nnx, nny !Size of nx and ny
    Integer :: ix(4), iy(4)  !Integer storage array for locations of 16 nearest-neighbors on the old grid
    Real*8 :: odx, ody !Grid distance (assumes uniform)
    Real*8 :: xc, yc !Current location on the new grid
    Real*8 :: mu, mu2, mu3 !Normalized current grid coordinate
    Real*8 ::  akm1, ak, akp1, akp2 !Storage for old data
    Real*8 :: a0, a1, a2, a3 !Interpolant coefficients
    Real*8 :: ry(4) !Array temporary after x interpolation
    Integer :: i1, i2, ij !loop variables
  
    nnx = Size(nx)
    nny = Size(ny)
    
    onx = Size(ox)
    ony = Size(oy)
    
    odx = ox(onx/2)-ox(onx/2-1)
    ody = oy(ony/2)-oy(ony/2-1)
    
    Do i2=1,nny
       yc = ny(i2)
       Do i1=1,nnx
          xc = nx(i1)        
          !For bicubic convolution interpolation, need nearest 2 points from old grid
          !  in each direction. Find them in O(1) time. (From Ben G. 2012)
          ix(2) = Mod(Floor((xc-ox(1))/odx),onx)+1
          iy(2) = Max(Floor((yc-oy(1))/ody), 1)
          ix(3) = Mod(ix(2)+1,onx)+1
          iy(3) = Min(iy(2)+1,ony)
          ix(1) = ix(2)-1
          If (ix(1) .lt. 1) ix(1) = onx+ix(1)
          iy(1) = Max(iy(2)-1,1)
          ix(4) = Mod(ix(2)+2,onx)+1
          iy(4) = Min(iy(2)+2,ony)
          
          mu = (xc-ox(ix(2)))/(ox(ix(3))-ox(ix(2)))
          mu2 = mu*mu
          mu3 = mu*mu2
          Do ij=1,4
             akm1 = arr(iy(ij),ix(1))
             ak   = arr(iy(ij),ix(2))
             akp1 = arr(iy(ij),ix(3))
             akp2 = arr(iy(ij),ix(4))
             
             a3 = ak 
             a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
             a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
             a2 = -0.5d0*akm1 + 0.5d0*akp1
             ry(ij) = a0*mu3+a1*mu2+a2*mu+a3           
          EndDo
          
          mu = (yc-oy(iy(2)))/(oy(iy(3))-oy(iy(2)))
          mu2 = mu*mu
          mu3 = mu*mu2
          
          akm1 = ry(1)
          ak   = ry(2)
          akp1 = ry(3)
          akp2 = ry(4)
          
          a3 = ak 
          a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
          a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
          a2 = -0.5d0*akm1 + 0.5d0*akp1
          newarr(i2,i1) = a0*mu3+a1*mu2+a2*mu+a3     
       EndDo
    EndDo
  End Subroutine Interp2d

End Module Interpolation
