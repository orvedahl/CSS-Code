Module Parallel

  !  This is where we keep all of the MPI calls.  In case we want to use a
  !  new message passing scheme, everything that needs to be changed is
  !  conveniently located right here!
  Use Constants
  Implicit None

  External Sum_Max_Op

  Include 'mpif.h'
 
  Type Boundary_Size
     Integer :: bnr,bnt,bnp,br1,br2,bt1,bt2,bp1,bp2
  End Type Boundary_Size

  Type Boundary_Info
     Type(Boundary_Size) :: bsize
     Real*8, Allocatable :: bdata(:,:,:,:)
  End Type Boundary_Info

  !Communicators
  Integer :: node_comm, top_comm, bot_comm, theta_top_comm, theta_bot_comm
  Integer :: check_comm, cw_comm, interp_comm, intread_comm, radial_comm, exchange_req(4)
  !Node Ranks
  Integer :: myrank, myradialrank, mythetarank, myphirank
  Integer :: mybot_rank, mytop_rank, mytheta_bot_rank, mytheta_top_rank
  Integer :: mycheck_rank, mycw_rank, myint_rank, myintr_rank, myrr_rank
  !MPI Cartesian Network Topology
  Integer :: cartnbrs(3,2), cartcoords(3), node_lup(3,2)
  !Number of CPUs and MPI err handle
  Integer :: nnodes, ierr, MPI_sum_max_op
  !Number and size of nearest-neighbor exchanges
  Integer :: nxchg(3), myxchgsize(3), myxchgsizesng(3), myxchgsize2(2), myxchgsizemag(3)
  !Integer subdomain sizes on node
  Integer :: mynr,myr1,myr2,mynth,myth1,myth2,mynphi,myphi1,myphi2

  !Boundary condition data and dimensions
  Type(Boundary_Info), Dimension(3,2) :: binfo, dbinfo, tbcinfo
  !Primary variable storage arrays
  Real*8, Target, Allocatable, Dimension(:,:,:,:) :: vars0, vars1, vars2, dvars1, dvars2, maga
  !Auxilliary arrays
  Real*8, Target, Allocatable, Dimension(:,:,:) :: temperature, invrho
  Real*8, Dimension(:,:), Allocatable :: exchange_snd, exchange_rcv

  !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: vars0, vars1, vars2, dvars1, dvars2, maga, temperature, invrho

  Interface assignment(=)
     Module Procedure Bsize_Assignment
  End Interface

  Interface Global_AllReduce
    Module Procedure GAR_0D,GAR_1D,GAR_2D,GAR_3D, GAR_INT_1D, GAR_INT_2D
  End Interface

  Interface Global_Reduce
    Module Procedure GR_0D,GR_1D,GR_2D,GR_3D
  End Interface

  Interface Send
    Module Procedure Send_0D, Send_1D, Send_2D, Send_3D, Send_4D, Send_Int_0D, Send_Int_1D,Send_Bool_1D
  End Interface

  Interface ISend
     Module Procedure ISend_3D, ISend_4D
  End Interface

  Interface Receive
    Module Procedure Receive_0D, Receive_1D, Receive_2D, Receive_3D, Receive_4D, Receive_Int_0D, Receive_Int_1D,Recv_Bool_1D
  End Interface

  Interface IReceive
     Module Procedure IReceive_3D, IReceive_4D
  End Interface

  Interface SendRecv
     Module Procedure SendRecv_4D, SendRecv_3D
  End Interface

Contains

  !------------------------------------------------------------------------

  Subroutine Initialize_Parallel_Part1()
    Call MPI_INIT(ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nnodes, ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr) !Main group is MPI_COMM_WORLD
    node_comm = MPI_COMM_WORLD
  End Subroutine Initialize_Parallel_Part1

  Subroutine Initialize_Parallel_Part2()
    Integer :: color, key
    Integer :: dummy, dims(3)
    Logical :: periodic(3)
    dims = (/nnt,nnp,nnr/)
    periodic = (/.False.,.True.,.False./)
    Call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periodic,1,node_comm,ierr)
    Call MPI_COMM_RANK(node_comm,myrank,ierr)
    Call MPI_CART_COORDS(node_comm,myrank,3,cartcoords,ierr)
    !Up in theta
    dummy = myrank
    Call MPI_CART_SHIFT(node_comm,0,-1,dummy,cartnbrs(2,1),ierr)
    !Down in theta
    dummy = myrank
    Call MPI_CART_SHIFT(node_comm,0,1,dummy,cartnbrs(2,2),ierr)
    !Down in phi
    dummy = myrank
    Call MPI_CART_SHIFT(node_comm,1,-1,dummy,cartnbrs(3,1),ierr)
    !Up in phi
    dummy = myrank
    Call MPI_CART_SHIFT(node_comm,1,1,dummy,cartnbrs(3,2),ierr)
    !Down in r
    dummy = myrank
    Call MPI_CART_SHIFT(node_comm,2,-1,dummy,cartnbrs(1,1),ierr)
    !Up in r
    dummy = myrank
    Call MPI_CART_SHIFT(node_comm,2,1,dummy,cartnbrs(1,2),ierr)
    
    Call Partition()
    Call Node_Partition()

    !Split off the bottom boundary communicator
    If (myr1 .eq. 1) Then
       color = 1
       key = mythetarank+nnt*myphirank
    Else
       color = MPI_UNDEFINED
       key = MPI_UNDEFINED
    EndIf
    Call MPI_COMM_SPLIT(node_comm,color,key,bot_comm,ierr)
    If (myr1 .eq. 1) Then
       Call MPI_COMM_RANK(bot_comm, mybot_rank, ierr)
    Else
       mybot_rank = -1
    EndIf

    !Split off the top boundary communicator
    If (myr2 .eq. nr) Then
       color = 2
       key = mythetarank+nnt*myphirank
    Else
       color = MPI_UNDEFINED
       key = MPI_UNDEFINED
    EndIf
    Call MPI_COMM_SPLIT(node_comm,color,key,top_comm,ierr)
    If (myr2 .eq. nr) Then
       Call MPI_COMM_RANK(top_comm, mytop_rank, ierr)
    Else
       mytop_rank = -1
    EndIf

    !Split off the theta bottom boundary communicator
    If (myth1 .eq. 1) Then
       color = 1
       key = myradialrank+nnr*myphirank
    Else
       color = MPI_UNDEFINED
       key = MPI_UNDEFINED
    EndIf
    Call MPI_COMM_SPLIT(node_comm,color,key,theta_bot_comm,ierr)
    If (myth1 .eq. 1) Then
       Call MPI_COMM_RANK(theta_bot_comm, mytheta_bot_rank, ierr)
    Else
       mytheta_bot_rank = -1
    EndIf

    !Split off the theta top boundary communicator
    If (myth2 .eq. nth) Then
       color = 2
       key = myradialrank+nnr*myphirank
    Else
       color = MPI_UNDEFINED
       key = MPI_UNDEFINED
    EndIf
    Call MPI_COMM_SPLIT(node_comm,color,key,theta_top_comm,ierr)
    If (myth2 .eq. nth) Then
       Call MPI_COMM_RANK(theta_top_comm, mytheta_top_rank, ierr)
    Else
       mytheta_top_rank = -1
    EndIf

    !Split off theta communicator for checkpointing
    color = mythetarank
    key = myradialrank + nnr*myphirank
    Call MPI_COMM_SPLIT(node_comm,color,key,check_comm,ierr)
    Call MPI_COMM_RANK(check_comm, mycheck_rank, ierr)

    !Split off radial communicator for checkpointing
    color = myradialrank
    key = mythetarank + nnt*myphirank
    Call MPI_COMM_SPLIT(node_comm,color,key,radial_comm,ierr)
    Call MPI_COMM_RANK(radial_comm, myrr_rank, ierr)

    !Split off communicator for writing checkpoint
    If (mycheck_rank .eq. 0) Then
       color = 3
       key = mythetarank
    Else
       color = MPI_UNDEFINED
       key = MPI_UNDEFINED
    EndIf
    Call MPI_COMM_SPLIT(node_comm,color,key,cw_comm,ierr)
    If (mycheck_rank .eq. 0) Then
       Call MPI_COMM_RANK(cw_comm, mycw_rank, ierr)
    Else
       mycw_rank = MPI_PROC_NULL
    EndIf

    !Split off radial communicator for interpolator
    color = myradialrank
    key = mythetarank + nnt*myphirank
    Call MPI_COMM_SPLIT(node_comm,color,key,interp_comm,ierr)
    Call MPI_COMM_RANK(interp_comm, myint_rank, ierr)

    !Split off communicator for reading checkpoint (interpolator)
    If (myint_rank .eq. 0) Then
       color = 4
       key = myradialrank
    Else
       color = MPI_UNDEFINED
       key = MPI_UNDEFINED
    EndIf
    Call MPI_COMM_SPLIT(node_comm,color,key,intread_comm,ierr)
    If (myint_rank .eq. 0) Then
       Call MPI_COMM_RANK(intread_comm, myintr_rank, ierr)
    Else
       myintr_rank = MPI_PROC_NULL
    EndIf

    Call MPI_Op_Create(Sum_Max_Op,.True.,MPI_sum_max_op,ierr)

  End Subroutine Initialize_Parallel_Part2

  !------------------------------------------------------------------------

  Subroutine Partition()
    Integer :: irem, ii
    Integer, Allocatable, Dimension(:) :: myd1, myd2, mynd

    !Determine node data bounds
    !myradialrank = Mod(myrank/nnt,nnr)
    !mythetarank = Mod(myrank,nnt)
    !myphirank = Mod(myrank/nnr/nnt,nnp)

    myradialrank = cartcoords(3)
    mythetarank = cartcoords(1)
    myphirank = cartcoords(2)
    
    mynr = nr/nnr
    If (mynr .lt. 9) Then 
       Print*, 'Too few interior points in r, number of radial nodes is too large.'
       Call Graceful_Exit()
    EndIf

    irem = Mod(nr,nnr) ! guaranteed to be less than nnr, so add one to each node less than nnr
    !Distrbute any remaining points over iremr nodes for load balancing
    If (irem .gt. 0) Then
       Allocate(myd1(nnr),myd2(nnr),mynd(nnr))
       myd1(1) = 1
       myd2(1) = mynr+1
       mynd(1) = mynr+1
       Do ii=2,nnr
          If (ii .le. irem) Then
             mynd(ii) = mynr+1
          Else
             mynd(ii) = mynr
          EndIf
          myd1(ii) = myd2(ii-1)+1
          myd2(ii) = myd1(ii)+mynd(ii)-1
       EndDo
       myr1 = myd1(myradialrank+1)
       myr2 = myd2(myradialrank+1)
       mynr = mynd(myradialrank+1)
       Deallocate(myd1,myd2,mynd)
    Else
       myr1 = myradialrank*mynr+1
       myr2 = myr1+mynr-1
    EndIf

    mynth = nth/nnt
    If (mynth .lt. 9) Then 
       Print*, 'Too few interior points in theta, number of theta nodes is too large.'
       Call Graceful_Exit()
    EndIf

    irem = Mod(nth,nnt) ! guaranteed to be less than nnr, so add one to each node less than nnr
    !Distrbute any remaining points over iremr nodes for load balancing
    If (irem .gt. 0) Then
       Allocate(myd1(nnt),myd2(nnt),mynd(nnt))
       myd1(1) = 1
       myd2(1) = mynth+1
       mynd(1) = mynth+1
       Do ii=2,nnt
          If (ii .le. irem) Then
             mynd(ii) = mynth+1
          Else
             mynd(ii) = mynth
          EndIf
          myd1(ii) = myd2(ii-1)+1
          myd2(ii) = myd1(ii)+mynd(ii)-1
       EndDo
       myth1 = myd1(mythetarank+1)
       myth2 = myd2(mythetarank+1)
       mynth = mynd(mythetarank+1)
       Deallocate(myd1,myd2,mynd)
    Else
       myth1 = mythetarank*mynth+1
       myth2 = myth1+mynth-1
    EndIf
    
    mynphi = nphi/nnp
    If (mynphi .lt. 9) Then 
       Print*, 'Too few interior points in phi, number of phi nodes is too large.'
       Call Graceful_Exit()
    EndIf

    irem = Mod(nphi,nnp) ! guaranteed to be less than nnr, so add one to each node less than nnr
    !Distrbute any remaining points over iremr nodes for load balancing
   
    If (irem .gt. 0) Then
       Allocate(myd1(nnp),myd2(nnp),mynd(nnp))
       myd1(1) = 1
       myd2(1) = mynphi+1
       mynd(1) = mynphi+1
       Do ii=2,nnp
          If (ii .le. irem) Then
             mynd(ii) = mynphi+1
          Else
             mynd(ii) = mynphi
          EndIf
          myd1(ii) = myd2(ii-1)+1
          myd2(ii) = myd1(ii)+mynd(ii)-1
       EndDo
       myphi1 = myd1(myphirank+1)
       myphi2 = myd2(myphirank+1)
       mynphi = mynd(myphirank+1)
       Deallocate(myd1,myd2,mynd)
    Else
       myphi1 = myphirank*mynphi+1
       myphi2 = myphi1+mynphi-1
    EndIf

  End Subroutine Partition

  Subroutine Node_Partition()
    Integer :: i,j
    Integer :: m1
    !Determine if a true boundary (on the edge of the domain)
    !Find nearest neighbors for boundary exchange

    If (myr1 .eq. 1) Then       
       Allocate(binfo(1,1)%bdata(mynth,mynphi,1,nv))       
       binfo(1,1)%bsize%bnr=1
       binfo(1,1)%bsize%bnt=mynth
       binfo(1,1)%bsize%bnp=mynphi
       binfo(1,1)%bsize%br1=1
       binfo(1,1)%bsize%br2=1
       binfo(1,1)%bsize%bt1=1
       binfo(1,1)%bsize%bt2=mynth
       binfo(1,1)%bsize%bp1=1
       binfo(1,1)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(1,1)%bdata(mynth,mynphi,1,1))
          tbcinfo(1,1)%bsize = binfo(1,1)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(1,1)%bdata(mynth,mynphi,1,5))
          dbinfo(1,1)%bsize = binfo(1,1)%bsize
       EndIf
    Else
       Allocate(binfo(1,1)%bdata(mynth,mynphi,4,nv))
       binfo(1,1)%bsize%bnr=4
       binfo(1,1)%bsize%bnt=mynth
       binfo(1,1)%bsize%bnp=mynphi
       binfo(1,1)%bsize%br1=1
       binfo(1,1)%bsize%br2=4
       binfo(1,1)%bsize%bt1=1
       binfo(1,1)%bsize%bt2=mynth
       binfo(1,1)%bsize%bp1=1
       binfo(1,1)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(1,1)%bdata(mynth,mynphi,4,1))
          tbcinfo(1,1)%bsize = binfo(1,1)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(1,1)%bdata(mynth,mynphi,4,5))
          dbinfo(1,1)%bsize = binfo(1,1)%bsize
       EndIf
    EndIf

    If (myr2 .eq. nr) Then
       Allocate(binfo(1,2)%bdata(mynth,mynphi,1,nv))
       binfo(1,2)%bsize%bnr=1
       binfo(1,2)%bsize%bnt=mynth
       binfo(1,2)%bsize%bnp=mynphi
       binfo(1,2)%bsize%br1=mynr
       binfo(1,2)%bsize%br2=mynr
       binfo(1,2)%bsize%bt1=1
       binfo(1,2)%bsize%bt2=mynth
       binfo(1,2)%bsize%bp1=1
       binfo(1,2)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(1,2)%bdata(mynth,mynphi,1,1))
          tbcinfo(1,2)%bsize = binfo(1,2)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(1,2)%bdata(mynth,mynphi,1,5))
          dbinfo(1,2)%bsize = binfo(1,2)%bsize
       EndIf
    Else
       Allocate(binfo(1,2)%bdata(mynth,mynphi,4,nv))
       binfo(1,2)%bsize%bnr=4
       binfo(1,2)%bsize%bnt=mynth
       binfo(1,2)%bsize%bnp=mynphi
       binfo(1,2)%bsize%br1=mynr-3
       binfo(1,2)%bsize%br2=mynr
       binfo(1,2)%bsize%bt1=1
       binfo(1,2)%bsize%bt2=mynth
       binfo(1,2)%bsize%bp1=1
       binfo(1,2)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(1,2)%bdata(mynth,mynphi,4,1))
          tbcinfo(1,2)%bsize = binfo(1,2)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(1,2)%bdata(mynth,mynphi,4,5))
          dbinfo(1,2)%bsize = binfo(1,2)%bsize
       EndIf
    EndIf

    If (myth1 .eq. 1) Then
       Allocate(binfo(2,1)%bdata(1,mynphi,mynr,nv))
       binfo(2,1)%bsize%bnr=mynr
       binfo(2,1)%bsize%bnt=1
       binfo(2,1)%bsize%bnp=mynphi
       binfo(2,1)%bsize%br1=1
       binfo(2,1)%bsize%br2=mynr
       binfo(2,1)%bsize%bt1=1
       binfo(2,1)%bsize%bt2=1
       binfo(2,1)%bsize%bp1=1
       binfo(2,1)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(2,1)%bdata(1,mynphi,mynr,1))
          tbcinfo(2,1)%bsize = binfo(2,1)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(2,1)%bdata(1,mynphi,mynr,5))
          dbinfo(2,1)%bsize = binfo(2,1)%bsize
       EndIf
    Else
       Allocate(binfo(2,1)%bdata(4,mynphi,mynr,nv))
       binfo(2,1)%bsize%bnr=mynr
       binfo(2,1)%bsize%bnt=4
       binfo(2,1)%bsize%bnp=mynphi
       binfo(2,1)%bsize%br1=1
       binfo(2,1)%bsize%br2=mynr
       binfo(2,1)%bsize%bt1=1
       binfo(2,1)%bsize%bt2=4
       binfo(2,1)%bsize%bp1=1
       binfo(2,1)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(2,1)%bdata(4,mynphi,mynr,1))
          tbcinfo(2,1)%bsize = binfo(2,1)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(2,1)%bdata(4,mynphi,mynr,5))
          dbinfo(2,1)%bsize = binfo(2,1)%bsize
       EndIf
    EndIf

    If (myth2 .eq. nth) Then
       If (Theta_Symmetric) Then
          Allocate(binfo(2,2)%bdata(4,mynphi,mynr,nv))
          binfo(2,2)%bsize%bnr=mynr
          binfo(2,2)%bsize%bnt=4
          binfo(2,2)%bsize%bnp=mynphi
          binfo(2,2)%bsize%br1=1
          binfo(2,2)%bsize%br2=mynr
          binfo(2,2)%bsize%bt1=mynth
          binfo(2,2)%bsize%bt2=mynth+4
          binfo(2,2)%bsize%bp1=1
          binfo(2,2)%bsize%bp2=mynphi
          If (do_radheat) Then 
             Allocate(tbcinfo(2,2)%bdata(4,mynphi,mynr,1))
             tbcinfo(2,2)%bsize = binfo(2,2)%bsize             
          EndIf
          If (magnetic) Then
             Allocate(dbinfo(2,2)%bdata(4,mynphi,mynr,5))
             dbinfo(2,2)%bsize = binfo(2,2)%bsize
          EndIf
       Else
          Allocate(binfo(2,2)%bdata(1,mynphi,mynr,nv))
          binfo(2,2)%bsize%bnr=mynr
          binfo(2,2)%bsize%bnt=1
          binfo(2,2)%bsize%bnp=mynphi
          binfo(2,2)%bsize%br1=1
          binfo(2,2)%bsize%br2=mynr
          binfo(2,2)%bsize%bt1=mynth
          binfo(2,2)%bsize%bt2=mynth
          binfo(2,2)%bsize%bp1=1
          binfo(2,2)%bsize%bp2=mynphi
          If (do_radheat) Then 
             Allocate(tbcinfo(2,2)%bdata(1,mynphi,mynr,1))
             tbcinfo(2,2)%bsize = binfo(2,2)%bsize             
          EndIf
          If (magnetic) Then
             Allocate(dbinfo(2,2)%bdata(1,mynphi,mynr,5))
             dbinfo(2,2)%bsize = binfo(2,2)%bsize
          EndIf
       EndIf
    Else
       Allocate(binfo(2,2)%bdata(4,mynphi,mynr,nv))
       binfo(2,2)%bsize%bnr=mynr
       binfo(2,2)%bsize%bnt=4
       binfo(2,2)%bsize%bnp=mynphi
       binfo(2,2)%bsize%br1=1
       binfo(2,2)%bsize%br2=mynr
       binfo(2,2)%bsize%bt1=mynth-3
       binfo(2,2)%bsize%bt2=mynth
       binfo(2,2)%bsize%bp1=1
       binfo(2,2)%bsize%bp2=mynphi
       If (do_radheat) Then 
          Allocate(tbcinfo(2,2)%bdata(4,mynphi,mynr,1))
          tbcinfo(2,2)%bsize = binfo(2,2)%bsize
       EndIf
       If (magnetic) Then
          Allocate(dbinfo(2,2)%bdata(4,mynphi,mynr,5))
          dbinfo(2,2)%bsize = binfo(2,2)%bsize
       EndIf
    EndIf

    Allocate(binfo(3,1)%bdata(mynth,4,mynr,nv))
    binfo(3,1)%bsize%bnr=mynr
    binfo(3,1)%bsize%bnt=mynth
    binfo(3,1)%bsize%bnp=4
    binfo(3,1)%bsize%br1=1
    binfo(3,1)%bsize%br2=mynr
    binfo(3,1)%bsize%bt1=1
    binfo(3,1)%bsize%bt2=mynth
    binfo(3,1)%bsize%bp1=1
    binfo(3,1)%bsize%bp2=4
    If (do_radheat) Then 
       Allocate(tbcinfo(3,1)%bdata(mynth,4,mynr,1))
       tbcinfo(3,1)%bsize = binfo(3,1)%bsize
    EndIf
    If (magnetic) Then
       Allocate(dbinfo(3,1)%bdata(mynth,4,mynr,5))
       dbinfo(3,1)%bsize = binfo(3,1)%bsize
    EndIf

    Allocate(binfo(3,2)%bdata(mynth,4,mynr,nv))
    binfo(3,2)%bsize%bnr=mynr
    binfo(3,2)%bsize%bnt=mynth
    binfo(3,2)%bsize%bnp=4
    binfo(3,2)%bsize%br1=1
    binfo(3,2)%bsize%br2=mynr
    binfo(3,2)%bsize%bt1=1
    binfo(3,2)%bsize%bt2=mynth
    binfo(3,2)%bsize%bp1=mynphi-3
    binfo(3,2)%bsize%bp2=mynphi
    If (do_radheat) Then 
       Allocate(tbcinfo(3,2)%bdata(mynth,4,mynr,1))
       tbcinfo(3,2)%bsize = binfo(3,2)%bsize
    EndIf
    If (magnetic) Then
       Allocate(dbinfo(3,2)%bdata(mynth,4,mynr,5))
       dbinfo(3,2)%bsize = binfo(3,2)%bsize
    EndIf
 
    node_lup = cartnbrs

    Do i=1,3
       Do j=1,2
          If (node_lup(i,j) .eq. MPI_PROC_NULL) Then
             node_lup(i,j) = -1
          EndIf
          binfo(i,j)%bdata=0d0
          If (do_radheat) Then
             tbcinfo(i,j)%bdata=0d0
             binfo(i,j)%bdata=0d0
          EndIf
          If (magnetic) dbinfo(i,j)%bdata=0d0
       EndDo
    EndDo
    
    !Determine number of exchanges
    nxchg = 0
    myxchgsize = 0
    myxchgsize2 = 0
    myxchgsizemag = 0
    myxchgsizesng = 0
    Do i=1,3
       Do j=1,2
          m1 = cartnbrs(i,j)
          If (m1 .ne. MPI_PROC_NULL) Then
             nxchg(i)=nxchg(i)+1
          EndIf
          myxchgsize(i) = max(myxchgsize(i),product(shape(binfo(i,j)%bdata)))
          If (magnetic) Then
             myxchgsizemag(i) = max(myxchgsizemag(i),product(shape(dbinfo(i,j)%bdata(:,:,:,1:3))))
             myxchgsizesng(i) = max(myxchgsizesng(i),product(shape(binfo(i,j)%bdata(:,:,:,1))))
          EndIf
          If (Laplacian) Then
             If (i .eq. 1) Then
                myxchgsize2(i) = max(myxchgsize2(i),product(shape(dbinfo(i,j)%bdata(:,:,:,1:4))))
             ElseIf (i.eq.2) Then
                myxchgsize2(i) = max(myxchgsize2(i),product(shape(dbinfo(i,j)%bdata(:,:,:,4:5))))
             EndIf
          EndIf
       EndDo
    EndDo

  End Subroutine Node_Partition

  !------------------------------------------------------------------------

  Subroutine Bsize_Assignment(b2,b1)
    Type(Boundary_Size), Intent(In) :: b1
    Type(Boundary_Size), Intent(Out) :: b2
    b2%bnr=b1%bnr
    b2%bnt=b1%bnt
    b2%bnp=b1%bnp
    b2%br1=b1%br1
    b2%br2=b1%br2
    b2%bt1=b1%bt1
    b2%bt2=b1%bt2
    b2%bp1=b1%bp1
    b2%bp2=b1%bp2
  End Subroutine Bsize_Assignment

  !------------------------------------------------------------------------

  Subroutine Boundary_Exchange_Init(direction)
    Integer, Intent(In) :: direction
    Integer :: i,j,k,shp(4)
    Integer :: dest_node, nsize, nsz(1)

    Allocate(exchange_snd(myxchgsize(direction),nxchg(direction)),exchange_rcv(myxchgsize(direction),nxchg(direction)))
    exchange_snd = 0d0
    exchange_rcv = 0d0
    exchange_req = MPI_REQUEST_NULL

    i = direction
    j = 1
    k = 1
    dest_node = cartnbrs(i,j)
    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(binfo(i,j)%bdata)
       nsz(1) = product(shp)
       nsize = nsz(1)
       Call MPI_IRECV(exchange_rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,exchange_req(1),ierr)
       exchange_snd(1:nsize,k) = reshape(vars1(binfo(i,j)%bsize%bt1:binfo(i,j)%bsize%bt2,binfo(i,j)%bsize%bp1:binfo(i,j)%bsize%bp2,binfo(i,j)%bsize%br1:binfo(i,j)%bsize%br2,:),nsz)
       Call MPI_ISEND(exchange_snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,exchange_req(2),ierr)
       k=k+1
    EndIf

    j = 2
    dest_node = cartnbrs(i,j)
    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(binfo(i,j)%bdata)
       nsz(1) = product(shp)
       nsize = nsz(1)
       Call MPI_IRECV(exchange_rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,exchange_req(3),ierr)
       exchange_snd(1:nsize,k) = reshape(vars1(binfo(i,j)%bsize%bt1:binfo(i,j)%bsize%bt2,binfo(i,j)%bsize%bp1:binfo(i,j)%bsize%bp2,binfo(i,j)%bsize%br1:binfo(i,j)%bsize%br2,:),nsz)
       Call MPI_ISEND(exchange_snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,exchange_req(4),ierr)
    EndIf

  End Subroutine Boundary_Exchange_Init

  Subroutine Boundary_Exchange_Finish(direction)
    Integer, Intent(In) :: direction
    Integer :: i,j,k,shp(4)
    Integer :: dest_node, nsize, status(MPI_STATUS_SIZE,4)

    i = direction
    k = 1
    j = 1
    dest_node = cartnbrs(i,j)
    Call MPI_WAITALL(4,exchange_req,status,ierr)
    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(binfo(i,j)%bdata)
       nsize = product(shp)
       binfo(i,j)%bdata = reshape(exchange_rcv(1:nsize,k),shp)
       k = k+1
    EndIf

    j = 2
    dest_node = cartnbrs(i,j)

    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(binfo(i,j)%bdata)
       nsize = product(shp)
       binfo(i,j)%bdata = reshape(exchange_rcv(1:nsize,k),shp)
    EndIf
    Deallocate(exchange_snd,exchange_rcv)
  End Subroutine Boundary_Exchange_Finish

  Subroutine Boundary_Exchange_Lapl(direction)
    Integer, Intent(In) :: direction
    Real*8, Allocatable, Dimension(:,:) :: snd,rcv
    Integer :: i,j,k,shp(4)
    Integer :: dest_node, req(4), nsize, nsz(1), status(MPI_STATUS_SIZE,4)
    req = MPI_REQUEST_NULL
    i = direction
    If (i.eq.1) Then
       Allocate(snd(myxchgsize2(i),nxchg(i)),rcv(myxchgsize2(i),nxchg(i)))
       snd = 0d0
       rcv = 0d0
       j = 1
       dest_node = cartnbrs(i,j)
       k = 1
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,1:4))
          nsz(1) = product(shp)
          nsize = nsz(1)          
          snd(1:nsize,k) = reshape(vars2(dbinfo(i,j)%bsize%bt1:dbinfo(i,j)%bsize%bt2,dbinfo(i,j)%bsize%bp1:dbinfo(i,j)%bsize%bp2,dbinfo(i,j)%bsize%br1:dbinfo(i,j)%bsize%br2,1:4),nsz)
          Call MPI_ISEND(snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,req(1),ierr)
          Call MPI_IRECV(rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,req(2),ierr)
          k=k+1
       EndIf

       j = 2
       dest_node = cartnbrs(i,j)
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,1:4))
          nsz(1) = product(shp)
          nsize = nsz(1)
          Call MPI_IRECV(rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,req(3),ierr)          
          snd(1:nsize,k) = reshape(vars2(dbinfo(i,j)%bsize%bt1:dbinfo(i,j)%bsize%bt2,dbinfo(i,j)%bsize%bp1:dbinfo(i,j)%bsize%bp2,dbinfo(i,j)%bsize%br1:dbinfo(i,j)%bsize%br2,1:4),nsz)          
          Call MPI_ISEND(snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,req(4),ierr)
       EndIf
       
       Call MPI_WAITALL(4,req,status,ierr)

       k = 1
       j = 1
       dest_node = cartnbrs(i,j)
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,1:4))
          nsize = product(shp)
          dbinfo(i,j)%bdata(:,:,:,1:4) = reshape(rcv(1:nsize,k),shp)
          k = k+1
       EndIf
       
       j = 2
       dest_node = cartnbrs(i,j)
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,1:4))
          nsize = product(shp)
          dbinfo(i,j)%bdata(:,:,:,1:4) = reshape(rcv(1:nsize,k),shp)
       EndIf
       Deallocate(snd,rcv)
    ElseIf (i.eq.2) Then
       Allocate(snd(myxchgsize2(i),nxchg(i)),rcv(myxchgsize2(i),nxchg(i)))
       snd = 0d0
       rcv = 0d0
       j = 1
       dest_node = cartnbrs(i,j)
       k = 1
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,4:5))
          nsz(1) = product(shp)
          nsize = nsz(1)
          snd(1:nsize,k) = reshape(vars2(dbinfo(i,j)%bsize%bt1:dbinfo(i,j)%bsize%bt2,dbinfo(i,j)%bsize%bp1:dbinfo(i,j)%bsize%bp2,dbinfo(i,j)%bsize%br1:dbinfo(i,j)%bsize%br2,4:5),nsz)
          Call MPI_ISEND(snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,req(1),ierr)
          Call MPI_IRECV(rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,req(2),ierr)
          k=k+1
       EndIf

       j = 2
       dest_node = cartnbrs(i,j)
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,4:5))
          nsz(1) = product(shp)
          nsize = nsz(1)
          Call MPI_IRECV(rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,req(3),ierr)
          
          snd(1:nsize,k) = reshape(vars2(dbinfo(i,j)%bsize%bt1:dbinfo(i,j)%bsize%bt2,dbinfo(i,j)%bsize%bp1:dbinfo(i,j)%bsize%bp2,dbinfo(i,j)%bsize%br1:dbinfo(i,j)%bsize%br2,4:5),nsz)
          
          Call MPI_ISEND(snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,req(4),ierr)
       EndIf
       
       Call MPI_WAITALL(4,req,status,ierr)

       k = 1
       j = 1
       dest_node = cartnbrs(i,j)
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,4:5))
          nsize = product(shp)
          dbinfo(i,j)%bdata(:,:,:,4:5) = reshape(rcv(1:nsize,k),shp)
          k = k+1
       EndIf
       
       j = 2
       dest_node = cartnbrs(i,j)
       If (dest_node .ne. MPI_PROC_NULL) Then
          shp = shape(dbinfo(i,j)%bdata(:,:,:,4:5))
          nsize = product(shp)
          dbinfo(i,j)%bdata(:,:,:,4:5) = reshape(rcv(1:nsize,k),shp)
       EndIf
       Deallocate(snd,rcv)
    EndIf
  End Subroutine Boundary_Exchange_Lapl

  !------------------------------------------------------------------------

  Subroutine Boundary_Exchange_Mag_Init(direction)
    Integer, Intent(In) :: direction
    Integer :: i,j,k,shp(4)
    Integer :: dest_node, nsize, nsz(1)

    Allocate(exchange_snd(myxchgsizemag(direction),nxchg(direction)),exchange_rcv(myxchgsizemag(direction),nxchg(direction)))
    exchange_snd = 0d0
    exchange_rcv = 0d0
    exchange_req = MPI_REQUEST_NULL

    i = direction
    j = 1
    k = 1
    dest_node = cartnbrs(i,j)
    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(dbinfo(i,j)%bdata(:,:,:,1:3))
       nsz(1) = product(shp)
       nsize = nsz(1)
       Call MPI_IRECV(exchange_rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,exchange_req(1),ierr)
       exchange_snd(1:nsize,k) = reshape(maga(dbinfo(i,j)%bsize%bt1:dbinfo(i,j)%bsize%bt2,dbinfo(i,j)%bsize%bp1:dbinfo(i,j)%bsize%bp2,dbinfo(i,j)%bsize%br1:dbinfo(i,j)%bsize%br2,1:3),nsz)
       Call MPI_ISEND(exchange_snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,exchange_req(2),ierr)
       k=k+1
    EndIf

    j = 2
    dest_node = cartnbrs(i,j)
    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(dbinfo(i,j)%bdata(:,:,:,1:3))
       nsz(1) = product(shp)
       nsize = nsz(1)
       Call MPI_IRECV(exchange_rcv(1,k),nsize,MPI_REAL8,dest_node,dest_node,node_comm,exchange_req(3),ierr)
       exchange_snd(1:nsize,k) = reshape(maga(dbinfo(i,j)%bsize%bt1:dbinfo(i,j)%bsize%bt2,dbinfo(i,j)%bsize%bp1:dbinfo(i,j)%bsize%bp2,dbinfo(i,j)%bsize%br1:dbinfo(i,j)%bsize%br2,1:3),nsz)
       Call MPI_ISEND(exchange_snd(1,k),nsize,MPI_REAL8,dest_node,myrank,node_comm,exchange_req(4),ierr)
    EndIf

  End Subroutine Boundary_Exchange_Mag_Init

  Subroutine Boundary_Exchange_Mag_Finish(direction)
    Integer, Intent(In) :: direction
    Integer :: i,j,k,shp(4)
    Integer :: dest_node, nsize, status(MPI_STATUS_SIZE,4)

    i = direction
    k = 1
    j = 1
    dest_node = cartnbrs(i,j)
    Call MPI_WAITALL(4,exchange_req,status,ierr)
    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(dbinfo(i,j)%bdata(:,:,:,1:3))
       nsize = product(shp)
       dbinfo(i,j)%bdata(:,:,:,1:3) = reshape(exchange_rcv(1:nsize,k),shp)
       k = k+1
    EndIf

    j = 2
    dest_node = cartnbrs(i,j)

    If (dest_node .ne. MPI_PROC_NULL) Then
       shp = shape(binfo(i,j)%bdata(:,:,:,1:3))
       nsize = product(shp)
       dbinfo(i,j)%bdata(:,:,:,1:3) = reshape(exchange_rcv(1:nsize,k),shp)
    EndIf
    Deallocate(exchange_snd,exchange_rcv)
  End Subroutine Boundary_Exchange_Mag_Finish

  !------------------------------------------------------------------------

  Subroutine Finalize_Parallel()
    Integer :: i,j

    Do i=1,3
       Do j=1,2
          If (Allocated(binfo(i,j)%bdata)) Deallocate(binfo(i,j)%bdata)
          If (Allocated(tbcinfo(i,j)%bdata)) Deallocate(tbcinfo(i,j)%bdata)
       EndDo
    EndDo
  End Subroutine Finalize_Parallel

  Subroutine Graceful_Exit()
    Call MPI_FINALIZE(ierr)
    Stop
  End Subroutine Graceful_Exit

  !  global reduce interface subroutines, can add more if needed
  Subroutine GAR_3D(send,recv,comm,operation)

    Real*8 :: send(:,:,:)  !  data to be reduced
    Real*8 :: recv(:,:,:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    !Call MPI_BARRIER(comm,ierr)
    Call MPI_ALLREDUCE(send(1,1,1),recv(1,1,1),n,MPI_REAL8,op,comm,ierr)    

  End Subroutine GAR_3D

  Subroutine GAR_2D(send,recv,comm,operation)

    Real*8 :: send(:,:)  !  data to be reduced
    Real*8 :: recv(:,:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    !Call MPI_BARRIER(comm,ierr)
    Call MPI_ALLREDUCE(send(1,1),recv(1,1),n,MPI_REAL8,op,comm,ierr)    

  End Subroutine GAR_2D

  Subroutine GAR_INT_1D(send,recv,comm,operation)

    Integer :: send(:)  !  data to be reduced
    Integer :: recv(:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    Call MPI_ALLREDUCE(send(1),recv(1),n,MPI_INTEGER,op,comm,ierr)    

  End Subroutine GAR_INT_1D

  Subroutine GAR_INT_2D(send,recv,comm,operation)

    Integer :: send(:,:)  !  data to be reduced
    Integer :: recv(:,:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    Call MPI_ALLREDUCE(send(1,1),recv(1,1),n,MPI_INTEGER,op,comm,ierr)    

  End Subroutine GAR_INT_2D

  Subroutine GAR_1D(send,recv,comm,operation)

    Real*8 :: send(:)  !  data to be reduced
    Real*8 :: recv(:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case ('csm') ; op=MPI_sum_max_op
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    !Call MPI_BARRIER(comm,ierr)
    Call MPI_ALLREDUCE(send(1),recv(1),n,MPI_REAL8,op,comm,ierr)    

  End Subroutine GAR_1D

  Subroutine GAR_0D(send,recv,comm,operation)

    Real*8, Intent(In) :: send  !  data to be reduced
    Real*8, Intent(Out) :: recv  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm
    Integer :: op

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    !Call MPI_BARRIER(comm,ierr)
    Call MPI_ALLREDUCE(send,recv,1,MPI_REAL8,op,comm,ierr)

  End Subroutine GAR_0D

  !  global reduce interface subroutines, can add more if needed
  Subroutine GR_3D(send,recv,root,comm,operation)

    Real*8 :: send(:,:,:)  !  data to be reduced
    Real*8 :: recv(:,:,:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer*4, Intent(In) :: comm,root
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    Call MPI_REDUCE(send(1,1,1),recv(1,1,1),n,MPI_REAL8,op,root,comm,ierr)

  End Subroutine GR_3D

  Subroutine GR_2D(send,recv,root,comm,operation)

    Real*8 :: send(:,:)  !  data to be reduced
    Real*8 :: recv(:,:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm,root
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    Call MPI_REDUCE(send(1,1),recv(1,1),n,MPI_REAL8,op,root,comm,ierr)

  End Subroutine GR_2D

  Subroutine GR_1D(send,recv,root,comm,operation)

    Real*8 :: send(:)  !  data to be reduced
    Real*8 :: recv(:)  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm,root
    Integer :: op,n

    n=Size(send)

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    Call MPI_REDUCE(send(1),recv(1),n,MPI_REAL8,op,root,comm,ierr)

  End Subroutine GR_1D

  Subroutine GR_0D(send,recv,root,comm,operation)

    Real*8, Intent(In) :: send  !  data to be reduced
    Real*8, Intent(Out) :: recv  !  result of global reduction
    Character*3, Intent(In) :: operation
    Integer, Intent(In) :: comm,root
    Integer :: op

    Select Case (operation)  !  note: mutual exclusivity is NOT checked for
      Case ('min')  ;  op=MPI_MIN
      Case ('max')  ;  op=MPI_MAX
      Case ('sum')  ;  op=MPI_SUM
      Case ('pro')  ;  op=MPI_PROD
      Case Default
        Print*,'ERROR in global reduce'
        Stop
    End Select

    Call MPI_REDUCE(send,recv,1,MPI_REAL8,op,root,comm,ierr)

  End Subroutine GR_0D

  Subroutine Send_Bool_1D(array,to,from,comm)
    Logical :: array(:)
    Integer, Intent(In) :: to, from
    Integer, Intent(In) :: comm
    Integer :: nelem
    nelem = Size(array)
    Call mpi_send(array(1),nelem,MPI_LOGICAL,to,from,comm,ierr)
  End Subroutine Send_Bool_1D

  Subroutine Send_Int_0D(array,to,from,comm)
     Integer :: array, to, from
     Integer, Intent(In) :: comm
     Call mpi_send(array,1,MPI_INTEGER,to,from,comm,ierr)
  End Subroutine Send_Int_0D 

  Subroutine Send_Int_1D(array,to,from,comm)
     Integer :: array(:),to,from
     Integer, Intent(In) :: comm
     Integer :: nelem
     nelem = Size(array)
     Call mpi_send(array(1),nelem,MPI_INTEGER,to,from,comm,ierr)
  End Subroutine Send_Int_1D

  Subroutine Send_0D(array,to,from,comm)
     Real*8 :: array
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm

     Call mpi_send(array, 1, MPI_REAL8, to, from, comm,ierr)

  End Subroutine Send_0D

  Subroutine Send_1D(array,to,from,comm)
     Real*8 :: array(:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_send(array(1), nelem, MPI_REAL8, to, from, comm,ierr)
  End Subroutine Send_1D

  Subroutine Send_2D(array,to,from,comm)
     Real*8 :: array(:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_send(array(1,1), nelem, MPI_REAL8, to, from, comm,ierr)
  End Subroutine Send_2D

  Subroutine Send_3D(array,to,from,comm)
     Real*8 :: array(:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_send(array(1,1,1), nelem, MPI_REAL8, to, from, comm,ierr)
  End Subroutine Send_3D

  Subroutine Send_4D(array,to,from,comm)
     Real*8 :: array(:,:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_send(array(1,1,1,1), nelem, MPI_REAL8, to, from, comm,ierr)
  End Subroutine Send_4D

  Subroutine ISend_3D(ain,to,from,comm,req)
     Real*8 :: ain(:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer, Intent(Out) :: req
     Integer :: nelem
     
     nelem = Size(ain)

     Call mpi_isend(ain(1,1,1), nelem, MPI_REAL8, to, from, comm, req, ierr)
  End Subroutine ISend_3D

  Subroutine ISend_4D(ain,to,from,comm,req)
     Real*8 :: ain(:,:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer, Intent(Out) :: req
     Integer :: nelem
     
     nelem = Size(ain)

     Call mpi_isend(ain(1,1,1,1), nelem, MPI_REAL8, to, from, comm, req, ierr)
  End Subroutine ISend_4D

  Subroutine SendRecv_3D(ain,aout,to,from,comm)
     Real*8 :: ain(:,:,:)
     Real*8 :: aout(:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: nelem
     
     nelem = Size(ain)

     Call mpi_sendrecv(ain(1,1,1), nelem, MPI_REAL8, to, from, aout(1,1,1), nelem, MPI_REAL8, from, to, comm, ierr)
  End Subroutine SendRecv_3D

  Subroutine SendRecv_4D(ain,aout,to,from,comm)
     Real*8 :: ain(:,:,:,:)
     Real*8 :: aout(:,:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: nelem
     
     nelem = Size(ain)

     Call mpi_sendrecv(ain(1,1,1,1), nelem, MPI_REAL8, to, from, aout(1,1,1,1), nelem, MPI_REAL8, from, to, comm, ierr)
  End Subroutine SendRecv_4D

  Subroutine IReceive_3D(aout,from,to,comm,req)
     Real*8 :: aout(:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer, Intent(Out) :: req
     Integer :: nelem
     
     nelem = Size(aout)

     Call mpi_irecv(aout(1,1,1), nelem, MPI_REAL8, from, to, comm, req, ierr)
  End Subroutine IReceive_3D

  Subroutine IReceive_4D(aout,from,to,comm,req)
     Real*8 :: aout(:,:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer, Intent(Out) :: req
     Integer :: nelem
     
     nelem = Size(aout)

     Call mpi_irecv(aout(1,1,1,1), nelem, MPI_REAL8, from, to, comm, req, ierr)
  End Subroutine IReceive_4D

  Subroutine Recv_Bool_1D(array,from,to,comm)
    Logical :: array(:)
    Integer, Intent(In) :: to, from
    Integer, Intent(In) :: comm
    Integer :: status_mpi(MPI_STATUS_SIZE)
    Integer :: nelem
    nelem = Size(array)
    Call mpi_recv(array(1),nelem,MPI_LOGICAL,from,to,comm, status_mpi, ierr)
  End Subroutine Recv_Bool_1D

  Subroutine Receive_Int_0D(array,from,to,comm)
     Integer :: array
     Integer, Intent(In) :: from, to
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)

     Call mpi_recv(array,1,MPI_INTEGER,from,to,comm, status_mpi,ierr)
  End Subroutine Receive_Int_0D 

  Subroutine Receive_Int_1D(array,from,to,comm)
     Integer :: array(:)
     Integer, Intent(In) :: from, to
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)
     Integer :: nelem
     nelem = Size(array)
     Call mpi_recv(array(1),nelem,MPI_INTEGER,from,to,comm, status_mpi, ierr)
  End Subroutine Receive_Int_1D

  Subroutine Receive_0D(array,from,to,comm)
     Real*8 :: array
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)

     Call mpi_recv(array, 1, MPI_REAL8, from, to, comm, status_mpi,ierr)
  End Subroutine Receive_0D

  Subroutine Receive_1D(array,from,to,comm)
     Real*8 :: array(:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_recv(array(1), nelem, MPI_REAL8, from, to, comm, status_mpi,ierr)
  End Subroutine Receive_1D

  Subroutine Receive_2D(array,from,to,comm)
     Real*8 :: array(:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_recv(array(1,1), nelem, MPI_REAL8, from, to, comm, status_mpi,ierr)
  End Subroutine Receive_2D

  Subroutine Receive_3D(array,from,to,comm)
    Real*8 :: array(:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_recv(array(1,1,1), nelem, MPI_REAL8, from, to, comm, status_mpi,ierr)
  End Subroutine Receive_3D

  Subroutine Receive_4D(array,from,to,comm)
     Real*8 :: array(:,:,:,:)
     Integer, Intent(In) :: to, from
     Integer, Intent(In) :: comm
     Integer :: status_mpi(MPI_STATUS_SIZE)
     Integer :: nelem
     
     nelem = Size(array)

     Call mpi_recv(array(1,1,1,1), nelem, MPI_REAL8, from, to, comm, status_mpi,ierr)
  End Subroutine Receive_4D

  Subroutine Barrier(grp)
    Integer :: grp
    Call MPI_BARRIER(grp,ierr)

  End Subroutine Barrier

  Subroutine MWait(req)
    Integer, Intent(In) :: req
    Integer :: status(MPI_STATUS_SIZE),ierr
    
    Call MPI_WAIT(req,status,ierr)
  End Subroutine MWait


  Subroutine MWaitAll(req)
    Integer, Intent(In) :: req(:)
    Integer, Allocatable, Dimension(:,:) :: status
    Integer :: ierr, n
    n = size(req)
    Allocate(status(n,MPI_STATUS_SIZE))
    Call MPI_WAITALL(n,req,status,ierr)
  End Subroutine MWaitAll

End Module Parallel
