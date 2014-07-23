!=============================================================================================!
!                                                                                             !
!   Module Derivatives contains routines that compute the sparse derivative matrices and      !
!   handle the matrix multiplication to evaluate the derivatives.  The derivatives are sixth  !
!   order finite differences everywhere in the interior. Boundary information is gathered     !
!   from "nearest-neighbor" cpus.                                                             !
!   Kyle Augustson JILA 2009                                                                  !
!                                                                                             !
!=============================================================================================!

!===================================================================================!
!                                                                                   !
!      IBC=0 uses neumann conditions on boundaries (Dydx=0)                         !
!      IBC=1 calculates value on boundaries                                         !
!      IBC=2 neumann on left end, Dirichelet on right                               ! 
!      IBC=3 vica versa                                                             !
!      IBC=4 uses neumann conditions on boundaries (Dydx<>0)                        !
!                                                                                   !
!    Note that the neuman BC must be set if either IBC = 2, 3 or 4                  !
!                                                                                   !
!===================================================================================!

Module Derivatives
  Use Parallel

  Implicit None

  Integer :: dtype1,dtype2,dtype3
  Integer, Allocatable, Dimension(:) :: ibcr,ibct,ibcp
  Real*8 :: dr,dth,dphi,dri,dti,dpi,dri2,dti2,dpi2
  Real*8 :: r1,r2,th1,th2,phi1,phi2,x_1,x_2
  Real*8, Allocatable, Dimension(:) :: sines,cosines,r_inv,dxdr,d2xdr2,r_ind,dr_ind
  Real*8, Allocatable, Dimension(:) :: x_ind,theta,phi
  Real*8, Allocatable, Dimension(:) :: gam1p, gam2p, gam3p,gam1p2, gam2p2, gam3p2
  Real*8, Allocatable, Dimension(:,:) :: gam1r, gam2r, gam3r, gam1t, gam2t, gam3t, gam1r2, gam2r2, gam3r2, gam1t2, gam2t2, gam3t2
  !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: sines, cosines, r_inv, dxdr, d2xdr2, r_ind, dr_ind, x_ind, theta, phi
Contains

  Subroutine Initialize_Derivatives()
    Implicit None
    Real*8 :: c1, c2, c3, alpha, a1, a2, g2, up1, up2, cn, b1
    Integer :: luptmp(2), ir, it, ip, iv
    If (magnetic) Then
       Allocate(ibcr(nv+4),ibct(nv+4),ibcp(nv+4))
    Else
       Allocate(ibcr(nv+1),ibct(nv+1),ibcp(nv+1))
    EndIf

    !Determine boundary types in r
    ibcr(1:6) = (/lbcr_rho,lbcr_w,lbcr_v,lbcr_u,lbcr_s,lbcr_t/)
    If (magnetic) ibcr(7:12) = (/lbcr_bt,lbcr_bp,lbcr_br,lbcr_at,lbcr_ap,lbcr_ar/)
    luptmp = node_lup(1,:)
    Call Determine_Dtype(luptmp,dtype1)

    !Determine boundary types in theta
    ibct(1:6) = (/lbct_rho,lbct_w,lbct_v,lbct_u,lbct_s,lbct_t/)
    If (magnetic) ibct(7:12) = (/lbct_bt,lbct_bp,lbct_br,lbct_at,lbct_ap,lbct_ar/)
    luptmp = node_lup(2,:)
    Call Determine_Dtype(luptmp,dtype2)

    If ((Theta_Symmetric).and.(myth2.eq.nth)) Then
       dtype2 = 3
    EndIf

    ibcp = -1 !All periodic in phi
    dtype3 = 3

    Allocate(x_ind(nr),r_ind(nr),r_inv(nr),theta(nth),phi(nphi),dr_ind(nr))

    If (non_uniform) Then
       dr = (x_2-x_1)/Dble(nr-1)
       x_ind = dr*(/(Dble(ir),ir=0,nr-1)/)+x_1
    Else      
       dr = (r2-r1)/Dble(nr-1)
       dr_ind = dr
       x_ind = dr*(/(Dble(ir),ir=0,nr-1)/)+r1
       r_ind = x_ind
       r_inv = 1d0/r_ind
    EndIf

    dth = (th2-th1)/Dble(nth-1)
    dphi = (phi2-phi1)/Dble(nphi)
    theta = dth*(/(Dble(it),it=0,nth-1)/)+th1
    phi = dphi*(/(Dble(ip),ip=0,nphi-1)/)+phi1

    !  replaces call to subroutine angles() in f77 version
    Allocate(sines(nth),cosines(nth))
    sines=Sin(theta)
    cosines=Cos(theta)

    dri = 1d0/(x_ind(myr2)-x_ind(myr1))
    dti = 1d0/(theta(myth2)-theta(myth1))
    dpi = 1d0/(phi(myphi2)-phi(myphi1))

    dri = Dble(mynr-1)*dri
    dti = Dble(mynth-1)*dti
    dpi = Dble(mynphi-1)*dpi

    dri2 = dri*dri
    dti2 = dti*dti
    dpi2 = dpi*dpi

    Allocate(dxdr(nr),d2xdr2(nr))

    !Allocate derivative coefs (radial)
    Allocate(gam1r(mynr,5),gam2r(mynr,5),gam3r(mynr,5))

    Do iv=0,4
       alpha = 1d0/3d0
       a2 = 0.75d0 
       g2 = 0.125d0 
       a1 = 6d0
       b1 = 3d0
       cn = 1d0/(a2-b1*g2) 
       up1 = (a1*a2-b1)*cn 
       If ((iv.ne.1).and.(iv.ne.2)) Then
          up2 = 0d0
       Else
          up2 = up1
       End If

       If (mod(iv,2).EQ.0) Then
          up1 = 0d0
       End If

       c1 = 1d0
       c2 = 1d0
       gam1r(1,iv+1) = c2
       gam2r(1,iv+1) = c2
       gam3r(1,iv+1) = c2
       If (myr1 .eq. 1) Then
          c3 = up1*c2
          gam1r(2,iv+1) = c3
          c2 = c1/(c1-g2*c3)
          gam2r(2,iv+1) = c2

          c3 = a2*c2
          gam1r(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r(3,iv+1) = c2

          gam3r(2,iv+1) = -g2*gam2r(2,iv+1)
       Else
          c3 = alpha*c2
          gam1r(2,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r(2,iv+1) = c2
          gam3r(2,iv+1) = -alpha*gam2r(2,iv+1)

          c3 = alpha*c2
          gam1r(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r(3,iv+1) = c2
       EndIf

       Do ir=4,mynr-2
          c3 = alpha*c2
          gam1r(ir,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r(ir,iv+1) = c2
       EndDo

       gam3r(3:mynr-2,iv+1) = -alpha*gam2r(3:mynr-2,iv+1)

       If (myr2 .eq. nr) Then
          c3 = alpha*c2
          gam1r(mynr-1,iv+1) = c3
          c2 = c1/(c1-a2*c3)
          gam2r(mynr-1,iv+1) = c2
          c3 = g2*c2
          gam1r(mynr,iv+1) = c3
          c2 = c1/(c1-up2*c3)
          gam2r(mynr,iv+1) = c2

          gam3r(mynr-1,iv+1) = -a2*gam2r(mynr-1,iv+1)
          gam3r(mynr,iv+1) = -up2*gam2r(mynr,iv+1)
       Else
          c3 = alpha*c2
          gam1r(mynr-1,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r(mynr-1,iv+1) = c2
          c3 = alpha*c2
          gam1r(mynr,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r(mynr,iv+1) = c2

          gam3r(mynr-1,iv+1) = -alpha*gam2r(mynr-1,iv+1)
          gam3r(mynr,iv+1) = -alpha*gam2r(mynr,iv+1)
       EndIf
    EndDo
    gam1r = -gam1r

    !Allocate 2nd derivative coefs (radial)
    Allocate(gam1r2(mynr,5),gam2r2(mynr,5),gam3r2(mynr,5))
    Do iv=0,4
       alpha = 2d0/11d0
       Select Case(iv)
       Case(0)
          up1 = 707d0/193d0
          up2 = 707d0/193d0
       Case(1)
          up1 = 137d0/13d0
          up2 = 137d0/13d0
       Case(2)
          up1 = 707d0/193d0
          up2 = 137d0/13d0
       Case(3)
          up1 = 137d0/13d0
          up2 = 707d0/193d0
       Case(4)
          up1 = 707d0/193d0
          up2 = 707d0/193d0
       EndSelect

       a2 = -0.35d0
       g2 = 0.1d0

       c1 = 1d0
       c2 = 1d0
       gam1r2(1,iv+1) = c2
       gam2r2(1,iv+1) = c2
       gam3r2(1,iv+1) = c2
       If (myr1 .eq. 1) Then
          c3 = up1
          gam1r2(2,iv+1) = c3
          c2 = c1/(c1-g2*c3)
          gam2r2(2,iv+1) = c2

          c3 = a2*c2
          gam1r2(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r2(3,iv+1) = c2

          gam3r2(2,iv+1) = -g2*gam2r2(2,iv+1)
       Else
          c3 = alpha*c2
          gam1r2(2,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r2(2,iv+1) = c2
          gam3r2(2,iv+1) = -alpha*gam2r2(2,iv+1)

          c3 = alpha*c2
          gam1r2(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r2(3,iv+1) = c2
       EndIf

       Do ir=4,mynr-2
          c3 = alpha*c2
          gam1r2(ir,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r2(ir,iv+1) = c2
       EndDo

       gam3r2(3:mynr-2,iv+1) = -alpha*gam2r2(3:mynr-2,iv+1)

       g2 = a2
       a2 = 0.1d0

       If (myr2 .eq. nr) Then
          c3 = alpha*c2
          gam1r2(mynr-1,iv+1) = c3
          c2 = c1/(c1-g2*c3)
          gam2r2(mynr-1,iv+1) = c2

          c3 = a2*c2
          gam1r2(mynr,iv+1) = c3
          c2 = c1/(c1-up2*c3)
          gam2r2(mynr,iv+1) = c2

          gam3r2(mynr-1,iv+1) = -g2*gam2r2(mynr-1,iv+1)
          gam3r2(mynr,iv+1) = -up2*gam2r2(mynr,iv+1)
       Else
          c3 = alpha*c2
          gam1r2(mynr-1,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r2(mynr-1,iv+1) = c2
          c3 = alpha*c2
          gam1r2(mynr,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2r2(mynr,iv+1) = c2

          gam3r2(mynr-1,iv+1) = -alpha*gam2r2(mynr-1,iv+1)
          gam3r2(mynr,iv+1) = -alpha*gam2r2(mynr,iv+1)
       EndIf
       gam1r2(:,iv+1) = -gam1r2(:,iv+1)
    EndDo

    !Allocate derivative coefs (theta)
    Allocate(gam1t(mynth,5),gam2t(mynth,5),gam3t(mynth,5))
    Do iv=0,4    
       alpha = 1d0/3d0
       a2 = 0.75d0 
       g2 = 0.125d0 
       a1 = 6d0
       b1 = 3d0
       cn = 1d0/(a2-b1*g2) 
       up1 = (a1*a2-b1)*cn 
       If ((iv.ne.1).and.(iv.ne.2)) Then
          up2 = 0d0
       Else
          up2 = up1
       End If

       If (mod(iv,2).EQ.0) Then
          up1 = 0d0
       End If

       c1 = 1d0
       c2 = 1d0
       gam1t(1,iv+1) = c1
       gam2t(1,iv+1) = c1
       gam3t(1,iv+1) = c1
       If (myth1 .eq. 1) Then
          c3 = up1*c2
          gam1t(2,iv+1) = c3
          c2 = c1/(c1-g2*c3)
          gam2t(2,iv+1) = c2

          c3 = a2*c2
          gam1t(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t(3,iv+1) = c2

          gam3t(2,iv+1) = -g2*gam2t(2,iv+1)
       Else
          c3 = alpha*c2
          gam1t(2,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t(2,iv+1) = c2
          gam3t(2,iv+1) = -alpha*gam2t(2,iv+1)

          c3 = alpha*c2
          gam1t(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t(3,iv+1) = c2
       EndIf

       Do ir=4,mynth-2
          c3 = alpha*c2
          gam1t(ir,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t(ir,iv+1) = c2
       EndDo

       If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
          c3 = alpha*c2
          gam1t(mynth-1,iv+1) = c3
          c2 = c1/(c1-a2*c3)
          gam2t(mynth-1,iv+1) = c2
          c3 = g2*c2
          gam1t(mynth,iv+1) = c3
          c2 = c1/(c1-up2*c3)
          gam2t(mynth,iv+1) = c2

          gam3t(mynth-1,iv+1) = -a2*gam2t(mynth-1,iv+1)
          gam3t(mynth,iv+1) = -up2*gam2t(mynth,iv+1)
       Else
          c3 = alpha*c2
          gam1t(mynth-1,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t(mynth-1,iv+1) = c2
          c3 = alpha*c2
          gam1t(mynth,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t(mynth,iv+1) = c2

          gam3t(mynth-1,iv+1) = -alpha*gam2t(mynth-1,iv+1)
          gam3t(mynth,iv+1) = -alpha*gam2t(mynth,iv+1)
       EndIf
       gam3t(3:mynth-2,iv+1) = -alpha*gam2t(3:mynth-2,iv+1) 
    EndDo
    gam1t = -gam1t

    !Allocate 2nd derivative coefs (theta)
    Allocate(gam1t2(mynth,5),gam2t2(mynth,5),gam3t2(mynth,5))
    Do iv=0,4
       alpha = 2d0/11d0
       Select Case(iv)
       Case(0)
          up1 = 707d0/193d0
          up2 = 707d0/193d0
       Case(1)
          up1 = 137d0/13d0
          up2 = 137d0/13d0
       Case(2)
          up1 = 707d0/193d0
          up2 = 137d0/13d0
       Case(3)
          up1 = 137d0/13d0
          up2 = 707d0/193d0
       Case(4)
          up1 = 707d0/193d0
          up2 = 707d0/193d0
       EndSelect

       a2 = -0.35d0
       g2 = 0.1d0

       c1 = 1d0
       c2 = 1d0
       gam1t2(1,iv+1) = c2
       gam2t2(1,iv+1) = c2
       gam3t2(1,iv+1) = c2
       If (myth1 .eq. 1) Then
          c3 = up1
          gam1t2(2,iv+1) = c3
          c2 = c1/(c1-g2*c3)
          gam2t2(2,iv+1) = c2

          c3 = a2*c2
          gam1t2(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t2(3,iv+1) = c2

          gam3t2(2,iv+1) = -g2*gam2t2(2,iv+1)
       Else
          c3 = alpha*c2
          gam1t2(2,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t2(2,iv+1) = c2
          gam3t2(2,iv+1) = -alpha*gam2t2(2,iv+1)

          c3 = alpha*c2
          gam1t2(3,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t2(3,iv+1) = c2
       EndIf

       Do ir=4,mynth-2
          c3 = alpha*c2
          gam1t2(ir,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t2(ir,iv+1) = c2
       EndDo

       gam3t2(3:mynth-2,iv+1) = -alpha*gam2t2(3:mynth-2,iv+1)

       g2 = a2
       a2 = 0.1d0

       If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
          c3 = alpha*c2
          gam1t2(mynth-1,iv+1) = c3
          c2 = c1/(c1-g2*c3)
          gam2t2(mynth-1,iv+1) = c2

          c3 = a2*c2
          gam1t2(mynth,iv+1) = c3
          c2 = c1/(c1-up2*c3)
          gam2t2(mynth,iv+1) = c2

          gam3t2(mynth-1,iv+1) = -g2*gam2t2(mynth-1,iv+1)
          gam3t2(mynth,iv+1) = -up2*gam2t2(mynth,iv+1)
       Else
          c3 = alpha*c2
          gam1t2(mynth-1,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t2(mynth-1,iv+1) = c2
          c3 = alpha*c2
          gam1t2(mynth,iv+1) = c3
          c2 = c1/(c1-alpha*c3)
          gam2t2(mynth,iv+1) = c2

          gam3t2(mynth-1,iv+1) = -alpha*gam2t2(mynth-1,iv+1)
          gam3t2(mynth,iv+1) = -alpha*gam2t2(mynth,iv+1)
       EndIf
       gam1t2(:,iv+1) = -gam1t2(:,iv+1)
    EndDo

    !Allocate derivative coefs (phi)
    Allocate(gam1p(mynphi),gam2p(mynphi),gam3p(mynphi))
    alpha = 1d0/3d0
    c1 = 1d0
    c2 = 1d0
    gam1p(1) = c1
    gam2p(1) = c1
    gam3p(1) = c1

    c3 = alpha*c2
    gam1p(2) = c3
    c2 = c1/(c1-alpha*c3)
    gam2p(2) = c2
    gam3p(2) = -alpha*gam2p(2)

    Do ir=3,mynphi
       c3 = alpha*c2
       gam1p(ir) = c3
       c2 = c1/(c1-alpha*c3)
       gam2p(ir) = c2
    EndDo
    gam3p(3:mynphi) = -alpha*gam2p(3:mynphi) 
    gam1p = -gam1p


    !Allocate 2nd derivative coefs (phi)
    Allocate(gam1p2(mynphi),gam2p2(mynphi),gam3p2(mynphi))
    alpha = 2d0/11d0
    c1 = 1d0
    c2 = 1d0
    gam1p2(1) = c1
    gam2p2(1) = c1
    gam3p2(1) = c1

    c3 = alpha*c2
    gam1p2(2) = c3
    c2 = c1/(c1-alpha*c3)
    gam2p2(2) = c2
    gam3p2(2) = -alpha*gam2p2(2)

    Do ir=3,mynphi
       c3 = alpha*c2
       gam1p2(ir) = c3
       c2 = c1/(c1-alpha*c3)
       gam2p2(ir) = c2
    EndDo
    gam3p2(3:mynphi) = -alpha*gam2p2(3:mynphi) 
    gam1p2 = -gam1p2

  End Subroutine Initialize_Derivatives

  Subroutine Finalize_Derivatives()
    Deallocate(r_ind,r_inv,theta,phi,sines,cosines,dr_ind,dxdr,d2xdr2)
    Deallocate(gam1r,gam2r,gam3r,gam1r2,gam2r2,gam3r2)
    Deallocate(gam1t,gam2t,gam3t,gam1t2,gam2t2,gam3t2)
    Deallocate(gam1p,gam2p,gam3p,gam1p2,gam2p2,gam3p2)
  End Subroutine Finalize_Derivatives

  Subroutine Determine_Dtype(lup,dtype)
    Implicit None
    Integer, Intent(Out) :: dtype
    Integer, Intent(InOut) :: lup(2)

    If ((lup(1) .ne. -1) .and. (lup(2) .ne. -1)) Then !internal boundary or periodic
       dtype = 3
    ElseIf ((lup(1) .eq. -1) .and. (lup(2) .ne. -1)) Then
       dtype = 2
    ElseIf ((lup(1) .ne. -1) .and. (lup(2) .eq. -1)) Then
       dtype = 1
    ElseIf ((lup(1) .eq. -1) .and. (lup(2) .eq. -1)) Then
       dtype = 0
    EndIf

  End Subroutine Determine_Dtype

  !#DEC$ ATTRIBUTES FORCEINLINE :: Dbydr
  Subroutine Dbydr(x,dxdy,b1,b2,var)
    Implicit None
    Real*8, Intent(InOut) :: x(:,:),b1(:,:),b2(:,:)
    Real*8, Intent(InOut) :: dxdy(:,:)
    Integer, Intent(In) :: var
    Real*8 :: dtmp(mynth,mynr),dtmp2(mynr)
    Real*8 :: c1,c2,c3,c4,c5,alpha2,gamma2,alpha1,beta1,const
    Integer :: j, bdry

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dtmp, dtmp2

    If (var .ge. 1) Then
       bdry = ibcr(var)
    Else
       bdry = -var
    EndIf

    Select Case (dtype1)
    Case (1)
       c1 = -1d0/180d0
       c2 = 9d0
       c3 = 45d0
       c4 = 7d0/9d0
       c5 = 1d0/36d0
       c2 = c1*c2
       c3 = c1*c3
       const = -(c3+c5)
       beta1 = c4-c2
       alpha1 = c1+c5
       c1 = -c1
       c4 = -c4
       !#DIR$ VECTOR ALIGNED
       dxdy(:,1) = c1*b1(:,1) + c2*b1(:,2) + const*b1(:,3) + c4*b1(:,4) + c3*x(:,1) + beta1*x(:,2) + alpha1*x(:,3)
       c4 = 7/9d0
       c5 = 1d0/36d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,2) = c4*(x(:,3)-x(:,1)) + c5*(x(:,4)-b1(:,4))

       c1 = 43d0/96d0
       c2 = 5d0/6d0
       c3 = -9d0/8d0
       c4 = -1d0/6d0
       c5 = 1d0/96d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynr-1) = c1*x(:,mynr) + c2*x(:,mynr-1) + c3*x(:,mynr-2) + c4*x(:,mynr-3) + c5*x(:,mynr-4)
       
       If ((bdry.eq. 3).or.(bdry.eq. 4)) Then
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynr) = b2(:,1)/dri
       Else
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynr)=0d0
       EndIf

       If ((bdry.eq.1).or.(bdry.eq.2))  Then
          ! here are the pentadiagonal and fifth order values for the boundary.
          ! precondition the matrix to make it tridiagonal
          alpha2 = 0.75d0
          gamma2 = 0.125d0
          alpha1 = 6d0
          beta1  = 3d0
          const = 1d0/(alpha2-beta1*gamma2)
          alpha2 = alpha2*const
          beta1 = -beta1*const
          ! fifth order coeff.
          c1 = 10d0/3d0
          c2 = 3d0
          c3 = -6d0
          c4 = -1d0/3d0
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynr) = c1*x(:,mynr) + c2*x(:,mynr-1) + c3*x(:,mynr-2) + c4*x(:,mynr-3)
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynr) = alpha2*dxdy(:,mynr) + beta1*dxdy(:,mynr-1)
       End If
       c4 = 7d0/9d0
       c5 = 1d0/36d0
    Case (2)
       If (mod(bdry,2).eq.0) Then
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1) = b1(:,1)/dri
       EndIf

       If (bdry.eq.0) Then
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1)=0d0
       EndIf
       
       c1 = -43d0/96d0
       c2 = -5d0/6d0
       c3 = 9d0/8d0
       c4 = 1d0/6d0
       c5 = -1d0/96d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,2) = c1*x(:,1) + c2*x(:,2) + c3*x(:,3) + c4*x(:,4) + c5*x(:,5)
       
       If ((bdry.eq.1).or.(bdry.eq.3)) Then
          ! here are the pentadiagonal and fifth order values for the boundary.
          ! precondition the matrix to make it tridiagonal
          alpha2 = 0.75d0
          gamma2 = 0.125d0
          alpha1 = 6d0
          beta1  = 3d0
          const  = 1d0/(alpha2-beta1*gamma2)
          alpha2 = alpha2*const
          beta1 = -beta1*const

          ! fifth order coeff.
          c1 = -10d0/3d0
          c2 = -3d0
          c3 = 6d0
          c4 = 1d0/3d0
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1) = c1*x(:,1) + c2*x(:,2) + c3*x(:,3) + c4*x(:,4)
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1) = alpha2*dxdy(:,1) + beta1*dxdy(:,2)
       End If
       
       c1 = 7d0/9d0
       c2 = 1d0/36d0
       c3 = -1d0/180d0
       c4 = 9d0
       c5 = 45d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynr-1) = c1*(x(:,mynr) - x(:,mynr-2)) + c2*(b2(:,1) - x(:,mynr-3))

       const = c2+c3*c5
       beta1 = c3*c4-c1
       alpha1 = -(c2+c3)
       c4 = -c3*c4
       c5 = -c3*c5
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynr) = c1*b2(:,1) + const*b2(:,2) + c4*b2(:,3) + c3*b2(:,4) + beta1*x(:,mynr-1) + alpha1*x(:,mynr-2) + c5*x(:,mynr)

       c4 = 7/9d0
       c5 = 1d0/36d0
    Case (3)

       c1 = -1d0/180d0
       c2 = 9d0
       c3 = 45d0
       c4 = 7d0/9d0
       c5 = 1d0/36d0
       c2 = c1*c2
       c3 = c1*c3
       const = -(c3+c5)
       beta1 = c4-c2
       alpha1 = c1+c5
       c1 = -c1
       c4 = -c4
       !#DIR$ VECTOR ALIGNED
       dxdy(:,1) = c1*b1(:,1) + c2*b1(:,2) + const*b1(:,3) + c4*b1(:,4) + c3*x(:,1) + beta1*x(:,2) + alpha1*x(:,3)
       c4 = 7/9d0
       c5 = 1d0/36d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,2) = c4*(x(:,3)-x(:,1)) + c5*(x(:,4)-b1(:,4))
 
       c1 = 7d0/9d0
       c2 = 1d0/36d0
       c3 = -1d0/180d0
       c4 = 9d0
       c5 = 45d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynr-1) = c1*(x(:,mynr) - x(:,mynr-2)) + c2*(b2(:,1) - x(:,mynr-3))

       const = c2+c3*c5
       beta1 = c3*c4-c1
       alpha1 = -(c2+c3)
       c4 = -c3*c4
       c5 = -c3*c5
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynr) = c1*b2(:,1) + const*b2(:,2) + c4*b2(:,3) + c3*b2(:,4) + beta1*x(:,mynr-1) + alpha1*x(:,mynr-2) + c5*x(:,mynr)

       c4 = 7/9d0
       c5 = 1d0/36d0
    EndSelect

    !#DEC$ LOOP COUNT MAX=252, MIN=8, AVG=12
    !#DIR$ VECTOR ALIGNED
    Do j=3,mynr-2
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j) = c4*(x(:,j+1)-x(:,j-1)) + c5*(x(:,j+2)-x(:,j-2))
    EndDo

    bdry = bdry+1
    alpha1 = -1d0/3d0
    c1 = gam3r(2,bdry)
    c2 = gam3r(mynr-1,bdry)
    c3 = gam3r(mynr,bdry)
    dtmp2 = gam2r(:,bdry)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,1) = dtmp2(1)*dxdy(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,2) = dtmp2(2)*dxdy(:,2)+c1*dtmp(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,3) = dtmp2(3)*(dxdy(:,3)+alpha1*dtmp(:,2))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,4) = dtmp2(4)*(dxdy(:,4)+alpha1*dtmp(:,3))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,5) = dtmp2(5)*(dxdy(:,5)+alpha1*dtmp(:,4))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,6) = dtmp2(6)*(dxdy(:,6)+alpha1*dtmp(:,5))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,7) = dtmp2(7)*(dxdy(:,7)+alpha1*dtmp(:,6))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,8) = dtmp2(8)*(dxdy(:,8)+alpha1*dtmp(:,7))

    !#DEC$ LOOP COUNT MAX=123, MIN=1, AVG=3
    !#DIR$ VECTOR ALIGNED
    Do j=9,mynr-3,2
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j) = dtmp2(j)*(dxdy(:,j)+alpha1*dtmp(:,j-1))
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j+1) = dtmp2(j+1)*(dxdy(:,j+1)+alpha1*dtmp(:,j))
    EndDo
    j = mynr-1
    !#DIR$ VECTOR ALIGNED
    dtmp(:,j) = dtmp2(j)*dxdy(:,j)+c2*dtmp(:,j-1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,mynr) = dtmp2(mynr)*dxdy(:,mynr)+c3*dtmp(:,j)

    !Backward substitution
    dtmp2 = gam1r(:,bdry)
    dxdy(:,mynr) = dtmp(:,mynr)
    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=4
    !#DIR$ VECTOR ALIGNED
    Do j=mynr-1,9,-2
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j) = dtmp(:,j)+dtmp2(j+1)*dxdy(:,j+1)
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j-1) = dtmp(:,j-1)+dtmp2(j)*dxdy(:,j)
    EndDo

    !#DIR$ VECTOR ALIGNED
    dxdy(:,7) = dtmp(:,7)+dtmp2(8)*dxdy(:,8)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,6) = dtmp(:,6)+dtmp2(7)*dxdy(:,7)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,5) = dtmp(:,5)+dtmp2(6)*dxdy(:,6)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,4) = dtmp(:,4)+dtmp2(5)*dxdy(:,5)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,3) = dtmp(:,3)+dtmp2(4)*dxdy(:,4)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,2) = dtmp(:,2)+dtmp2(3)*dxdy(:,3)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,1) = dtmp(:,1)+dtmp2(2)*dxdy(:,2)

    !#DIR$ VECTOR ALIGNED
    dxdy = dri*dxdy
  End Subroutine Dbydr

  !#DEC$ ATTRIBUTES FORCEINLINE :: D2bydr2
  Subroutine D2bydr2(x,d2xdy2,b1,b2,var)
    Implicit None
    Real*8, Intent(InOut) :: x(:,:),b1(:,:),b2(:,:)
    Real*8, Intent(InOut) :: d2xdy2(:,:)
    Integer, Intent(In) :: var
    Real*8 :: dtmp(mynth,mynr),dtmp2(mynr)
    Real*8 :: c1,c2,c3,c4,c5,c6,c7,alpha1
    Integer :: j, bdry

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dtmp, dtmp2

    If (var .ge. 1) Then
       bdry = ibcr(var)
    Else
       bdry = -var
    EndIf

    Select Case (dtype1)
    Case (1)
       c1 = 12d0/11d0
       c2 = 3d0/44d0
       c3 = -1d0/990d0  
       c4 = 2d0
       c5 = -27d0
       c6 = 270d0
       c7 = -490d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = c3*(c4*(b1(:,1)+x(:,3))+c5*(b1(:,2)+x(:,2))+c6*(b1(:,3)+x(:,1))+c7*b1(:,4))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = d2xdy2(:,1)+c1*(x(:,2)-2d0*x(:,1)+b1(:,4))+c2*(b1(:,3)-2d0*x(:,1)+x(:,3))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,2) = c1*(x(:,3)-2d0*x(:,2)+x(:,1))+c2*(b1(:,4)-2d0*x(:,2)+x(:,4))

       c1=1.2375d0
       c2=-3d0
       c3=2.325d0
       c4=-0.6d0
       c5=0.0375d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr-1)=c1*x(:,mynr)+c2*x(:,mynr-1)+c3*x(:,mynr-2)+c4*x(:,mynr-3)+c5*x(:,mynr-4)

       If ((bdry.eq.1).or.(bdry.eq.2)) Then
          !Fifth order boundary
          c1 = 1955d0/156d0
          c2 = -4057d0/156d0
          c3 = 1117d0/78d0
          c4 = -55d0/78d0
          c5 = -29d0/156d0
          c6 = 7d0/156d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,mynr)=c1*x(:,mynr)+c2*x(:,mynr-1)+c3*x(:,mynr-2)+c4*x(:,mynr-3)+c5*x(:,mynr-4)+c6*x(:,mynr-5)
       Else
          !Fifth Order, with 1st order derivative specified
          c1 = 575d0/193d0/dri !Note dri here
          c2 = -5827d0/2316d0
          c3 = 1987d0/1158d0
          c4 = 465d0/386d0
          c5 = -547d0/1158d0
          c6 = 157d0/2316d0
          
          If (bdry .eq. 0) c1 = 0d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,mynr)=c1*b2(:,1)+c2*x(:,mynr-1)+c3*x(:,mynr-2)+c4*x(:,mynr-3)+c5*x(:,mynr-4)+c6*x(:,mynr-5)
       EndIf
       c1 = 12d0/11d0
       c2 = 3d0/44d0
    Case (2)
       If (mod(bdry,2).eq.0) Then
          !Fifth Order, with 1st order derivative specified
          c1 = -575d0/193d0/dri !Note dri here
          c2 = -5827d0/2316d0
          c3 = 1987d0/1158d0
          c4 = 465d0/386d0
          c5 = -547d0/1158d0
          c6 = 157d0/2316d0

          If (bdry .eq. 0) c1 = 0d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,1)=c1*b1(:,1)+c2*x(:,2)+c3*x(:,3)+c4*x(:,4)+c5*x(:,5)+c6*x(:,6)
       Else
          !Fifth order boundary
          c1 = 1955d0/156d0
          c2 = -4057d0/156d0
          c3 = 1117d0/78d0
          c4 = -55d0/78d0
          c5 = -29d0/156d0
          c6 = 7d0/156d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,1)=c1*x(:,1)+c2*x(:,2)+c3*x(:,3)+c4*x(:,4)+c5*x(:,5)+c6*x(:,6)
       EndIf

       c1=1.2375d0
       c2=-3d0
       c3=2.325d0
       c4=-0.6d0
       c5=0.0375d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,2)=c1*x(:,1)+c2*x(:,2)+c3*x(:,3)+c4*x(:,4)+c5*x(:,5)

       c1 = 12d0/11d0
       c2 = 3d0/44d0
       c3 = -1d0/990d0
       c4 = 2d0
       c5 = -27d0
       c6 = 270d0
       c7 = -490d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr-1) = c1*(x(:,mynr)-2d0*x(:,mynr-1)+x(:,mynr-2))+c2*(x(:,mynr-3)-2d0*x(:,mynr-1)+b2(:,1))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr) = c1*(b2(:,1)-2d0*x(:,mynr)+x(:,mynr-1))+c2*(x(:,mynr-2)-2d0*x(:,mynr)+b2(:,2))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr) = d2xdy2(:,mynr)+c3*(c4*(x(:,mynr-2)+b2(:,4))+c5*(x(:,mynr-1)+b2(:,3))+c6*(x(:,mynr)+b2(:,2))+c7*b2(:,1))
    Case (3) 
       c1 = 12d0/11d0
       c2 = 3d0/44d0
       c3 = -1d0/990d0
       c4 = 2d0
       c5 = -27d0
       c6 = 270d0
       c7 = -490d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = c3*(c4*(b1(:,1)+x(:,3))+c5*(b1(:,2)+x(:,2))+c6*(b1(:,3)+x(:,1))+c7*b1(:,4))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = d2xdy2(:,1)+c1*(x(:,2)-2d0*x(:,1)+b1(:,4))+c2*(b1(:,3)-2d0*x(:,1)+x(:,3))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,2) = c1*(x(:,3)-2d0*x(:,2)+x(:,1))+c2*(b1(:,4)-2d0*x(:,2)+x(:,4))
   
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr-1) = c1*(x(:,mynr)-2d0*x(:,mynr-1)+x(:,mynr-2))+c2*(x(:,mynr-3)-2d0*x(:,mynr-1)+b2(:,1))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr) = c1*(b2(:,1)-2d0*x(:,mynr)+x(:,mynr-1))+c2*(x(:,mynr-2)-2d0*x(:,mynr)+b2(:,2))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynr) = d2xdy2(:,mynr)+c3*(c4*(x(:,mynr-2)+b2(:,4))+c5*(x(:,mynr-1)+b2(:,3))+c6*(x(:,mynr)+b2(:,2))+c7*b2(:,1))
    EndSelect

    !#DEC$ LOOP COUNT MAX=252, MIN=8, AVG=12
    !#DIR$ VECTOR ALIGNED
    Do j=3,mynr-2
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j) = c1*(x(:,j+1)-2d0*x(:,j)+x(:,j-1))+c2*(x(:,j-2)-2d0*x(:,j)+x(:,j+2))
    End Do

    bdry = bdry+1
    alpha1 = -2d0/11d0
    c1 = gam3r2(2,bdry)
    c2 = gam3r2(mynr-1,bdry)
    c3 = gam3r2(mynr,bdry)
    dtmp2 = gam2r2(:,bdry)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,1) = dtmp2(1)*d2xdy2(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,2) = dtmp2(2)*d2xdy2(:,2)+c1*dtmp(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,3) = dtmp2(3)*(d2xdy2(:,3)+alpha1*dtmp(:,2))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,4) = dtmp2(4)*(d2xdy2(:,4)+alpha1*dtmp(:,3))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,5) = dtmp2(5)*(d2xdy2(:,5)+alpha1*dtmp(:,4))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,6) = dtmp2(6)*(d2xdy2(:,6)+alpha1*dtmp(:,5))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,7) = dtmp2(7)*(d2xdy2(:,7)+alpha1*dtmp(:,6))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,8) = dtmp2(8)*(d2xdy2(:,8)+alpha1*dtmp(:,7))

    !#DEC$ LOOP COUNT MAX=123, MIN=1, AVG=3
    !#DIR$ VECTOR ALIGNED
    Do j=9,mynr-3,2
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j) = dtmp2(j)*(d2xdy2(:,j)+alpha1*dtmp(:,j-1))
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j+1) = dtmp2(j+1)*(d2xdy2(:,j+1)+alpha1*dtmp(:,j))
    EndDo !2*(mynr-16) + 6*mynr+58 = 8*mynr+26
    j = mynr-1
    !#DIR$ VECTOR ALIGNED
    dtmp(:,j) = dtmp2(j)*d2xdy2(:,j)+c2*dtmp(:,j-1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,mynr) = dtmp2(mynr)*d2xdy2(:,mynr)+c3*dtmp(:,j)

    !Backward substitution
    dtmp2 = gam1r2(:,bdry)
    d2xdy2(:,mynr) = dtmp(:,mynr)
    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=4
    !#DIR$ VECTOR ALIGNED
    Do j=mynr-1,9,-2
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j) = dtmp(:,j)+dtmp2(j+1)*d2xdy2(:,j+1)
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j-1) = dtmp(:,j-1)+dtmp2(j)*d2xdy2(:,j)
    EndDo !8*mynr+26 + 2*(mynr-16) = 10*mynr-6

    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,7) = dtmp(:,7)+dtmp2(8)*d2xdy2(:,8)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,6) = dtmp(:,6)+dtmp2(7)*d2xdy2(:,7)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,5) = dtmp(:,5)+dtmp2(6)*d2xdy2(:,6)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,4) = dtmp(:,4)+dtmp2(5)*d2xdy2(:,5)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,3) = dtmp(:,3)+dtmp2(4)*d2xdy2(:,4)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,2) = dtmp(:,2)+dtmp2(3)*d2xdy2(:,3)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,1) = dtmp(:,1)+dtmp2(2)*d2xdy2(:,2)

    !#DIR$ VECTOR ALIGNED
    d2xdy2 = dri2*d2xdy2
  End Subroutine D2bydr2

  !#DEC$ ATTRIBUTES FORCEINLINE :: Dbydt
  Subroutine Dbydt(xt,dxdyt,b1t,b2t,var)
    Implicit None
    Real*8, Intent(InOut) :: xt(:,:), dxdyt(:,:), b1t(:,:), b2t(:,:)
    Integer, Intent(In) :: var
    Real*8 :: b1(mynphi,binfo(2,1)%bsize%bnt), b2(mynphi,binfo(2,2)%bsize%bnt)
    Real*8 :: dtmp(mynphi,mynth), dtmp2(mynth), x(mynphi,mynth), dxdy(mynphi,mynth)
    Real*8 :: c1,c2,c3,c4,c5,alpha2,gamma2,alpha1,beta1,const
    Integer :: j, bdry

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dtmp, dtmp2, x, dxdy, b1, b2

    If (var .ge. 1) Then
       bdry = ibct(var)
    Else
       bdry = -var
    EndIf

    b1 = transpose(b1t)
    b2 = transpose(b2t)
    x = transpose(xt)

    Select Case (dtype2)
    Case (1)
       c1 = -1d0/180d0
       c2 = 9d0
       c3 = 45d0
       c4 = 7d0/9d0
       c5 = 1d0/36d0
       c2 = c1*c2
       c3 = c1*c3
       const = -(c3+c5)
       beta1 = c4-c2
       alpha1 = c1+c5
       c1 = -c1
       c4 = -c4
       !#DIR$ VECTOR ALIGNED
       dxdy(:,1) = c1*b1(:,1) + c2*b1(:,2) + const*b1(:,3) + c4*b1(:,4) + c3*x(:,1) + beta1*x(:,2) + alpha1*x(:,3)
       c4 = 7/9d0
       c5 = 1d0/36d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,2) = c4*(x(:,3)-x(:,1)) + c5*(x(:,4)-b1(:,4))

       c1 = 43d0/96d0
       c2 = 5d0/6d0
       c3 = -9d0/8d0
       c4 = -1d0/6d0
       c5 = 1d0/96d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynth-1) = c1*x(:,mynth) + c2*x(:,mynth-1) + c3*x(:,mynth-2) + c4*x(:,mynth-3) + c5*x(:,mynth-4)
       
       If ((bdry.eq. 3).or.(bdry.eq. 4)) Then
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynth) = b2(:,1)/dti
       Else
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynth)=0d0
       EndIf

       If ((bdry.eq.1).or.(bdry.eq.2))  Then
          ! here are the pentadiagonal and fifth order values for the boundary.
          ! precondition the matrix to make it tridiagonal
          alpha2 = 0.75d0
          gamma2 = 0.125d0
          alpha1 = 6d0
          beta1  = 3d0
          const = 1d0/(alpha2-beta1*gamma2)
          alpha2 = alpha2*const
          beta1 = -beta1*const
          ! fifth order coeff.
          c1 = 10d0/3d0
          c2 = 3d0
          c3 = -6d0
          c4 = -1d0/3d0
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynth) = c1*x(:,mynth) + c2*x(:,mynth-1) + c3*x(:,mynth-2) + c4*x(:,mynth-3)
          !#DIR$ VECTOR ALIGNED
          dxdy(:,mynth) = alpha2*dxdy(:,mynth) + beta1*dxdy(:,mynth-1)
       End If
       c4 = 7d0/9d0
       c5 = 1d0/36d0
    Case (2)
       If (mod(bdry,2).eq.0) Then
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1) = b1(:,1)/dti
       EndIf

       If (bdry.eq.0) Then
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1)=0d0
       EndIf
       
       c1 = -43d0/96d0
       c2 = -5d0/6d0
       c3 = 9d0/8d0
       c4 = 1d0/6d0
       c5 = -1d0/96d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,2) = c1*x(:,1) + c2*x(:,2) + c3*x(:,3) + c4*x(:,4) + c5*x(:,5)
       
       If ((bdry.eq.1).or.(bdry.eq.3)) Then
          ! here are the pentadiagonal and fifth order values for the boundary.
          ! precondition the matrix to make it tridiagonal
          alpha2 = 0.75d0
          gamma2 = 0.125d0
          alpha1 = 6d0
          beta1  = 3d0
          const  = 1d0/(alpha2-beta1*gamma2)
          alpha2 = alpha2*const
          beta1 = -beta1*const

          ! fifth order coeff.
          c1 = -10d0/3d0
          c2 = -3d0
          c3 = 6d0
          c4 = 1d0/3d0
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1) = c1*x(:,1) + c2*x(:,2) + c3*x(:,3) + c4*x(:,4)
          !#DIR$ VECTOR ALIGNED
          dxdy(:,1) = alpha2*dxdy(:,1) + beta1*dxdy(:,2)
       End If
       
       c1 = 7d0/9d0
       c2 = 1d0/36d0
       c3 = -1d0/180d0
       c4 = 9d0
       c5 = 45d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynth-1) = c1*(x(:,mynth) - x(:,mynth-2)) + c2*(b2(:,1) - x(:,mynth-3))

       const = c2+c3*c5
       beta1 = c3*c4-c1
       alpha1 = -(c2+c3)
       c4 = -c3*c4
       c5 = -c3*c5
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynth) = c1*b2(:,1) + const*b2(:,2) + c4*b2(:,3) + c3*b2(:,4) + beta1*x(:,mynth-1) + alpha1*x(:,mynth-2) + c5*x(:,mynth)

       c4 = 7/9d0
       c5 = 1d0/36d0
    Case (3)
       c1 = -1d0/180d0
       c2 = 9d0
       c3 = 45d0
       c4 = 7d0/9d0
       c5 = 1d0/36d0
       c2 = c1*c2
       c3 = c1*c3
       const = -(c3+c5)
       beta1 = c4-c2
       alpha1 = c1+c5
       c1 = -c1
       c4 = -c4
       !#DIR$ VECTOR ALIGNED
       dxdy(:,1) = c1*b1(:,1) + c2*b1(:,2) + const*b1(:,3) + c4*b1(:,4) + c3*x(:,1) + beta1*x(:,2) + alpha1*x(:,3)
       c4 = 7/9d0
       c5 = 1d0/36d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,2) = c4*(x(:,3)-x(:,1)) + c5*(x(:,4)-b1(:,4))
 
       c1 = 7d0/9d0
       c2 = 1d0/36d0
       c3 = -1d0/180d0
       c4 = 9d0
       c5 = 45d0
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynth-1) = c1*(x(:,mynth) - x(:,mynth-2)) + c2*(b2(:,1) - x(:,mynth-3))

       const = c2+c3*c5
       beta1 = c3*c4-c1
       alpha1 = -(c2+c3)
       c4 = -c3*c4
       c5 = -c3*c5
       !#DIR$ VECTOR ALIGNED
       dxdy(:,mynth) = c1*b2(:,1) + const*b2(:,2) + c4*b2(:,3) + c3*b2(:,4) + beta1*x(:,mynth-1) + alpha1*x(:,mynth-2) + c5*x(:,mynth)

       c4 = 7/9d0
       c5 = 1d0/36d0
    EndSelect

    !#DEC$ LOOP COUNT MAX=252, MIN=8, AVG=28
    !#DIR$ VECTOR ALIGNED
    Do j=3,mynth-2
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j) = c4*(x(:,j+1)-x(:,j-1)) + c5*(x(:,j+2)-x(:,j-2))
    EndDo

    bdry = bdry+1
    alpha1 = -1d0/3d0
    c1 = gam3t(2,bdry)
    c2 = gam3t(mynth-1,bdry)
    c3 = gam3t(mynth,bdry)
    dtmp2 = gam2t(:,bdry)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,1) = dtmp2(1)*dxdy(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,2) = dtmp2(2)*dxdy(:,2)+c1*dtmp(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,3) = dtmp2(3)*(dxdy(:,3)+alpha1*dtmp(:,2))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,4) = dtmp2(4)*(dxdy(:,4)+alpha1*dtmp(:,3))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,5) = dtmp2(5)*(dxdy(:,5)+alpha1*dtmp(:,4))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,6) = dtmp2(6)*(dxdy(:,6)+alpha1*dtmp(:,5))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,7) = dtmp2(7)*(dxdy(:,7)+alpha1*dtmp(:,6))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,8) = dtmp2(8)*(dxdy(:,8)+alpha1*dtmp(:,7))

    !#DEC$ LOOP COUNT MAX=123, MIN=1, AVG=11
    !#DIR$ VECTOR ALIGNED
    Do j=9,mynth-3,2
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j) = dtmp2(j)*(dxdy(:,j)+alpha1*dtmp(:,j-1))
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j+1) = dtmp2(j+1)*(dxdy(:,j+1)+alpha1*dtmp(:,j))
    EndDo
    j = mynth-1
    !#DIR$ VECTOR ALIGNED
    dtmp(:,j) = dtmp2(j)*dxdy(:,j)+c2*dtmp(:,j-1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,mynth) = dtmp2(mynth)*dxdy(:,mynth)+c3*dtmp(:,j)

    !Backward substitution
    dtmp2 = gam1t(:,bdry)
    dxdy(:,mynth) = dtmp(:,mynth)
    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=12
    !#DIR$ VECTOR ALIGNED
    Do j=mynth-1,9,-2
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j) = dtmp(:,j)+dtmp2(j+1)*dxdy(:,j+1)
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j-1) = dtmp(:,j-1)+dtmp2(j)*dxdy(:,j)
    EndDo

    !#DIR$ VECTOR ALIGNED
    dxdy(:,7) = dtmp(:,7)+dtmp2(8)*dxdy(:,8)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,6) = dtmp(:,6)+dtmp2(7)*dxdy(:,7)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,5) = dtmp(:,5)+dtmp2(6)*dxdy(:,6)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,4) = dtmp(:,4)+dtmp2(5)*dxdy(:,5)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,3) = dtmp(:,3)+dtmp2(4)*dxdy(:,4)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,2) = dtmp(:,2)+dtmp2(3)*dxdy(:,3)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,1) = dtmp(:,1)+dtmp2(2)*dxdy(:,2)

    !#DIR$ VECTOR ALIGNED
    dxdyt = dti*transpose(dxdy)
  End Subroutine Dbydt !(11*mynth + 19)*mynphi


  Subroutine D2bydt2(xt,d2xdy2t,b1t,b2t,var)
    Implicit None
    Real*8, Intent(InOut) :: xt(:,:), d2xdy2t(:,:), b1t(:,:), b2t(:,:)
    Integer, Intent(In) :: var
    Real*8 :: b1(mynphi,binfo(2,1)%bsize%bnt), b2(mynphi,binfo(2,2)%bsize%bnt)
    Real*8 :: dtmp(mynphi,mynth), dtmp2(mynth), x(mynphi,mynth), d2xdy2(mynphi,mynth)
    Real*8 :: c1,c2,c3,c4,c5,c6,c7,alpha1
    Integer :: j, bdry

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dtmp, dtmp2, x, d2xdy2, b1, b2

    If (var .ge. 1) Then
       bdry = ibct(var)
    Else
       bdry = -var
    EndIf

    b1 = transpose(b1t)
    b2 = transpose(b2t)
    x = transpose(xt)

    Select Case (dtype2)
    Case (1)
       c1 = 12d0/11d0
       c2 = 3d0/44d0
       c3 = -1d0/990d0  
       c4 = 2d0
       c5 = -27d0
       c6 = 270d0
       c7 = -490d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = c3*(c4*(b1(:,1)+x(:,3))+c5*(b1(:,2)+x(:,2))+c6*(b1(:,3)+x(:,1))+c7*b1(:,4))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = d2xdy2(:,1)+c1*(x(:,2)-2d0*x(:,1)+b1(:,4))+c2*(b1(:,3)-2d0*x(:,1)+x(:,3))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,2) = c1*(x(:,3)-2d0*x(:,2)+x(:,1))+c2*(b1(:,4)-2d0*x(:,2)+x(:,4))

       c1=1.2375d0
       c2=-3d0
       c3=2.325d0
       c4=-0.6d0
       c5=0.0375d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth-1)=c1*x(:,mynth)+c2*x(:,mynth-1)+c3*x(:,mynth-2)+c4*x(:,mynth-3)+c5*x(:,mynth-4)

       If ((bdry.eq.1).or.(bdry.eq.2)) Then
          !Fifth order boundary
          c1 = 1955d0/156d0
          c2 = -4057d0/156d0
          c3 = 1117d0/78d0
          c4 = -55d0/78d0
          c5 = -29d0/156d0
          c6 = 7d0/156d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,mynth)=c1*x(:,mynth)+c2*x(:,mynth-1)+c3*x(:,mynth-2)+c4*x(:,mynth-3)+c5*x(:,mynth-4)+c6*x(:,mynth-5)
       Else
          !Fifth Order, with 1st order derivative specified
          c1 = 575d0/193d0/dri !Note dri here
          c2 = -5827d0/2316d0
          c3 = 1987d0/1158d0
          c4 = 465d0/386d0
          c5 = -547d0/1158d0
          c6 = 157d0/2316d0
          
          If (bdry .eq. 0) c1 = 0d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,mynth)=c1*b2(:,1)+c2*x(:,mynth-1)+c3*x(:,mynth-2)+c4*x(:,mynth-3)+c5*x(:,mynth-4)+c6*x(:,mynth-5)
       EndIf
       c1 = 12d0/11d0
       c2 = 3d0/44d0
    Case (2)
       If (mod(bdry,2).eq.0) Then
          !Fifth Order, with 1st order derivative specified
          c1 = -575d0/193d0/dri !Note dri here
          c2 = -5827d0/2316d0
          c3 = 1987d0/1158d0
          c4 = 465d0/386d0
          c5 = -547d0/1158d0
          c6 = 157d0/2316d0

          If (bdry .eq. 0) c1 = 0d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,1)=c1*b1(:,1)+c2*x(:,2)+c3*x(:,3)+c4*x(:,4)+c5*x(:,5)+c6*x(:,6)
       Else
          !Fifth order boundary
          c1 = 1955d0/156d0
          c2 = -4057d0/156d0
          c3 = 1117d0/78d0
          c4 = -55d0/78d0
          c5 = -29d0/156d0
          c6 = 7d0/156d0
          !#DIR$ VECTOR ALIGNED
          d2xdy2(:,1)=c1*x(:,1)+c2*x(:,2)+c3*x(:,3)+c4*x(:,4)+c5*x(:,5)+c6*x(:,6)
       EndIf

       c1=1.2375d0
       c2=-3d0
       c3=2.325d0
       c4=-0.6d0
       c5=0.0375d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,2)=c1*x(:,1)+c2*x(:,2)+c3*x(:,3)+c4*x(:,4)+c5*x(:,5)

       c1 = 12d0/11d0
       c2 = 3d0/44d0
       c3 = -1d0/990d0
       c4 = 2d0
       c5 = -27d0
       c6 = 270d0
       c7 = -490d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth-1) = c1*(x(:,mynth)-2d0*x(:,mynth-1)+x(:,mynth-2))+c2*(x(:,mynth-3)-2d0*x(:,mynth-1)+b2(:,1))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth) = c1*(b2(:,1)-2d0*x(:,mynth)+x(:,mynth-1))+c2*(x(:,mynth-2)-2d0*x(:,mynth)+b2(:,2))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth) = d2xdy2(:,mynth)+c3*(c4*(x(:,mynth-2)+b2(:,4))+c5*(x(:,mynth-1)+b2(:,3))+c6*(x(:,mynth)+b2(:,2))+c7*b2(:,1))
    Case (3) 
       c1 = 12d0/11d0
       c2 = 3d0/44d0
       c3 = -1d0/990d0
       c4 = 2d0
       c5 = -27d0
       c6 = 270d0
       c7 = -490d0
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = c3*(c4*(b1(:,1)+x(:,3))+c5*(b1(:,2)+x(:,2))+c6*(b1(:,3)+x(:,1))+c7*b1(:,4))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,1) = d2xdy2(:,1)+c1*(x(:,2)-2d0*x(:,1)+b1(:,4))+c2*(b1(:,3)-2d0*x(:,1)+x(:,3))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,2) = c1*(x(:,3)-2d0*x(:,2)+x(:,1))+c2*(b1(:,4)-2d0*x(:,2)+x(:,4))
   
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth-1) = c1*(x(:,mynth)-2d0*x(:,mynth-1)+x(:,mynth-2))+c2*(x(:,mynth-3)-2d0*x(:,mynth-1)+b2(:,1))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth) = c1*(b2(:,1)-2d0*x(:,mynth)+x(:,mynth-1))+c2*(x(:,mynth-2)-2d0*x(:,mynth)+b2(:,2))
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,mynth) = d2xdy2(:,mynth)+c3*(c4*(x(:,mynth-2)+b2(:,4))+c5*(x(:,mynth-1)+b2(:,3))+c6*(x(:,mynth)+b2(:,2))+c7*b2(:,1))
    EndSelect

    !#DEC$ LOOP COUNT MAX=252, MIN=8, AVG=12
    !#DIR$ VECTOR ALIGNED
    Do j=3,mynth-2
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j) = c1*(x(:,j+1)-2d0*x(:,j)+x(:,j-1))+c2*(x(:,j-2)-2d0*x(:,j)+x(:,j+2))
    End Do

    bdry = bdry+1
    alpha1 = -2d0/11d0
    c1 = gam3t2(2,bdry)
    c2 = gam3t2(mynth-1,bdry)
    c3 = gam3t2(mynth,bdry)
    dtmp2 = gam2t2(:,bdry)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,1) = dtmp2(1)*d2xdy2(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,2) = dtmp2(2)*d2xdy2(:,2)+c1*dtmp(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,3) = dtmp2(3)*(d2xdy2(:,3)+alpha1*dtmp(:,2))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,4) = dtmp2(4)*(d2xdy2(:,4)+alpha1*dtmp(:,3))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,5) = dtmp2(5)*(d2xdy2(:,5)+alpha1*dtmp(:,4))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,6) = dtmp2(6)*(d2xdy2(:,6)+alpha1*dtmp(:,5))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,7) = dtmp2(7)*(d2xdy2(:,7)+alpha1*dtmp(:,6))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,8) = dtmp2(8)*(d2xdy2(:,8)+alpha1*dtmp(:,7))

    !#DEC$ LOOP COUNT MAX=123, MIN=1, AVG=3
    !#DIR$ VECTOR ALIGNED
    Do j=9,mynth-3,2
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j) = dtmp2(j)*(d2xdy2(:,j)+alpha1*dtmp(:,j-1))
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j+1) = dtmp2(j+1)*(d2xdy2(:,j+1)+alpha1*dtmp(:,j))
    EndDo !2*(mynth-16) + 6*mynth+58 = 8*mynth+26
    j = mynth-1
    !#DIR$ VECTOR ALIGNED
    dtmp(:,j) = dtmp2(j)*d2xdy2(:,j)+c2*dtmp(:,j-1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,mynth) = dtmp2(mynth)*d2xdy2(:,mynth)+c3*dtmp(:,j)

    !Backward substitution
    dtmp2 = gam1t2(:,bdry)
    d2xdy2(:,mynth) = dtmp(:,mynth)
    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=4
    !#DIR$ VECTOR ALIGNED
    Do j=mynth-1,9,-2
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j) = dtmp(:,j)+dtmp2(j+1)*d2xdy2(:,j+1)
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j-1) = dtmp(:,j-1)+dtmp2(j)*d2xdy2(:,j)
    EndDo !8*mynth+26 + 2*(mynth-16) = 10*mynth-6

    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,7) = dtmp(:,7)+dtmp2(8)*d2xdy2(:,8)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,6) = dtmp(:,6)+dtmp2(7)*d2xdy2(:,7)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,5) = dtmp(:,5)+dtmp2(6)*d2xdy2(:,6)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,4) = dtmp(:,4)+dtmp2(5)*d2xdy2(:,5)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,3) = dtmp(:,3)+dtmp2(4)*d2xdy2(:,4)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,2) = dtmp(:,2)+dtmp2(3)*d2xdy2(:,3)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,1) = dtmp(:,1)+dtmp2(2)*d2xdy2(:,2)

    !#DIR$ VECTOR ALIGNED
    d2xdy2t = dti2*transpose(d2xdy2)
  End Subroutine D2bydt2

  !#DEC$ ATTRIBUTES FORCEINLINE :: Dbydp
  Subroutine Dbydp(x,dxdy,b1,b2,var)
    Implicit None
    Real*8, Intent(InOut) :: x(:,:),b1(:,:),b2(:,:)
    Real*8, Intent(InOut) :: dxdy(:,:)
    Integer, Intent(In) :: var
    Real*8 :: dtmp(mynth,mynphi)
    Real*8 :: c1,c2,c3,c4,c5,c6,c7,c8
    Integer :: j

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dtmp

    c1 = -1d0/180d0
    c2 = 9d0
    c3 = 45d0
    c4 = 7d0/9d0
    c5 = 1d0/36d0
    c2 = c1*c2
    c3 = c1*c3
    c6 = -(c3+c5)
    c7 = c4-c2
    c8 = c1+c5
    c1 = -c1
    c4 = -c4
    !#DIR$ VECTOR ALIGNED
    dxdy(:,1) = c1*b1(:,1) + c2*b1(:,2) + c6*b1(:,3) + c4*b1(:,4) + c3*x(:,1) + c7*x(:,2) + c8*x(:,3)
    c4 = 7/9d0
    c5 = 1d0/36d0

    !#DIR$ VECTOR ALIGNED
    dxdy(:,2) = c4*(x(:,3)-x(:,1)) + c5*(x(:,4)-b1(:,4))

    !#DEC$ LOOP COUNT MAX=252, MIN=8, AVG=28
    !#DIR$ VECTOR ALIGNED
    Do j=3,mynphi-2
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j) = c4*(x(:,j+1)-x(:,j-1))+c5*(x(:,j+2)-x(:,j-2))
    EndDo

    c1 = 7d0/9d0
    c2 = 1d0/36d0
    c3 = -1d0/180d0
    c4 = 9d0
    c5 = 45d0
    !#DIR$ VECTOR ALIGNED
    dxdy(:,mynphi-1) = c1*(x(:,mynphi) - x(:,mynphi-2)) + c2*(b2(:,1) - x(:,mynphi-3))

    c6 = c2+c3*c5
    c7 = c3*c4-c1
    c8 = -(c2+c3)
    c4 = -c3*c4
    c5 = -c3*c5
    !#DIR$ VECTOR ALIGNED
    dxdy(:,mynphi) = c1*b2(:,1) + c6*b2(:,2) + c4*b2(:,3) + c3*b2(:,4) + c7*x(:,mynphi-1) + c8*x(:,mynphi-2) + c5*x(:,mynphi)

    c1 = -1d0/3d0
    !#DIR$ VECTOR ALIGNED
    dtmp(:,1) = gam2p(1)*dxdy(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,2) = gam2p(2)*(dxdy(:,2)+c1*dtmp(:,1))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,3) = gam2p(3)*(dxdy(:,3)+c1*dtmp(:,2))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,4) = gam2p(4)*(dxdy(:,4)+c1*dtmp(:,3))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,5) = gam2p(5)*(dxdy(:,5)+c1*dtmp(:,4))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,6) = gam2p(6)*(dxdy(:,6)+c1*dtmp(:,5))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,7) = gam2p(7)*(dxdy(:,7)+c1*dtmp(:,6))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,8) = gam2p(8)*(dxdy(:,8)+c1*dtmp(:,7))

    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=12
    !#DIR$ VECTOR ALIGNED
    Do j=9,mynphi-1,2
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j) = gam2p(j)*(dxdy(:,j)+c1*dtmp(:,j-1))
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j+1) = gam2p(j+1)*(dxdy(:,j+1)+c1*dtmp(:,j))
    EndDo

    j = mynphi
    dxdy(:,j) = dtmp(:,j)

    !Backward substitution
    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=12
    !#DIR$ VECTOR ALIGNED
    Do j=mynphi-1,9,-2
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j) = dtmp(:,j)+gam1p(j+1)*dxdy(:,j+1)
       !#DIR$ VECTOR ALIGNED
       dxdy(:,j-1) = dtmp(:,j-1)+gam1p(j)*dxdy(:,j)
    EndDo

    !#DIR$ VECTOR ALIGNED
    dxdy(:,7) = dtmp(:,7)+gam1p(8)*dxdy(:,8)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,6) = dtmp(:,6)+gam1p(7)*dxdy(:,7)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,5) = dtmp(:,5)+gam1p(6)*dxdy(:,6)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,4) = dtmp(:,4)+gam1p(5)*dxdy(:,5)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,3) = dtmp(:,3)+gam1p(4)*dxdy(:,4)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,2) = dtmp(:,2)+gam1p(3)*dxdy(:,3)
    !#DIR$ VECTOR ALIGNED
    dxdy(:,1) = dtmp(:,1)+gam1p(2)*dxdy(:,2)

    !#DIR$ VECTOR ALIGNED
    dxdy = dpi*dxdy
  End Subroutine Dbydp !(11*mynphi + 20)*mynth

  !#DEC$ ATTRIBUTES FORCEINLINE :: D2bydp2
  Subroutine D2bydp2(x,d2xdy2,b1,b2,var)
    Implicit None
    Real*8, Intent(InOut) :: x(:,:),b1(:,:),b2(:,:)
    Real*8, Intent(InOut) :: d2xdy2(:,:)
    Integer, Intent(In) :: var
    Real*8 :: dtmp(mynth,mynphi)
    Real*8 :: c1,c2,c3,c4,c5,c6,c7
    Integer :: j

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: dtmp

    c1 = 12d0/11d0
    c2 = 3d0/44d0
    c3 = -1d0/990d0
    c4 = 2d0
    c5 = -27d0
    c6 = 270d0
    c7 = -490d0
    d2xdy2(:,1) = c3*(c4*(b1(:,1)+x(:,3))+c5*(b1(:,2)+x(:,2))+c6*(b1(:,3)+x(:,1))+c7*b1(:,4))
    d2xdy2(:,1) = d2xdy2(:,1)+c1*(x(:,2)-2d0*x(:,1)+b1(:,4))+c2*(b1(:,3)-2d0*x(:,1)+x(:,3))
    d2xdy2(:,2) = c1*(x(:,3)-2d0*x(:,2)+x(:,1))+c2*(b1(:,4)-2d0*x(:,2)+x(:,4))

    !#DEC$ LOOP COUNT MAX=252, MIN=8, AVG=28
    Do j=3,mynphi-2
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j) = c1*(x(:,j+1)-2d0*x(:,j)+x(:,j-1))+c2*(x(:,j-2)-2d0*x(:,j)+x(:,j+2))
    End Do

    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,mynphi-1) = c1*(x(:,mynphi)-2d0*x(:,mynphi-1)+x(:,mynphi-2))+c2*(x(:,mynphi-3)-2d0*x(:,mynphi-1)+b2(:,1))
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,mynphi) = c1*(b2(:,1)-2d0*x(:,mynphi)+x(:,mynphi-1))+c2*(x(:,mynphi-2)-2d0*x(:,mynphi)+b2(:,2))
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,mynphi) = d2xdy2(:,mynphi)+c3*(c4*(x(:,mynphi-2)+b2(:,4))+c5*(x(:,mynphi-1)+b2(:,3))+c6*(x(:,mynphi)+b2(:,2))+c7*b2(:,1))

    c1 = -2d0/11d0
    !#DIR$ VECTOR ALIGNED
    dtmp(:,1) = gam2p2(1)*d2xdy2(:,1)
    !#DIR$ VECTOR ALIGNED
    dtmp(:,2) = gam2p2(2)*(d2xdy2(:,2)+c1*dtmp(:,1))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,3) = gam2p2(3)*(d2xdy2(:,3)+c1*dtmp(:,2))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,4) = gam2p2(4)*(d2xdy2(:,4)+c1*dtmp(:,3))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,5) = gam2p2(5)*(d2xdy2(:,5)+c1*dtmp(:,4))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,6) = gam2p2(6)*(d2xdy2(:,6)+c1*dtmp(:,5))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,7) = gam2p2(7)*(d2xdy2(:,7)+c1*dtmp(:,6))
    !#DIR$ VECTOR ALIGNED
    dtmp(:,8) = gam2p2(8)*(d2xdy2(:,8)+c1*dtmp(:,7))

    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=12
    Do j=9,mynphi-1,2
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j) = gam2p2(j)*(d2xdy2(:,j)+c1*dtmp(:,j-1))
       !#DIR$ VECTOR ALIGNED
       dtmp(:,j+1) = gam2p2(j+1)*(d2xdy2(:,j+1)+c1*dtmp(:,j))
    EndDo !2*(mynphi-16) + 6*mynphi+58 = 8*mynphi+26

    !Backward substitution
    j = mynphi
    d2xdy2(:,j) = dtmp(:,j)
    !#DEC$ LOOP COUNT MAX=124, MIN=2, AVG=12
    Do j=mynphi-1,9,-2
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j) = dtmp(:,j)+gam1p2(j+1)*d2xdy2(:,j+1)
       !#DIR$ VECTOR ALIGNED
       d2xdy2(:,j-1) = dtmp(:,j-1)+gam1p2(j)*d2xdy2(:,j)
    EndDo !8*mynphi+26 + 2*(mynphi-16) = 10*mynphi-6

    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,7) = dtmp(:,7)+gam1p2(8)*d2xdy2(:,8)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,6) = dtmp(:,6)+gam1p2(7)*d2xdy2(:,7)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,5) = dtmp(:,5)+gam1p2(6)*d2xdy2(:,6)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,4) = dtmp(:,4)+gam1p2(5)*d2xdy2(:,5)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,3) = dtmp(:,3)+gam1p2(4)*d2xdy2(:,4)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,2) = dtmp(:,2)+gam1p2(3)*d2xdy2(:,3)
    !#DIR$ VECTOR ALIGNED
    d2xdy2(:,1) = dtmp(:,1)+gam1p2(2)*d2xdy2(:,2)

    !#DIR$ VECTOR ALIGNED
    d2xdy2 = dpi2*d2xdy2 !11*mynphi+24

  End Subroutine D2bydp2

  Subroutine VectorDivergence()
    Implicit None
    Real*8 :: cscth(mynth,mynphi), sinth(mynth,mynphi), rads(mynth,mynr), irads(mynth,mynr), work1(mynth,mynr)
    Real*8 :: radrb1(mynth,binfo(1,1)%bsize%bnr), radrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: workt(mynth,mynphi), workp(mynth,mynphi), workr(mynth,mynr)
    Real*8 :: dworkt(mynth,mynphi), dworkp(mynth,mynphi), dworkr(mynth,mynr)
    Real*8 :: sintb1(binfo(2,1)%bsize%bnt,mynphi), worktb1(binfo(2,1)%bsize%bnt,mynphi), workpb1(mynth,binfo(3,1)%bsize%bnp), workrb1(mynth,binfo(1,1)%bsize%bnr)
    Real*8 :: sintb2(binfo(2,2)%bsize%bnt,mynphi), worktb2(binfo(2,2)%bsize%bnt,mynphi), workpb2(mynth,binfo(3,2)%bsize%bnp), workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: invr
    Integer :: ir, it, ip, irr, itt, ipp, bdry

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: cscth, sinth, workt, dworkt, sintb1, sintb2, worktb1, worktb2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: rads, irads, work1, radrb1, radrb2, workr, dworkr, workrb1, workrb2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: workp, dworkp, workpb1, workpb2

    Do ip=1,mynphi
       sinth(:,ip) = sines(myth1:myth2)
       !#DIR$ VECTOR ALIGNED
       cscth(:,ip) = 1d0/sines(myth1:myth2)
    EndDo

    If (myth1.eq.1) Then
       sintb1(1,:) = sines(1)
    Else
       Do ipp=1,mynphi
          sintb1(:,ipp) = sines(myth1-4:myth1-1)
       EndDo
    EndIf

    If (myth2.eq.nth) Then
       If (Theta_Symmetric) Then
          Do ipp=1,mynphi
             sintb2(:,ipp) = sin(th2+dth*(/(Dble(it),it=1,4)/))
          EndDo
       Else
          sintb2(1,:) = sines(nth)
       EndIf
    Else
       Do ipp=1,mynphi
          sintb2(:,ipp) = sines(myth2+1:myth2+4)
       EndDo
    EndIf

    bdry = -1

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do irr=1,mynr
       invr = r_inv(myr1+irr-1)          
       !dBt/dt
       !#DIR$ VECTOR ALIGNED
       workt = sinth*vars1(:,:,irr,6)
       !#DIR$ VECTOR ALIGNED
       worktb1 = sintb1*binfo(2,1)%bdata(:,:,irr,6)
       !#DIR$ VECTOR ALIGNED
       worktb2 = sintb2*binfo(2,2)%bdata(:,:,irr,6)
       Call dbydt(workt,dworkt,worktb1,worktb2,bdry)
       !#DIR$ VECTOR ALIGNED
       invrho(:,:,irr) = invr*dworkt*cscth
    EndDo

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do irr=1,mynr
       invr = r_inv(myr1+irr-1)
       !dBp/dp
       workp = vars1(:,:,irr,7)
       workpb1 = binfo(3,1)%bdata(:,:,irr,7)
       workpb2 = binfo(3,2)%bdata(:,:,irr,7)
       Call dbydp(workp,dworkp,workpb1,workpb2,0)
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ipp=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             invrho(itt,ipp,irr) = invrho(itt,ipp,irr) + invr*dworkp(itt,ipp)*cscth(itt,ipp)
          EndDo
       EndDo
    EndDo

    If (non_uniform) Then
       Do ipp=1,mynth
          work1(ipp,:)=dxdr(myr1:myr2)
       EndDo
    EndIf

    Do ipp=1,mynth
       !#DIR$ VECTOR ALIGNED
       rads(ipp,:) = (r_ind(myr1:myr2)/r2)**2
       !#DIR$ VECTOR ALIGNED
       irads(ipp,:) = (r2*r_inv(myr1:myr2))**2
    EndDo

    If (myr1 .eq. 1) Then
       !#DIR$ VECTOR ALIGNED
       radrb1(:,1) = (r1/r2)**2
    Else
       !#DIR$ VECTOR ALIGNED
       Do ipp=1,mynth
          !#DIR$ VECTOR ALIGNED
          radrb1(ipp,:) = (r_ind(myr1-4:myr1-1)/r2)**2
       EndDo
    EndIf

    If (myr2 .eq. nr) Then
       radrb2(:,1) = 1d0
    Else
       !#DIR$ VECTOR ALIGNED
       Do ipp=1,mynth
          !#DIR$ VECTOR ALIGNED
          radrb2(ipp,:) = (r_ind(myr2+1:myr2+4)/r2)**2
       EndDo
    EndIf

    bdry = -1

    !dBr/dr
    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    !#DIR$ VECTOR ALIGNED
    Do ipp=1,mynphi
       !#DIR$ VECTOR ALIGNED
       workr = rads*vars1(:,ipp,:,8)
       !#DIR$ VECTOR ALIGNED
       workrb1 = radrb1*binfo(1,1)%bdata(:,ipp,:,8)
       !#DIR$ VECTOR ALIGNED
       workrb2 = radrb2*binfo(1,2)%bdata(:,ipp,:,8)
       Call dbydr(workr,dworkr,workrb1,workrb2,bdry)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          dworkr=dworkr*work1
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do irr=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             invrho(itt,ipp,irr) = invrho(itt,ipp,irr) + irads(itt,irr)*dworkr(itt,irr)
          EndDo
       EndDo
    EndDo
  End Subroutine VectorDivergence

  Subroutine Compute_B(ind)
    Implicit None
    Integer, Intent(In) :: ind
    
    If (ind .eq. 0) Then
       maga = vars0(:,:,:,6:8)
    Else
       maga = vars1(:,:,:,6:8)
    EndIf
    
    Call VectorCurlA()
    
    !Boundary Conditions
    !  zero normal mag field on latitudinal boundaries
    If (Theta_Perf_Conductor) Then
       If (myth1 .eq. 1) Then
          dvars1(1,:,:,6)=0d0
       EndIf
       If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
          dvars1(mynth,:,:,6)=0d0
       EndIf
    Else
       !  zero tangential field on theta boundaries
       If (myth1 .eq. 1) Then
          dvars1(1,:,:,7)=0d0
          dvars1(1,:,:,8)=0d0
       EndIf
       If ((myth2 .eq. nth).and.(.not.Theta_Symmetric)) Then
          dvars1(mynth,:,:,7)=0d0
          dvars1(mynth,:,:,8)=0d0
       EndIf
    EndIf

    !  zero tangential field on radial boundaries
    If (myr1.Eq.1) Then
       If (.not. bot_inout) Then
          If (Bottom_Radial_Field) Then
             dvars1(:,:,1,6)=0d0
             If (.not. Unstratified) dvars1(:,:,1,7)=0d0
          ElseIf (Bottom_Perfect_Cond) Then
             dvars1(:,:,1,8) = 0d0
          EndIf
       EndIf
    EndIf
    If (myr2 .eq. nr) Then
       If (.not. top_inout) then
          If (Top_Radial_Field) Then
             dvars1(:,:,mynr,6)=0d0
             If (.not. Unstratified) dvars1(:,:,mynr,7)=0d0
          ElseIf (Top_Perfect_Cond) Then
             dvars1(:,:,mynr,8) = 0d0
          EndIf
       EndIf
    EndIf

    If (ind .eq. 0) Then
       vars0(:,:,:,6:8) = dvars1(:,:,:,6:8)
    Else
       vars1(:,:,:,6:8) = dvars1(:,:,:,6:8)
    EndIf

    dvars1(:,:,:,6:8) = 0d0
  End Subroutine Compute_B

  Subroutine VectorCurl()
    Implicit None
    Real*8 :: cscth(mynth,mynphi), sinth(mynth,mynphi), work1(mynth,mynr), finvr(mynth,mynr), indr(mynth,mynr)
    Real*8 :: radrb1(mynth,binfo(1,1)%bsize%bnr), radrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: workt(mynth,mynphi), workp(mynth,mynphi), workr(mynth,mynr)
    Real*8 :: dworkt(mynth,mynphi), dworkp(mynth,mynphi), dworkr(mynth,mynr)
    Real*8 :: sintb1(binfo(2,1)%bsize%bnt,mynphi), worktb1(binfo(2,1)%bsize%bnt,mynphi), workpb1(mynth,binfo(3,1)%bsize%bnp), workrb1(mynth,binfo(1,1)%bsize%bnr)
    Real*8 :: sintb2(binfo(2,2)%bsize%bnt,mynphi), worktb2(binfo(2,2)%bsize%bnt,mynphi), workpb2(mynth,binfo(3,2)%bsize%bnp), workrb2(mynth,binfo(1,2)%bsize%bnr)
    Real*8 :: invr
    Integer :: it, irr, itt, ipp, bdry

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: cscth, sinth, workt, dworkt, sintb1, sintb2, worktb1, worktb2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: indr, finvr, work1, radrb1, radrb2, workr, dworkr, workrb1, workrb2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: workp, dworkp, workpb1, workpb2

    Do ipp=1,mynphi
       sinth(:,ipp) = sines(myth1:myth2)
       !#DIR$ VECTOR ALIGNED
       cscth(:,ipp) = 1d0/sines(myth1:myth2)
    EndDo

    If (myth1.eq.1) Then
       sintb1(1,:) = sines(1)
    Else
       Do ipp=1,mynphi
          sintb1(:,ipp) = sines(myth1-4:myth1-1)
       EndDo
    EndIf

    If (myth2.eq.nth) Then
       If (Theta_Symmetric) Then
          Do ipp=1,mynphi
             sintb2(:,ipp) = sin(th2+dth*(/(Dble(it),it=1,4)/))
          EndDo
       Else
          sintb2(1,:) = sines(nth)
       EndIf
    Else
       Do ipp=1,mynphi
          sintb2(:,ipp) = sines(myth2+1:myth2+4)
       EndDo
    EndIf

    If (Theta_Normal_Field) Then
       bdry = -1
    Else
       bdry = 0
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do irr=1,mynr
       invr = r_inv(myr1+irr-1)
       !dBp/dt
       !#DIR$ VECTOR ALIGNED
       workt = sinth*vars1(:,:,irr,7)
       !#DIR$ VECTOR ALIGNED
       worktb1 = sintb1*binfo(2,1)%bdata(:,:,irr,7)
       !#DIR$ VECTOR ALIGNED
       worktb2 = sintb2*binfo(2,2)%bdata(:,:,irr,7)
       Call dbydt(workt,dworkt,worktb1,worktb2,bdry)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,irr,8) = dworkt*cscth*invr

       !dBr/dt
       workt = vars1(:,:,irr,8)
       worktb1 = binfo(2,1)%bdata(:,:,irr,8)
       worktb2 = binfo(2,2)%bdata(:,:,irr,8)
       Call dbydt(workt,dworkt,worktb1,worktb2,bdry)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,irr,7) = -dworkt*invr
    EndDo

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do irr=1,mynr
       invr = r_inv(myr1+irr-1)
       !dBr/dp
       workp = vars1(:,:,irr,8)
       workpb1 = binfo(3,1)%bdata(:,:,irr,8)
       workpb2 = binfo(3,2)%bdata(:,:,irr,8)
       Call dbydp(workp,dworkp,workpb1,workpb2,0)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,irr,6) = invr*dworkp*cscth

       !dBt/dp
       workp = vars1(:,:,irr,6)
       workpb1 = binfo(3,1)%bdata(:,:,irr,6)
       workpb2 = binfo(3,2)%bdata(:,:,irr,6)
       Call dbydp(workp,dworkp,workpb1,workpb2,0)
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ipp=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             dvars1(itt,ipp,irr,8) = dvars1(itt,ipp,irr,8) - invr*dworkp(itt,ipp)*cscth(itt,ipp)
          EndDo
       EndDo
    EndDo

    If (non_uniform) Then
       Do ipp=1,mynth
          work1(ipp,:)=dxdr(myr1:myr2)
       EndDo
    EndIf
    Do ipp=1,mynth
       finvr(ipp,:) = r_inv(myr1:myr2)
       indr(ipp,:) = r_ind(myr1:myr2)
    EndDo

    If (myr1 .eq. 1) Then
       radrb1(:,1) = r1
    Else
       Do ipp=1,mynth
          radrb1(ipp,:) = r_ind(myr1-4:myr1-1)
       EndDo
    EndIf

    If (myr2 .eq. nr) Then
       radrb2(:,1) = r2
    Else
       Do ipp=1,mynth
          radrb2(ipp,:) = r_ind(myr2+1:myr2+4)
       EndDo
    EndIf

    If ((Bottom_Radial_Field).and.(Top_Radial_Field)) Then
       bdry = -1
    EndIf
       
    If ((.not. Bottom_Radial_Field).and.(Top_Radial_Field)) Then
       bdry = -2
       If (myr1.eq.1) Then
          binfo(1,1)%bdata(:,:,:,6) = 0d0
          binfo(1,1)%bdata(:,:,:,7) = 0d0
       EndIf
    EndIf

    If ((Bottom_Radial_Field).and.(.not.Top_Radial_Field)) Then
       bdry = -3
       If (myr2.eq.nr) Then
          binfo(1,2)%bdata(:,:,:,6) = 0d0
          binfo(1,2)%bdata(:,:,:,7) = 0d0
       EndIf
    EndIf

    If ((.not.Bottom_Radial_Field).and.(.not.Top_Radial_Field)) Then
       bdry = 0
    EndIf

    If (Unstratified) Then
       bdry=-4
       If (myr1.eq.1) Then
          binfo(1,1)%bdata(:,:,:,6) = 0d0
          binfo(1,1)%bdata(:,:,1,7) = -vars1(:,:,1,7)/r1
       EndIf
       If (myr2.eq.nr) Then
          binfo(1,2)%bdata(:,:,:,6) = 0d0
          binfo(1,2)%bdata(:,:,1,7) = -vars1(:,:,mynr,7)/r2
       EndIf
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    !#DIR$ VECTOR ALIGNED
    Do ipp=1,mynphi       
       !dBt/dr
       !#DIR$ VECTOR ALIGNED
       workr = indr*vars1(:,ipp,:,6)
       !#DIR$ VECTOR ALIGNED
       workrb1 = radrb1*binfo(1,1)%bdata(:,ipp,:,6)
       !#DIR$ VECTOR ALIGNED
       workrb2 = radrb2*binfo(1,2)%bdata(:,ipp,:,6)
       Call dbydr(workr,dworkr,workrb1,workrb2,7)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          dworkr=dworkr*work1
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do irr=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             dvars1(itt,ipp,irr,7) = dvars1(itt,ipp,irr,7) + dworkr(itt,irr)*finvr(itt,irr)
          EndDo
       EndDo

       !dBp/dr
       !#DIR$ VECTOR ALIGNED
       workr = indr*vars1(:,ipp,:,7)
       !#DIR$ VECTOR ALIGNED
       workrb1 = radrb1*binfo(1,1)%bdata(:,ipp,:,7)
       !#DIR$ VECTOR ALIGNED
       workrb2 = radrb2*binfo(1,2)%bdata(:,ipp,:,7)
       Call dbydr(workr,dworkr,workrb1,workrb2,8)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          dworkr=dworkr*work1
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do irr=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             dvars1(itt,ipp,irr,6) = dvars1(itt,ipp,irr,6) - dworkr(itt,irr)*finvr(itt,irr)
          EndDo
       EndDo
    EndDo
  End Subroutine VectorCurl

  Subroutine VectorCurlA()
    Implicit None
    Real*8 :: cscth(mynth,mynphi), sinth(mynth,mynphi), work1(mynth,mynr), finvr(mynth,mynr), indr(mynth,mynr)
    Real*8 :: radrb1(mynth,dbinfo(1,1)%bsize%bnr), radrb2(mynth,dbinfo(1,2)%bsize%bnr)
    Real*8 :: workt(mynth,mynphi), workp(mynth,mynphi), workr(mynth,mynr)
    Real*8 :: dworkt(mynth,mynphi), dworkp(mynth,mynphi), dworkr(mynth,mynr)
    Real*8 :: sintb1(dbinfo(2,1)%bsize%bnt,mynphi), worktb1(dbinfo(2,1)%bsize%bnt,mynphi), workpb1(mynth,dbinfo(3,1)%bsize%bnp), workrb1(mynth,dbinfo(1,1)%bsize%bnr)
    Real*8 :: sintb2(dbinfo(2,2)%bsize%bnt,mynphi), worktb2(dbinfo(2,2)%bsize%bnt,mynphi), workpb2(mynth,dbinfo(3,2)%bsize%bnp), workrb2(mynth,dbinfo(1,2)%bsize%bnr)
    Real*8 :: invr
    Integer :: irr, itt, ipp, it

    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: cscth, sinth, workt, dworkt, sintb1, sintb2, worktb1, worktb2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: indr, finvr, work1, radrb1, radrb2, workr, dworkr, workrb1, workrb2
    !#DIR$ATTRIBUTES ALIGN : VECWIDTH :: workp, dworkp, workpb1, workpb2

    Call Boundary_Exchange_Mag_Init(2)

    Do ipp=1,mynphi
       sinth(:,ipp) = sines(myth1:myth2)
       !#DIR$ VECTOR ALIGNED
       cscth(:,ipp) = 1d0/sines(myth1:myth2)
    EndDo

    If (myth1.eq.1) Then
       sintb1(1,:) = sines(1)
    Else
       Do ipp=1,mynphi
          sintb1(:,ipp) = sines(myth1-4:myth1-1)
       EndDo
    EndIf

    If (myth2.eq.nth) Then
       If (Theta_Symmetric) Then
          Do ipp=1,mynphi
             sintb2(:,ipp) = sin(th2+dth*(/(Dble(it),it=1,4)/))
          EndDo
       Else
          sintb2(1,:) = sines(nth)
       EndIf
    Else
       Do ipp=1,mynphi
          sintb2(:,ipp) = sines(myth2+1:myth2+4)
       EndDo
    EndIf

    Call Boundary_Exchange_Mag_Finish(2)
    Call Boundary_Exchange_Mag_Init(3)

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    !#DIR$ VECTOR ALIGNED
    Do irr=1,mynr
       invr = r_inv(myr1+irr-1)

       !dBp/dt
       !#DIR$ VECTOR ALIGNED
       workt = sinth*maga(:,:,irr,2)
       !#DIR$ VECTOR ALIGNED
       worktb1 = sintb1*dbinfo(2,1)%bdata(:,:,irr,2)
       !#DIR$ VECTOR ALIGNED
       worktb2 = sintb2*dbinfo(2,2)%bdata(:,:,irr,2)
       Call dbydt(workt,dworkt,worktb1,worktb2,11)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,irr,8) = dworkt*cscth*invr

       !dBr/dt
       workt = maga(:,:,irr,3)
       worktb1 = dbinfo(2,1)%bdata(:,:,irr,3)
       worktb2 = dbinfo(2,2)%bdata(:,:,irr,3)
       Call dbydt(workt,dworkt,worktb1,worktb2,12)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,irr,7) = -dworkt*invr
    EndDo

    Call Boundary_Exchange_Mag_Finish(3)
    Call Boundary_Exchange_Mag_Init(1)

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=16
    Do irr=1,mynr
       invr = r_inv(myr1+irr-1)
       !dBr/dp
       workp = maga(:,:,irr,3)
       workpb1 = dbinfo(3,1)%bdata(:,:,irr,3)
       workpb2 = dbinfo(3,2)%bdata(:,:,irr,3)
       Call dbydp(workp,dworkp,workpb1,workpb2,0)
       !#DIR$ VECTOR ALIGNED
       dvars1(:,:,irr,6) = invr*dworkp*cscth

       !dBt/dp
       workp = maga(:,:,irr,1)
       workpb1 = dbinfo(3,1)%bdata(:,:,irr,1)
       workpb2 = dbinfo(3,2)%bdata(:,:,irr,1)
       Call dbydp(workp,dworkp,workpb1,workpb2,0)
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do ipp=1,mynphi
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             dvars1(itt,ipp,irr,8) = dvars1(itt,ipp,irr,8)-invr*dworkp(itt,ipp)*cscth(itt,ipp)
          EndDo
       EndDo
    EndDo

    If (non_uniform) Then
       Do ipp=1,mynth
          work1(ipp,:)=dxdr(myr1:myr2)
       EndDo
    EndIf
    Do ipp=1,mynth
       finvr(ipp,:) = r_inv(myr1:myr2)
       indr(ipp,:) = r_ind(myr1:myr2)
    EndDo

    If (myr1 .eq. 1) Then
       radrb1(:,1) = r1
    Else
       Do ipp=1,mynth
          radrb1(ipp,:) = r_ind(myr1-4:myr1-1)
       EndDo
    EndIf

    If (myr2 .eq. nr) Then
       radrb2(:,1) = r2
    Else
       Do ipp=1,mynth
          radrb2(ipp,:) = r_ind(myr2+1:myr2+4)
       EndDo
    EndIf

    Call Boundary_Exchange_Mag_Finish(1)

    If ((Bottom_Radial_Field) .and. (myr1.eq.1)) Then
       dbinfo(1,1)%bdata(:,:,:,1) = 0d0
       dbinfo(1,1)%bdata(:,:,:,2) = 0d0
    EndIf

    If ((Top_Radial_Field).and.(myr2.eq.nr)) Then
       dbinfo(1,2)%bdata(:,:,:,1) = 0d0
       dbinfo(1,2)%bdata(:,:,:,2) = 0d0
    EndIf

    !#DEC$ LOOP COUNT MAX=256, MIN=16, AVG=32
    !#DIR$ VECTOR ALIGNED
    Do ipp=1,mynphi       
       !dBt/dr
       !#DIR$ VECTOR ALIGNED
       workr = indr*maga(:,ipp,:,1)
       !#DIR$ VECTOR ALIGNED
       workrb1 = radrb1*dbinfo(1,1)%bdata(:,ipp,:,1)
       !#DIR$ VECTOR ALIGNED
       workrb2 = radrb2*dbinfo(1,2)%bdata(:,ipp,:,1)
       Call dbydr(workr,dworkr,workrb1,workrb2,10)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          dworkr=dworkr*work1
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do irr=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             dvars1(itt,ipp,irr,7) = dvars1(itt,ipp,irr,7) + dworkr(itt,irr)*finvr(itt,irr)
          EndDo
       EndDo

       !dBp/dr
       !#DIR$ VECTOR ALIGNED
       workr = indr*maga(:,ipp,:,2)
       !#DIR$ VECTOR ALIGNED
       workrb1 = radrb1*dbinfo(1,1)%bdata(:,ipp,:,2)
       !#DIR$ VECTOR ALIGNED
       workrb2 = radrb2*dbinfo(1,2)%bdata(:,ipp,:,2)
       Call dbydr(workr,dworkr,workrb1,workrb2,11)
       If (non_uniform) Then
          !#DIR$ VECTOR ALIGNED
          dworkr=dworkr*work1
       EndIf
       !#DIR$ VECTOR ALIGNED
       !#DIR$ IVDEP
       Do irr=1,mynr
          !#DIR$ VECTOR ALIGNED
          !#DIR$ IVDEP
          Do itt=1,mynth
             dvars1(itt,ipp,irr,6) = dvars1(itt,ipp,irr,6) - dworkr(itt,irr)*finvr(itt,irr)
          EndDo
       EndDo
    EndDo
  End Subroutine VectorCurlA

  Subroutine Dbyd1_dt0(x,dxdy,bdat1,bdat2,n1,ibc)
    Implicit None
    Integer, Intent(In) :: n1,ibc
    Real*8, Intent(InOut) :: dxdy(n1)
    Real*8, Intent(InOut) :: x(n1), bdat1, bdat2
    Real*8 :: dtmp(n1),gam(n1)
    Real*8 :: c1,c2,c3,c4,c5,alpha,alpha1,alpha2,beta1,gamma2,const,up1,up2
    Integer :: j

    !-------------------------------------------------------------
    c1 = 7d0/9d0 !nops = 3
    c2 = 1d0/36d0 !nops = 4

    If (ibc .eq. 2) Then
       dxdy(1) = bdat1
    ElseIf (ibc .eq. 3) Then
       dxdy(n1) = bdat2
    ElseIf (ibc .eq. 4) Then
       dxdy(1)  = bdat1
       dxdy(n1) = bdat2
    EndIf

    Do j=3,n1-2
       dxdy(j) = c1*(x(j+1)-x(j-1))+c2*(x(j+2)-x(j-2)) !nops = 5*(n1-5)+24 = 5*n1-1
    End Do

    alpha=1d0/3d0
    ! here are the pentadiagonal and fifth order values for the boundary.
    alpha2=3d0/4d0
    gamma2=1d0/8d0
    alpha1=6d0
    beta1=3d0
    ! precondition the matrix to make it tridiagonal
    const=1d0/(alpha2-beta1*gamma2)
    up1=(alpha1*alpha2-beta1)*const

    c1 = -43d0/96d0
    c2 = -5d0/6d0
    c3 = 9d0/8d0
    c4 = 1d0/6d0
    c5 = -1d0/96d0
    dxdy(2) = c1*x(1) + c2*x(2) + c3*x(3) + c4*x(4) + c5*x(5) 

    ! fifth order coeff.
    c1 = -10d0/3d0
    c2 = -3d0
    c3 = 6d0
    c4 = 1d0/3d0
    If (ibc.eq.0) Then
       dxdy(1)=0d0
       ! if ibc = 2,4 Then dxdy(*,1) must be set by the calling routine
    Else If ((ibc.eq.1).or.(ibc.eq.3)) Then
       dxdy(1)=c1*x(1)+c2*x(2)+c3*x(3)+c4*x(4)
    End If

    If (mod(ibc,2).ne.0) Then
       dxdy(1)=(alpha2*dxdy(1)-beta1*dxdy(2))*const
    End If

    If (mod(ibc,2).EQ.0) Then
       up1 = 0d0
    End If

    c1 = -43d0/96d0
    c2 = -5d0/6d0
    c3 = 9d0/8d0
    c4 = 1d0/6d0
    c5 = -1d0/96d0
    dxdy(n1-1) =-(c1*x(n1) + c2*x(n1-1) + c3*x(n1-2) + c4*x(n1-3) + c5*x(n1-4))

    ! fifth order coeff.
    c1 = -10d0/3d0
    c2 = -3d0
    c3 = 6d0
    c4 = 1d0/3d0
    If (ibc.eq.0) Then
       dxdy(n1)=0d0
       ! if ibc = 3,4 Then dxdy(*,n2) must be set by the calling routine
    Else If ((ibc.eq.1).or.(ibc.eq.2))  Then
       dxdy(n1)=-(c1*x(n1)+c2*x(n1-1)+c3*x(n1-2)+c4*x(n1-3))
    End If

    If ((ibc.eq.1).or.(ibc.eq.2)) Then   
       dxdy(n1)=(alpha2*dxdy(n1)-beta1*dxdy(n1-1))*const
    End If

    If (ibc.ge.0) Then
       If ((ibc.ne.1).and.(ibc.ne.2)) Then
          up2 = 0d0
       Else
          up2 = up1
       End If
    End If

    c1 = 1d0
    c2 = 1d0
    dtmp(1) = dxdy(1)

    c3 = up1*c2 !1 
    gam(2) = c3
    c2 = c1/(c1-gamma2*c3) !3 
    dtmp(2) = (dxdy(2)-gamma2*dtmp(1))*c2 !3 

    c3 = alpha2*c2 !1 
    gam(3) = c3
    c2 = c1/(c1-alpha*c3) !3
    dtmp(3) = (dxdy(3)-alpha*dtmp(2))*c2 !3  

    c3 = alpha*c2 !1 
    gam(4) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(4) = (dxdy(4)-alpha*dtmp(3))*c2 !3  

    c3 = alpha*c2 !1 
    gam(5) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(5) = (dxdy(5)-alpha*dtmp(4))*c2 !3  

    c3 = alpha*c2 !1 
    gam(6) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(6) = (dxdy(6)-alpha*dtmp(5))*c2 !3  

    c3 = alpha*c2 !1 
    gam(7) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(7) = (dxdy(7)-alpha*dtmp(6))*c2 !3  

    c3 = alpha*c2 !1 
    gam(8) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(8) = (dxdy(8)-alpha*dtmp(7))*c2 !3  

    c3 = alpha*c2 !1 
    gam(9) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(9) = (dxdy(9)-alpha*dtmp(8))*c2 !3  

    c3 = alpha*c2 !1 
    gam(10) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(10) = (dxdy(10)-alpha*dtmp(9))*c2 !3  

    c3 = alpha*c2 !1 
    gam(11) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(11) = (dxdy(11)-alpha*dtmp(10))*c2 !3  

    c3 = alpha*c2 !1 
    gam(12) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(12) = (dxdy(12)-alpha*dtmp(11))*c2 !3  

    c3 = alpha*c2 !1 
    gam(13) = c3
    c2 = c1/(c1-alpha*c3) !3 
    dtmp(13) = (dxdy(13)-alpha*dtmp(12))*c2 !3  

    Do j=14,n1-2
       c3 = alpha*c2 !1 
       gam(j) = c3
       c2 = c1/(c1-alpha*c3) !3 
       dtmp(j) = (dxdy(j)-alpha*dtmp(j-1))*c2 !3 
    EndDo

    c3 = alpha*c2 !1 
    gam(n1-1) = c3
    c2 = c1/(c1-alpha2*c3) !3 
    dtmp(n1-1) = (dxdy(n1-1)-alpha2*dtmp(n1-2))*c2 !3 

    c3 = gamma2*c2 !1 
    gam(n1) = c3
    c2 = c1/(c1-up2*c3) !3 
    dtmp(n1) = (dxdy(n1)-up2*dtmp(n1-1))*c2 !3

    dxdy(n1) = dtmp(n1)
    Do j=n1-1,1,-1
       dxdy(j) = dtmp(j)-gam(j+1)*dxdy(j+1)
    End Do

  End Subroutine Dbyd1_dt0

End Module Derivatives
