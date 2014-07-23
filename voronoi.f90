!*****************************************************************************
!
!! VORONOI module
!
!  Discussion:
!
!    VORONOI uses GEOMPACK to get some Voronoi diagram information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Adapted by:
!
!    Kyle Augustson
!
!    11 October 2012
!
!  Original Author:
!
!    John Burkardt
!

Module Voronoi

  Implicit None

Contains

  Subroutine Compute_Voronoi(g_x, g_y, g_degree, g_start, g_face, v_num, v_xy, ierr)

    !*****************************************************************************
    !
    !! Compute_Voronoi returns data defining the Voronoi diagram.
    !
    !  Discussion:
    !
    !    The routine first determines the Delaunay triangulation.
    !
    !    The Voronoi diagram is then determined from this information.
    !
    !    In particular, the circumcenter of each Delaunay triangle
    !    is a vertex of a Voronoi polygon.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 February 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) G_XY(2,G_NUM), the point coordinates.
    !
    !    Output, integer ( kind = 4 ) G_DEGREE(G_NUM), the degree of each 
    !    Voronoi cell.
    !
    !    Output, integer ( kind = 4 ) G_START(G_NUM), the index in G_FACE of the 
    !    first vertex at which to begin a traversal of the boundary of the 
    !    cell associated with each point.
    !
    !    Output, integer ( kind = 4 ) G_FACE(6*G_NUM), the sequence of vertices to 
    !    be used in a traversal of the boundary of the cell associated with each 
    !    point.
    !
    !    Output, integer ( kind = 4 ) V_NUM, the number of vertices of the Voronoi 
    !    diagram.
    !
    !    Output, real ( kind = 8 ) V_XY(2,V_NUM), the coordinates of the vertices
    !    of the Voronoi diagram.
    !
    Implicit None

    Integer(kind=4), Intent(InOut) :: ierr
    Integer(kind=4), parameter :: dim_num = 2
    Real(kind=8), Intent(In), Dimension(:) :: g_x, g_y
    Real(kind=8), Intent(InOut), Dimension(:,:) :: v_xy

    Integer(kind=4), Intent(InOut), Dimension(:) :: g_degree, g_face, g_start
    Integer(kind=4), Allocatable, Dimension(:,:) :: tnbr, nodtri
    Integer(kind=4) :: g_num, count, g, g_next, i, i1, i2, i3, j, k, s, sp1, s_next, v, v_num, v_next, v_old, v_save

    Real(kind=8), Allocatable, Dimension(:,:) :: g_xy, tmptri

    g_num = size(g_x)

    Allocate(g_xy(2,g_num), tmptri(dim_num,3))
    Allocate(tnbr(3,2*g_num), nodtri(3,2*g_num))

    g_xy(1,:) = g_x
    g_xy(2,:) = g_y

    !
    !  Compute the Delaunay triangulation.
    !
    call dtris2(g_num, g_xy, v_num, nodtri, tnbr, ierr)

    !
    !  Determine the degree of each cell.
    !
    g_degree(1:g_num) = 0

    do j = 1, v_num
       do i = 1, 3
          k = nodtri(i,j)
          if ( 0 < k ) then
             g_degree(k) = g_degree(k) + 1
          end if
       end do
    end do

    !
    !  Each (finite) Delaunay triangle contains a vertex of the Voronoi polygon,
    !  at the triangle's circumcenter.
    !
    do v = 1, v_num

       i1 = nodtri(1,v)
       i2 = nodtri(2,v)
       i3 = nodtri(3,v)

       tmptri(:,1) = g_xy(:,i1)
       tmptri(:,2) = g_xy(:,i2)
       tmptri(:,3) = g_xy(:,i3)

       call triangle_circumcenter_2d(tmptri, v_xy(:,v))

    end do

    !
    !  For each generator G:
    !    Determine if its region is infinite.
    !      Find a Delaunay triangle containing G.
    !      Seek another triangle containing the next node in that triangle.
    !
    count = 0
    g_start = 0

    do g = 1, g_num
       s_next = 0
       v_next = 0
       do v = 1, v_num
          do s = 1, 3
             if ( nodtri(s,v) == g ) then
                v_next = v
                s_next = s
                exit
             end if
          end do
          if ( v_next /= 0 ) then
             exit
          end if
       end do

       v_save = v_next

       do

          s_next = i4_wrap(s_next + 1, 1, 3)
          g_next = nodtri(s_next,v_next)

          if ( g_next == g ) then
             s_next = i4_wrap(s_next + 1, 1, 3)
             g_next = nodtri(s_next,v_next)
          end if

          v_old = v_next
          v_next = 0

          do v = 1, v_num

             if ( v == v_old ) then
                cycle
             end if

             do s = 1, 3

                if ( nodtri(s,v) == g ) then
                   sp1 = i4_wrap(s + 1, 1, 3)
                   if ( nodtri(sp1,v) == g_next ) then
                      v_next = v
                      s_next = sp1
                      exit
                   end if
                   sp1 = i4_wrap( s + 2, 1, 3)
                   if ( nodtri(sp1,v) == g_next ) then
                      v_next = v
                      s_next = sp1
                      exit
                   end if
                end if

             end do
             if ( v_next /= 0 ) then
                exit
             end if
          end do

          if ( v_next == v_save ) then
             exit
          end if

          if ( v_next == 0 ) then
             v_next = v_old
             exit
          end if

       end do
       !
       !  Now, starting in the current triangle, V_NEXT, cycle again,
       !  and copy the list of nodes into the array.
       !
       v_save = v_next

       count = count + 1
       g_start(g) = count
       g_face(count) = v_next

       do

          s_next = i4_wrap(s_next + 1, 1, 3)
          g_next = nodtri(s_next,v_next)

          if ( g_next == g ) then
             s_next = i4_wrap(s_next + 1, 1, 3)
             g_next = nodtri(s_next,v_next)
          end if

          v_old = v_next
          v_next = 0
          do v = 1, v_num
             if ( v == v_old ) then
                cycle
             end if

             do s = 1, 3

                if ( nodtri(s,v) == g ) then
                   sp1 = i4_wrap(s + 1, 1, 3)
                   if ( nodtri(sp1,v) == g_next ) then
                      v_next = v
                      s_next = sp1
                      exit
                   end if
                   sp1 = i4_wrap(s + 2, 1, 3)
                   if ( nodtri(sp1,v) == g_next ) then
                      v_next = v
                      s_next = sp1
                      exit
                   end if
                end if

             end do
             if ( v_next /= 0 ) then
                exit
             end if
          end do

          if ( v_next == v_save ) then
             exit
          end if

          if ( v_next == 0 ) then
             exit
          end if

          count = count + 1
          g_face(count) = v_next

       end do
    end do

    Deallocate(g_xy,tmptri,nodtri,tnbr)
  End Subroutine Compute_Voronoi

  Function diaedg(x0, y0, x1, y1, x2, y2, x3, y3)

    !*****************************************************************************
    !
    !! DIAEDG chooses a diagonal edge.
    !
    !  Discussion:
    !
    !    The routine determines whether 0--2 or 1--3 is the diagonal edge
    !    that should be chosen, based on the circumcircle criterion, where
    !    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
    !    quadrilateral in counterclockwise order.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    19 February 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
    !    coordinates of the vertices of a quadrilateral, given in
    !    counter clockwise order.
    !
    !    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
    !    +1, if diagonal edge 02 is chosen;
    !    -1, if diagonal edge 13 is chosen;
    !     0, if the four vertices are cocircular.
    !
    Implicit None

    Real(kind=8), Intent(In) :: x0, x1, x2, x3, y0, y1, y2, y3
    Real(kind=8) :: ca, cb, dx10, dx12, dx30, dx32, dy10, dy12, dy30, dy32, s, tol, tola, tolb
    Integer(kind=4) :: diaedg

    tol = 100.0D+00 * epsilon(tol)

    dx10 = x1 - x0
    dy10 = y1 - y0
    dx12 = x1 - x2
    dy12 = y1 - y2
    dx30 = x3 - x0
    dy30 = y3 - y0
    dx32 = x3 - x2
    dy32 = y3 - y2

    tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
    tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

    ca = dx10 * dx30 + dy10 * dy30
    cb = dx12 * dx32 + dy12 * dy32

    if ( tola < ca .and. tolb < cb ) then

       diaedg = -1

    else if ( ca < -tola .and. cb < -tolb ) then

       diaedg = 1

    else

       tola = max ( tola, tolb )
       s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

       if ( tola < s ) then
          diaedg = -1
       else if ( s < -tola ) then
          diaedg = 1
       else
          diaedg = 0
       end if

    end if

    return
  End Function diaedg

  Subroutine dtris2(point_num, point_xy, tri_num, tri_vert, tri_nabe, ierr)

    !*****************************************************************************
    !
    !! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
    !
    !  Discussion:
    !
    !    The routine constructs the Delaunay triangulation of a set of 2D vertices
    !    using an incremental approach and diagonal edge swaps.  Vertices are
    !    first sorted in lexicographically increasing (X,Y) order, and
    !    then are inserted one at a time from outside the convex hull.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    15 January 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of vertices.
    !
    !    Input/output, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the vertices.
    !    On output, the vertices have been sorted into dictionary order.
    !
    !    Output, integer ( kind = 4 ) TRI_NUM, the number of triangles in the 
    !    triangulation;  TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the 
    !    number of boundary vertices.
    !
    !    Output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make up 
    !    each triangle.  The elements are indices of POINT_XY.  The vertices of the 
    !    triangles are in counter clockwise order.
    !
    !    Output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
    !    list.  Positive elements are indices of TIL; negative elements are used 
    !    for links of a counter clockwise linked list of boundary edges; 
    !    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
    !    to the neighbor along edge from vertex J to J+1 (mod 3).
    !
    Implicit None

    Integer(kind=4), parameter :: dim_num = 2
    Integer(kind=4), Intent(InOut) :: ierr, tri_num
    Integer(kind=4), Intent(In) :: point_num

    Integer(kind=4), Intent(InOut), Dimension(:,:) :: tri_nabe, tri_vert
    Real(kind=8), Intent(InOut), Dimension(:,:) :: point_xy

    Real(kind=8) :: cmax, tol
    Integer(kind=4), Allocatable, Dimension(:) :: indx, stack
    Integer(kind=4) :: e, i, j, k, l, ledg, lr, ltri, m, m1, m2, n, redg, rtri, t, top


    Allocate(indx(point_num), stack(point_num))

    tol = 100.0D+00 * epsilon(tol)

    ierr = 0
    !
    !  Sort the vertices by increasing (x,y).
    !
    call r82vec_sort_heap_index_a(point_num, point_xy, indx)

    call r82vec_permute(point_num, point_xy, indx)
    !
    !  Make sure that the data points are "reasonably" distinct.
    !
    m1 = 1

    do i = 2, point_num

       m = m1
       m1 = i

       k = 0

       do j = 1, 2

          cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

          if ( tol * ( cmax + 1.0D+00 ) &
               < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
             k = j
             exit
          end if

       end do

       if ( k == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a,i8)' ) '  Fails for point number I = ', i
          write ( *, '(a,i8)' ) '  M = ', m
          write ( *, '(a,i8)' ) '  M1 = ', m1
          write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
          write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
          ierr = 224
          return
       end if

    end do
    !
    !  Starting from points M1 and M2, search for a third point M that
    !  makes a "healthy" triangle (M1,M2,M)
    !
    m1 = 1
    m2 = 2
    j = 3

    do

       if ( point_num < j ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          ierr = 225
          return
       end if

       m = j

       lr = lrline(point_xy(1,m), point_xy(2,m), point_xy(1,m1), point_xy(2,m1), point_xy(1,m2), point_xy(2,m2))

       if ( lr /= 0 ) then
          exit
       end if

       j = j + 1

    end do
    !
    !  Set up the triangle information for (M1,M2,M), and for any other
    !  triangles you created because points were collinear with M1, M2.
    !
    tri_num = j - 2

    if ( lr == -1 ) then

       tri_vert(1,1) = m1
       tri_vert(2,1) = m2
       tri_vert(3,1) = m
       tri_nabe(3,1) = -3

       do i = 2, tri_num

          m1 = m2
          m2 = i+1
          tri_vert(1,i) = m1
          tri_vert(2,i) = m2
          tri_vert(3,i) = m
          tri_nabe(1,i-1) = -3 * i
          tri_nabe(2,i-1) = i
          tri_nabe(3,i) = i - 1

       end do

       tri_nabe(1,tri_num) = -3 * tri_num - 1
       tri_nabe(2,tri_num) = -5
       ledg = 2
       ltri = tri_num

    else

       tri_vert(1,1) = m2
       tri_vert(2,1) = m1
       tri_vert(3,1) = m
       tri_nabe(1,1) = -4

       do i = 2, tri_num
          m1 = m2
          m2 = i+1
          tri_vert(1,i) = m2
          tri_vert(2,i) = m1
          tri_vert(3,i) = m
          tri_nabe(3,i-1) = i
          tri_nabe(1,i) = -3 * i - 3
          tri_nabe(2,i) = i - 1
       end do

       tri_nabe(3,tri_num) = -3 * tri_num
       tri_nabe(2,1) = -3 * tri_num - 2
       ledg = 2
       ltri = 1

    end if
    !
    !  Insert the vertices one at a time from outside the convex hull,
    !  determine visible boundary edges, and apply diagonal edge swaps until
    !  Delaunay triangulation of vertices (so far) is obtained.
    !
    top = 0

    do i = j+1, point_num

       m = i
       m1 = tri_vert(ledg,ltri)

       if ( ledg <= 2 ) then
          m2 = tri_vert(ledg+1,ltri)
       else
          m2 = tri_vert(1,ltri)
       end if

       lr = lrline(point_xy(1,m), point_xy(2,m), point_xy(1,m1), point_xy(2,m1), point_xy(1,m2), point_xy(2,m2))

       if ( 0 < lr ) then
          rtri = ltri
          redg = ledg
          ltri = 0
       else
          l = -tri_nabe(ledg,ltri)
          rtri = l / 3
          redg = mod(l,3) + 1
       end if

       call vbedg(point_xy(1,m), point_xy(2,m), point_xy, tri_vert, tri_nabe, ltri, ledg, rtri, redg)

       n = tri_num + 1
       l = -tri_nabe(ledg,ltri)

       do

          t = l / 3
          e = mod ( l, 3 ) + 1
          l = -tri_nabe(e,t)
          m2 = tri_vert(e,t)

          if ( e <= 2 ) then
             m1 = tri_vert(e+1,t)
          else
             m1 = tri_vert(1,t)
          end if

          tri_num = tri_num + 1
          tri_nabe(e,t) = tri_num
          tri_vert(1,tri_num) = m1
          tri_vert(2,tri_num) = m2
          tri_vert(3,tri_num) = m
          tri_nabe(1,tri_num) = t
          tri_nabe(2,tri_num) = tri_num - 1
          tri_nabe(3,tri_num) = tri_num + 1
          top = top + 1

          if ( point_num < top ) then
             ierr = 8
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
             write ( *, '(a)' ) '  Stack overflow.'
             return
          end if

          stack(top) = tri_num

          if ( t == rtri .and. e == redg ) then
             exit
          end if

       end do

       tri_nabe(ledg,ltri) = -3 * n - 1
       tri_nabe(2,n) = -3 * tri_num - 2
       tri_nabe(3,tri_num) = -l
       ltri = n
       ledg = 2

       call swapec(m, top, ltri, ledg, point_num, point_xy, tri_vert, tri_nabe, stack, ierr)

       if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a)' ) '  Error return from SWAPEC.'
          return
       end if

    end do
    !
    !  Now account for the sorting that we did.
    !
    do i = 1, 3
       do j = 1, tri_num
          tri_vert(i,j) = indx ( tri_vert(i,j) )
       end do
    end do

    call perm_inv(point_num, indx)

    call r82vec_permute(point_num, point_xy, indx)

    Deallocate(stack, indx)

  End Subroutine dtris2

  Function i4_modp(i, j)

    !*****************************************************************************
    !
    !! I4_MODP returns the nonnegative remainder of I4 division.
    !
    !  Discussion:
    !
    !    If
    !      NREM = I4_MODP ( I, J )
    !      NMULT = ( I - NREM ) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus, suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do, if A was positive, but if A
    !    was negative, your result would be between -360 and 0.
    !
    !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
    !
    !  Example:
    !
    !        I     J     MOD  I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number to be divided.
    !
    !    Input, integer ( kind = 4 ) J, the number that divides I.
    !
    !    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
    !    divided by J.
    !
    Implicit None

    Integer(kind=4), Intent(In) :: i, j
    Integer(kind=4) :: i4_modp

    if ( j == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'I4_MODP - Fatal error!'
       write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
       stop
    end if

    i4_modp = mod(i, j)

    if ( i4_modp < 0 ) then
       i4_modp = i4_modp + abs(j)
    end if

    return
  End Function i4_modp

  Function i4_wrap(ival, ilo, ihi)

    !*****************************************************************************
    !
    !! I4_WRAP forces an I4 to lie between given limits by wrapping.
    !
    !  Example:
    !
    !    ILO = 4, IHI = 8
    !
    !    I  I4_WRAP
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    15 July 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) IVAL, an integer value.
    !
    !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
    !
    !    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
    !
    Implicit None

    Integer(kind=4) :: i4_wrap, wide
    Integer(kind=4), Intent(In) :: ihi, ilo, ival

    wide = ihi + 1 - ilo

    if ( wide == 0 ) then
       i4_wrap = ilo
    else
       i4_wrap = ilo + i4_modp(ival - ilo, wide)
    end if

    return
  End Function i4_wrap

  Subroutine i4vec_indicator(n, a)

    !*****************************************************************************
    !
    !! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    09 November 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of elements of A.
    !
    !    Output, integer ( kind = 4 ) A(N), the array to be initialized.
    !
    Implicit None

    Integer(kind=4), Intent(In) :: n
    Integer(kind=4), Intent(InOut), Dimension(:) ::  a
    Integer(kind=4) :: i

    do i = 1, n
       a(i) = i
    end do
  End Subroutine i4vec_indicator

  Subroutine line_exp_normal_2d(p1, p2, normal)

    !*****************************************************************************
    !
    !! LINE_EXP_NORMAL_2D computes the unit normal vector to a line in 2D.
    !
    !  Discussion:
    !
    !    The explicit form of a line in 2D is:
    !
    !      ( P1, P2 ) = ( (X1,Y1), (X2,Y2) ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    01 January 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the line.
    !
    !    Output, real ( kind = 8 ) NORMAL(2), a unit normal vector to the line.
    !
    Implicit None

    Integer(kind=4), Parameter :: dim_num = 2
    Real(kind=8), Intent(In) :: p1(dim_num), p2(dim_num)
    Real(kind=8), Intent(InOut) :: normal(dim_num)
    Real(kind=8) :: norm

    norm = sqrt((p2(1)-p1(1))**2 + (p2(2)-p1(2))**2)

    if (norm == 0.0D+00) then
       normal(1:dim_num) = 0.0D+00
       return
    end if

    normal(1) =  (p2(2) - p1(2))/norm
    normal(2) = -(p2(1) - p1(1))/norm

    return
  End Subroutine line_exp_normal_2d

  Function lrline(xu, yu, xv1, yv1, xv2, yv2)

    !*****************************************************************************
    !
    !! LRLINE determines if a point is left of, right or, or on a directed line.
    !
    !  Discussion:
    !
    !    The directed line is parallel to, and at a signed distance DV from
    !    a directed base line from (XV1,YV1) to (XV2,YV2).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    14 July 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
    !    position relative to the directed line is to be determined.
    !
    !    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
    !    that determine the directed base line.
    !
    !    Output, integer ( kind = 4 ) LRLINE, the result:
    !    +1, the point is to the right of the directed line;
    !     0, the point is on the directed line;
    !    -1, the point is to the left of the directed line.
    !
    Implicit None

    Integer(kind=4) :: lrline
    Real(kind=8), Intent(In) :: xu, xv1, xv2, yu, yv1, yv2
    Real(kind=8) :: dx, dxu, dy, dyu, t, tol, tolabs

    tol = 100.0D+00 * epsilon(tol)

    dx = xv2 - xv1
    dy = yv2 - yv1
    dxu = xu - xv1
    dyu = yu - yv1

    tolabs = tol * max(abs(dx), abs(dy), abs(dxu), abs(dyu))

    t = dy * dxu - dx * dyu

    if ( tolabs < t ) then
       lrline = 1
    else if ( -tolabs <= t ) then
       lrline = 0
    else
       lrline = -1
    end if

    return
  End Function lrline

  Subroutine perm_inv(n, p)

    !*****************************************************************************
    !
    !! PERM_INV inverts a permutation "in place".
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    25 July 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of objects being permuted.
    !
    !    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard 
    !    index form.  On output, P describes the inverse permutation
    !
    Implicit None

    Integer(kind=4), Intent(In) :: n
    Integer(kind=4), Intent(InOut), Dimension(:) :: p
    Integer(kind=4) :: i, i0, i1, i2, is

    if ( n <= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PERM_INV - Fatal error!'
       write ( *, '(a,i8)' ) '  Input value of N = ', n
       stop
    end if

    is = 1

    do i = 1, n

       i1 = p(i)

       do while ( i < i1 )
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
       end do

       is = -sign ( 1, p(i) )
       p(i) = sign ( p(i), is )

    end do

    do i = 1, n

       i1 = -p(i)

       if ( 0 <= i1 ) then

          i0 = i

          do

             i2 = p(i1)
             p(i1) = i0

             if ( i2 < 0 ) then
                exit
             end if

             i0 = i1
             i1 = i2

          end do

       end if

    end do

    return
  End Subroutine perm_inv

  Subroutine r82vec_permute(n, a, p)

    !*****************************************************************************
    !
    !! R82VEC_PERMUTE permutes an R82VEC in place.
    !
    !  Discussion:
    !
    !    This routine permutes an array of real "objects", but the same
    !    logic can be used to permute an array of objects of any arithmetic
    !    type, or an array of objects of any complexity.  The only temporary
    !    storage required is enough to store a single object.  The number
    !    of data movements made is N + the number of cycles of order 2 or more,
    !    which is never more than N + N/2.
    !
    !  Example:
    !
    !    Input:
    !
    !      N = 5
    !      P = (   2,    4,    5,    1,    3 )
    !      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
    !          (11.0, 22.0, 33.0, 44.0, 55.0 )
    !
    !    Output:
    !
    !      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
    !             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    11 January 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of objects.
    !
    !    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
    !
    !    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
    !    that the I-th element of the output array should be the J-th
    !    element of the input array.  P must be a legal permutation
    !    of the integers from 1 to N, otherwise the algorithm will
    !    fail catastrophically.
    !
    Implicit None

    Integer(kind=4), Intent(In) :: n
    Integer(kind=4), Intent(InOut), Dimension(:) :: p
    Real(kind=8), Intent(InOut), Dimension(:,:) :: a
    Real(kind=8) :: a_temp(2)
    Integer(kind=4) :: iget, iput, istart

    !
    !  Search for the next element of the permutation that has not been used.
    !
    do istart = 1, n

       if ( p(istart) < 0 ) then

          cycle

       else if ( p(istart) == istart ) then

          p(istart) = - p(istart)
          cycle

       else

          a_temp(1:2) = a(1:2,istart)
          iget = istart
          !
          !  Copy the new value into the vacated entry.
          !
          do

             iput = iget
             iget = p(iget)

             p(iput) = - p(iput)

             if ( iget < 1 .or. n < iget ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
                stop
             end if

             if ( iget == istart ) then
                a(1:2,iput) = a_temp(1:2)
                exit
             end if

             a(1:2,iput) = a(1:2,iget)

          end do

       end if

    end do
    !
    !  Restore the signs of the entries.
    !
    p(1:n) = -p(1:n)

    return
  End Subroutine r82vec_permute

  Subroutine r82vec_sort_heap_index_a(n, a, indx)

    !*****************************************************************************
    !
    !! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
    !
    !  Discussion:
    !
    !    The sorting is not actually carried out.  Rather an index array is
    !    created which defines the sorting.  This array may be used to sort
    !    or index the array, or to sort or index related arrays keyed on the
    !    original array.
    !
    !    Once the index array is computed, the sorting can be carried out
    !    "implicitly:
    !
    !      A(1:2,INDX(I)), I = 1 to N is sorted,
    !
    !    or explicitly, by the call
    !
    !      call R82VEC_PERMUTE ( N, A, INDX )
    !
    !    after which A(1:2,I), I = 1 to N is sorted.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    11 January 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries in the array.
    !
    !    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
    !
    !    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
    !    I-th element of the sorted array is A(1:2,INDX(I)).
    !
    Implicit None

    Integer(kind=4), Intent(In) :: n
    Real(kind=8), Intent(In), Dimension(:,:) :: a
    integer ( kind = 4 ), Intent(InOut), Dimension(:) ::  indx(n)

    Real(kind=8) :: aval(2)
    Integer(kind=4) :: i, indxt, ir, j, l

    if ( n < 1 ) then
       return
    end if

    if ( n == 1 ) then
       indx(1) = 1
       return
    end if

    call i4vec_indicator(n, indx)

    l = n / 2 + 1
    ir = n

    do

       if ( 1 < l ) then

          l = l - 1
          indxt = indx(l)
          aval(1:2) = a(1:2,indxt)

       else

          indxt = indx(ir)
          aval(1:2) = a(1:2,indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir == 1 ) then
             indx(1) = indxt
             exit
          end if

       end if

       i = l
       j = l + l

       do while ( j <= ir )

          if ( j < ir ) then
             if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
                  ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
                  a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
                j = j + 1
             end if
          end if

          if (   aval(1) <  a(1,indx(j)) .or. &
               ( aval(1) == a(1,indx(j)) .and. &
               aval(2) <  a(2,indx(j)) ) ) then
             indx(i) = indx(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if

       end do

       indx(i) = indxt

    end do

    return
  End Subroutine r82vec_sort_heap_index_a

  Subroutine swapec(i, top, btri, bedg, point_num, point_xy, tri_vert, tri_nabe, stack, ierr)

    !*****************************************************************************
    !
    !! SWAPEC swaps diagonal edges until all triangles are Delaunay.
    !
    !  Discussion:
    !
    !    The routine swaps diagonal edges in a 2D triangulation, based on
    !    the empty circumcircle criterion, until all triangles are Delaunay,
    !    given that I is the index of the new vertex added to the triangulation.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    14 July 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the index of the new vertex.
    !
    !    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
    !    On output, TOP is zero.
    !
    !    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are 
    !    the triangle and edge indices of a boundary edge whose updated indices
    !    must be recorded.  On output, these may be updated because of swaps.
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
    !
    !    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
    !    of the points.
    !
    !    Input/output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle 
    !    incidence list.  May be updated on output because of swaps.
    !
    !    Input/output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle 
    !    neighbor list; negative values are used for links of the counter-clockwise 
    !    linked list of boundary edges;  May be updated on output because of swaps.
    !
    !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    !
    !    Workspace, integer ( kind = 4 ) STACK(MAXST); on input, entries 1 through 
    !    TOP contain the indices of initial triangles (involving vertex I)
    !    put in stack; the edges opposite I should be in interior;  entries
    !    TOP+1 through MAXST are used as a stack.
    !
    !    Output, integer IERR is set to 8 for abnormal return.
    !
    Implicit None

    Integer(kind=4), Intent(In) :: i, point_num
    Real(kind=8), Intent(In), Dimension(:,:) :: point_xy
    Integer(kind=4), Intent(InOut) :: ierr, top, btri, bedg, tri_vert(:,:), tri_nabe(:,:), stack(:)

    Integer(kind=4) :: a, b, c, e, ee, em1, ep1, f, fm1, fp1, l, r, s, swap, t, tt, u
    Real(kind=8) :: x, y

    !
    !  Determine whether triangles in stack are Delaunay, and swap
    !  diagonal edge of convex quadrilateral if not.
    !
    x = point_xy(1,i)
    y = point_xy(2,i)

    do

       if ( top <= 0 ) then
          exit
       end if

       t = stack(top)
       top = top - 1

       if ( tri_vert(1,t) == i ) then
          e = 2
          b = tri_vert(3,t)
       else if ( tri_vert(2,t) == i ) then
          e = 3
          b = tri_vert(1,t)
       else
          e = 1
          b = tri_vert(2,t)
       end if

       a = tri_vert(e,t)
       u = tri_nabe(e,t)

       if ( tri_nabe(1,u) == t ) then
          f = 1
          c = tri_vert(3,u)
       else if ( tri_nabe(2,u) == t ) then
          f = 2
          c = tri_vert(1,u)
       else
          f = 3
          c = tri_vert(2,u)
       end if

       swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
            point_xy(2,c), point_xy(1,b), point_xy(2,b) )

       if ( swap == 1 ) then

          em1 = i4_wrap ( e - 1, 1, 3 )
          ep1 = i4_wrap ( e + 1, 1, 3 )
          fm1 = i4_wrap ( f - 1, 1, 3 )
          fp1 = i4_wrap ( f + 1, 1, 3 )

          tri_vert(ep1,t) = c
          tri_vert(fp1,u) = i
          r = tri_nabe(ep1,t)
          s = tri_nabe(fp1,u)
          tri_nabe(ep1,t) = u
          tri_nabe(fp1,u) = t
          tri_nabe(e,t) = s
          tri_nabe(f,u) = r

          if ( 0 < tri_nabe(fm1,u) ) then
             top = top + 1
             stack(top) = u
          end if

          if ( 0 < s ) then

             if ( tri_nabe(1,s) == u ) then
                tri_nabe(1,s) = t
             else if ( tri_nabe(2,s) == u ) then
                tri_nabe(2,s) = t
             else
                tri_nabe(3,s) = t
             end if

             top = top + 1

             if ( point_num < top ) then
                ierr = 8
                return
             end if

             stack(top) = t

          else

             if ( u == btri .and. fp1 == bedg ) then
                btri = t
                bedg = e
             end if

             l = - ( 3 * t + e - 1 )
             tt = t
             ee = em1

             do while ( 0 < tri_nabe(ee,tt) )

                tt = tri_nabe(ee,tt)

                if ( tri_vert(1,tt) == a ) then
                   ee = 3
                else if ( tri_vert(2,tt) == a ) then
                   ee = 1
                else
                   ee = 2
                end if

             end do

             tri_nabe(ee,tt) = l

          end if

          if ( 0 < r ) then

             if ( tri_nabe(1,r) == t ) then
                tri_nabe(1,r) = u
             else if ( tri_nabe(2,r) == t ) then
                tri_nabe(2,r) = u
             else
                tri_nabe(3,r) = u
             end if

          else

             if ( t == btri .and. ep1 == bedg ) then
                btri = u
                bedg = f
             end if

             l = - ( 3 * u + f - 1 )
             tt = u
             ee = fm1

             do while ( 0 < tri_nabe(ee,tt) )

                tt = tri_nabe(ee,tt)

                if ( tri_vert(1,tt) == b ) then
                   ee = 3
                else if ( tri_vert(2,tt) == b ) then
                   ee = 1
                else
                   ee = 2
                end if

             end do

             tri_nabe(ee,tt) = l

          end if

       end if

    end do

    return
  End Subroutine swapec

  Subroutine triangle_circumcenter_2d(t, center)

    !*****************************************************************************
    !
    !! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
    !
    !  Discussion:
    !
    !    The circumcenter of a triangle is the center of the circumcircle, the
    !    circle that passes through the three vertices of the triangle.
    !
    !    The circumcircle contains the triangle, but it is not necessarily the
    !    smallest triangle to do so.
    !
    !    If all angles of the triangle are no greater than 90 degrees, then
    !    the center of the circumscribed circle will lie inside the triangle.
    !    Otherwise, the center will lie outside the triangle.
    !
    !    The circumcenter is the intersection of the perpendicular bisectors
    !    of the sides of the triangle.
    !
    !    In geometry, the circumcenter of a triangle is often symbolized by "O".
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    18 February 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) CENTER(2), the circumcenter of the triangle.
    !
    Implicit None

    Integer(kind=4), Parameter :: dim_num = 2
    Real(kind=8), Intent(In), Dimension(:,:) :: t
    Real(kind=8), Intent(InOut), Dimension(:) :: center
    Real(kind=8) :: asq, bot, csq, top(dim_num)

    asq = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
    csq = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2

    top(1) =    ( t(2,2) - t(2,1) ) * csq - ( t(2,3) - t(2,1) ) * asq
    top(2) =  - ( t(1,2) - t(1,1) ) * csq + ( t(1,3) - t(1,1) ) * asq

    bot  =  ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) &
         - ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )

    center(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / bot

    return
  End Subroutine triangle_circumcenter_2d

  Subroutine vbedg(x, y, point_xy, tri_vert, tri_nabe, ltri, ledg, rtri, redg)

    !*****************************************************************************
    !
    !! VBEDG determines which boundary edges are visible to a point.
    !
    !  Discussion:
    !
    !    The point (X,Y) is assumed to be outside the convex hull of the
    !    region covered by the 2D triangulation.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    25 August 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, Y, a point outside the convex hull
    !    of the current triangulation.
    !
    !    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the vertices.
    !
    !    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle 
    !    incidence list.
    !
    !    Input, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
    !    list; negative values are used for links of a counter clockwise linked list
    !    of boundary edges;
    !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    !
    !    Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these 
    !    values are assumed to be already computed and are not changed, else they 
    !    are updated.  On output, LTRI is the index of boundary triangle to the left 
    !    of the leftmost boundary triangle visible from (X,Y), and LEDG is the 
    !    boundary edge of triangle LTRI to the left of the leftmost boundary edge 
    !    visible from (X,Y).  1 <= LEDG <= 3.
    !
    !    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the 
    !    boundary triangle to begin the search at.  On output, the index of the 
    !    rightmost boundary triangle visible from (X,Y).
    !
    !    Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that is 
    !    visible from (X,Y).  1 <= REDG <= 3.
    !
    Implicit None

    Real(kind=8), Intent(In) :: x, y, point_xy(:,:)
    Integer(kind=4), Intent(In) :: tri_vert(:,:), tri_nabe(:,:)
    Integer(kind=4), Intent(InOut) :: ltri, ledg, rtri, redg
    Integer(kind=4) :: a, b, e, l, lr, t
    Logical :: ldone
    !
    !  Find the rightmost visible boundary edge using links, then possibly
    !  leftmost visible boundary edge using triangle neighbor information.
    !
    if ( ltri == 0 ) then
       ldone = .false.
       ltri = rtri
       ledg = redg
    else
       ldone = .true.
    end if

    do

       l = -tri_nabe(redg,rtri)
       t = l / 3
       e = mod ( l, 3 ) + 1
       a = tri_vert(e,t)

       if ( e <= 2 ) then
          b = tri_vert(e+1,t)
       else
          b = tri_vert(1,t)
       end if

       lr = lrline(x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), point_xy(2,b))

       if ( lr <= 0 ) then
          exit
       end if

       rtri = t
       redg = e

    end do

    if ( ldone ) then
       return
    end if

    t = ltri
    e = ledg

    do

       b = tri_vert(e,t)
       e = i4_wrap(e-1, 1, 3)

       do while ( 0 < tri_nabe(e,t) )

          t = tri_nabe(e,t)

          if ( tri_vert(1,t) == b ) then
             e = 3
          else if ( tri_vert(2,t) == b ) then
             e = 1
          else
             e = 2
          end if

       end do

       a = tri_vert(e,t)

       lr = lrline(x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), point_xy(2,b))

       if ( lr <= 0 ) then
          exit
       end if

    end do

    ltri = t
    ledg = e

    return
  End Subroutine vbedg

End Module Voronoi
