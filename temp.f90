    module kinds
      integer, parameter::I4 = selected_int_kind(9)
      integer, parameter::R4 = kind(1.0)
      integer, parameter::R8 = kind(1.d0)
!
! Maximum number of particles (can be changed by the user)
      integer(I4), parameter::NMAX = 100

    end module kinds

!==============================================================================
    module constants
      use kinds
      real(R8), parameter::ZERO  = 0.0_R8
      real(R8), parameter::ONE   = 1.0_R8
      real(R8), parameter::TWO   = 2.0_R8
      real(R8), parameter::THREE = 3.0_R8
      real(R8), parameter::FOUR  = 4.0_R8
      real(R8), parameter::HALF  = 0.5_R8
      real(R8), parameter::THIRD = 1.0_R8 / 3.0_R8
      real(R8), parameter::SIXTH = 1.0_R8 / 6.0_R8
!
      real(R8), parameter::PI     = 3.141592653589793_R8
      real(R8), parameter::TWOPI  = PI * TWO
      real(R8), parameter::PIBY2  = PI * HALF
      real(R8), parameter::ROOT2  = 1.41421356237309_R8
      real(R8), parameter::ROOT10 = 3.16227766016838_R8
!
! Conversion factor from degrees to radians
      real(R8), parameter::DR    = PI / 180.0_R8
!
! A very large positive number
      real(R8), parameter::BIG_NUMBER   = huge(ONE)
!
! A very small positive number (bigger than machine precision)
      real(R8), parameter::SMALL_NUMBER = 1.0e-30_R8
!
! Universal constants
      real(R8), parameter::G = 6.67428e-8_R8
!
! Solar System related constants
      real(R8), parameter::AU    = 1.49597870700e13_R8
      real(R8), parameter::GMSUN = 1.32712440041e26_R8
      real(R8), parameter::MSUN  = GMSUN / G
      real(R8), parameter::MEARTH   = MSUN / 332946.0487_R8
      real(R8), parameter::MJUPITER = MSUN / 1047.348644_R8
      real(R8), parameter::DAY  = 86400.0_R8
      real(R8), parameter::YEAR = DAY * 365.25_R8
!
! Unit conversion factors for energy and angular momentum
      real(R8), parameter::CONVERT_EN  = MSUN * AU * AU / (DAY * DAY)
      real(R8), parameter::CONVERT_ANG = MSUN * AU * AU / DAY
!
! Maximum number of unique particle indices (slightly less than 224 ^ 3)
      real(R8), parameter::INDEX_MAX = 11239423.99_R8
    end module constants

    module globals
      use kinds
      integer(I4)::npair, ipair(NMAX), jpair(NMAX)
      real(R8)::denergy
      integer(I4)::nfrag
      real(R8)::mfrag_min
    end module globals


    module functions

    interface
      subroutine fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
        rce_hill,rad,rcrit,status,index,name)
      use kinds
      integer(I4), intent(in)::itarg,iproj
      real(R8),    intent(inout)::m1,m2
      integer(I4), intent(inout)::n,nbig,index(:)
      real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
      real(R8), intent(inout)::rho(:),rce_hill(:),rad(:)
      real(R8), intent(inout)::rcrit(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine fragment_bodies
    end interface

    interface
      subroutine other_fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
        rce_hill,rad,rcrit,status,index,name)
      use kinds
      integer(I4), intent(in)::itarg,iproj
      real(R8),    intent(inout)::m1,m2
      integer(I4), intent(inout)::n,nbig,index(:)
      real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
      real(R8), intent(inout)::rho(:),rce_hill(:),rad(:)
      real(R8), intent(inout)::rcrit(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine other_fragment_bodies
    end interface

    interface
      function cross_product (a,b)
      use kinds
      real(R8), intent(in)::a(3),b(3)
      real(R8)::cross_product(3)
      end function cross_product
    end interface

    interface
      subroutine calc_relative_coords (m,x,v,i,j,xrel,vrel,xcom,vcom)
      use kinds
      integer(I4), intent(in)::i,j
      real(R8), intent(in)::m(:),x(:,:),v(:,:)
      real(R8), intent(out)::xrel(:),vrel(:),xcom(:),vcom(:)
      end subroutine calc_relative_coords
    end interface


    end module functions

!==============================================================================
! Calculates the cross product of two vectors.
!
    function cross_product (a,b)
    use constants
    implicit none
!
    real(R8), intent(in)::a(3),b(3)
    real(R8)::cross_product(3)
!------------------------------------------------------------------------------
    cross_product(1) = a(2) * b(3)  -  a(3) * b(2)
    cross_product(2) = a(3) * b(1)  -  a(1) * b(3)
    cross_product(3) = a(1) * b(2)  -  a(2) * b(1)
!
    end function cross_product

!==============================================================================
! Resolves a collision between target ITARG and projectile IPROJ
! that results in the formation of a single large remnant of mass M1 
! plus some equal-mass fragments.
! Fragments are assumed to be ejected at slightly above escape velocity in
! uniform directions within a plane centred on the centre of mass.
!
    subroutine fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
      rce_hill,rad,rcrit,status,index,name)
!
    use constants;    use globals;   
    use functions, only: cross_product, calc_relative_coords
    implicit none
    integer(I4), intent(in)::itarg,iproj
    real(R8),    intent(inout)::m1,m2
    integer(I4), intent(inout)::n,nbig,index(:)
    real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
    real(R8), intent(inout)::rho(:),rce_hill(:),rad(:)
    real(R8), intent(inout)::rcrit(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::ifrag,jfrag,nnew,nold,ntheta,final_value
    real(R8)::msum,rsum,m_ejecta,mfrag,rej,v2esc,v2ej,vej,theta
    real(R8)::xcom(3),vcom(3),xrel(3),vrel(3),en0,en1,dx(3),dv(3),mredu
    real(R8)::mxsum(3),mvsum(3),xoff(3),voff(3),l(3),p(3),z(3)
!------------------------------------------------------------------------------
! Save initial number of bodies
    nold = n
!
    msum = m(itarg)    +  m(iproj)
    rsum = rad(itarg)  +  rad(iproj)
!
! Calculate centre of mass of the target and projectile
    call calc_relative_coords (m,x,v,itarg,iproj,xrel,vrel,xcom,vcom)
!
! Reduced mass of target and projectile
    mredu = m(itarg) * m(iproj) / msum
!
! Initial energy of the target and projectile (center of mass frame)
    en0 = HALF * mredu * dot_product(vrel, vrel) &
        - G * m(itarg) * m(iproj) / sqrt(dot_product(xrel, xrel))
!
! Calculate number of fragments NNEW and individual fragment mass MFRAG
    m_ejecta = msum  -  (m1  +  m2)
    nnew  = int(m_ejecta / mfrag_min)
    write(*,*) "nnew: ",nnew
!
    if (nnew == 0) then
      mfrag = ZERO
      if (m2 > ZERO) then
        m2 = msum  -  m1
      else
        m1 = msum
      end if
    else
      mfrag = m_ejecta / dble(nnew)
    end if
!
! Number NTHETA of particles apart from the largest remnant
    ntheta = nnew
    if (m2 > ZERO) ntheta = nnew  +  1
!
! Distance REJ and velocity VEJ of fragments with respect to centre of mass
    rej = rsum * FOUR
    v2esc = TWO * G * msum / rsum
    v2ej = 1.1_R8 * v2esc  -  TWO * G * msum * (ONE / rsum  -  ONE / rej)
    vej  = sqrt(v2ej)
!
! Set up unit vectors in the collision reference frame
!  l = parallel to target velocity
!  p = perpendicular to target velocity in collision plane
!  z = normal to the collision plane
    l(:) = vrel(:)
    z(:) = cross_product (vrel, xrel)
    p(:) = cross_product (z, vrel)
!
    l = l / sqrt(dot_product(l, l))
    z = z / sqrt(dot_product(z, z))
    p = p / sqrt(dot_product(p, p))
!
    write (*,'(3a,es11.4)') '   Remnant:  ',name(itarg), '  m=', m1 / MSUN
    if (m2.gt.ZERO) then
       write (*,'(3a,es11.4)') '   Remnant:  ',name(iproj), '  m=', m2 / MSUN
    endif
!
    mxsum = ZERO;      mvsum = ZERO
!
! Create each new fragment

! First, create the first fragment with the same name and space in the array as m2/iproj
    if (m2.eq.ZERO.and.nnew.gt.0) then

       m(iproj)        = mfrag
       s(:,iproj)      = ZERO
!      rho(iproj)      = rho(iproj)  ! The data are already there
!      rce_hill(iproj) = rce_hill(iproj)  ! The data are already there
!      rcrit(iproj)    = rcrit(iproj)  ! The data are already there
       rad(iproj)      = (THREE * m(iproj) / (FOUR * PI * rho(iproj)))**THIRD
       !index(iproj)    = 0   !Don't do this, body already exists!
!      status(iproj)   = 'big  '  ! The data are already there

       theta = TWOPI * ONE / dble(ntheta)
       x(:,iproj) = xcom(:)  -  rej * cos(theta) * l(:)  - rej * sin(theta) * p(:)
       write(*,*) x(:,iproj)
       v(:,iproj) = vcom(:)  -  vej * cos(theta) * l(:)  - vej * sin(theta) * p(:)
       mxsum(:) = mxsum(:)  +  m(iproj) * x(:,iproj)
       mvsum(:) = mvsum(:)  +  m(iproj) * v(:,iproj)
       write (*,'(3a,es11.4)') '   Remnant:  ',name(iproj), '  m=', m(iproj) / MSUN
    endif

!To give the right value to loop to in what follows
    final_value = nold+nnew
    if (m2.eq.ZERO) then
          final_value = nold+nnew - 1
    endif

! Then add in the rest of the fragments at the end of the array, if more than one fragment

    if ( (m2.gt.ZERO.and.nnew.gt.0).or.(nnew.gt.1)) then

       do ifrag = nold + 1, final_value
          n = n  +  1
          nbig = nbig  +  1
!
          ! Set up mass, spin, density, radius etc
          m(ifrag)        = mfrag
          s(:,ifrag)      = ZERO
          rho(ifrag)      = rho(iproj)
          rce_hill(ifrag) = rce_hill(iproj)
          rcrit(ifrag)    = rcrit(iproj)
          rad(ifrag)      = (THREE * m(n) / (FOUR * PI * rho(n)))**THIRD
          index(ifrag)    = 0
          status(ifrag)   = 'big  '
!
          if (m2.gt.ZERO) then
             theta = TWOPI * dble(ifrag-nold) / dble(ntheta)
          else
             theta = TWOPI * dble(ifrag+1-nold) / dble(ntheta)
          endif
          x(:,ifrag) = xcom(:)  -  rej * cos(theta) * l(:)  - rej * sin(theta) * p(:)
          v(:,ifrag) = vcom(:)  -  vej * cos(theta) * l(:)  - vej * sin(theta) * p(:)
          mxsum(:) = mxsum(:)  +  m(ifrag) * x(:,ifrag)
          mvsum(:) = mvsum(:)  +  m(ifrag) * v(:,ifrag)
!
! Choose a unique name for the fragment
          nfrag = nfrag  +  1
          name(ifrag) = 'FRAG0000'
!
          if (nfrag < 10) then
             write (name(n)(8:8),'(i1)') nfrag
          else if (nfrag < 100) then
             write (name(n)(7:8),'(i2)') nfrag
          else if (nfrag < 1000) then
             write (name(n)(6:8),'(i3)') nfrag
          else if (nfrag < 10000) then
             write (name(n)(5:8),'(i4)') nfrag
          else if (nfrag < 100000) then
             write (name(n)(4:8),'(i5)') nfrag
          else if (nfrag < 1000000) then
             write (name(n)(3:8),'(i6)') nfrag
          else if (nfrag < 10000000) then
             write (name(n)(2:8),'(i7)') nfrag
          else
             write (*,*) "There are far too many fragments in the sim, over 1e7!"
          end if
          write (*,'(3a,es11.4)') '   Fragment: ',name(n), '  m=', m(n) / MSUN
       end do
    endif
!
! Masses and velocities of the remnant(s)
    m(itarg) = m1
    x(:,itarg) = xcom(:)
    v(:,itarg) = vcom(:)
    s(:,itarg) = s(:,itarg)  +  mredu * cross_product(xrel,vrel)
    if (m2 == ZERO) s(:,itarg) = s(:,itarg)  +  s(:,iproj)
!
    if (m2 > ZERO ) then
      m(iproj) = m2
      x(:,iproj) = xcom(:)  -  rej * l(:)
      v(:,iproj) = vcom(:)  -  vej * l(:)
!    else  ! This was removed when I decided to have iproj become one of the fragments, to decrease the number of fragment names
!      m(iproj) = ZERO;             s(:,iproj) = ZERO
!      x(:,iproj) = -x(:,iproj)
!      v(:,iproj) = -v(:,iproj)
!      status(iproj) = 'dead '
    elseif (nnew == 0) then
       m(iproj) = ZERO;             s(:,iproj) = ZERO
       x(:,iproj) = -x(:,iproj)
       v(:,iproj) = -v(:,iproj)
       status(iproj) = 'dead '
    endif
  
!
    if (m2.gt.ZERO) then
       mxsum(:) = mxsum(:)  +  m(itarg) * x(:,itarg)  +  m(iproj) * x(:,iproj)
       mvsum(:) = mvsum(:)  +  m(itarg) * v(:,itarg)  +  m(iproj) * v(:,iproj)
    else
       mxsum(:) = mxsum(:)  +  m(itarg) * x(:,itarg)  !+  m(iproj) * x(:,iproj)
       mvsum(:) = mvsum(:)  +  m(itarg) * v(:,itarg)  !+  m(iproj) * v(:,iproj)
    endif
!
! Adjust the velocities to conserve momentum
    xoff(:) = xcom(:)  -  mxsum(:) / msum
    voff(:) = vcom(:)  -  mvsum(:) / msum
!
    x(:,itarg) = x(:,itarg)  +  xoff(:)
    v(:,itarg) = v(:,itarg)  +  voff(:)
    x(:,iproj) = x(:,iproj)  +  xoff(:) !This needs to happen whether m2 is ZERO or not
    v(:,iproj) = v(:,iproj)  +  voff(:) !This too
!
    do ifrag = nold + 1, final_value
       write(*,*) "inside loop"
       write(*,*) ifrag
      x(:,ifrag) = x(:,ifrag)  +  xoff(:)
      v(:,ifrag) = v(:,ifrag)  +  voff(:)
    end do
!
! Final energy of remnants and fragments
    dv(:) = v(:,itarg)  -  vcom(:)
    en1 = HALF * m(itarg) * dot_product(dv, dv)
    dv(:) = v(:,iproj)  -  vcom(:)
    en1 = en1  +  HALF * m(iproj) * dot_product(dv, dv)
!
    dx(:) = x(:,itarg)  -  x(:,iproj)
    en1 = en1  -  G * m(itarg) * m(iproj) / sqrt(dot_product(dx, dx))
!
    do ifrag = nold + 1, final_value
      dv(:) = v(:,ifrag)  -  vcom(:)
      en1 = en1  +  HALF * mfrag * dot_product(dv, dv)
!
      dx(:) = x(:,ifrag)  -  x(:,itarg)
      en1 = en1  -  G * m(itarg) * m(ifrag) / sqrt(dot_product(dx, dx))
      dx(:) = x(:,ifrag)  -  x(:,iproj)
      en1 = en1  -  G * m(iproj) * m(ifrag) / sqrt(dot_product(dx, dx))
!
      do jfrag = ifrag + 1, final_value
        dx(:) = x(:,ifrag)  -  x(:,jfrag)
        en1 = en1  -  G * mfrag * mfrag / sqrt(dot_product(dx, dx))
      end do
    end do
!
! Calculate energy loss due to the collision
    denergy = denergy  +  (en0  -  en1)
    write(*,*) "denergy, ",  denergy

! Add entries in the list of particle pairs within critical distance
    if ( (m2.gt.ZERO.and.nnew.gt.0).or.(nnew.gt.1)) then
       !write (*,*)
!
       do ifrag = nold + 1, final_value 
          npair = npair  +  1
          ipair(npair) = itarg !pair up the new fragments with the target
          jpair(npair) = ifrag
!          write (*,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
!
          !if (m2 /= ZERO) then   !This conditional removed with the change to m2 stuff that I've done, now this should always run
             npair = npair  +  1
             ipair(npair) = iproj
             jpair(npair) = ifrag
!             write (*,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
          !end if
          !
          do jfrag = nold + 1, ifrag - 1
             npair = npair  +  1
             ipair(npair) = ifrag
             jpair(npair) = jfrag
!             write (*,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
          end do
       end do
    endif
!
    end subroutine fragment_bodies



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============================================================================
! Resolves a collision between target ITARG and projectile IPROJ
! that results in the formation of a single large remnant of mass M1 
! plus some equal-mass fragments.
! Fragments are assumed to be ejected at slightly above escape velocity in
! uniform directions within a plane centred on the centre of mass.
!
    subroutine other_fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
      rce_hill,rad,rcrit,status,index,name)
!
    use constants;    use globals;   
    use functions, only: cross_product, calc_relative_coords
    implicit none
    integer(I4), intent(in)::itarg,iproj
    real(R8),    intent(inout)::m1,m2
    integer(I4), intent(inout)::n,nbig,index(:)
    real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
    real(R8), intent(inout)::rho(:),rce_hill(:),rad(:)
    real(R8), intent(inout)::rcrit(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::ifrag,jfrag,nnew,nold,ntheta,final_value
    real(R8)::msum,rsum,m_ejecta,mfrag,rej,v2esc,v2ej,vej,theta
    real(R8)::xcom(3),vcom(3),xrel(3),vrel(3),en0,en1,dx(3),dv(3),mredu
    real(R8)::mxsum(3),mvsum(3),xoff(3),voff(3),l(3),p(3),z(3)
!------------------------------------------------------------------------------
! Save initial number of bodies
    nold = n
!
    msum = m(itarg)    +  m(iproj)
    rsum = rad(itarg)  +  rad(iproj)
!
! Calculate centre of mass of the target and projectile
    call calc_relative_coords (m,x,v,itarg,iproj,xrel,vrel,xcom,vcom)
!
! Reduced mass of target and projectile
    mredu = m(itarg) * m(iproj) / msum
!
! Initial energy of the target and projectile (center of mass frame)
    en0 = HALF * mredu * dot_product(vrel, vrel) &
        - G * m(itarg) * m(iproj) / sqrt(dot_product(xrel, xrel))
!
! Calculate number of fragments NNEW and individual fragment mass MFRAG
    m_ejecta = msum  -  (m1  +  m2)
    nnew  = int(m_ejecta / mfrag_min)
    write(*,*) "nnew: ", nnew
!
    if (nnew == 0) then
      mfrag = ZERO
      if (m2 > ZERO) then
        m2 = msum  -  m1
      else
        m1 = msum
      end if
    else
      mfrag = m_ejecta / dble(nnew)
    end if
!
! Number NTHETA of particles apart from the largest remnant
    ntheta = nnew
    if (m2 > ZERO) ntheta = nnew  +  1
!
! Distance REJ and velocity VEJ of fragments with respect to centre of mass
    rej = rsum * FOUR
    v2esc = TWO * G * msum / rsum
    v2ej = 1.1_R8 * v2esc  -  TWO * G * msum * (ONE / rsum  -  ONE / rej)
    vej  = sqrt(v2ej)
!
! Set up unit vectors in the collision reference frame
!  l = parallel to target velocity
!  p = perpendicular to target velocity in collision plane
!  z = normal to the collision plane
    l(:) = vrel(:)
    z(:) = cross_product (vrel, xrel)
    p(:) = cross_product (z, vrel)
!
    l = l / sqrt(dot_product(l, l))
    z = z / sqrt(dot_product(z, z))
    p = p / sqrt(dot_product(p, p))
!
    write (*,'(3a,es11.4)') '   Remnant:  ',name(itarg), '  m=', m1 / MSUN
    write (*,'(3a,es11.4)') '   Remnant:  ',name(iproj), '  m=', m2 / MSUN
!
    mxsum = ZERO;      mvsum = ZERO
!
! Create each new fragment
    do ifrag = nold + 1, nold + nnew
      n = n  +  1
      nbig = nbig  +  1
!
! Set up mass, spin, density, radius etc
      m(ifrag)        = mfrag
      s(:,ifrag)      = ZERO
      rho(ifrag)      = rho(iproj)
      rce_hill(ifrag) = rce_hill(iproj)
      rcrit(ifrag)    = rcrit(iproj)
      rad(ifrag)      = (THREE * m(n) / (FOUR * PI * rho(n)))**THIRD
      index(ifrag)    = 0
      status(ifrag)   = 'big  '
!
      theta = TWOPI * dble(ifrag-nold) / dble(ntheta)
      x(:,ifrag) = xcom(:)  -  rej * cos(theta) * l(:)  - rej * sin(theta) * p(:)
      v(:,ifrag) = vcom(:)  -  vej * cos(theta) * l(:)  - vej * sin(theta) * p(:)
      mxsum(:) = mxsum(:)  +  m(ifrag) * x(:,ifrag)
      mvsum(:) = mvsum(:)  +  m(ifrag) * v(:,ifrag)
!
! Choose a unique name for the fragment
      nfrag = nfrag  +  1
      name(ifrag) = 'FRAG0000'
!
      if (nfrag < 10) then
        write (name(n)(8:8),'(i1)') nfrag
      else if (nfrag < 100) then
        write (name(n)(7:8),'(i2)') nfrag
      else if (nfrag < 1000) then
        write (name(n)(6:8),'(i3)') nfrag
      else if (nfrag < 10000) then
        write (name(n)(5:8),'(i4)') nfrag
      else if (nfrag < 100000) then
        write (name(n)(4:8),'(i5)') nfrag
      else if (nfrag < 1000000) then
        write (name(n)(3:8),'(i6)') nfrag
      else if (nfrag < 10000000) then
        write (name(n)(2:8),'(i7)') nfrag
      else
         write (*,*) "There are far too many fragments in the sim, over 1e7!"
      end if
      write (*,'(3a,es11.4)') '   Fragment: ',name(n), '  m=', m(n) / MSUN
    end do

!
! Masses and velocities of the remnant(s)
    m(itarg) = m1
    x(:,itarg) = xcom(:)
    v(:,itarg) = vcom(:)
    s(:,itarg) = s(:,itarg)  +  mredu * cross_product(xrel,vrel)
    if (m2 == ZERO) s(:,itarg) = s(:,itarg)  +  s(:,iproj)
!
    if (m2 > ZERO) then
      m(iproj) = m2
      x(:,iproj) = xcom(:)  -  rej * l(:)
      v(:,iproj) = vcom(:)  -  vej * l(:)
    else
      m(iproj) = ZERO;             s(:,iproj) = ZERO
      x(:,iproj) = -x(:,iproj)
      v(:,iproj) = -v(:,iproj)
      status(iproj) = 'dead '
    end if    
!
    mxsum(:) = mxsum(:)  +  m(itarg) * x(:,itarg)  +  m(iproj) * x(:,iproj)
    mvsum(:) = mvsum(:)  +  m(itarg) * v(:,itarg)  +  m(iproj) * v(:,iproj)
!
! Adjust the velocities to conserve momentum
    xoff(:) = xcom(:)  -  mxsum(:) / msum
    voff(:) = vcom(:)  -  mvsum(:) / msum
!
    x(:,itarg) = x(:,itarg)  +  xoff(:)
    x(:,iproj) = x(:,iproj)  +  xoff(:)
    v(:,itarg) = v(:,itarg)  +  voff(:)
    v(:,iproj) = v(:,iproj)  +  voff(:)
!
    do ifrag = nold + 1, nold + nnew
      x(:,ifrag) = x(:,ifrag)  +  xoff(:)
      v(:,ifrag) = v(:,ifrag)  +  voff(:)
    end do
!
! Final energy of remnants and fragments
    dv(:) = v(:,itarg)  -  vcom(:)
    en1 = HALF * m(itarg) * dot_product(dv, dv)
    dv(:) = v(:,iproj)  -  vcom(:)
    en1 = en1  +  HALF * m(iproj) * dot_product(dv, dv)
!
    dx(:) = x(:,itarg)  -  x(:,iproj)
    en1 = en1  -  G * m(itarg) * m(iproj) / sqrt(dot_product(dx, dx))
!
    do ifrag = nold + 1, nold + nnew
      dv(:) = v(:,ifrag)  -  vcom(:)
      en1 = en1  +  HALF * mfrag * dot_product(dv, dv)
!
      dx(:) = x(:,ifrag)  -  x(:,itarg)
      en1 = en1  -  G * m(itarg) * m(ifrag) / sqrt(dot_product(dx, dx))
      dx(:) = x(:,ifrag)  -  x(:,iproj)
      en1 = en1  -  G * m(iproj) * m(ifrag) / sqrt(dot_product(dx, dx))
!
      do jfrag = ifrag + 1, nold + nnew
        dx(:) = x(:,ifrag)  -  x(:,jfrag)
        en1 = en1  -  G * mfrag * mfrag / sqrt(dot_product(dx, dx))
      end do
    end do
!
! Calculate energy loss due to the collision
    denergy = denergy  +  (en0  -  en1)
    write(*,*) "denergy: ",denergy

!
! Add entries in the list of particle pairs within critical distance
!    if (nnew > 0) write (*,*)
!
    do ifrag = nold + 1, nold + nnew
      npair = npair  +  1
      ipair(npair) = itarg
      jpair(npair) = ifrag
!      write (*,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
!
      if (m2 /= ZERO) then
        npair = npair  +  1
        ipair(npair) = iproj
        jpair(npair) = ifrag
!        write (*,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
      end if
!
      do jfrag = nold + 1, ifrag - 1
        npair = npair  +  1
        ipair(npair) = ifrag
        jpair(npair) = jfrag
!        write (*,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
      end do
    end do
!
    end subroutine other_fragment_bodies



!==============================================================================
! Calculates the centre of mass, and relative coordinates and velocities, of
! particles I and J given their absolute coordinates and velocities.
!
    subroutine calc_relative_coords (m,x,v,i,j,xrel,vrel,xcom,vcom)
!
    use constants;    use globals
    implicit none
    integer(I4), intent(in)::i,j
    real(R8), intent(in)::m(:),x(:,:),v(:,:)
    real(R8), intent(out)::xrel(:),vrel(:),xcom(:),vcom(:)
!
    real(R8)::msum
!------------------------------------------------------------------------------
    msum = m(i)  +  m(j)
    xcom(:) = (m(i) * x(:,i)  +  m(j) * x(:,j)) / msum
    vcom(:) = (m(i) * v(:,i)  +  m(j) * v(:,j)) / msum
!
    xrel(:) = x(:,i)  -  x(:,j)
    vrel(:) = v(:,i)  -  v(:,j)
!
    end subroutine calc_relative_coords



    program  hey_there

      use constants;    use globals
      use functions

      implicit none

      integer(I4):: itarg,iproj,n,nbig,index(NMAX)
      real(R8):: x0(3,NMAX),v0(3,NMAX),m(NMAX),s(3,NMAX),ngf(3,NMAX)
      real(R8):: x(3,NMAX),v(3,NMAX)
      real(R8):: x_save(3,NMAX),v_save(3,NMAX)
      real(R8):: rho(NMAX),rce_hill(NMAX),rad(NMAX),rcrit(NMAX)
      character(8)::name(NMAX)
      character(5)::status(NMAX),status_orig(NMAX),status_save(NMAX)
      
      real(R8)::m1,m2,m1_orig,m2_orig,m_array_1_orig,m_array_2_orig
      
      s = ZERO
      ngf = ZERO
      rho = THREE
      rce_hill = ZERO
      rad = ONE
      rcrit = 10.0_R8
      x0 = ZERO
      v0 = ZERO
      
     
      name(1) = '     one'
      name(2) = '     two'
      name(3) = '   three'
      name(4) = '    four'
      name(5) = '    five'
      n = 5
      nbig = n
      index(1) = 1
      index(2) = 2
      index(3) = 3
      index(4) = 4
      index(5) = 5
      
      x0(1,1) = ONE
      x0(2,1) = ZERO
      x0(3,1) = ZERO
      x0(1,2) = 1.1_R8
      x0(2,2) = ZERO
      x0(3,2) = ZERO
      v0(1,1) = ZERO
      v0(2,1) = ONE
      v0(3,1) = ZERO
      v0(1,1) = ZERO
      v0(2,1) = -1.0_R8
      v0(3,1) = ZERO
      itarg = 1
      iproj = 2
      m(1:5) = ONE
      m_array_1_orig = 5.015
      m_array_2_orig = 5.015
      m(1) =  m_array_1_orig
      m(2) =  m_array_2_orig

      mfrag_min = 0.6_R8
      
      m1_orig = 8.9220
      m2_orig = ZERO
      m1 = m1_orig
      m2 = m2_orig

      status_orig = 'big  '
      status(:) = status_orig
      
      x(:,:) = x0
      v(:,:) = v0

      !write(*,*) x(1:2,1)
      !write(*,*) x(1:2,2)
      !write(*,*) x(1:2,3)
      !write(*,*) npair, ipair(1:10), jpair(1:10),denergy,nfrag,mfrag_min

      !write(*,*)

      call fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
           rce_hill,rad,rcrit,status,index,name)


      write(*,*)  "Masses, then positions, of new code"
      write(*,*) m(1:12)

      write(*,*) x(1:2,1)
      write(*,*) x(1:2,2)
      write(*,*) x(1:2,3)
      write(*,*) x(1:2,4)
      write(*,*) x(1:2,5)
      write(*,*) x(1:2,6)
      write(*,*) x(1:2,7)
      write(*,*) x(1:2,8)
      write(*,*) x(1:2,9)
      write(*,*) x(1:2,10)
      write(*,*) x(1:2,11)
      write(*,*) x(1:2,12)
!      write(*,*) x(1:2,13)
!      write(*,*) x(1:2,14)
!      write(*,*) x(1:2,15)
!      write(*,*) x(1:2,16)

      
      s = ZERO
      ngf = ZERO
      rho = THREE
      rce_hill = ZERO
      rad = ONE
      rcrit = 10.0_R8
      x0 = ZERO
      v0 = ZERO
      

      
      name(1) = '     one'
      name(2) = '     two'
      name(3) = '   three'
      name(4) = '    four'
      name(5) = '    five'
      n = 5
      nbig = n
      index(1) = 1
      index(2) = 2
      index(3) = 3
      index(4) = 4
      index(5) = 5
      
      x0(1,1) = ONE
      x0(2,1) = ZERO
      x0(3,1) = ZERO
      x0(1,2) = 1.1_R8
      x0(2,2) = ZERO
      x0(3,2) = ZERO
      v0(1,1) = ZERO
      v0(2,1) = ONE
      v0(3,1) = ZERO
      v0(1,1) = ZERO
      v0(2,1) = -1.0_R8
      v0(3,1) = ZERO
      itarg = 1
      iproj = 2
      m(1:5) = ONE
      m(1) =  m_array_1_orig
      m(2) =  m_array_2_orig
      
      m1 = m1_orig
      m2 = m2_orig
      
      x_save(:,:) = x
      v_save(:,:) = v
      status_save(:) = status
      x(:,:) = x0
      v(:,:) = v0
      status(:) = status_orig
      npair = 0
      ipair = 0
      jpair = 0
      denergy = ZERO
      nfrag = 0
      !write(*,*) npair, ipair(1:10), jpair(1:10),denergy,nfrag,mfrag_min
      write(*,*) '---------------------------------------'
      call other_fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
           rce_hill,rad,rcrit,status,index,name)

      write(*,*)  "Masses, then positions, of old code"
      write(*,*) m(1:12)

      write(*,*) x(1:2,1)
      write(*,*) x(1:2,2)
      write(*,*) x(1:2,3)
      write(*,*) x(1:2,4)
      write(*,*) x(1:2,5)
      write(*,*) x(1:2,6)
      write(*,*) x(1:2,7)
      write(*,*) x(1:2,8)
      write(*,*) x(1:2,9)
      write(*,*) x(1:2,10)
      write(*,*) x(1:2,11)
      write(*,*) x(1:2,12)
!      write(*,*) x(1:2,13)
!      write(*,*) x(1:2,14)
!      write(*,*) x(1:2,15)
!      write(*,*) x(1:2,16)
!      write(*,*) x(1:2,17)
!      write(*,*) x(1:2,18)
!      write(*,*) x(1:2,19)
!      write(*,*) x(1:2,20)
!      write(*,*) x(1:2,21)
!      write(*,*) x(1:2,22)
      
      

      if (all(x_save(:,:).eq.x(:,:).and.all(v_save(:,:).eq.v(:,:)))) then
         write(*,*) "The arrays are purely equal, which makes sense ONLY if m2 is not zero OR body two is being removed"
         write(*,*) "m2 is: ", m2_orig
         write(*,*) "m2 status is: ", status_save(2), "  ", status(2)
         
      else

         IF(ALL(x_save(:,1).EQ.x(:,1)).and.all(x_save(:,8:99).EQ.x(:,9:))) THEN
            write(*,*) "x is equal"
         ELSE
            write(*,*) "x is ****unequal"
         ENDIF
         
         IF(ALL(v_save(:,1).EQ.v(:,1)).and.all(v_save(:,8:99).EQ.v(:,9:))) THEN
            write(*,*) "v is equal"
         ELSE
            write(*,*) "v is ****unequal"
         ENDIF
         
         if (.not. (all(x_save(:,1).EQ.x(:,1))) ) then
            write(*,*) "first elements *****unequal"
            write(*,*) x_save(:,1)
            write(*,*) x(:,1)
         endif
         
         if (.not. (all(x_save(:,2).EQ.x(:,6))) ) then
            write(*,*) "second elements ****unequal"
            write(*,*) x_save(:,2)
            write(*,*) x(:,6)
         endif
         
         if (.not. (all(x_save(:,6:99).EQ.x(:,7:))) ) then
            write(*,*) "remaining elements *******unequal"
         else
            write(*,*) "remaining elements equal"
         endif

      endif

      if (all(status_save(:).eq.status(:))) then
         write(*,*) "The status arrays are purely equal, which makes sense ONLY if m2 is not zero OR body two is getting removed"
         write(*,*) "m2 is: ", m2_orig
         write(*,*) "m2 status is: ", status_save(2), "  ", status(2)
         
      else


         IF((status_save(1).eq.status(1)).and.(status_save(8).eq.status(9)).and.(status_save(9).eq.status(10))) THEN !:99).eq.status(9:))) THEN
            write(*,*) "status is equal"
         ELSE
            write(*,*) "status is ****unequal"
         ENDIF
         
         
         if (.not. (status_save(1).EQ.status(1)) ) then
            write(*,*) "first status *****unequal"
            write(*,*) status_save(1)
            write(*,*) status(1)
         endif
         
         if (.not. (status_save(2).EQ.status(6)) ) then
            write(*,*) "second status ****unequal"
            write(*,*) status_save(2)
            write(*,*) status(6)
         endif
         
         if (.not. ( (status_save(6).eq.status(7)).and.(status_save(7).eq.status(8))) ) then !status_save(6:99).EQ.status(7:)) ) then
            write(*,*) "remaining status *******unequal"
            write(*,*) status_save(6)
            write(*,*) status(7)
            write(*,*) status_save(7)
            write(*,*) status(8)
         else
            write(*,*) "remaining status equal"
         endif
         

      endif
      
    end program hey_there
