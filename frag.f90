!
! Version of Mercury7 that includes fragmentation during collisions
! based on results by Leinhardt and Stewart (2011).
!
! *** POSSIBLE CHANGE TO WAY OF STORING PARTICLES DURING HYBRID BS STEP:
!     INSTEAD OF MAKING NEW LOCAL ARRAYS XBS ETC, MOVE ALL CRITICAL PARTICLES
!     TO THE START OF THE GLOBAL ARRAYS. ONLY CONSIDER THE FIRST NCRIT PARTICLES
!     WHEN DOING THE BS INTEGRATION
!
! *** NEED NEW SUBROUTINE AFTER BS STEP THAT MAKES SURE ALL SMALL PARTICLES
!     ARE AT THE END OF THE GLOBAL PARTICLE ARRAYS
!
! ** THINK MORE CAREFULLY ABOUT WHEN FLAG_ACCEL NEEDS TO BE SET ***
!
! ** COLLIDE CENTRAL: DO H2B THEN COLLIDE THEN B2H ***
!
! Units: input and messages use AU, solar masses, days.
! Internal calculations and compressed output use CGS units.
!
! Last updated: 23 March 2012
!
!==============================================================================
    module kinds
      integer, parameter::I4 = selected_int_kind(9)
      integer, parameter::R4 = kind(1.0)
      integer, parameter::R8 = kind(1.d0)
!
! Maximum number of particles (can be changed by the user)
      integer(I4), parameter::NMAX = 2048
!
! Encounter data type:
!   I, J     = indices of particles involved in an encounter
!   D        = distance of closest apporoach
!   T        = time of closest approach
!   IXV, JXV = coords and velocities of particles at closest approach
      type::encounter
        integer(I4)::i,j
        real(R8)::d,t,im,jm,ix(3),iv(3),jx(3),jv(3)
      end type encounter
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
!==============================================================================
    module globals
      use kinds
      integer(I4)::algor, ndump, nfun, nindex, npair, ipair(NMAX), jpair(NMAX)
      integer(I4)::opt_time, icrit(NMAX)
      real(R8)::time,tstart,tstop,dt0,tout,tdump,tfun,tlog,dtfun,dtout,dtdump
      real(R8)::energy0, angmom0, denergy, dangmom
      real(R8)::mcen, scen(3), j2, j4, j6
      real(R8)::rcen, rmax, tol, cefac
      logical::flag_accel, flag_ejections, flag_ngf, flag_oblate, flag_stop
      logical::opt_user_force, opt_no_encounters, opt_collisions
      logical::opt_fragment
      character(80)::outfile(3), dumpfile(4)
      character(5)::goal
      character(6)::time_string
    end module globals
!==============================================================================
    module swift
      use kinds
! Maximum array size
     integer(I4), parameter::NPLMAX = 2048  ! max number of planets, including the Sun 
     integer(I4), parameter::NTPMAX = 2048 ! max number of test particles
! Size of the test particle status flag
      integer(I4), parameter::NSTAT = 3
! convergence criteria for danby
      real(R8), parameter::DANBYAC= 1.0d-14
      real(R8), parameter::DANBYB = 1.0d-13
! loop limits in the Laguerre attempts
      integer(I4), parameter::NLAG1 = 50
      integer(I4), parameter::NLAG2 = 400
! A small number
      real(R8), parameter::TINY=4.D-15
! trig stuff
      real(R8), parameter::PI = 3.141592653589793D0
      real(R8), parameter::TWOPI = 2.0_R8 * PI
      real(R8), parameter::PIBY2 = PI / 2.0_R8
      real(R8), parameter::DEGRAD = 180.0D0 / PI
    end module swift
!==============================================================================
    module frag_globals
      use kinds
      integer(I4)::nfrag
      real(R8)::mfrag_min
    end module frag_globals
!==============================================================================
    module frag_interfaces
!
    interface
      subroutine calc_absolute_coords (xrel,vrel,xcom,vcom,i,j,m,x,v)
      use kinds
      integer(I4), intent(in)::i,j
      real(R8), intent(in)::xrel(:),vrel(:),xcom(:),vcom(:),m(:)
      real(R8), intent(out)::x(:,:),v(:,:)
      end subroutine calc_absolute_coords
    end interface
!
    interface
      subroutine bounce_bodies (i,j,m,x,v,s,rad,name)
      use kinds
      integer(I4), intent(in)::i,j
      real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),rad(:)
      character(8), intent(inout)::name(:)
      end subroutine bounce_bodies
    end interface
!
    interface
      subroutine calc_cartesian_velocities (xrel,vrel,vr,vtheta,vphi)
      use kinds
      real(R8), intent(in)::xrel(:),vr,vtheta,vphi
      real(R8), intent(out)::vrel(:)
      end subroutine calc_cartesian_velocities
    end interface
!
    interface
      subroutine calc_impact_geometry (xrel,vrel,msum,rsum,b,v2imp)
      use kinds
      real(R8), intent(in)::xrel(:),vrel(:),msum,rsum
      real(R8), intent(out)::b,v2imp
      end subroutine calc_impact_geometry
    end interface
!
    interface
      function calc_largest_remnant (mtarg,mproj,rtarg,rproj,b,v2imp)
      use kinds
      real(R8), intent(in)::mtarg,mproj,rtarg,rproj,b,v2imp
      real(R8)::calc_largest_remnant
      end function calc_largest_remnant
    end interface
!
    interface
      function calc_qstar (mtarg,mproj,alpha)
      use kinds
      real(R8), intent(in)::mtarg,mproj,alpha
      real(R8)::calc_qstar
      end function calc_qstar
    end interface
!
    interface
      subroutine calc_relative_coords (m,x,v,i,j,xrel,vrel,xcom,vcom)
      use kinds
      integer(I4), intent(in)::i,j
      real(R8), intent(in)::m(:),x(:,:),v(:,:)
      real(R8), intent(out)::xrel(:),vrel(:),xcom(:),vcom(:)
      end subroutine calc_relative_coords
    end interface
!
    interface
      function calc_remnant_mass (mtarg,mproj,q,qstar)
      use kinds
      real(R8), intent(in)::mtarg,mproj,q,qstar
      real(R8)::calc_remnant_mass
      end function calc_remnant_mass
    end interface
!
    interface
      function calc_second_remnant (mtarg,mproj,rtarg,rproj,b,v2imp)
      use kinds
      real(R8), intent(in)::mtarg,mproj,rtarg,rproj,b,v2imp
      real(R8)::calc_second_remnant
      end function calc_second_remnant
    end interface
!
    interface
      subroutine calc_spherical_velocities (xrel,vrel,vr,vtheta,vphi)
      use kinds
      real(R8), intent(in)::xrel(:),vrel(:)
      real(R8), intent(out)::vr,vtheta,vphi
      end subroutine calc_spherical_velocities
    end interface
!
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
!
    interface
      subroutine write_coll_text (t,tname,pname,text,ltext)
      use kinds
      real(R8),      intent(in)::t
      integer(I4),   intent(in)::ltext
      character(8),  intent(in)::tname,pname
      character(25), intent(in)::text
      end subroutine write_coll_text
    end interface
!
    end module frag_interfaces
!==============================================================================
    module interfaces
!
    interface
      function accel_all (n,nbig,m,x,v,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:)
      real(R8)::accel_all(3,size(m))
      end function accel_all
    end interface
!
    interface
      function accel_bs2 (n,nbig,m,x,v,ngf,rcrit,flag)
      use kinds
      integer(I4), intent(in)::n,nbig,flag
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:),rcrit(:)
      real(R8)::accel_bs2(3,size(m))
      end function accel_bs2
    end interface
!
    interface
      function accel_hkce (n,nbig,m,x,v,ngf,rcrit)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:),rcrit(:)
      real(R8)::accel_hkce(3,size(m))
      end function accel_hkce
    end interface
!
    interface
      function accel_hybrid (n,nbig,m,x,v,ngf,rcrit)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:),rcrit(:)
      real(R8)::accel_hybrid(3,size(m))
      end function accel_hybrid
    end interface
!
    interface
      function accel_mvs (n,nbig,m,xh,vh,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),xh(:,:),vh(:,:),ngf(:,:)
      real(R8)::accel_mvs(3,size(m))
      end function accel_mvs
    end interface
!
    interface
      function accel_ngf (n,nbig,m,x,v,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:)
      real(R8)::accel_ngf(3,size(m))
      end function accel_ngf
    end interface
!
    interface
      function accel_oblate (n,nbig,m,x,v)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:)
      real(R8)::accel_oblate(3,size(m))
      end function accel_oblate
    end interface
!
    interface
      function accel_small (n,nbig,m,x,v,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:)
      real(R8)::accel_small(3,size(m))
      end function accel_small
    end interface
!
    interface
      function accel_user (n,nbig,m,x,v)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:)
      real(R8)::accel_user(3,size(m))
      end function accel_user
    end interface
!
    interface
      function acosh (x)
      use kinds
      real(R8), intent(in)::x
      real(R8)::acosh
      end function acosh
    end interface
!
    interface
      subroutine bs1 (dt,dt_nextcall,n,nbig,m,x0,v0,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),ngf(:,:)
      real(R8),    intent(inout)::dt,dt_nextcall,x0(:,:),v0(:,:)
      end subroutine bs1
    end interface
!
    interface
      subroutine bs2 (dt,dt_nextcall,n,nbig,m,x0,v0,ngf,rcrit,flag)
      use kinds
      integer(I4), intent(in)::n,nbig,flag
      real(R8),    intent(in)::m(:),ngf(:,:),rcrit(:)
      real(R8),    intent(inout)::dt,dt_nextcall,x0(:,:),v0(:,:)
      end subroutine bs2
    end interface
!
    interface
      subroutine calc_barycentric_coords (n,m,x,v,xb,vb,xcen,vcen)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::m(:),x(:,:),v(:,:)
      real(R8),    intent(out)::xb(:,:),vb(:,:),xcen(3),vcen(3)
      end subroutine calc_barycentric_coords
    end interface
!
    interface
      subroutine calc_bounding_box (dt,n,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::dt,x0(:,:),v0(:,:),x1(:,:),v1(:,:)
      real(R8),    intent(out)::xmin(:),xmax(:),ymin(:),ymax(:)
      end subroutine calc_bounding_box
    end interface
!
    interface
      subroutine calc_cartesian_coords (gmsum,q,e,incl,long,node,mean,x,v)
      use kinds
      real(R8), intent(in)::gmsum,q,e,incl,long,node,mean
      real(R8), intent(out)::x(3),v(3)
      end subroutine calc_cartesian_coords
    end interface
!
    interface
      subroutine calc_date (jd1,year1,month1,day1)
      use kinds
      real(R8),    intent(in)::jd1
      integer(I4), intent(out)::year1,month1
      real(R8),    intent(out)::day1
      end subroutine calc_date
    end interface
!
    interface
      subroutine calc_energy (n,nbig,m,x,v,s,energy,angmom)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),s(:,:)
      real(R8),    intent(out)::energy,angmom
      end subroutine calc_energy
    end interface
!
    interface
      character(8) function calc_float_string (x)
      use kinds
      real(R8), intent(in)::x
      end function calc_float_string
    end interface
!
    interface
      subroutine calc_jacobi_coords (n,nbig,m,x,xj)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:)
      real(R8),    intent(out)::xj(:,:)
      end subroutine calc_jacobi_coords
    end interface
!
    interface
      function calc_kepler (ecc,mean0)
      use kinds
      real(R8), intent(in)::ecc,mean0
      real(R8)::calc_kepler
      end function calc_kepler
    end interface
!
    interface
      subroutine calc_minimum (dt,s0,s1,ds0,ds1,smin,tmin)
      use kinds
      real(R8), intent(in)::dt,s0,s1,ds0,ds1
      real(R8), intent(out)::smin,tmin
      end subroutine calc_minimum
    end interface
!
    interface
      function calc_radius (n,m,rho)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::m(:),rho(:)
      real(R8)::calc_radius(size(m))
      end function calc_radius
    end interface
!
    interface
      function calc_rce (n,m,x,v,rce_hill)
      use kinds
      integer(I4), intent(in)::n
      real(R8), intent(in)::m(:),x(:,:),v(:,:),rce_hill(:)
      real(R8)::calc_rce(size(m))
      end function calc_rce
    end interface
!
    interface
      function calc_rhill (n,m,x,v)
      use kinds
      integer(I4), intent(in)::n
      real(R8), intent(in)::m(:),x(:,:),v(:,:)
      real(R8)::calc_rhill(size(m))
      end function calc_rhill
    end interface
!
    interface
      function calc_rcrit (dt,n,m,x,v)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::dt,m(:),x(:,:),v(:,:)
      real(R8)::calc_rcrit(size(m))
      end function calc_rcrit
    end interface
!
    interface
      character(8) function calc_real_string (x,xmin,xmax)
      use kinds
      real(R8), intent(in)::x,xmin,xmax
      end function calc_real_string
    end interface
!
    interface
      subroutine calc_sin_cos (x0,sx,cx)
      use kinds
      real(R8), intent(in)::x0
      real(R8), intent(out)::sx,cx
      end subroutine calc_sin_cos
    end interface
!
    interface
      subroutine calc_sinh_cosh (x,sx,cx)
      use kinds
      real(R8), intent(in)::x
      real(R8), intent(out)::sx,cx
      end subroutine calc_sinh_cosh
    end interface
!
    interface
      subroutine calc_spherical_polars (x,mag,theta,phi)
      use kinds
      real(R8), intent(in)::x(3)
      real(R8), intent(out)::mag,theta,phi
      end subroutine calc_spherical_polars
    end interface
!
    interface
      subroutine calc_substring (string,nsub,delimit)
      use kinds
      character(*), intent(in)::string
      integer(I4), intent(out)::nsub,delimit(:,:)
      end subroutine calc_substring
    end interface
!
    interface
      subroutine check_central (dt,n,nbig,m,x0,v0,x1,v1,nhit,hit)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::dt,m(:),x0(:,:),v0(:,:),x1(:,:),v1(:,:)
      integer(I4), intent(out)::nhit
      type(encounter), intent(out)::hit(:)
      end subroutine check_central
    end interface
!
    interface
      subroutine check_critical (dt,n,nbig,x0,v0,x1,v1,rcrit,ncrit,ncrit_big, &
        name)
      use kinds
      real(R8),     intent(in)::dt,x0(:,:),v0(:,:),x1(:,:),v1(:,:),rcrit(:)
      integer(I4),  intent(in)::n,nbig
      character(8), intent(in)::name(:)
      integer(I4),  intent(out)::ncrit,ncrit_big
      end subroutine check_critical
    end interface
!
    interface
      subroutine check_encounters (t,dt,n,nbig,m,x0,v0,x1,v1,rad,rce_hill, &
        name,nclo,clo,nhit,hit)
      use kinds
      integer(I4),  intent(in)::n,nbig
      real(R8),     intent(in)::t,dt
      real(R8),     intent(in)::m(:),x0(:,:),v0(:,:),x1(:,:),v1(:,:),rad(:),rce_hill(:)
      character(8), intent(in)::name(:)
      type(encounter), intent(inout)::clo(:),hit(:)
      end subroutine check_encounters
    end interface
!
    interface
      subroutine collide_bodies (t,i,j,n,nbig,m,x0,v0,x,v,s,ngf,rho,rce_hill,rad, &
        rcrit,status,index,name)
      use kinds
      integer(I4),  intent(in)::i,j
      real(R8),     intent(in)::t
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(in)::x0(:,:),v0(:,:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:)
      real(R8),     intent(inout)::rce_hill(:),rad(:),rcrit(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine collide_bodies
    end interface
!
    interface
      subroutine collide_central (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index, &
        name,nhit,hit)
      use kinds
      integer(I4),  intent(in)::nhit
      type(encounter), intent(in)::hit(:)
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:)
      real(R8),     intent(inout)::ngf(:,:),rho(:),rce_hill(:),rad(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine collide_central
    end interface
!
    interface
      subroutine convert_helio_to_jacobi (n,nbig,m,x,v)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:)
      real(R8),    intent(inout)::x(:,:),v(:,:)
      end subroutine convert_helio_to_jacobi
    end interface
!
    interface
      subroutine convert_jacobi_to_helio (n,nbig,m,x,v)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:)
      real(R8),    intent(inout)::x(:,:),v(:,:)
      end subroutine convert_jacobi_to_helio
    end interface
!
    interface
      subroutine convert_vbary_to_helio (n,m,v)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::m(:)
      real(R8),    intent(inout)::v(:,:)
      end subroutine convert_vbary_to_helio
    end interface
!
    interface
      subroutine convert_vhelio_to_bary (n,m,v)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::m(:)
      real(R8),    intent(inout)::v(:,:)
      end subroutine convert_vhelio_to_bary
    end interface
!
    interface
      subroutine corrector (dt,n,nbig,m,x,v,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::dt,m(:),ngf(:,:)
      real(R8),    intent(inout)::x(:,:),v(:,:)
      end subroutine corrector
    end interface
!
    interface
      function cross_product (a,b)
      use kinds
      real(R8), intent(in)::a(3),b(3)
      real(R8)::cross_product(3)
      end function cross_product
    end interface
!
    interface
      subroutine drift (dt,gmsum,x,v)
      use kinds
      real(R8), intent(in)::dt,gmsum
      real(R8), intent(inout)::x(3),v(3)
      end subroutine drift
    end interface
!
    interface
      subroutine driver_hybrid (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
      use kinds
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine driver_hybrid
    end interface
!
    interface
      subroutine driver_mvs (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
      use kinds
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine driver_mvs
    end interface
!
    interface
      subroutine driver_variable_step (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
      use kinds
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine driver_variable_step
    end interface
!
    interface
      subroutine ejections (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
      use kinds
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:)
      real(R8),     intent(inout)::ngf(:,:),rho(:),rce_hill(:),rad(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine ejections
    end interface
!
    interface
      subroutine finish (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
      use kinds
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine finish
    end interface
!
    interface
      subroutine inverse_corrector (dt,n,nbig,m,x,v,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::dt,m(:),ngf(:,:)
      real(R8),    intent(inout)::x(:,:),v(:,:)
      end subroutine inverse_corrector
    end interface
!
    interface
      subroutine jump (dt,n,m,x,v)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(in)::dt,m(:),v(:,:)
      real(R8),    intent(inout)::x(:,:)
      end subroutine jump
    end interface
!
    interface
      subroutine merge_bodies (itarg,iproj,m,x,v,s,status,name)
      use kinds
      integer(I4),  intent(in)::itarg,iproj
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine merge_bodies
    end interface
!
    interface
      subroutine output_coords (n,nbig,m,x,v,s,rho,status,index,name)
      use kinds
      integer(I4),  intent(in)::n,nbig
      real(R8),     intent(in)::m(:),x(:,:),v(:,:),s(:,:),rho(:)
      character(8), intent(in)::name(:)
      character(5), intent(in)::status(:)
      integer(I4),  intent(inout)::index(:)
      end subroutine output_coords
    end interface
!
    interface
      subroutine output_codes (n,nbig,status,index,name,m)
      use kinds
      integer(I4),  intent(in)::n,nbig
      character(8), intent(in)::name(:)
      character(5), intent(in)::status(:)
      integer(I4),  intent(inout)::index(:)
      real(R8),     intent(in)::m(:)
      end subroutine output_codes
    end interface
!
    interface
      subroutine output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
      use kinds
      integer(I4),  intent(in)::n,nbig,index(:)
      real(R8),     intent(in)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
      character(8), intent(in)::name(:)
      end subroutine output_datadump
    end interface
!
    interface
      subroutine output_encounters (n,nbig,rho,status,index,name,nclo,clo,m)
      use kinds
      integer(I4) , intent(in)::n,nbig,nclo
      real(R8),     intent(in)::rho(:)
      character(8), intent(in)::name(:)
      character(5), intent(in)::status(:)
      type(encounter), intent(in)::clo(:)
      integer(I4),  intent(inout)::index(:)
      real(R8),     intent(in)::m(:)
      end subroutine output_encounters
    end interface
!
    interface
      subroutine output_progress (n,nbig,m,x,v,s)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),x(:,:),v(:,:),s(:,:)
      end subroutine output_progress
    end interface
!
    interface
      subroutine ra15 (dt,dt_nextcall,n,nbig,m,x0,v0,ngf)
      use kinds
      integer(I4), intent(in)::n,nbig
      real(R8),    intent(in)::m(:),ngf(:,:)
      real(R8),    intent(inout)::dt,dt_nextcall,x0(:,:),v0(:,:)
      end subroutine ra15
    end interface
!
    interface
      subroutine remove_dead_bodies (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status, &
        index,name)
      use kinds
      integer(I4),  intent(inout)::n,nbig,index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:)
      real(R8),     intent(inout)::rce_hill(:),rad(:)
      character(8), intent(inout)::name(:)
      character(5), intent(inout)::status(:)
      end subroutine remove_dead_bodies
    end interface
!
    interface
      subroutine setup (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
      use kinds
      integer(I4),  intent(out)::n,nbig,index(:)
      real(R8),     intent(out)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
      character(8), intent(out)::name(:)
      character(5), intent(out)::status(:)
      end subroutine setup
    end interface
!
    interface
      subroutine sort_array (n,x,order)
      use kinds
      integer(I4), intent(in)::n
      real(R8),    intent(inout)::x(n)
      integer(I4), intent(out)::order(n)
      end subroutine
    end interface
!
    interface
      subroutine synch_epochs (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name,epoch)
      use kinds
      integer(I4),  intent(in)::n,nbig
      integer(I4),  intent(inout)::index(:)
      real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
      real(R8),     intent(inout)::rho(:),rce_hill(:),epoch(:)
      character(8), intent(inout)::name(:)
      end subroutine synch_epochs
    end interface
!

!!!!!!!! Added by Joshua Wallace
    interface
      subroutine mco_x2ov (m,x,v,fr,fv)
      use kinds
      real(R8),  intent(in)::m
      real(R8),  intent(out)::fr,fv
      real(R8),  intent(in)::x(3),v(3)
      end subroutine mco_x2ov
   end interface

    end module interfaces
!==============================================================================
    program mercury7_0
!
    use constants;    use globals
    use interfaces, only: driver_hybrid, driver_mvs, driver_variable_step, &
      finish, setup
    implicit none
!
    integer(I4)::n,nbig,index(NMAX)
    real(R8)::m(NMAX),x(3,NMAX),v(3,NMAX),s(3,NMAX),ngf(3,NMAX),rho(NMAX),rce_hill(NMAX)
    character(8)::name(NMAX)
    character(5)::status(NMAX)
    real(R8)::time_start, time_end
    integer(I4)::days, hours, minutes, seconds
!------------------------------------------------------------------------------
! Start timing
    call cpu_time(time_start)
! Set up the initial conditions
    call setup (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
    write(*,*) "Number of big bodies at time 0 years is: ", nbig
!
! Do the integration
    if (algor == 10) &
      call driver_hybrid (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    if (algor == 1.or.algor == 9) &
      call driver_mvs (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    if (algor == 2.or.algor == 3.or.algor == 4) &
      call driver_variable_step (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
! Finish
    call finish (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
! End timing
    call cpu_time(time_end)
    write(*,*) "CPU time: ", time_end-time_start, " seconds"
    days = int( (time_end-time_start)/86400.0_R8) 
    hours = int( mod((time_end-time_start)/3600.0_R8, 24.0_R8) )
    minutes = int( mod((time_end-time_start)/60.0_R8, 60.0_R8) )
    seconds = int( mod((time_end-time_start), 60.0_R8) )
    write(*,*) "CPU time, (D,H,M,S): ",days, hours, minutes, seconds

    end program mercury7_0
!==============================================================================
! Calculates the accelerations on all particles using a variable timestep
! algorithm.
!
    function accel_all (n,nbig,m,x,v,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_small
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:)
    real(R8)::accel_all(3,size(m))
!
    integer(I4)::i,j
    real(R8)::p(3),dx(3),temp1,temp2,d_1,d2,d_3,r_1,r2,r_3(NMAX)
!------------------------------------------------------------------------------
! Small forces
    accel_all(:,1:n) = accel_small (n,nbig,m,x,v,ngf)
!
! Direct forces between pairs of bodies
    do i = 1, nbig
      do j = i + 1, n
        dx = x(:,j)  -  x(:,i)
        d2 = dot_product (dx, dx)
        d_1 = ONE / sqrt(d2)
        d_3 = d_1 * d_1 * d_1
        temp1 = G * m(i) * d_3
        temp2 = G * m(j) * d_3
        accel_all(:,j) = accel_all(:,j)  -  temp1 * dx(:)
        accel_all(:,i) = accel_all(:,i)  +  temp2 * dx(:)
      end do
    end do
!
! Indirect terms (add these on last to reduce roundoff error)
    p = ZERO
    do i = 1, n
      r2 = dot_product (x(:,i), x(:,i))
      r_1 = ONE / sqrt(r2)
      r_3(i) = r_1 * r_1 * r_1
      temp1 = G * m(i) * r_3(i)
      p(:) = p(:)  -  temp1 * x(:,i)
    end do
!
    do i = 1, n
      temp1 = G * mcen * r_3(i)
      accel_all(:,i) = accel_all(:,i)  +  p(:)  -  temp1 * x(:,i)
    end do
!
    end function accel_all
!==============================================================================
! Calculates the accelerations on all particles when using the BS2 subroutine.
!
    function accel_bs2 (n,nbig,m,x,v,ngf,rcrit,flag)
    use constants;    use globals
    use interfaces, only: accel_all, accel_hkce
    implicit none
!
    integer(I4), intent(in)::n,nbig,flag
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:),rcrit(:)
    real(R8)::accel_bs2(3,size(m))
!------------------------------------------------------------------------------
    if (flag == 1) then
      accel_bs2(:,1:n) = accel_hkce (n,nbig,m,x,v,ngf,rcrit)
    else
      accel_bs2(:,1:n) = accel_all  (n,nbig,m,x,v,ngf)
    end if
!
    end function accel_bs2
!==============================================================================
! Calculates the accelerations on particles undergoing close encounters
! using a hybrid symplectic algorithm.
!
    function accel_hkce (n,nbig,m,x,v,ngf,rcrit)
!
    use constants;    use globals
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:),rcrit(:)
    real(R8)::accel_hkce(3,size(m))
!
    integer(I4)::i,j,k
    real(R8)::temp,fac,faci,facj,d_1,d2,q,q2,rcrit1,rcrit2,dx(3)
    logical::flag
!------------------------------------------------------------------------------
! Initialize accelerations
    accel_hkce(:,1:n) = ZERO
!
! Direct terms
    do k = 1, npair
      i = ipair(k);      j = jpair(k)
      flag = .true.
      dx = x(:,j)  -  x(:,i)
      d2 = dot_product(dx, dx)
      rcrit1 = max (rcrit(i), rcrit(j))
      rcrit2 = rcrit1 * rcrit1
!
      if (d2 < rcrit2) then
        d_1 = ONE / sqrt(d2)
        fac = d_1 * d_1 * d_1
!
        if (d2 > 0.01 * rcrit2) then
          temp = rcrit1 * d_1
          q = (ONE  -  0.1_R8 * temp) / (0.9_R8 * temp)
          q2 = q * q
          fac = fac * (ONE  -  q * q2 * (10.0_R8  -  15.0_R8 * q  +  6.0_R8 * q2))
        end if
!
        faci = G * m(i) * fac
        facj = G * m(j) * fac
        accel_hkce(:,j) = accel_hkce(:,j)  -  faci * dx
        accel_hkce(:,i) = accel_hkce(:,i)  +  facj * dx
      end if
    end do
!
! Solar terms
    do i = 1, n
      d2 = dot_product (x(:,i), x(:,i))
      d_1 = ONE / sqrt(d2)
      fac = G * mcen * d_1 * d_1 * d_1
      accel_hkce(:,i) = accel_hkce(:,i)  -  fac * x(:,i)
    end do
!
    end function accel_hkce
!==============================================================================
! Calculates the accelerations on all particles for the hybrid symplectic
! algorithm.
!
    function accel_hybrid (n,nbig,m,x,v,ngf,rcrit)
!
    use constants;    use globals
    use interfaces, only: accel_small
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:),rcrit(:)
    real(R8)::accel_hybrid(3,size(m))
!
    integer(I4)::i,j
    real(R8)::dx(3),rcrit1,rcrit2,d2,d_1,temp,fac,faci,facj,q,q2
!------------------------------------------------------------------------------
! Small forces
    accel_hybrid(:,1:n) = accel_small (n,nbig,m,x,v,ngf)
!
! Direct forces between pairs of bodies
    do i = 1, nbig
      do j = i + 1, n
        dx(:) = x(:,j)  -  x(:,i)
        d2 = dot_product(dx, dx)
        rcrit1 = max(rcrit(i), rcrit(j))
        rcrit2 = rcrit1 * rcrit1
!
        fac = ZERO
        if (d2 > 0.01_R8 * rcrit2) then
          d_1 = ONE / sqrt(d2)
          fac = d_1 * d_1 * d_1
!
          if (d2 < rcrit2) then
            temp = rcrit1 * d_1
            q = (ONE  -  0.1_R8 * temp) / (0.9_R8 * temp)
            q2 = q * q
            fac = fac * q * q2 * (10.0_R8  -  15.0_R8 * q  +  6.0_R8 * q2)
          end if
        end if
!
        faci = G * m(i) * fac
        facj = G * m(j) * fac
        accel_hybrid(:,j) = accel_hybrid(:,j)  -  faci * dx
        accel_hybrid(:,i) = accel_hybrid(:,i)  +  facj * dx
      end do
    end do
!
    end function accel_hybrid
!==============================================================================
! Calculates the accelerations on all particles using a mixed-variable
! symplectic algorithm.
!
    function accel_mvs (n,nbig,m,xh,vh,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_small, calc_jacobi_coords
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),xh(:,:),vh(:,:),ngf(:,:)
    real(R8)::accel_mvs(3,size(m))
!
    integer(I4)::i,j,k
    real(R8)::xj(3,size(m)),a1(3,size(m)),a2(3,size(m)),a3(3,size(m))
    real(R8)::gmcen,fac0,fac2,minside,d_1,d2,d_3,faci,facj
    real(R8)::rh_1,rh2,rh_3,rj_1,rj2,rj_3,dx(3),a0(3),a0tp(3)
!------------------------------------------------------------------------------
! Initialize accelerations to zero
    a0 = ZERO;               a1(:,1:n) = ZERO
    a2(:,1:n) = ZERO;        a3(:,1:n) = ZERO
    minside = ZERO;             gmcen = G * mcen
!
! Calculate Jacobi coordinates
    call calc_jacobi_coords (n,nbig,m,xh,xj)
!
! Calculate indirect acceleration terms (zero for innermost body)
    do k = 2, nbig
      minside = minside  +  m(k-1)
      rh2  = dot_product(xh(:,k), xh(:,k))
      rj2  = dot_product(xj(:,k), xj(:,k))
      rh_1 = ONE / sqrt(rh2)
      rj_1 = ONE / sqrt(rj2)
      rh_3 = rh_1 * rh_1 * rh_1
      rj_3 = rj_1 * rj_1 * rj_1
!
! Add to A0 term
      fac0  = G * m(k) * rh_3
      a0(:) = a0(:)  -  fac0 * xh(:,k)
!
! Calculate A1 and A2 for this body
      a1(:,k) = gmcen * (xj(:,k) * rj_3  -  xh(:,k) * rh_3)
      fac2  = gmcen * m(k) * rj_3 / (minside + mcen)
      a2(:,k) = a2(:,k-1)  +  fac2 * xj(:,k)
    end do
!
! Calculate A0 for test particles
    rh2  = dot_product(xh(:,1), xh(:,1))
    rh_1 = ONE / sqrt(rh2)
    rh_3 = rh_1 * rh_1 * rh_1
    fac0 = m(1) * rh_3
    a0tp(:) = a0(:)  -  fac0 * xh(:,1)
!
! Calculate A3 (direct terms)
    do i = 1, nbig - 1
      do j = i + 1, nbig
        dx(:) = xh(:,j)  -  xh(:,i)
        d2 = dot_product(dx, dx)
        d_1 = ONE / sqrt(d2)
        d_3 = d_1 * d_1 * d_1
        faci = G * m(i) * d_3
        facj = G * m(j) * d_3
        a3(:,j) = a3(:,j)  -  faci * dx(:)
        a3(:,i) = a3(:,i)  +  facj * dx(:)
      end do
    end do
!
    do i = 1, nbig
      do j = nbig + 1, n
        dx(:) = xh(:,j)  -  xh(:,i)
        d2 = dot_product(dx, dx)
        d_1 = ONE / sqrt(d2)
        d_3 = d_1 * d_1 * d_1
        faci = G * m(i) * d_3
        a3(:,j) = a3(:,j)  -  faci * dx(:)
      end do
    end do
!
! Big-body accelerations
    accel_mvs(1,1:nbig) = a0(1)  +  a1(1,1:nbig)  +  a2(1,1:nbig)  +  a3(1,1:nbig)
    accel_mvs(2,1:nbig) = a0(2)  +  a1(2,1:nbig)  +  a2(2,1:nbig)  +  a3(2,1:nbig)
    accel_mvs(3,1:nbig) = a0(3)  +  a1(3,1:nbig)  +  a2(3,1:nbig)  +  a3(3,1:nbig)
!
! Small-body accelerations
    accel_mvs(1,nbig+1:n) = a0tp(1)  +  a3(1,nbig+1:n)
    accel_mvs(2,nbig+1:n) = a0tp(2)  +  a3(2,nbig+1:n)
    accel_mvs(3,nbig+1:n) = a0tp(3)  +  a3(3,nbig+1:n)
!
! Small forces
    accel_mvs(:,1:n) = accel_mvs(:,1:n)  +  accel_small (n,nbig,m,xh,vh,ngf)
!
    end function accel_mvs
!==============================================================================
! Calculates the accelerations on all particles due to non-gravitational
! cometary jet forces. The number of active particles is N.
!
    function accel_ngf (n,nbig,m,x,v,ngf)
    use constants;    use globals
    use interfaces, only: cross_product
    implicit none
    real(R8), parameter::AA =  0.111262_R8
    real(R8), parameter::NN =  5.093_R8
    real(R8), parameter::KK = -4.6142_R8
    real(R8), parameter::MM = -2.15_R8
    real(R8), parameter::R0 =  2.808_R8
!
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:)
    real(R8)::accel_ngf(3,size(m))
    integer(I4)::i
    real(R8)::tran(3),norm(3),r2(NMAX),r,rv,q,gfac,a1,a2,a3
!------------------------------------------------------------------------------
    accel_ngf(:,1:n) = ZERO
!
    do i = 1, n
      r2(i) = dot_product(x(:,i), x(:,i))
      if (r2(i) > 88_R8.and.abs(ngf(1,i)) < 1.e-7_R8.and. &
        abs(ngf(2,i)) < 1.e-7_R8.and.abs(ngf(3,i)) < 1.e-7_R8) cycle
!
! If body is < 9.36 from central body or non-grav forces are very large...
      r = sqrt(r2(i))
      rv = dot_product(x(:,i), v(:,i))
!
! Calculate Q = R / R0 and GFAC (see Marsden et al. 1973)
      q = r / R0
      gfac = AA * q**MM * (ONE  +  q**NN)**KK
!
! Within-orbital-plane transverse vector components
      tran(:) = r2(i) * v(:,i)  -  rv * x(:,i)
!
! Orbit-normal vector components
      norm = cross_product (x(:,i), v(:,i))
!
! Multiplication factors
      a1 = ngf(1,i) * gfac / r
      a2 = ngf(2,i) * gfac / sqrt(dot_product(tran, tran))
      a3 = ngf(3,i) * gfac / sqrt(dot_product(norm, norm))
!
! X,Y and Z components of non-gravitational acceleration
      accel_ngf(:,i) = a1 * x(:,i)  +  a2 * tran(:)  +  a3 * norm(:)
    end do
!
    end function accel_ngf
!==============================================================================
! Calculates the acceleration on all particles due to the oblateness of the
! central body. The number of active particles is N.
!
    function accel_oblate (n,nbig,m,x,v)
    use constants;    use globals
    implicit none
!
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:)
    real(R8)::accel_oblate(3,size(m))
    integer(I4)::i
    real(R8)::acen(3)
    real(R8)::jr2,jr4,jr6,r2,r_1,r_2,u2,u4,u6,tmp1,tmp2,tmp3,tmp4
!------------------------------------------------------------------------------
    accel_oblate(:,1:n) = ZERO;    acen = ZERO
!
    do i = 1, n
!
! Calculate barycentric accelerations on the bodies
      r2 = dot_product(x(:,i), x(:,i))
      r_1 = ONE / sqrt(r2)
      r_2 = r_1 * r_1
      jr2 = j2 * r_2
      jr4 = j4 * r_2 * r_2
      jr6 = j6 * r_2 * r_2 * r_2
      u2 = x(3,i) * x(3,i) * r_2
      u4 = u2 * u2
      u6 = u4 * u2
!
      tmp1 = G * mcen * r_2 * r_1
      tmp2 = jr2 * (7.5_R8 * u2  -  1.5_R8) &
           + jr4 * (39.375_R8 * u4  -  26.25_R8 * u2  +  1.875_R8) &
           + jr6 * (187.6875_R8 * u6  -  216.5625_R8 * u4  +  59.0625_R8 * u2  &
           - 2.1875_R8)
      tmp3 = jr2 * THREE &
           + jr4 * (17.5_R8 * u2  -  7.5_R8) &
           + jr6 * (86.625_R8 * u4  -  78.75_R8 * u2  +  13.125_R8)
!
      accel_oblate(1,i) = x(1,i) * tmp1 * tmp2
      accel_oblate(2,i) = x(2,i) * tmp1 * tmp2
      accel_oblate(3,i) = x(3,i) * tmp1 * (tmp2  -  tmp3)
!
! Calculate barycentric accelerations on the central body
      tmp4 = m(i) / mcen
      acen(:) = acen(:)  -  tmp4 * accel_oblate(:,i)
    end do
!
! Convert to accelerations with respect to the central body
    accel_oblate(1,1:n) = accel_oblate(1,1:n)  -  acen(1)
    accel_oblate(2,1:n) = accel_oblate(2,1:n)  -  acen(2)
    accel_oblate(3,1:n) = accel_oblate(3,1:n)  -  acen(3)
!
    end function accel_oblate
!==============================================================================
! Calculates the small accelerations on all the particles.
! The number of active particles is N.
!
    function accel_small (n,nbig,m,x,v,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_ngf, accel_oblate, accel_user
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),ngf(:,:)
    real(R8)::accel_small(3,size(m))
!------------------------------------------------------------------------------
    accel_small(:,1:n) = ZERO
!
! Oblateness of central body
    if (flag_oblate) accel_small(:,1:n) = accel_small(:,1:n) &
                                        +  accel_oblate (n,nbig,m,x,v)
!
! Cometary non-gravitational acceleration
    if (flag_ngf)    accel_small(:,1:n) = accel_small(:,1:n) &
                                        +  accel_ngf (n,nbig,m,x,v,ngf)
!
! User-defined force
    if (opt_user_force) accel_small(:,1:n) = accel_small(:,1:n) &
                                           +  accel_user (n,nbig,m,x,v)
!
    end function accel_small
!==============================================================================
! Calculates the acceleration on all particles due to a user-defined force.
! The number of active particles is N.
    function accel_user (n,nbig,m,x,v)
    use constants;    use globals
    implicit none
!
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:)
    real(R8)::accel_user(3,size(m))
!------------------------------------------------------------------------------
    accel_user(:,1:n) = ZERO
!
    end function accel_user
!==============================================================================
! Calculates the inverse hyperbolic cosine of an angle X.
!
    function acosh (x)
    use constants
    implicit none
!
    real(R8), intent(in)::x
    real(R8)::acosh
!------------------------------------------------------------------------------
    if (x >= ONE) then
      acosh = log (x  +  sqrt(x * x  -  ONE))
    else
      acosh = ZERO
    end if
!
    end function acosh
!==============================================================================
! Advances all particles for a time DT using a general purpose Bulirsch-Stoer
! algorithm.
! On entry: DT is the desired timestep
! On exit:  DT is the actual timestep used
!           DT_NEXTCALL is the recommend timestep next time the routine is called
!
    subroutine bs1 (dt,dt_nextcall,n,nbig,m,x0,v0,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_all
    implicit none
    real(R8),parameter::SHRINK = 0.55_R8
    real(R8),parameter::GROW   = 1.3_R8
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),ngf(:,:)
    real(R8),    intent(inout)::dt,dt_nextcall,x0(:,:),v0(:,:)
!
    integer(I4)::i,j,j1,k,ns
    real(R8)::tmp0,tmp1,tmp2,errmax,h,hx2,h2(8)
    real(R8)::xscal(NMAX),vscal(NMAX),x(3,NMAX),v(3,NMAX)
    real(R8)::xend(3,NMAX),vend(3,NMAX),a(3,NMAX),a0(3,NMAX),dx(3,NMAX,8),dv(3,NMAX,8)
!------------------------------------------------------------------------------
! Calculate accelerations at the start of the step
    a0(:,1:n) = accel_all (n,nbig,m,x0,v0,ngf)
!
! Arrays used to scale the relative error (1 / absolute distance or velocity)
    do i = 1, n
      xscal(i) = ONE / dot_product (x0(:,i), x0(:,i))
      vscal(i) = ONE / dot_product (x0(:,i), x0(:,i))
    end do
!
! Try to do the step using current stepsize
    do
!
! For each value of NS, do a modified-midpoint integration with 2N substeps
      do ns = 1, 8
        h = dt / (TWO * dble(ns))
        h2(ns) = 0.25_R8 / (ns * ns)
        hx2 = h * TWO
        x(:,1:n) = x0(:,1:n)  +  h * v0(:,1:n)
        v(:,1:n) = v0(:,1:n)  +  h * a0(:,1:n)
!
        a(:,1:n) = accel_all (n,nbig,m,x,v,ngf)
        xend(:,1:n) = x0(:,1:n)  +  hx2 * v(:,1:n)
        vend(:,1:n) = v0(:,1:n)  +  hx2 * a(:,1:n)
!
        do j = 2, ns
          a(:,1:n) = accel_all (n,nbig,m,xend,vend,ngf)
          x(:,1:n) = x(:,1:n)  +  hx2 * vend(:,1:n)
          v(:,1:n) = v(:,1:n)  +  hx2 * a(:,1:n)
!
          a(:,1:n) = accel_all (n,nbig,m,x,v,ngf)
          xend(:,1:n) = xend(:,1:n)  +  hx2 * v(:,1:n)
          vend(:,1:n) = vend(:,1:n)  +  hx2 * a(:,1:n)
        end do
!
        a(:,1:n) = accel_all (n,nbig,m,xend,vend,ngf)
        dx(:,1:n,ns) = HALF * (xend(:,1:n)  +  x(:,1:n) + h * vend(:,1:n))
        dv(:,1:n,ns) = HALF * (vend(:,1:n)  +  v(:,1:n) + h * a(:,1:n))
!
! Update the DX and DV arrays used for polynomial extrapolation
        do j = ns - 1, 1, -1
          j1 = j + 1
          tmp0 = ONE / (h2(j) - h2(ns))
          tmp1 = tmp0 * h2(j1)
          tmp2 = tmp0 * h2(ns)
          dx(:,1:n,j) = tmp1 * dx(:,1:n,j1)  -  tmp2 * dx(:,1:n,j)
          dv(:,1:n,j) = tmp1 * dv(:,1:n,j1)  -  tmp2 * dv(:,1:n,j)
        end do
!
! After several integrations, test the relative error on extrapolated values
        if (ns > 3) then
!
! Maximum relative position and velocity errors (last D terms added)
          errmax = ZERO
          do k = 1, n
            tmp1 = dot_product (dx(:,k,1), dx(:,k,1))
            tmp2 = dot_product (dv(:,k,1), dv(:,k,1))
            errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
          end do
!
! If error is smaller than TOL, update position and velocity arrays, and exit
          if (errmax <= tol * tol) then
            x0(:,1:n) = ZERO;      v0(:,1:n) = ZERO
            do j = 1, ns
              x0(:,1:n) = x0(:,1:n)  +  dx(:,1:n, j)
              v0(:,1:n) = v0(:,1:n)  +  dv(:,1:n, j)
            end do
!
! Recommend a stepsize for the next call to this subroutine
            if (ns >= 8) dt_nextcall = dt * SHRINK
            if (ns < 7)  dt_nextcall = dt * GROW
            return
          end if
        end if
!
      end do
!
! If errors were too large, redo the step with half the previous step size.
      dt = dt * HALF
    end do
!
    end subroutine bs1
!==============================================================================
! Advances all particles for a time DT using a Bulirsch-Stoer algorithm designed
! for forces that depend on coordinates only (not velocities).
! On entry: DT is the desired timestep
! On exit:  DT is the actual timestep used
!           DT_NEXTCALL is the recommend timestep next time the routine is called
!
    subroutine bs2 (dt,dt_nextcall,n,nbig,m,x0,v0,ngf,rcrit,flag)
!
    use constants;    use globals
    use interfaces, only: accel_bs2
    implicit none
    real(R8), parameter::SHRINK = 0.55_R8
    real(R8), parameter::GROW   = 1.3_R8
    integer(I4), intent(in)::n,nbig,flag
    real(R8),    intent(in)::m(:),ngf(:,:),rcrit(:)
    real(R8),    intent(inout)::dt,dt_nextcall,x0(:,:),v0(:,:)
!
    integer(I4)::i,j,j1,k,ns
    real(R8)::tmp0,tmp1,tmp2,errmax,h,h2(12),hby2,h2by2
    real(R8)::xend(3,NMAX),xscal(NMAX),vscal(NMAX)
    real(R8)::dx(3,NMAX,12),dv(3,NMAX,12),a(3,NMAX),a0(3,NMAX),b(3,NMAX),c(3,NMAX)
!------------------------------------------------------------------------------
! Calculate accelerations at the start of the step
    a0(:,1:n) = accel_bs2 (n,nbig,m,x0,v0,ngf,rcrit,flag)
!
! Arrays used to scale the relative error (1 / absolute distance or velocity)
    do i = 1, n
      xscal(i) = ONE / dot_product (x0(:,i), x0(:,i))
      vscal(i) = ONE / dot_product (x0(:,i), x0(:,i))
    end do
!
! Try to do the step using current stepsize
    do
!
! For each value of NS, do a modified-midpoint integration with NS substeps
      do ns = 1, 12
        h = dt / (dble(ns))
        hby2  = HALF * h
        h2(ns) = h * h
        h2by2 = HALF * h2(ns)
!
        b(:,1:n) = HALF * a0(:,1:n)
        c(:,1:n) = ZERO
        xend(:,1:n) = h2by2 * a0(:,1:n)  +  h * v0(:,1:n)  +  x0(:,1:n)
!
        do j = 2, ns
          a(:,1:n) = accel_bs2 (n,nbig,m,xend,v0,ngf,rcrit,flag)
          tmp0 = h * dble(j)
          b(:,1:n) = b(:,1:n)  +  a(:,1:n)
          c(:,1:n) = c(:,1:n)  +  b(:,1:n)
          xend(:,1:n) = h2(ns) * c(:,1:n)  +  h2by2 * a0(:,1:n) &
                               +  tmp0 * v0(:,1:n) +  x0(:,1:n)
        end do
!
        a(:,1:n) = accel_bs2 (n,nbig,m,xend,v0,ngf,rcrit,flag)
        dx(:,1:n,ns) = xend(:,1:n)
        dv(:,1:n,ns) = h * b(:,1:n)  +  hby2 * a(:,1:n)  +  v0(:,1:n)
!
! Update the DX and DV arrays used for polynomial extrapolation
        do j = ns - 1, 1, -1
          j1 = j + 1
          tmp0 = ONE / (h2(j)  -  h2(ns))
          tmp1 = tmp0 * h2(j1)
          tmp2 = tmp0 * h2(ns)
          dx(:,1:n,j) = tmp1 * dx(:,1:n,j1)  -  tmp2 * dx(:,1:n,j)
          dv(:,1:n,j) = tmp1 * dv(:,1:n,j1)  -  tmp2 * dv(:,1:n,j)
        end do
!
! After several integrations, test the relative error on extrapolated values
        if (ns > 3) then
!
! Maximum relative position and velocity errors (last D term added)
          errmax = ZERO
          do k = 1, n
            tmp1 = dot_product (dx(:,k,1), dx(:,k,1))
            tmp2 = dot_product (dv(:,k,1), dv(:,k,1))
            errmax = max( errmax, tmp1*xscal(k), tmp2*vscal(k) )
          end do
!
! If error is smaller than TOL, update position and velocity arrays and exit
          if (errmax <= tol * tol) then
            x0(:,1:n) = ZERO;      v0(:,1:n) = ZERO
            do j = 1, ns
              x0(:,1:n) = x0(:,1:n)  +  dx(:,1:n, j)
              v0(:,1:n) = v0(:,1:n)  +  dv(:,1:n, j)
            end do
!
! Recommend a stepsize for the next call to this subroutine
            if (ns >= 8) dt_nextcall = dt * SHRINK
            if (ns < 7)  dt_nextcall = dt * GROW
            return
          end if
        end if
!
      end do
!
! If errors were too large, redo the step with half the previous step size.
      dt = dt * HALF
    end do
!
    end subroutine bs2
!==============================================================================
! Calculates barycentric coordinates XB and velocities VB for a set of particles
! given their coordinates X and velocities V with respect to the central body
! and their masses M.
! The number of active particles is N.
!
    subroutine calc_barycentric_coords (n,m,x,v,xb,vb,xcen,vcen)
!
    use constants;      use globals
    implicit none
    integer(I4), intent(in)::n
    real(R8),    intent(in)::m(:),x(:,:),v(:,:)
    real(R8),    intent(out)::xb(:,:),vb(:,:),xcen(3),vcen(3)
!
    real(R8)::mtot
!------------------------------------------------------------------------------
! Calculate coordinates and velocities of the central body
    xcen(1) = sum(m(1:n) * x(1,1:n))
    xcen(2) = sum(m(1:n) * x(2,1:n))
    xcen(3) = sum(m(1:n) * x(3,1:n))
    vcen(1) = sum(m(1:n) * v(1,1:n))
    vcen(2) = sum(m(1:n) * v(2,1:n))
    vcen(3) = sum(m(1:n) * v(3,1:n))
!
    mtot = sum(m(1:n))  +  mcen
    xcen(:) = -xcen(:) / mtot
    vcen(:) = -vcen(:) / mtot
!
! Calculate the barycentric coordinates and velocities
    xb(1,1:n) = x(1,1:n)  +  xcen(1)
    xb(2,1:n) = x(2,1:n)  +  xcen(2)
    xb(3,1:n) = x(3,1:n)  +  xcen(3)
    vb(1,1:n) = v(1,1:n)  +  vcen(1)
    vb(2,1:n) = v(2,1:n)  +  vcen(2)
    vb(3,1:n) = v(3,1:n)  +  vcen(3)
!
    end subroutine calc_barycentric_coords
!==============================================================================
! Calculates rectangular X-Y boxes encompassing the trajectory of each particle
! during a time interval DT given the particle's initial and final coordinates
! and velocities. The number of active particles is N.
! The box's coordinates are XMIN, XMAX, YMIN and YMAX.
!
    subroutine calc_bounding_box (dt,n,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
    use constants
    implicit none
!
    integer(I4), intent(in)::n
    real(R8),    intent(in)::dt,x0(:,:),v0(:,:),x1(:,:),v1(:,:)
    real(R8),    intent(out)::xmin(:),xmax(:),ymin(:),ymax(:)
    real(R8)::temp(size(x0,2))
!------------------------------------------------------------------------------
    xmin(1:n) = min ( x0(1,1:n), x1(1,1:n) )
    xmax(1:n) = max ( x0(1,1:n), x1(1,1:n) )
    ymin(1:n) = min ( x0(2,1:n), x1(2,1:n) )
    ymax(1:n) = max ( x0(2,1:n), x1(2,1:n) )
!
! If X velocity changes sign during interval, find extreme value using interpolation
    where (v0(1,1:n) * v1(1,1:n) < ZERO)
      temp(1:n) = (v0(1,1:n) * x1(1,1:n)  -  v1(1,1:n) * x0(1,1:n)  &
           -  HALF * dt * v0(1,1:n) * v1(1,1:n)) / (v0(1,1:n)  -  v1(1,1:n))
!
      xmin(1:n) = min (xmin(1:n), temp(1:n))
      xmax(1:n) = max (xmax(1:n), temp(1:n))
    end where
!
! If Y velocity changes sign during interval, find extreme value using interpolation
    where (v0(2,1:n) * v1(2,1:n) < ZERO)
      temp(1:n) = (v0(2,1:n) * x1(2,1:n)  -  v1(2,1:n) * x0(2,1:n)  &
           -  HALF * dt * v0(2,1:n) * v1(2,1:n)) / (v0(2,1:n)  -  v1(2,1:n))
!
      ymin(1:n) = min (ymin(1:n), temp(1:n))
      ymax(1:n) = max (ymax(1:n), temp(1:n))
    end where
!
    end subroutine calc_bounding_box
!==============================================================================
! Calculates the Cartesian coordinates X and velocities V of a particle given its
! orbital elements. GMSUN is the product of the gravitational constant and the
! sum of the particle mass plus the central mass.
! INCL is the orbital inclination.
! LONG is the longitude of pericentre.
! NODE is the nodal longitude.
! MEAN is the mean anomaly.
!
    subroutine calc_cartesian_coords (gmsum,q,e,incl,long,node,mean,x,v)
!
    use constants
    use interfaces, only: calc_kepler, calc_sin_cos, calc_sinh_cosh
    implicit none
    real(R8), intent(in)::gmsum,q,e,incl,long,node,mean
    real(R8), intent(out)::x(3),v(3)
!
    real(R8)::d1(3),d2(3),a,peri,orbel_fhybrid,orbel_zget
    real(R8)::ci,si,cn,sn,cw,sw,ce,se,romes,temp,z1,z2,z3,z4
!------------------------------------------------------------------------------
! Change from longitude of perihelion to argument of perihelion
    peri = long  -  node
!
! Rotation factors
    call calc_sin_cos (incl,si,ci)
    call calc_sin_cos (peri,sw,cw)
    call calc_sin_cos (node,sn,cn)
    z1 = cw * cn;      z2 = cw * sn
    z3 = sw * cn;      z4 = sw * sn
    d1 = (/ z1  -  z4 * ci,  z2  +  z3 * ci, sw * si/)
    d2 = (/-z3  -  z2 * ci, -z4  +  z1 * ci, cw * si/)
!
! Semi-major axis
    a = q / (ONE  -  e)
!
! Ellipse
    if (e < ONE) then
      romes = sqrt(ONE  -  e * e)
      temp = calc_kepler (e,mean)
      call calc_sin_cos (temp,se,ce)
      z1 = a * (ce  -  e)
      z2 = a * romes * se
      temp = sqrt(gmsum / a) / (ONE  -  e * ce)
      z3 = -se * temp
      z4 = romes * ce * temp
! Parabola
    else if (e == ONE) then
      ce = orbel_zget(mean)
      z1 = q * (ONE  -  ce * ce)
      z2 = TWO * q * ce
      z4 = sqrt(TWO * gmsum / q) / (ONE  +  ce * ce)
      z3 = -ce * z4
! Hyperbola
    else
      romes = sqrt(e * e  -  ONE)
      temp = orbel_fhybrid(e,mean)
      call calc_sinh_cosh (temp,se,ce)
      z1 = a * (ce  -  e)
      z2 = -a * romes * se
      temp = sqrt(gmsum / abs(a)) / (e * ce  -  ONE)
      z3 = -se * temp
      z4 = romes * ce * temp
    endif
!
    x(:) = d1(:) * z1  +  d2(:) * z2
    v(:) = d1(:) * z3  +  d2(:) * z4
!
    end subroutine calc_cartesian_coords
!==============================================================================
! Calculates a calendar date (year, month, day) given a Julian day number.
!
    subroutine calc_date (jd1,year1,month1,day1)
    use constants
    implicit none
    integer(I4), parameter::FOUR_YEAR    = 1461
    integer(I4), parameter::FOUR_CENTURY = 146097
!
    real(R8),    intent(in)::jd1
    integer(I4), intent(out)::year1,month1
    real(R8),    intent(out)::day1
    integer(I4)::jd,a,b,c,d,e,f
    real(R8)::frac
!------------------------------------------------------------------------------
    jd = floor(jd1  +  HALF)
    frac = (jd1  +  HALF)  -  dble(jd)
!
! Gregorian calendar
    if (jd >= 2299160.5_R8) then
      a = jd   +  32044
      b = (4 * a  +  3) / FOUR_CENTURY
      c = a  -  b * FOUR_CENTURY / 4
! Julian calendar
    else
      b = 0
      c = jd  +  32082
    end if
!
    d = (4 * c  +  3) / FOUR_YEAR
    e = c  -  FOUR_YEAR * d / 4
    f = (5 * e  +  2) / 153
!
    day1   = dble(e  -  (153 * f  +  2) / 5  +  1)  +  frac
    month1 = f  +  3  -  12 * (f / 10)
    year1  = b * 100  +  d  -  4800  +  f / 10
!
    end subroutine calc_date
!==============================================================================
! Calculates the total energy and angular momentum of a system of particles given
! their masses M, coordinates X and velocities V with respect to the central body
! and their spin angular momenta S. The number of active particles is N.
! Interactions between particles with indices > NBIG are ignored.
!
    subroutine calc_energy (n,nbig,m,x,v,s,energy,angmom)
!
    use constants;      use globals
    use interfaces, only: cross_product, calc_barycentric_coords
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),s(:,:)
    real(R8),    intent(out)::energy,angmom
!
    integer(I4)::i,j
    real(R8)::xb(3,size(m)),vb(3,size(m)),ang(3),dx(3),xcen(3),vcen(3)
    real(R8)::ke,pe,r2,r_1,r_2,r_4,r_6,u2,u4,u6
!------------------------------------------------------------------------------
    ke = ZERO;      pe = ZERO;      ang = ZERO
!
! Calculate barycentric coordinates and velocities of all bodies
    call calc_barycentric_coords (n,m,x,v,xb,vb,xcen,vcen)
!
! Angular momentum and kinetic energy of central body
    ang(:) = mcen * cross_product (xcen, vcen)  +  scen(:)
    ke = mcen * dot_product (vcen, vcen)
!
    do i = 1, n
      ang(:) = ang(:)  +  m(i) * cross_product (xb(:,i), vb(:,i))  +  s(:,i)
      ke  = ke  +  m(i) * dot_product (vb(:,i), vb(:,i))
    end do
!
! Potential energy terms due to pairs of bodies
    do i = 1, nbig
      do j = i + 1, n
        dx(:) = xb(:,j)  -  xb(:,i)
        r2 = dot_product(dx, dx)
        if (r2 /= 0) pe = pe  -  m(i) * m(j) / sqrt(r2)
      end do
    end do
!
! Potential energy terms involving the central body
    do i = 1, n
      dx(:) = xb(:,i)  -  xcen(:)
      r2 = dot_product(dx, dx)
      if (r2 /= 0) pe = pe  -  mcen * m(i) / sqrt(r2)
    end do
!
! Corrections for oblateness of the central body
    if (flag_oblate) then
      do i = 1, n
        dx(:) = xb(:,i)  -  xcen(:)
        r2 = dot_product(dx(:), dx(:))
        r_1 = ONE / sqrt(r2)
        r_2 = r_1 * r_1
        r_4 = r_2 * r_2
        r_6 = r_4 * r_2
        u2 = dx(3) * dx(3) * r_2
        u4 = u2 * u2
        u6 = u4 * u2
        pe = pe  +  mcen * m(i) * r_1 &
           * (j2 * r_2 * (1.5_R8 * u2  -  HALF) &
           +  j4 * r_4 * (4.375_R8 * u4  -  3.75_R8 * u2  +  .375_R8) &
           +  j6 * r_6 * (14.4375_R8 * u6 - 19.6875_R8 * u4 + 6.5625_R8 * u2 - .3125_R8))
      end do
    end if
!
    energy = HALF * ke  +  G * pe
    angmom = sqrt(dot_product(ang, ang))
!
    end subroutine calc_energy
!==============================================================================
! Converts a floating point number X into a compressed character string.
!
    character(8) function calc_float_string (x)
    use constants
    use interfaces, only: calc_real_string
    implicit none
!
    real(R8), intent(in)::x
    integer(I4)::ex
    real(R8)::ax,y
!------------------------------------------------------------------------------
    if (x == 0) then
      y = HALF
    else
      ax = abs(x)
      ex = int(log10(ax))
      if (ax >= 1) ex = ex + 1
      y = ax * (10.0_R8**(-ex))
      if (y == 1) then
        y = y * 0.1_R8
        ex = ex + 1
      end if
      y = sign(y,x) * HALF  +  HALF
    end if
!
    calc_float_string(1:8) = calc_real_string (y, ZERO, ONE)
    ex = ex  +  112
    if (ex > 223) ex = 223
    if (ex < 0) ex = 0
    calc_float_string(8:8) = char(ex + 32)
!
    end function calc_float_string
!==============================================================================
! For a set of particles, calculates Jacobi coordinates given the coordinates
! with respect to the central body.
!
    subroutine calc_jacobi_coords (n,nbig,m,x,xj)
!
    use constants;    use globals
    implicit none
    integer(I4), intent(in)::nbig,n
    real(R8),    intent(in)::m(:),x(:,:)
    real(R8),    intent(out)::xj(:,:)
!
    integer(I4)::i
    real(R8)::mx(3),temp,mtot
!------------------------------------------------------------------------------
    mtot = ZERO;      mx = ZERO
!
! Big bodies
    do i = 1, nbig
      temp = ONE / (mtot  +  mcen)
      mtot = mtot  +  m(i)
      xj(:,i) = x(:,i)  -  temp * mx(:)
      mx(:) = mx(:)  +  m(i) * x(:,i)
    end do
!
! Small bodies
    xj(:,nbig+1:n) = x(:,nbig+1:n)
!
    end subroutine calc_jacobi_coords
!==============================================================================
! Solves Kepler's equation for an ellipse, returning the eccentric anomaly given 
! the eccentricity and mean anomaly.
!
    function calc_kepler (ecc,mean0)
!
    use constants
    implicit none
    real(R8), parameter::PI2 = PI * PI
    real(R8), parameter::TWOTHIRD = 2.0_R8 / 3.0_R8
    real(R8), parameter::I24TH    = 1.0_R8 / 24.0_R8
    real(R8), parameter::CONST1   = 3.0_R8 * PI2 / (PI2  -  6.0_R8)
    real(R8), parameter::CONST2   = 1.6_R8 * PI  / (PI2  -  6.0_R8)
    real(R8), intent(in)::ecc,mean0
    real(R8)::calc_kepler
!
    real(R8)::sgn,ome,mean,mean2,alpha,d,q,q2,r,w
    real(R8)::bige,ese,ece,f0,f1,f2,f3,f4,d3,d4,d4_2,d5
!------------------------------------------------------------------------------
! Reduce mean anomaly to range 0 to PI
    mean = modulo (mean0, TWOPI)
    sgn = ONE
    if (mean > PI) then
      mean = TWOPI  -  mean;      sgn = -ONE
    end if
!
! Initial approximation
    ome = ONE  -  ecc
    mean2 = mean * mean
!
    alpha = CONST1  +  CONST2 * (PI  -  mean) / (ONE  +  ecc)
    d = THREE * ome  +  alpha * ecc
    q = TWO * alpha * d * ome  -  mean2
    q2 = q * q
    r = mean * (THREE * alpha * d * (d  -  ome)  +  mean2)
    w = (r  +  sqrt(q * q2  +  r * r))**TWOTHIRD
    bige = (TWO * r * w / (w * w  +  w * q  +  q2)  +  mean) / d
!
! Refine using high-order Newton's method
    ese = ecc * sin(bige)
    ece = ecc * cos(bige)
    f0 =  bige  -  ese  -  mean
    f1 =  ONE  -  ece
    f2 =  ese
    f3 =  ece
    f4 = -ese
!
    d3 = -f0 * f1 / (f1 * f1  -  HALF * f0 * f2)
    d4 = -f0 / (SIXTH * d3 * d3 * f3  +  HALF * d3 * f2  +  f1)
    d4_2 = d4 * d4
    d5 = -f0 / (I24TH * d4 * d4_2 * f4  +  SIXTH * d4_2 * f3  &
         +  HALF * d4 * f2  +  f1)
!
    calc_kepler = (bige  +  d5) * sgn
!
    end function calc_kepler
!==============================================================================
! Estimates the minimum value of a variable in a time interval DT given the
! initial and final values and their derivatives, using cubic interpolation.
! The estimated minimum is SMIN, while the time of minimum is TMIN.
!
    subroutine calc_minimum (dt,s0,s1,ds0,ds1,smin,tmin)
    use constants
    implicit none
!
    real(R8), intent(in)::dt,s0,s1,ds0,ds1
    real(R8), intent(out)::smin,tmin
    real(R8)::a,b,c,temp,tau
!------------------------------------------------------------------------------
! If a minimum occurs during the time interval...
    if (ds0 * dt < 0.and.ds1 * dt > 0) then
      temp = 6.0_R8 * (s0  -  s1)
      a = temp  +  THREE * dt * (ds0  +  ds1)
      b = temp  +  TWO   * dt * (ds0  +  ds1 * TWO)
      c = dt * ds1
!
      tau = ZERO
      temp = -HALF * (b  +  sign (sqrt(max(b * b  -  FOUR * a * c, ZERO)), b) )
      if (temp /= 0) tau = c / temp
!
! Make sure TAU falls in the interval -1 < TAU < 0
      tau = min (tau,  ZERO)
      tau = max (tau, -ONE)
!
! Calculate minimum value SMIN, and time of minimum TMIN
      tmin = tau * dt
      temp = ONE  +  tau
      smin = tau  * tau  * ((THREE  +  TWO * tau) * s0  +  temp * dt * ds0) &
           + temp * temp * ((ONE    -  TWO * tau) * s1  +  tau  * dt * ds1)
!
! Make sure minimum value is not negative
      smin = max(smin, ZERO)
!
! If there is no minimum during the time interval...
    else
      if (s0 <= s1) then
        smin = s0;        tmin = -dt
      else
        smin = s1;        tmin = ZERO
      end if
    end if
!
    end subroutine calc_minimum
!==============================================================================
! Calculates the radius of all particles given their mass M and density RHO.
! The number of active particles is N.
!
    function calc_radius (n,m,rho)
    use constants
    implicit none
    real(R8), parameter::FACTOR = THREE / (FOUR * PI)
!
    integer(I4), intent(in)::n
    real(R8),    intent(in)::m(:),rho(:)
    real(R8)::calc_radius(size(m))
!------------------------------------------------------------------------------
    calc_radius(1:n) = (FACTOR * m(1:n) / rho(1:n))**THIRD
!
    end function calc_radius
!==============================================================================
! Converts a real number X between XMIN and XMAX into a compressed character
! string.
!
    character(8) function calc_real_string (x,xmin,xmax)
    use constants
    implicit none
!
    real(R8), intent(in)::x,xmin,xmax
    integer(I4)::j
    real(R8)::y,z
!------------------------------------------------------------------------------
    calc_real_string(1:8) = '        '
    y = (x - xmin) / (xmax - xmin)
!
    if (y >= 1) then
      do j = 1, 8
        calc_real_string(j:j) = char(255)
      end do
    else if (y > 0) then
      z = y
      do j = 1, 8
        z = mod(z, ONE) * 224.0_R8
        calc_real_string(j:j) = char(int(z) + 32)
      end do
    end if
!
    end function calc_real_string
!==============================================================================
! Calculates the Hill radius for all particles
!   RHILL = a * (M / 3 Mcen) ^ (1/3)
! When the semi-major axis A is negative, the distance R is used instead.
!
    function calc_rhill (n,m,x,v)
!
    use constants;      use globals
    implicit none
    integer(I4), intent(in)::n
    real(R8), intent(in)::m(:),x(:,:),v(:,:)
    real(R8)::calc_rhill(size(m))
!
    integer::i
    real(R8)::r,v2,gmsum,a
!------------------------------------------------------------------------------
    do i = 1, n
      r  = sqrt(dot_product(x(:,i), x(:,i)))
      v2 =      dot_product(v(:,i), v(:,i))
      gmsum = G * (mcen  +  m(i))
!
      a = gmsum * r / (2.0_R8 * gmsum  -  r * v2)
      if (a <= ZERO) a = r
!
      calc_rhill(i) = a * (THIRD * m(i) / mcen)**THIRD
    end do
!
    end function calc_rhill
!==============================================================================
! Calculates the distances at which close encounters will be recorded given the
! close encounter limits in Hill radii.
!
    function calc_rce (n,m,x,v,rce_hill)
!
    use constants;      use globals
    use interfaces, only: calc_rhill
    implicit none
    integer(I4), intent(in)::n
    real(R8), intent(in)::m(:),x(:,:),v(:,:),rce_hill(:)
    real(R8)::calc_rce(size(m))
!
    real(R8)::rhill(NMAX)
!------------------------------------------------------------------------------
    rhill(1:n) = calc_rhill(n,m,x,v)
!
    calc_rce(1:n) = rce_hill(1:n) * rhill(1:n)
!
    end function calc_rce
!==============================================================================
! Calculates the critical distance for using the Bulirsch-Stoer integrator for
! all particles. The number of active particles is N.
!
    function calc_rcrit (dt,n,m,x,v)
!
    use constants;      use globals
    implicit none
    real(R8), parameter::N2 = 0.4_R8
    integer(I4), intent(in)::n
    real(R8),    intent(in)::dt,m(:),x(:,:),v(:,:)
    real(R8)::calc_rcrit(size(m))
!
    integer(I4)::i
    real(R8)::a(NMAX),rhill(NMAX),temp,amin,vmax,r,v2,gmsum
!------------------------------------------------------------------------------
! Calculate Hill radii
    do i = 1, n
      r  = sqrt(dot_product(x(:,i), x(:,i)))
      v2 =      dot_product(v(:,i), v(:,i))
      gmsum = G * (mcen  +  m(i))
      a(i) = gmsum * r / (2.0_R8 * gmsum  -  r * v2)
!
      temp = a(i)
      if (a(i) <= ZERO) temp = r
      rhill(i) = temp * (THIRD * m(i) / mcen)**THIRD
    end do
!
    amin = minval (a(1:n))
    vmax = sqrt (G * mcen / amin)
    temp = N2 * abs(dt) * vmax
    calc_rcrit(1:n) = max(rhill(1:n) * cefac, temp)
!
    end function calc_rcrit
!==============================================================================
! Calculates the sine and cosine of an angle X0.
!
    subroutine calc_sin_cos (x0,sx,cx)
!
    use constants
    implicit none
    real(R8), intent(in)::x0
    real(R8), intent(out)::sx,cx
!
    real(R8)::x
!------------------------------------------------------------------------------
    x = modulo (x0, TWOPI)
    cx = cos(x)
    sx = sqrt(ONE  -  cx * cx)
    if (x > PI) sx = -sx
!
    end subroutine calc_sin_cos
!==============================================================================
! Calculates the hyperbolic sine and hyperbolic cosine of an angle X.
!
    subroutine calc_sinh_cosh (x,sx,cx)
!
    use constants
    implicit none
    real(R8), intent(in)::x
    real(R8), intent(out)::sx,cx
!------------------------------------------------------------------------------
    sx = sinh(x)
    cx = sqrt (ONE  +  sx * sx)
!
    end subroutine calc_sinh_cosh
!==============================================================================
! Calculates spherical polar coordinates given Cartesian coordinates.
! The angle THETA is measured with respect to the X-Y plane.
!
    subroutine calc_spherical_polars (x,r,theta,phi)
!
    use constants
    implicit none
    real(R8), intent(in)::x(3)
    real(R8), intent(out)::r,theta,phi
!    real(R8)::p
!------------------------------------------------------------------------------
!    p = sqrt(x(1) * x(1)  +  x(2) * x(2))
!
    r     = sqrt(dot_product(x, x))
    theta = acos (x(3) / r)
    phi   = modulo (atan2 (x(2), x(1)), TWOPI)
!
    end subroutine calc_spherical_polars
!==============================================================================
! Given a string STRING, calculates the start and end point of each substring,
! where substrings are delimited by spaces of `=' characters.
!
    subroutine calc_substring (string,nsub,delimit)
!
    use constants
    implicit none
    character(*), intent(in)::string
    integer(I4), intent(out)::nsub,delimit(:,:)
!
    integer(I4)::j,k,len1
    character(150)::c
!------------------------------------------------------------------------------
    len1 = len(string)
    nsub = 0;      j = 0;      c = ' '
    delimit(1,1) = -1
!
! Find the start of substring
 10 j = j + 1
    if (j > len1) return
    c = string(j:j)
    if (c == ' '.or.c == '=') goto 10
!
! Find the end of substring
    k = j
 20 k = k + 1
    if (k > len1) goto 30
    c = string(k:k)
    if (c /= ' '.and.c /= '=') goto 20
!
! Store details for this substring
 30 nsub = nsub + 1
    delimit(1,nsub) = j
    delimit(2,nsub) = k - 1
!
    if (k < len1) then
      j = k
      goto 10
    end if
!
    end subroutine calc_substring
!==============================================================================
! Checks all particles to see if they approach within RCEN of the central
! body during a time interval DT. The initial and final coordinates are
! X0 and X1. The initial and final velocities are V0 and V1.
!
    subroutine check_central (dt,n,nbig,m,x0,v0,x1,v1,nhit,hit)
!
    use constants;    use globals
    use interfaces, only: cross_product, acosh
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::dt,m(:),x0(:,:),v0(:,:),x1(:,:),v1(:,:)
    integer(I4), intent(out)::nhit
    type(encounter), intent(out)::hit(:)
!
    integer(I4)::j
    real(R8)::rcen2,a,q,u0,uhit,m0,mhit,mean_motion,r0
    real(R8)::ang(3),ang2,p,r20,r21,rv0,rv1,temp,e,v2,r2min,gmsum,vesc
!------------------------------------------------------------------------------
    nhit = 0;      rcen2 = rcen * rcen
!
! Escape velocity of central body
    vesc = sqrt(2.0_R8 * G * mcen / rcen)
!
! Maximum distance from which object could hit central body during time DT
    temp = rcen  +  dt * vesc
    r2min = temp * temp
!
! Check for collisions with the central body
    do j = 1, n
      r20 = dot_product(x0(:,j), x0(:,j))
      r21 = dot_product(x1(:,j), x1(:,j))
      rv0 = dot_product(x0(:,j), v0(:,j))
      rv1 = dot_product(x1(:,j), v1(:,j))
!
      if (min(r20,r21) > rcen2.and.(rv0 * dt > 0.or.rv1 * dt < 0)) cycle
!
! If inside the central body, or passing through pericentre, check if
! pericentric distance Q is less than radius of central body
      ang(:) = cross_product (x0(:,j), v0(:,j))
      v2     = dot_product   (v0(:,j), v0(:,j))
      gmsum = G * (mcen  +  m(j))
      ang2 = dot_product(ang, ang)
      p = ang2 / gmsum
      r0 = sqrt(r20)
      temp = ONE  +  p * (v2 / gmsum  -  2.0_R8 / r0)
      e = sqrt( max(temp, ZERO) )
      q = p / (ONE  +  e)
      if (q > rcen) cycle
!
! If the body hit the central body
      nhit = nhit  +  1
      hit(nhit) % i = 0;      hit(nhit) % j = j
      hit(nhit) % d = rcen
!
! Time of impact relative to the end of the timestep
! (u = eccentric anomaly, m = mean anomaly, n = mean motion)
      if (e < 1) then
        a = q / (ONE  -  e)
        uhit = sign (acos((ONE  -  rcen / a) / e), -dt)
        u0   = sign (acos((ONE  -  r0   / a) / e), rv0)
        mhit = mod (uhit  -  e * sin(uhit) + PI, TWOPI)  -  PI
        m0   = mod (u0    -  e * sin(u0)   + PI, TWOPI)  -  PI
      else
        a = q / (e  -  ONE)
        uhit = sign (acosh((ONE  -  rcen / a) / e), -dt)
        u0   = sign (acosh((ONE  -  r0   / a) / e), rv0)
        mhit = mod (uhit  -  e * sinh(uhit) + PI, TWOPI)  -  PI
        m0   = mod (u0    -  e * sinh(u0)   + PI, TWOPI)  -  PI
      end if
!
      mean_motion = sqrt(gmsum / (a * a * a))
      hit(nhit) % t = (mhit  -  m0) / mean_motion  +  time
    end do
!
    end subroutine check_central
!==============================================================================
! Checks each pair of interacting particles to see if they pass within 
! RCRIT of each other during a time interval DT. The initial and final
! coordinates are X0 and X1. The initial and final velocities are V0
! and V1.
! The routine returns a list of NPAIR pairs of particles that have close
! encounters within the time interval. The indices of each pair are stored
! in IPAIR and JPAIR, while the indices of all particles involved in close
! encounters are stored in ICRIT.
!
    subroutine check_critical (dt,n,nbig,x0,v0,x1,v1,rcrit,ncrit,ncrit_big, &
      name)
!
    use constants;    use globals
    use interfaces, only: calc_bounding_box, calc_minimum
    implicit none
    real(R8),     intent(in)::dt,x0(:,:),v0(:,:),x1(:,:),v1(:,:),rcrit(:)
    integer(I4),  intent(in)::n,nbig
    character(8), intent(in)::name(:)
    integer(I4),  intent(out)::ncrit,ncrit_big
!
    integer(I4)::i,j,k,back(NMAX),crit(NMAX)
    real(R8)::xmin(NMAX),xmax(NMAX),ymin(NMAX),ymax(NMAX)
    real(R8)::d0,d1,d0t,d1t,d2min,temp,tmin,rcrit2,dx0(3),dv0(3),dx1(3),dv1(3)
    integer(I4)::itmp,jtmp
!------------------------------------------------------------------------------
    ncrit = 0;      npair = 0;      ncrit_big = 0
    crit = 0;      back = 0
!
! Calculate box of maximum and minimum values of x and y coordinates
    call calc_bounding_box (dt,n,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
!
! Adjust box values for the Big bodies by symplectic critical distance
    xmin(1:nbig) = xmin(1:nbig)  -  rcrit(1:nbig)
    xmax(1:nbig) = xmax(1:nbig)  +  rcrit(1:nbig)
    ymin(1:nbig) = ymin(1:nbig)  -  rcrit(1:nbig)
    ymax(1:nbig) = ymax(1:nbig)  +  rcrit(1:nbig)
!
! Identify pairs whose x-y boxes overlap, and calculate minimum separation
    do i = 1, nbig
      do j = i + 1, n
        if (xmax(i) < xmin(j).or.xmax(j) < xmin(i).or. &
            ymax(i) < ymin(j).or.ymax(j) < ymin(i)) cycle
!
! Calculate initial and final separations and derivatives
        dx0(:) = x0(:,i)  -  x0(:,j)
        dx1(:) = x1(:,i)  -  x1(:,j)
        dv0(:) = v0(:,i)  -  v0(:,j)
        dv1(:) = v1(:,i)  -  v1(:,j)
!
        d0 = dot_product(dx0, dx0)
        d1 = dot_product(dx1, dx1)
        d0t = TWO * dot_product(dx0, dv0)
        d1t = TWO * dot_product(dx1, dv1)
!
! If separation derivative changes sign, find the minimum separation
        d2min = BIG_NUMBER
        if (d0t*dt <= 0.and.d1t*dt >= 0) &
          call calc_minimum (dt,d0,d1,d0t,d1t,d2min,tmin)
!
! Determine the maximum separation that would qualify as an encounter
        temp = max(rcrit(i), rcrit(j))
        rcrit2 = temp * temp
!
! If minimum separation is small enough, flag this as a critical pair
        d2min = min (d0,d1,d2min)
        if (d2min <= rcrit2) then
          npair = npair  +  1
          ipair(npair) = i;      jpair(npair) = j
          crit(i) = 2;           crit(j) = 2
        end if
      end do
    end do
!
! Create array containing all bodies that are within critical distance
    do i = 1, n
      if (crit(i) /= 0) then
        ncrit = ncrit  +  1
        icrit(ncrit) = i
        back(i) = ncrit
        if (i <= nbig) ncrit_big = ncrit_big  +  1
      end if
    end do
!
! Convert indices for critical pairs to local indices used by ACCEL_HKCE
    do k = 1, npair
      itmp = ipair(k);    jtmp = jpair(k)
      ipair(k) = back(ipair(k))
      jpair(k) = back(jpair(k))
    end do
!
    end subroutine check_critical
!==============================================================================
! Checks all pairs of interacting particles to see if they are (i) within 
! RCE_HILL Hill radii of each other, and (ii) going through a close encounter
! minimum, during a time interval DT.
! The initial and final coordinatates are X0 and X1. The initial and final
! velocities are V0 and V1.
! If encounters are not allowed, the routine sets a flag to stop the
! integration.
!
    subroutine check_encounters (t,dt,n,nbig,m,x0,v0,x1,v1,rad,rce_hill, &
      name,nclo,clo,nhit,hit)
!
    use constants;    use globals
    use interfaces, only: calc_bounding_box, calc_date, calc_minimum, &
      calc_rce
    implicit none
    integer(I4),  intent(in)::n,nbig
    real(R8),     intent(in)::t,dt
    real(R8),     intent(in)::m(:),x0(:,:),v0(:,:),x1(:,:),v1(:,:),rad(:),rce_hill(:)
    character(8), intent(in)::name(:)
    type(encounter), intent(inout)::clo(:),hit(:)
!
    real(R8)::dx0(3),dv0(3),dx1(3),dv1(3),rce(NMAX)
    integer(I4)::i,j,itarg,iproj,nclo,nhit,year1,month1
    real(R8)::xmin(NMAX),xmax(NMAX),ymin(NMAX),ymax(NMAX)
    real(R8)::d0,d1,d0t,d1t,tmp0,tmp1,d2min,d2ce,d2hit,temp,tmin,t1
    logical::flag
!------------------------------------------------------------------------------
    nclo = 0;      nhit = 0;      flag = .false.
!
! Convert close encounter limits from Hill radii to absolute distances
    rce(1:n) = calc_rce (n,m,x0,v0,rce_hill)
!
! Calculate box of maximum and minimum values of x and y coordinates
    call calc_bounding_box (dt,n,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
!
! Adjust values by the maximum close-encounter radius plus a fudge factor
    do i = 1, n
      temp = rce(i) * 1.2_R8
      xmin(i) = xmin(i)  -  temp
      xmax(i) = xmax(i)  +  temp
      ymin(i) = ymin(i)  -  temp
      ymax(i) = ymax(i)  +  temp
    end do
!
! Identify pairs whose x-y boxes overlap, and calculate minimum separation
    do i = 1, nbig
      do j = i + 1, n
        if (xmax(i) < xmin(j).or.xmax(j) < xmin(i).or. &
            ymax(i) < ymin(j).or.ymax(j) < ymin(i)) cycle
!
! More massive body is the target, less massive body is the projectile
        itarg = i;      iproj = j
        if (m(j) > m(i).and.j <= nbig) then
          itarg = j;      iproj = i
        end if
!
! Calculate initial and final separations and derivatives
        dx0(:) = x0(:,i)  -  x0(:,j)
        dx1(:) = x1(:,i)  -  x1(:,j)
        dv0(:) = v0(:,i)  -  v0(:,j)
        dv1(:) = v1(:,i)  -  v1(:,j)
!
        d0 = dot_product(dx0, dx0)
        d1 = dot_product(dx1, dx1)
        d0t = TWO * dot_product(dx0, dv0)
        d1t = TWO * dot_product(dx1, dv1)
!
! Estimate minimum separation during the time interval, using interpolation
        call calc_minimum (dt,d0,d1,d0t,d1t,d2min,tmin)
        d2min = max(ZERO, d2min)
!
        temp  = max (rce(i), rce(j))
        d2ce  = temp * temp
        temp  = rad(i)  +  rad(j)
        d2hit = temp * temp
!
! Flag any close encounters
        if (d2min <= d2ce) then
          flag = .true.
        end if
!
! If minimum separation qualifies as an encounter or a collision occurs, store details
        if ((d2min <= d2ce.and.d0t*dt <= 0.and.d1t*dt >= 0).or.(d2min <= d2hit)) then
          nclo = nclo  +  1
          if (nclo > NMAX) then
 230        open (23,file=outfile(3),status='old',access='append',err=230)
            write (23,'(/,2a,/,a)') ' WARNING: ', &
              'Total number of current close encounters exceeds CMAX.'
            close (23)
          else
            clo(nclo) % i = itarg
            clo(nclo) % j = iproj
            clo(nclo) % t = tmin  +  t
            clo(nclo) % d = sqrt(d2min)
            clo(nclo) % im= m(i)
            clo(nclo) % jm= m(j)
!
! Make linear interpolation to estimate coordinates at time of closest approach
            tmp1 = -tmin / dt
            tmp0 = ONE  -  tmp1
            clo(nclo) % ix(:) = tmp0 * x0(:,i)  +  tmp1 * x1(:,i)
            clo(nclo) % iv(:) = tmp0 * v0(:,i)  +  tmp1 * v1(:,i)
            clo(nclo) % jx(:) = tmp0 * x0(:,j)  +  tmp1 * x1(:,j)
            clo(nclo) % jv(:) = tmp0 * v0(:,j)  +  tmp1 * v1(:,j)
            ! Store details of any collisions
            if (d2min <= d2hit) then
               nhit = nhit  +  1
               hit(nhit) % i = itarg
               hit(nhit) % j = iproj
               hit(nhit) % t = tmin  +  t
               hit(nhit) % d = sqrt(d2min)
               hit(nhit) % im= m(i)
               hit(nhit) % jm= m(j)
               hit(nhit) % ix(:) = tmp0 * x0(:,i)  +  tmp1 * x1(:,i)
               hit(nhit) % iv(:) = tmp0 * v0(:,i)  +  tmp1 * v1(:,i)
               hit(nhit) % jx(:) = tmp0 * x0(:,j)  +  tmp1 * x1(:,j)
               hit(nhit) % jv(:) = tmp0 * v0(:,j)  +  tmp1 * v1(:,j)
!               write(*,*) hit(nhit) % ix(:)
!               write(*,*) hit(nhit) % iv(:)
!               write(*,*) hit(nhit) % jx(:)
!               write(*,*) hit(nhit) % jv(:)
            end if
          end if
        end if
!

!
! Move on to the next pair of bodies
      end do
    end do
!
! If encounters are not allowed, set a flag to stop the integration
    if (opt_no_encounters.and.flag) then
      flag_stop = .true.
      i = clo(1) % i;      j = clo(1) % j
!
  20  open (23, file=outfile(3), status='old', access='append',err=20)
! If time style is Gregorian date then...
      if (opt_time == 1) then
        call calc_date (time/DAY,year1,month1,t1)
        write (23,'(5a,/,9x,a,i10,1x,i2,1x,f4.1)') ' WARNING: ', &
          'Stopping integration due to an encounter between ', &
          name(i),',',name(j),' at ',year1,month1,t1
! Otherwise...
      else
        if (opt_time == 0) t1 = time / DAY
        if (opt_time == 2) t1 = (time  -  tstart) / DAY
        if (opt_time == 3) t1 = (time  -  tstart) / YEAR
        write (23,'(5a,/,9x,a,f16.3,a)') ' WARNING: ', &
          'Stopping integration due to an encounter between ', &
          name(i),',',name(j),' at ',t1,time_string
      end if
      close(23)
    end if
!
    end subroutine check_encounters
!==============================================================================
! For a set of particles, converts coordinates and velocities with respect to
! the central body to Jacobi coordinates and velocities.
!
    subroutine convert_helio_to_jacobi (n,nbig,m,x,v)
!
    use constants;    use globals
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:)
    real(R8),    intent(inout)::x(:,:),v(:,:)
!
    integer(I4)::i
    real(R8)::xj(3,size(m)),vj(3,size(m)),mx(3),mv(3),temp,mtot
!------------------------------------------------------------------------------
    mtot = ZERO;      mx = ZERO;      mv = ZERO
!
    do i = 1, nbig
      temp = ONE / (mtot  +  mcen)
      mtot = mtot  +  m(i)
      xj(:,i) = x(:,i)  -  temp * mx(:)
      vj(:,i) = v(:,i)  -  temp * mv(:)
      mx(:) = mx(:)  +  m(i) * x(:,i)
      mv(:) = mv(:)  +  m(i) * v(:,i)
    end do
!
    x(:,1:nbig) = xj(:,1:nbig)
    v(:,1:nbig) = vj(:,1:nbig)
!
    end subroutine convert_helio_to_jacobi
!==============================================================================
! For a set of particles, converts Jacobi coordinates and velocities to
! coordinates and velocities with respect to the central body.
!
    subroutine convert_jacobi_to_helio (n,nbig,m,x,v)
!
    use constants;    use globals
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:)
    real(R8),    intent(inout)::x(:,:),v(:,:)
!
    integer(I4)::i
    real(R8)::xh(3,size(m)),vh(3,size(m)),mx(3),mv(3),temp,mtot
!------------------------------------------------------------------------------
    mtot = ZERO;      mx = ZERO;      mv = ZERO
!
    do i = 1, nbig
      xh(:,i) = x(:,i)  +  mx(:)
      vh(:,i) = v(:,i)  +  mv(:)
      mtot = mtot  +  m(i)
      temp = m(i) / (mtot  +  mcen)
      mx(:) = mx(:)  +  temp * x(:,i)
      mv(:) = mv(:)  +  temp * v(:,i)
    enddo
!
    x(:,1:nbig) = xh(:,1:nbig)
    v(:,1:nbig) = vh(:,1:nbig)
!
    end subroutine convert_jacobi_to_helio
!==============================================================================
! For a set of particles, converts barycentric velocities to velocities with 
! respect to the central body, given the masses M.
! The number of active particles is N.
!
    subroutine convert_vbary_to_helio (n,m,v)
!
    use constants;      use globals
    implicit none
    integer(I4), intent(in)::n
    real(R8),    intent(in)::m(:)
    real(R8),    intent(inout)::v(:,:)
!
    real(R8)::mvsum(3)
!------------------------------------------------------------------------------
    mvsum(1) = sum ( m(1:n) * v(1,1:n) ) 
    mvsum(2) = sum ( m(1:n) * v(2,1:n) )
    mvsum(3) = sum ( m(1:n) * v(3,1:n) )
!
    mvsum = mvsum / mcen
!
    v(1,1:n) = v(1,1:n)  +  mvsum(1)
    v(2,1:n) = v(2,1:n)  +  mvsum(2)
    v(3,1:n) = v(3,1:n)  +  mvsum(3)
!
    end subroutine convert_vbary_to_helio
!==============================================================================
! For a set of particles, converts velocities with respect to the central body
! to barycentric velocities, given the masses M.
! The number of active particles is N.
!
    subroutine convert_vhelio_to_bary (n,m,v)
!
    use constants;      use globals
    implicit none
    integer(I4), intent(in)::n
    real(R8),    intent(in)::m(:)
    real(R8),    intent(inout)::v(:,:)
!
    real(R8)::mvsum(3),mtot
!------------------------------------------------------------------------------
    mvsum(1) = sum ( m(1:n) * v(1,1:n) )
    mvsum(2) = sum ( m(1:n) * v(2,1:n) )
    mvsum(3) = sum ( m(1:n) * v(3,1:n) )
!
    mtot = sum(m(1:n))  +  mcen
    mvsum = mvsum / mtot
!
    v(1,1:n) = v(1,1:n)  -  mvsum(1)
    v(2,1:n) = v(2,1:n)  -  mvsum(2)
    v(3,1:n) = v(3,1:n)  -  mvsum(3)
!
    end subroutine convert_vhelio_to_bary
!==============================================================================
! Writes a message to the info file describing a collision that has
! taken place.
    subroutine write_coll_text (t,tname,pname,text,ltext)
!
    use constants;    use globals
    use frag_globals
    use interfaces, only: calc_date
    real(R8),      intent(in)::t
    integer(I4),   intent(in)::ltext
    character(8),  intent(in)::tname,pname
    character(25), intent(in)::text
!
    integer(I4)::month1,year1
    real(R8)::t1
!------------------------------------------------------------------------------
      write (23,*)
!
! Convert time to appropriate output format
    if (opt_time == 1) then
      call calc_date (t/DAY,year1,month1,t1)
      write (23,'(1x,a8,a,a8,a,i10,1x,i2,1x,f4.1)') pname, &
        text(1:ltext),tname,' at ',year1,month1,t1
    else
      if (opt_time == 0) t1 = t / DAY
      if (opt_time == 2) t1 = (t  -  tstart) / DAY
      if (opt_time == 3) t1 = (t  -  tstart) / YEAR
      write (23,'(1x,a8,a,a8,a,1x,f16.3,a)') pname, &
        text(1:ltext),tname,' at ',t1,time_string
    end if
!
    end subroutine write_coll_text
!==============================================================================
! Resolves a collision between particles I and J. Also writes a message
! to the information file.
!
    subroutine collide_bodies (t,i,j,n,nbig,m,x0,v0,x,v,s,ngf,rho,rce_hill,rad, &
      rcrit,status,index,name)
!
    use constants;    use globals
    use frag_globals
    use interfaces, only: calc_radius, merge_bodies
    use frag_interfaces, only: calc_relative_coords, calc_impact_geometry, &
      calc_largest_remnant, calc_second_remnant, fragment_bodies, &
      write_coll_text
    implicit none
    integer(I4),  intent(in)::i,j
    real(R8),     intent(in)::t
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(in)::x0(:,:),v0(:,:) !These are for the impact geometry calculation
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:)
    real(R8),     intent(inout)::rce_hill(:),rad(:),rcrit(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    real(R8), parameter::C1 =  2.43_R8
    real(R8), parameter::C2 = -0.0408_R8
    real(R8), parameter::C3 =  1.86_R8
    real(R8), parameter::C4 =  1.08_R8
!
    integer(I4)::itarg,iproj
    real(R8)::xrel(3),vrel(3),xcom(3),vcom(3),rsum,msum
    real(R8)::b,v2imp,m1,m2,v2esc,v2gm,zeta,fac
    character(25)::text
!------------------------------------------------------------------------------
! Write message to info file
 10 open (23, file=outfile(3), status='old', access='append', err=10)
!
! Check that the less massive particle is treated as the projectile
    itarg = i;    iproj = j
    if (m(j) > m(i).and.j <= nbig) then
      itarg = j;    iproj = i
    end if
!
! If colliding bodies always merge
    if (.not.opt_fragment) then
      text = ' merged with             '
      call write_coll_text (t,name(itarg),name(iproj),text,13)
      call merge_bodies (itarg,iproj,m,x,v,s,status,name)
!
! More realistic collision model based on Leinhardt + Stewart (2012)
    else
      write (23,*)
      write (23,'(2a)') '-------------------------------------------------', &
         '--------------------'
!
      msum = m(i)    +  m(j)
      rsum = rad(i)  +  rad(j)
      v2esc = TWO * G * msum / rsum
!
!..impact parameter, impact velocity and mass of largest remnant
      call calc_relative_coords (m,x0,v0,itarg,iproj,xrel,vrel,xcom,vcom)
      call calc_impact_geometry (xrel,vrel,msum,rsum,b,v2imp)
      m1 = calc_largest_remnant (m(itarg),m(iproj),rad(itarg),rad(iproj),b,v2imp)
!
!..boundary between graze & merge and hit & run (Genda et al. 2012)
      zeta = ((m(itarg)  -  m(iproj)) / msum)**2
      fac = (ONE  -  b / rsum)**2.5_R8
      v2gm = v2esc * (C1 * zeta * fac  +  C2 * zeta  +  C3 * fac  +  C4)
!
      write (23,'(a,f9.4)')   '  Mp / Mt:          ', m(iproj) / m(itarg)
      write (23,'(a,f9.4)')   '  b  / Rtarg:       ', b / rad(itarg)
      write (23,'(a,f9.4)')   '  Vimp / Vesc:      ', sqrt(v2imp / v2esc)
      write (23,'(a,f9.4)')   '  Vgm  / Vesc:      ', sqrt(v2gm  / v2esc)
      write (23,*)
      write (23,'(a,f9.4)')   '  M1 / Msum:        ', m1 / msum
      write (23,'(a,f9.4)')   '  M1 / Mtarg:       ', m1 / m(itarg)
      write (23,'(a,f9.4)')   '  Mfrag / Mfragmin: ', (msum - m1) / mfrag_min
!
! Simple merger
      if (v2imp <= v2esc) then
        text = ' simply merged with      '
        call write_coll_text (t,name(itarg),name(iproj),text,25)
        call merge_bodies (itarg,iproj,m,x,v,s,status,name)
      else
!
! NON-GRAZING REGIME
        if (b < rad(itarg)) then
!
!..effective merger
          if ((msum - m1) < mfrag_min) then
            text = ' effectively merged with '
            call write_coll_text (t,name(itarg),name(iproj),text,25)
            call merge_bodies (itarg,iproj,m,x,v,s,status,name)
!
!..merger with fragmentation
          else
            text = ' head-on smashed into    '
            call write_coll_text (t,name(itarg),name(iproj),text,25)
            m1 = max(m1, mfrag_min)
            m2 = ZERO
            call fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
              rce_hill,rad,rcrit,status,index,name)
          end if
!
! GRAZING REGIME
        else
!
!..graze and merge
          if (v2imp <= v2gm) then
            text = ' grazed and merged with  '
            call write_coll_text (t,name(itarg),name(iproj),text,25)
            call merge_bodies (itarg,iproj,m,x,v,s,status,name)
!
!..hit and run collision
          else if (m1 >= m(itarg)) then
            m1 = m(itarg)
            m2 = calc_second_remnant (m(itarg),m(iproj),rad(itarg),rad(iproj), &
              b,v2imp)
            m2 = max(m2, mfrag_min)
            if ((m(iproj) - m2) < mfrag_min) m2 = m(iproj)
            write (23,'(a,f9.4)')   '  M2 / Mproj:       ', m2 / m(iproj)
            text = ' hit and run with        '
            call write_coll_text (t,name(itarg),name(iproj),text,25)
            call fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
              rce_hill,rad,rcrit,status,index,name)
!
!..erosion with fragments
          else
            m1 = max(m1, mfrag_min)
            m2 = ZERO
            write (23,'(a,f9.4)')   '  M2 / Mproj:       ', m2 / m(iproj)
            text = ' grazing smashed into    '
            call write_coll_text (t,name(itarg),name(iproj),text,25)
            call fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
              rce_hill,rad,rcrit,status,index,name)
          end if
        end if
      end if
    end if
!
! Recompute physical radii
    rad(1:n) = calc_radius (n,m,rho)
!
    close (23)
!
    end subroutine collide_bodies
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
    use constants;    use globals;    use frag_globals
    use interfaces, only: cross_product
    use frag_interfaces, only: calc_relative_coords, &
      calc_largest_remnant
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
    integer(I4)::ifrag,jfrag,nnew,nold,ntheta
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
    write (23,'(3a,es11.4)') '   Remnant:  ',name(itarg), '  m=', m1 / MSUN
    write (23,'(3a,es11.4)') '   Remnant:  ',name(iproj), '  m=', m2 / MSUN
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
      else
        write (name(n)(5:8),'(i4)') nfrag
      end if
      write (23,'(3a,es11.4)') '   Fragment: ',name(n), '  m=', m(n) / MSUN
    end do
!
    write (23,*)
    write (23,*) 'Initial bodies'
    write (23,123) name(itarg), (x(1,itarg) - xcom(1)) / AU, &
      (x(2,itarg) - xcom(2)) / AU, &
      (x(3,itarg) - xcom(3)) / AU, &
      (v(1,itarg) - vcom(1)) / AU * DAY, &
      (v(2,itarg) - vcom(2)) / AU * DAY, &
      (v(3,itarg) - vcom(3)) / AU * DAY, &
      m(itarg) / MSUN * 1d7
    write (23,123) name(iproj), (x(1,iproj) - xcom(1)) / AU, &
      (x(2,iproj) - xcom(2)) / AU, &
      (x(3,iproj) - xcom(3)) / AU, &
      (v(1,iproj) - vcom(1)) / AU * DAY, &
      (v(2,iproj) - vcom(2)) / AU * DAY, &
      (v(3,iproj) - vcom(3)) / AU * DAY, &
      m(iproj) / MSUN * 1d7
 123 format (2x,a,2x,3(1x,f10.7),2x,3(1x,f8.6),3x,f8.4)
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
!
    write (23,*)
    write (23,*) 'Final bodies'
    write (23,123) name(itarg), (x(1,itarg) - xcom(1)) / AU, &
      (x(2,itarg) - xcom(2)) / AU, &
      (x(3,itarg) - xcom(3)) / AU, &
      (v(1,itarg) - vcom(1)) / AU * DAY, &
      (v(2,itarg) - vcom(2)) / AU * DAY, &
      (v(3,itarg) - vcom(3)) / AU * DAY, &
      m(itarg) / MSUN * 1d7
    write (23,123) name(iproj), (x(1,iproj) - xcom(1)) / AU, &
      (x(2,iproj) - xcom(2)) / AU, &
      (x(3,iproj) - xcom(3)) / AU, &
      (v(1,iproj) - vcom(1)) / AU * DAY, &
      (v(2,iproj) - vcom(2)) / AU * DAY, &
      (v(3,iproj) - vcom(3)) / AU * DAY, &
      m(iproj) / MSUN * 1d7
    do ifrag = nold + 1, nold + nnew
      write (23,123) name(ifrag), (x(1,ifrag) - xcom(1)) / AU, &
        (x(2,ifrag) - xcom(2)) / AU, &
        (x(3,ifrag) - xcom(3)) / AU, &
        (v(1,ifrag) - vcom(1)) / AU * DAY, &
        (v(2,ifrag) - vcom(2)) / AU * DAY, &
        (v(3,ifrag) - vcom(3)) / AU * DAY, &
        m(ifrag) / MSUN * 1d7
    end do
!
! Add entries in the list of particle pairs within critical distance
    if (nnew > 0) write (23,*)
!
    do ifrag = nold + 1, nold + nnew
      npair = npair  +  1
      ipair(npair) = itarg
      jpair(npair) = ifrag
      write (23,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
!
      if (m2 /= ZERO) then
        npair = npair  +  1
        ipair(npair) = iproj
        jpair(npair) = ifrag
        write (23,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
      end if
!
      do jfrag = nold + 1, ifrag - 1
        npair = npair  +  1
        ipair(npair) = ifrag
        jpair(npair) = jfrag
        write (23,*) '    New Pair: ',name(ipair(npair)),'  ',name(jpair(npair))
      end do
    end do
!
    end subroutine fragment_bodies
!==============================================================================
! Calculates the impact velocity and impact angle of a collision given the
! relative coordinates and velocities at some time before the collision, 
! assuming the two particles move solely due to their mutual gravity.
!
    subroutine calc_impact_geometry (xrel,vrel,msum,rsum,b,v2imp)
!
    use constants;    use globals
    use interfaces, only: cross_product
    implicit none
    real(R8), intent(in)::xrel(:),vrel(:),msum,rsum
    real(R8), intent(out)::b,v2imp
!
    real(R8)::v2rel,x2rel,x1rel,h(3),h2
!------------------------------------------------------------------------------
! Separation and velocity
    x2rel = dot_product (xrel, xrel)
    v2rel = dot_product (vrel, vrel)
    x1rel = sqrt(x2rel)
!
! Angular momentum
    h = cross_product(xrel, vrel)
    h2 = dot_product(h, h)
!
! Relative velocity at impact, using vis-viva equation
    v2imp = v2rel  +  TWO * G * msum * (ONE / rsum  -  ONE / x1rel)
!
! Impact parameter
    b = sqrt(h2 / v2imp) !This is really B = b*R
!
    end subroutine calc_impact_geometry
!==============================================================================
! Calculates the mass of the largest remnant from a collision between a
! target of mass MTARG and a projectile of mass MPROJ, given the impact 
! parameter B and impact velocity squared V2IMP.
! The formula is taken from Leinhardt and Stewart 2011.
!
    function calc_largest_remnant (mtarg,mproj,rtarg,rproj,b,v2imp)
!
    use constants;    use globals;    use frag_globals
    use frag_interfaces, only: calc_qstar, calc_remnant_mass
    implicit none
    real(R8), intent(in)::mtarg,mproj,rtarg,rproj,b,v2imp
    real(R8)::calc_largest_remnant
!
    real(R8)::l,alpha,msum,q,qstar
!------------------------------------------------------------------------------
! Fraction of projectile that intersects the target
    l = rtarg  +  rproj  -  b
    l = min(l, TWO * rproj)
    alpha = l * l * (THREE * rproj  -  l) / (FOUR * rproj * rproj * rproj)
    alpha = min(ONE, alpha)
!
! Impact energy per unit mass
    msum = mtarg  +  mproj
    q = HALF * v2imp * mtarg * mproj / (msum * msum)
!
! Critical impact energy per unit mass
    qstar = calc_qstar (mtarg, mproj, alpha)
!
! Mass of escaping projectile
    calc_largest_remnant = calc_remnant_mass (mtarg,mproj,q,qstar)
!
    end function calc_largest_remnant
!==============================================================================
! Calculates the mass of the second remnant (escaping projectile) for a hit
! and run collision between a target of mass MTARG and a projectile of mass
! MPROJ, given the impact parameter B and impact velocity squared V2IMP.
! The formula is taken from Leinhardt and Stewart 2011.
!
    function calc_second_remnant (mtarg,mproj,rtarg,rproj,b,v2imp)
!
    use constants;    use globals;    use frag_globals
    use frag_interfaces, only: calc_qstar, calc_remnant_mass
    implicit none
    real(R8), intent(in)::mtarg,mproj,rtarg,rproj,b,v2imp
    real(R8)::calc_second_remnant
!
    real(R8)::l,phi,biga,bigl,mint,temp,msum,q,qstar
!------------------------------------------------------------------------------
! Mass MINT of target that interacts with the projectile
    l = rtarg  +  rproj  -  b
    phi = TWO * acos((l - rproj) / rproj)
    biga = rproj * rproj * (PI  -  HALF * (phi - sin(phi)))
    bigl = TWO * sqrt(rtarg * rtarg  -  (rtarg  - HALF * l)**2)
    temp = biga * bigl / (FOUR * THIRD * PI * rtarg * rtarg * rtarg)
    mint = mtarg * min(temp, ONE)
!
! Impact energy per unit mass
    msum = mint  +  mproj
    q = HALF * v2imp * mint * mproj / (msum * msum)
!
! Critical impact energy per unit mass
    qstar = calc_qstar (mproj, mint, ONE)
!
! Mass of escaping projectile
    calc_second_remnant = calc_remnant_mass (mproj,mint,q,qstar)
    calc_second_remnant = min(mproj, calc_second_remnant)
!
    end function calc_second_remnant
!==============================================================================
! Calculates the critical impact energy per unit mass Qstar for a
! collision between a target MTARG and a projectile MPROJ given the
! fraction ALPHA of projectile mass that overlaps the target, following
! Leinhardt and Stewart (2012).
!
    function calc_qstar (mtarg,mproj,alpha)
!
    use constants;    use globals
    implicit none
    real(R8), parameter::CSTAR = 1.8_R8
    real(R8), intent(in)::mtarg,mproj,alpha
    real(R8)::calc_qstar
!
    real(R8)::gamma,msum,mu,mu_alpha,rc1,gamma1,temp,fac
!------------------------------------------------------------------------------
    if (alpha == ZERO) then
      calc_qstar = BIG_NUMBER
      return
    end if
!
! Total mass and mass ratio
    msum = mtarg  +  mproj
    gamma = mproj / mtarg
!
! Reduced mass and reduced interacting mass
    mu = mtarg * mproj / msum
    mu_alpha = alpha * mtarg * mproj / (mtarg  +  alpha * mproj)
!
! Radius of total mass at density = 1 g/cm^3
    rc1 = (THREE * msum / (FOUR * PI)) ** THIRD
!
! Modification factor for interacting mass fraction and mass ratio
    gamma1 = gamma  +  ONE
    temp = mu / mu_alpha
    fac = temp * sqrt(temp) * 0.25_R8 * gamma1 * gamma1 / gamma
!
! Qstar
    calc_qstar = fac * CSTAR * 0.8_R8 * PI * G * rc1 * rc1
!
    end function calc_qstar
!==============================================================================
! Calculates the mass of the largest remnant of a collision between a target
! mass MTARG and a projectile mass MPROJ, given the impact energy per unit
! mass Q and its critical value QSTAR.
! The routine can also be used to calculate the mass of the second remnant
! (escaping projectile) for a hit and run collision.
!
    function calc_remnant_mass (mtarg,mproj,q,qstar)
!
    use constants;    use globals
    real(R8), intent(in)::mtarg,mproj,q,qstar
    real(R8)::calc_remnant_mass
!
    real(R8)::qratio,msum,fac
!------------------------------------------------------------------------------
    msum = mtarg  +  mproj
    qratio = q / qstar
!
    if (qratio < 1.8_R8) then
      calc_remnant_mass = msum * (ONE  -  HALF * qratio)
    else
      fac = qratio / 1.8_R8
      calc_remnant_mass = msum * 0.1_R8 * fac**(-1.5_R8)
    end if
!
    end function calc_remnant_mass
!==============================================================================
! Does an elastic bounce between a target ITARG and a projectile IPROJ.
! The coordinates remain the same but the radial velocity changes sign.
!
    subroutine bounce_bodies (itarg,iproj,m,x,v,s,rad,name)
!
    use constants;    use globals
    use frag_interfaces, only: calc_relative_coords, calc_absolute_coords, &
      calc_spherical_velocities, calc_cartesian_velocities
    use interfaces, only: calc_date
    implicit none
    integer(I4), intent(in)::itarg,iproj
    real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),rad(:)
    character(8), intent(inout)::name(:)
!
    real(R8)::xrel(3),vrel(3),xcom(3),vcom(3),vr,vtheta,vphi,gmsum,dt
!------------------------------------------------------------------------------
! Calculate centre of mass and relative coordinates
    call calc_relative_coords (m,x,v,itarg,iproj,xrel,vrel,xcom,vcom)
!
! Move back in time before the collision
    gmsum = G * (m(itarg)  +  m(iproj))
    dt = TWO * sqrt(dot_product(xrel, xrel) / dot_product(vrel, vrel))
    call drift (-dt,gmsum,xrel,vrel)
!
! Calculate change in relative velocity due to bounce
    call calc_spherical_velocities (xrel,vrel,vr,vtheta,vphi)
    vr = -vr
    call calc_cartesian_velocities (xrel,vrel,vr,vtheta,vphi)
!
! Move forward in time again
    call drift (dt,gmsum,xrel,vrel)
!
! Calculate new absolute velocities
    call calc_absolute_coords (xrel,vrel,xcom,vcom,itarg,iproj,m,x,v)
!
    end subroutine bounce_bodies
!==============================================================================
! Calculates spherical velocity components given the Cartesian coordinates 
! and velocity components.
!
    subroutine calc_spherical_velocities (xrel,vrel,vr,vtheta,vphi)
!
    use constants;    use globals
    implicit none
    real(R8), intent(in)::xrel(:),vrel(:)
    real(R8), intent(out)::vr,vtheta,vphi
!
    real(R8)::r,p
!------------------------------------------------------------------------------
    r = sqrt(dot_product(xrel, xrel))
    p = sqrt(xrel(1) * xrel(1)  +  xrel(2) * xrel(2))
!
    vr     = dot_product(xrel, vrel) / r
    vtheta = p * vrel(3) / r &
           - xrel(3) * (xrel(1) * vrel(1)  +  xrel(2) * vrel(2)) / (r * p)
    vphi   = (xrel(1) * vrel(2)  -  xrel(2) * vrel(1)) / p
!
    end subroutine calc_spherical_velocities
!==============================================================================
! Calculates Cartesian velocity components given the Cartesian coordinates 
! and spherical velocity components.
!
    subroutine calc_cartesian_velocities (xrel,vrel,vr,vtheta,vphi)
!
    use constants;    use globals
    implicit none
    real(R8), intent(in)::xrel(:),vr,vtheta,vphi
    real(R8), intent(out)::vrel(:)
!
    real(R8)::r,p,temp
!------------------------------------------------------------------------------
    r = sqrt(dot_product(xrel, xrel))
    p = sqrt(xrel(1) * xrel(1)  +  xrel(2) * xrel(2))
!
    temp = p * vr  -  xrel(3) * vtheta
!
    vrel(1) = xrel(1) * temp / (r * p)  -  xrel(2) * vphi / p
    vrel(2) = xrel(2) * temp / (r * p)  +  xrel(1) * vphi / p
    vrel(3) = (xrel(3) * vr  +  p * vtheta) / r
!
    end subroutine calc_cartesian_velocities
!==============================================================================
! Calculates the absolute coordinates and velocities of particles I and J
! given their centre of mass and relative coordinates and velocities.
!
    subroutine calc_absolute_coords (xrel,vrel,xcom,vcom,i,j,m,x,v)
!
    use constants;    use globals
    implicit none
    integer(I4), intent(in)::i,j
    real(R8), intent(in)::xrel(:),vrel(:),xcom(:),vcom(:),m(:)
    real(R8), intent(out)::x(:,:),v(:,:)
!
    real(R8)::msum
!------------------------------------------------------------------------------
    msum = m(i)  +  m(j)
!
   x(:,i) = xcom(:)  +  m(j) * xrel(:) / msum
   x(:,j) = xcom(:)  -  m(i) * xrel(:) / msum
!
   v(:,i) = vcom(:)  +  m(j) * vrel(:) / msum
   v(:,j) = vcom(:)  -  m(i) * vrel(:) / msum
!
    end subroutine calc_absolute_coords
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
!==============================================================================
! Resolves a collision between target ITARG and projectile IPROJ
! that results in the formation of a single large remnant of mass M1 
! plus some equal-mass fragments.
! Fragments are assumed to be ejected at slightly above escape velocity in
! uniform directions within a plane centred on the centre of mass.
!
    subroutine old_fragment_bodies (itarg,iproj,m1,m2,n,nbig,m,x,v,s,ngf,rho, &
      rce_hill,rad,rcrit,status,index,name)
!
    use constants;    use globals;    use frag_globals
    use interfaces, only: cross_product
    use frag_interfaces, only: calc_relative_coords, &
      calc_largest_remnant
    implicit none
    integer(I4), intent(in)::itarg,iproj
    real(R8),    intent(in)::m1,m2
    integer(I4), intent(inout)::n,nbig,index(:)
    real(R8), intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
    real(R8), intent(inout)::rho(:),rce_hill(:),rad(:)
    real(R8), intent(inout)::rcrit(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::ifrag,jfrag,nnew,nold,ntheta
    real(R8)::msum,rsum,ssum(3),m_ejecta,mfrag,r,v2esc,v2,v1,theta
    real(R8)::xcom(3),vcom(3),xrel(3),vrel(3),en0,en1,dx(3),dv(3),mredu
    real(R8)::mxsum(3),mvsum(3),xoff(3),voff(3)
!------------------------------------------------------------------------------
! Save initial number of bodies
    nold = n
!
    msum = m(itarg)    +  m(iproj)
    rsum = rad(itarg)  +  rad(iproj)
    ssum(:) = s(:,itarg)  +  s(:,iproj)
!
    write (23,'(3a,es11.4)') '   Remnant:  ',name(itarg), '  m=', m1 / MSUN
    write (23,'(3a,es11.4)') '   Remnant:  ',name(iproj), '  m=', m2 / MSUN
!
! Calculate centre of mass of the target and projectile
    call calc_relative_coords (m,x,v,itarg,iproj,xrel,vrel,xcom,vcom)
!
! Reduced mass of target and projectile
    mredu = m(itarg) * m(iproj) / msum
!
! Initial energy of the target and projectile (center of mass frame)
    en0 = HALF * m(itarg) * m(iproj) * dot_product(vrel, vrel) / msum &
        - G * m(itarg) * m(iproj) / sqrt(dot_product(xrel, xrel))
!
! Calculate number of fragments NNEW and individual fragment mass MFRAG
    m_ejecta = msum  -  (m1  +  m2)
    nnew  = int(m_ejecta / mfrag_min)
    mfrag = m_ejecta / dble(nnew)
!
! Distance R and velocity V1 of fragments with respect to centre of mass
    r = rsum * FOUR
    v2esc = TWO * G * msum / rsum
    v2 = 1.1_R8 * v2esc  -  TWO * G * msum * (ONE / rsum  -  ONE / r)
    v1 = sqrt(v2)
!
! Number NTHETA of particles apart from the largest remnant
    ntheta = nnew
    if (m2 > ZERO) ntheta = nnew  +  1
    mxsum = ZERO;      mvsum = ZERO
!
! Create each new fragment
    do ifrag = nold + 1, nold + nnew
      n = n  +  1
      nbig = nbig  +  1
      nfrag = nfrag  +  1
!
! Set up mass, spin, density, radius etc
      m(n)        = mfrag
      s(:,n)      = ZERO
      rho(n)      = rho(iproj)
      rce_hill(n) = rce_hill(iproj)
      rcrit(n)    = rcrit(iproj)
      rad(n)      = (THREE * m(n) / (FOUR * PI * rho(n)))**THIRD
      index(n)    = 0
      status(n)   = 'big  '
!
! Choose a unique name for the fragment
      name(n) = 'FRAG0000'
      if (nfrag < 10) then
        write (name(n)(8:8),'(i1)') nfrag
      else if (nfrag < 100) then
        write (name(n)(7:8),'(i2)') nfrag
      else if (nfrag < 1000) then
        write (name(n)(6:8),'(i3)') nfrag
      else
        write (name(n)(5:8),'(i4)') nfrag
      end if
      write (23,'(3a,es11.4)') '   Fragment: ',name(n), '  m=', m(n) / MSUN
!
! Arrange the fragments in a plane symmetrically about the centre of mass
      theta = TWOPI * dble(ifrag) / dble(ntheta)
      x(1,n) = xcom(1)  +  r * cos(theta)
      x(2,n) = xcom(2)  +  r * sin(theta)
      x(3,n) = xcom(3)
      mxsum(:) = mxsum(:)  +  m(n) * x(:,n)
!
      v(1,n) = vcom(1)  +  v1 * cos(theta)
      v(2,n) = vcom(2)  +  v1 * sin(theta)
      v(3,n) = vcom(3)
      mvsum(:) = mvsum(:)  +  m(n) * v(:,n)
!
! Add an entry in the list of particle pairs within critical distance
      npair = npair  +  1
      ipair(npair) = itarg
      jpair(npair) = n
!
      do jfrag = nold + 1, ifrag - 1
        npair = npair  +  1
        ipair(npair) = ifrag
        jpair(npair) = jfrag
      end do
    end do
!
    write (23,*)
    write (23,*) 'Initial bodies'
    write (23,123) name(itarg), x(1,itarg)/AU, x(2,itarg)/AU, &
      v(1,itarg)/AU*DAY, v(2,itarg)/AU*DAY, m(itarg)/MSUN*1d7
    write (23,123) name(iproj), x(1,iproj)/AU, x(2,iproj)/AU, &
      v(1,iproj)/AU*DAY, v(2,iproj)/AU*DAY, m(iproj)/MSUN*1d7
 123 format (2x,a,2x,2(1x,f10.7),2x,2(1x,f8.6),3x,f8.4)
!
! New mass and spin of the target. Initially, place at centre of mass.
    x(:,itarg) = xcom(:)
    v(:,itarg) = vcom(:)
    s(:,itarg) = s(:,itarg)  +  mredu * cross_product(xrel,vrel)
    if (m2 <= ZERO) s(:,itarg) = s(:,itarg)  +  s(:,iproj)
    m(itarg) = m1
    mxsum(:) = mxsum(:)  +  m1 * x(:,itarg)
    mvsum(:) = mvsum(:)  +  m1 * v(:,itarg)
!
! Update the projectile (if necessary, flag it for removal and move away
! from the collision site.
    if (m2 > ZERO) then
      m(iproj) = m2
      x(1,iproj) = xcom(1)  +  r
      x(2,iproj) = xcom(2)
      x(3,iproj) = xcom(3)
!
      v(1,iproj) = vcom(1)  +  v1
      v(2,iproj) = vcom(2)
      v(3,iproj) = vcom(3)
    else
      m(iproj) = ZERO;             x(:,iproj) = -x(:,iproj)
      v(:,iproj) = -v(:,iproj);    s(:,iproj) = ZERO
      status(iproj) = 'dead '
    end if
    mxsum(:) = mxsum(:)  +  m2 * x(:,iproj)
    mvsum(:) = mvsum(:)  +  m2 * v(:,iproj)
!
! Adjust all coords and velocities to conserve centre of mass and momentum
    xoff(:) = xcom - mxsum(:) / msum
    voff(:) = vcom - mvsum(:) / msum
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
      do jfrag = ifrag + 1, n
        dx(:) = x(:,ifrag)  -  x(:,jfrag)
        en1 = en1  -  G * mfrag * mfrag / sqrt(dot_product(dx, dx))
      end do
    end do
!
! Calculate energy loss due to the collision
    denergy = denergy  +  (en0  -  en1)
!
    write (23,*)
    write (23,*) 'xcom, vcom:'
    write (23,123) 'xvcom   ', xcom(1)/AU, xcom(2)/AU, &
      vcom(1)/AU*DAY, vcom(2)/AU*DAY
!
    write (23,*)
    write (23,*) 'Final bodies'
    write (23,123) name(itarg), x(1,itarg)/AU, x(2,itarg)/AU, &
      v(1,itarg)/AU*DAY, v(2,itarg)/AU*DAY, m(itarg)/MSUN*1d7
    write (23,123) name(iproj), x(1,iproj)/AU, x(2,iproj)/AU, &
      v(1,iproj)/AU*DAY, v(2,iproj)/AU*DAY, m(iproj)/MSUN*1d7
    do ifrag = nold + 1, nold + nnew
      write (23,123) name(ifrag), x(1,ifrag)/AU, x(2,ifrag)/AU, &
        v(1,ifrag)/AU*DAY, v(2,ifrag)/AU*DAY, m(ifrag)/MSUN*1d7
    end do
!
    end subroutine old_fragment_bodies
!==============================================================================
! Merges all particles that collide with the central body into the central
! body. Also removes the lost particles from the global arrays and calculates
! the energy change of the system.
!
    subroutine collide_central (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index, &
      name,nhit,hit)
!
    use constants;    use globals
    use interfaces, only: calc_date, calc_energy, cross_product, remove_dead_bodies
    implicit none
    integer(I4),  intent(in)::nhit
    type(encounter), intent(in)::hit(:)
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:)
    real(R8),     intent(inout)::ngf(:,:),rho(:),rce_hill(:),rad(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::i,k,yr,month
    real(R8)::t,t1,dxcen(3),dvcen(3),xrel(3),vrel(3),msum,mredu,energy1,energy2,angmom
!------------------------------------------------------------------------------
    if (nhit <= 0) return
!
! Calculate initial energy
    call calc_energy (n,nbig,m,x,v,s,energy1,angmom)
!
! Open info file
 10 open (23, file=outfile(3), status='old', access='append', err=10)
!
    do k = 1, nhit
      i = hit(k) % j;      t = hit(k) % t
!
! Write message to info file
      write (23,*)
      if (opt_time == 1) then
        call calc_date (t/DAY,yr,month,t1)
        write (23,'(1x,a8,a,i10,1x,i2,1x,f8.5)') name(i), &
          ' collided with the central body at ',yr,month,t1
      else
        if (opt_time == 0) t1 = t / DAY
        if (opt_time == 2) t1 = (t  -  tstart) / DAY
        if (opt_time == 3) t1 = (t  -  tstart) / YEAR
        write (23,'(1x,a8,a,f20.7,a)') name(i), &
          ' collided with the central body at ',t1,time_string
      end if
!
! Merge the body with the central body
      msum = mcen  +  m(i)
      mredu = mcen * m(i) / msum
      xrel(:) = x(:,i)
      vrel(:) = v(:,i)
!
! Shift the heliocentric coordinates and velocities of all bodies
      dxcen(:) = m(i) * x(:,i) / msum
      dvcen(:) = m(i) * v(:,i) / msum
      x(1,1:n) = x(1,1:n)  -  dxcen(1)
      x(2,1:n) = x(2,1:n)  -  dxcen(2)
      x(3,1:n) = x(3,1:n)  -  dxcen(3)
      v(1,1:n) = v(1,1:n)  -  dvcen(1)
      v(2,1:n) = v(2,1:n)  -  dvcen(2)
      v(3,1:n) = v(3,1:n)  -  dvcen(3)
!
! New mass and spin angular momentum of the central body
      mcen = msum
      scen(:) = scen(:)  +  s(:,i)  +  mredu * cross_product(xrel, vrel)
!
! Flag the lost body for removal
      m(i) = ZERO
      s(:,i) = ZERO
      status(i) = 'dead '
    end do
!
    close (23)
!
! Remove ejected objects from global arrays, update energy and angular momentum
    if (nhit > 0) then
      call remove_dead_bodies (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
      call calc_energy (n,nbig,m,x,v,s,energy2,angmom)
      denergy = denergy  +  (energy1  -  energy2)
    end if
!
    end subroutine collide_central
!==============================================================================
! Applies a symplectic corrector for a time DT to a set of particles.
!
    subroutine corrector (dt,n,nbig,m,x,v,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_mvs, drift, &
      convert_helio_to_jacobi, convert_jacobi_to_helio
    implicit none
    real(R8),parameter::HA(2) = (/-ROOT10 /  5.0_R8, -ROOT10 * 3.0_R8 / 10.0_R8/)
    real(R8),parameter::HB(2) = (/-ROOT10 / 24.0_R8,  ROOT10 / 72.0_R8/)
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::dt,m(:),ngf(:,:)
    real(R8),    intent(inout)::x(:,:),v(:,:)
!
    integer(I4)::i,k
    real(R8)::gm(NMAX),a(3,NMAX),minside,msofar,temp
!------------------------------------------------------------------------------
! Calculate effective central masses for Kepler drifts
    minside = mcen
    do i = 1, nbig
      msofar = minside  +  m(i)
      gm(i) = G * mcen * msofar / minside
      minside = msofar
    end do
    gm(nbig+1:n) = G * mcen
!
    do k = 1, 2
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call convert_helio_to_jacobi (n,nbig,m,x,v)
      temp = dt * HA(k)
      do i = 1, n
        call drift (temp,gm(i),x(:,i),v(:,i))
      end do
      call convert_jacobi_to_helio (n,nbig,m,x,v)
!
! Advance Interaction Hamiltonian
      a(:,1:n) = accel_mvs (n,nbig,m,x,v,ngf)
      temp = dt * HB(k)
      v(:,1:n) = v(:,1:n)  +  temp * a(:,1:n)
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call convert_helio_to_jacobi (n,nbig,m,x,v)
      temp = -2.0_R8 * dt * HA(k)
      do i = 1, n
        call drift (temp,gm(i),x(:,i),v(:,i))
      end do
      call convert_jacobi_to_helio (n,nbig,m,x,v)
!
! Advance Interaction Hamiltonian
      a(:,1:n) = accel_mvs (n,nbig,m,x,v,ngf)
      temp = -dt * HB(k)
      v(:,1:n) = v(:,1:n)  +  temp * a(:,1:n)
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call convert_helio_to_jacobi (n,nbig,m,x,v)
      temp = dt * HA(k)
      do i = 1, n
        call drift (temp,gm(i),x(:,i),v(:,i))
      end do
      call convert_jacobi_to_helio (n,nbig,m,x,v)
!
    end do
!
    end subroutine corrector
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
! Advances the Cartesian coordinates and velocities of a single particle along
! a Keplerian orbit for a time DT. The gravitational constant times the central
! mass plus the particle mass is GMSUM.
!
    subroutine drift (dt,gmsum,x,v)
    use constants
    implicit none
!
    real(R8), intent(in)::dt,gmsum
    real(R8), intent(inout)::x(3),v(3)
    integer(I4)::k,flag
    real(R8)::dttmp
!------------------------------------------------------------------------------
    call drift_dan (gmsum,x(1),x(2),x(3),v(1),v(2),v(3),dt,flag)
!
! If error is too large, try again with a smaller stepsize
    if(flag /= 0) then
      dttmp = dt * 0.1_R8
      do k = 1, 10
        call drift_dan (gmsum,x(1),x(2),x(3),v(1),v(2),v(3),dttmp,flag)
        if(flag /= 0) return
      end do
    endif
!
    end subroutine drift
!==============================================================================
! Does an integration using the hybrid symplectic algorithm.
!
    subroutine driver_hybrid (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    use constants;    use globals
    use interfaces, only: accel_hybrid, bs2, calc_radius, calc_rcrit, calc_rce, &
      check_central, check_critical, check_encounters, collide_bodies, &
      collide_central, convert_vhelio_to_bary, convert_vbary_to_helio, &
      drift, ejections, jump, output_coords, output_datadump, output_encounters, &
      output_progress, remove_dead_bodies, output_codes
    implicit none
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::i,j,k,nclo,nhit,ncrit,ncrit_big,ncrit_old
    real(R8)::dt,dtby2,gmcen,t,temp
    real(R8)::a(3,NMAX),rcrit(NMAX),rad(NMAX),x0(3,NMAX),v0(3,NMAX)
    type(encounter)::clo(NMAX),hit(NMAX)
    logical::flag_collision
!
! Variables used by Bulirsch-Stoer integrator
    integer(I4)::indexbs(NMAX)
    real(R8)::time_bs,dt_bs,dt_bs_nextcall
    real(R8)::mbs(NMAX),xbs(3,NMAX),vbs(3,NMAX),sbs(3,NMAX),ngfbs(3,NMAX)
    real(R8)::rhobs(NMAX),rcebs(NMAX),radbs(NMAX),rcritbs(NMAX)
    character(8)::namebs(NMAX)
    character(5)::statusbs(NMAX)
!------------------------------------------------------------------------------
! Initialize variables
    flag_accel = .true.
    tdump = time;      tfun  = time;      tlog = time
!
! Make sure the stepsize has the correct sign
    if (goal == 'start') dt = sign(dt0, tstart - time)
    if (goal == 'stop ') dt = sign(dt0, tstop  - time)
    dtby2 = dt * HALF
    dt_bs_nextcall = dt
!
! Calculate physical radii, close encounter, and critical encounter distances
    rad(1:n) = calc_radius (n,m,rho)
    rcrit(1:n) = calc_rcrit  (dt,n,m,x,v)
!
    do
      flag_stop = .false.
!
! See if we have reached the start time for the integration
      if (goal == 'start'.and.abs(time - tstart) <= HALF * dt0) then
        goal ='stop '
        dt = sign(dt0, tstop - time)
        dtby2 = dt * HALF
        dt_bs_nextcall = dt
        tout = tstart  -  sign(dtout, dt)
        flag_accel = .true.
      end if
!
! Output the coords and velocities if it is time
      if (goal == 'stop '.and.abs(time - tout) >= dtout - HALF * dt0) then
        call output_coords   (n,nbig,m,x,v,s,rho,status,index,name)
        call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
        tout  = tout   +  sign(dtout,  dt)
        tdump = tdump  +  sign(dtdump, dt)
      end if
!
! See if the integration has finished
      if (goal == 'stop '.and.abs(time - tstop) <= HALF * dt0) exit
!
!------------------------------------------------------------------------------
!  ADVANCE  ONE  TIMESTEP
!
! If necessary, recalculate accelerations
      if (flag_accel) then
        a(:,1:n) = accel_hybrid (n,nbig,m,x,v,ngf,rcrit)
        flag_accel = .false.
      end if
!
! Advance interaction Hamiltonian for half a timestep
      v(:,1:n) = v(:,1:n)  +  dtby2 * a(:,1:n)
!
! Advance jump Hamiltonian for half a timestep
      call convert_vhelio_to_bary (n,m,v)
      call jump (dtby2,n,m,x,v)
!
! Save the current coordinates and velocities
      do
        x0(:,1:n) = x(:,1:n)
        v0(:,1:n) = v(:,1:n)
!
! Advance Keplerian Hamiltonian for a full timestep
        gmcen = G * mcen
        do i = 1, n
          call drift (dt,gmcen,x(:,i),v(:,i))
        end do
!
! Check for collisions with the central body
        call check_central (dt,n,nbig,m,x0,v0,x,v,nhit,hit)
        if (nhit == 0) exit
!
! If particles hit the central body, restore initial coords and merge
        x(:,1:n) = x0(:,1:n)
        v(:,1:n) = v0(:,1:n)
        call convert_vbary_to_helio (n,m,v)
        call collide_central (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index, &
          name,nhit,hit)
        call convert_vhelio_to_bary (n,m,v)
        flag_accel = .true.
      end do
!
!------------------------------------------------------------------------------
!  ADVANCE  SOME  BODIES  USING  BULIRSCH-STOER
!
! Check whether any bodies passed within RCRIT of one another during Keplerian step
      call check_critical (dt,n,nbig,x0,v0,x,v,rcrit,ncrit,ncrit_big,name)
      flag_collision = .false.
!
! If some bodies came within RCRIT, advance these bodies using Bulirsch-Stoer
      ncrit_old = ncrit
      if (ncrit > 0) then
!
! Put data for critical bodies into local arrays for use with BS routine
        do k = 1, ncrit
          i = icrit(k)
          mbs(k)     = m(i);          sbs(:,k)   = s(:,i)
          xbs(:,k)   = x0(:,i);       vbs(:,k)   = v0(:,i)
          ngfbs(:,k) = ngf(:,i);      rcebs(k)   = rce_hill(i)
          radbs(k)   = rad(i);        rcritbs(k) = rcrit(i)
          rhobs(k)   = rho(i);        statusbs(k)= status(i)
          namebs(k)  = name(i);       indexbs(k) = index(i)
        end do
        time_bs = ZERO
!
! Do the Bulirsch-Stoer integration
        do
          dt_bs = sign(dt_bs_nextcall, dt)
          if (abs(time_bs + dt_bs) > abs(dt)) dt_bs = sign(abs(dt) - abs(time_bs), dt)
!
! Save the current coordinates and velocities
          x0(:,1:ncrit) = xbs(:,1:ncrit)
          v0(:,1:ncrit) = vbs(:,1:ncrit)
!
! Advance all critical bodies
          call bs2 (dt_bs,dt_bs_nextcall,ncrit,ncrit_big,mbs,xbs,vbs,ngfbs,rcritbs,1)
          time_bs = time_bs  +  dt_bs
!
! Check for close-encounter minima and collisions
          temp = time  +  time_bs
          call check_encounters (temp,dt_bs,ncrit,ncrit_big,mbs,x0,v0,xbs,vbs,radbs, &
            rcebs,namebs,nclo,clo,nhit,hit)
          if (nhit > 0) call output_encounters (ncrit,ncrit_big,rhobs,statusbs, &
              indexbs,namebs,nhit,hit,m)
!
! Resolve any collisions that occurred
          if (nhit > 0.and.opt_collisions) then
            do k = 1, nhit
              i = hit(k) % i;      j = hit(k) % j;      t = hit(k) % t
              call collide_bodies (t,i,j,ncrit,ncrit_big,mbs,x0,v0,xbs,vbs,sbs,ngfbs, &
                rhobs,rcebs,radbs,rcritbs,statusbs,indexbs,namebs)
              flag_collision = .true.
              write(*,*) "Number of big bodies at time ", t/YEAR, " years is: ", nbig
            end do
          end if
!
! If we have integrated for long enough, stop the Bulirsch-Stoer integration
          if (abs(time_bs) >= abs(dt)) exit
        end do
!
! Return data for the critical bodies to global arrays
        do k = 1, ncrit
          if (k <= ncrit_old) then
            i = icrit(k)
          else
            n = n  +  1
            nbig = nbig  +  1
            i = n
          end if
!
          m(i)     = mbs(k);          s(:,i)      = sbs(:,k)
          x(:,i)   = xbs(:,k);        v(:,i)      = vbs(:,k)
          ngf(:,i) = ngfbs(:,k);      rce_hill(i) = rcebs(k)
          rad(i)   = radbs(k);        rcrit(i)    = rcritbs(k)
          rho(i)   = rhobs(k);        status(i)   = statusbs(k)
          name(i)  = namebs(k);       index(i)    = indexbs(k)
        end do
        do k = 1, nclo
          clo(k) % i = icrit(clo(k) % i)
          clo(k) % j = icrit(clo(k) % j)
        end do
      end if
!------------------------------------------------------------------------------
!  CONTINUE  THE  NORMAL  TIME  STEP
!
! Advance jump Hamiltonian for half a timestep
      call jump (dtby2,n,m,x,v)
!
! Convert to heliocentric velocities
      call convert_vbary_to_helio (n,m,v)
!
! Advance interaction Hamiltonian for half a timestep
      a(:,1:n) = accel_hybrid (n,nbig,m,x,v,ngf,rcrit)
      v(:,1:n) = v(:,1:n)  +  dtby2 * a(:,1:n)
      time = time  +  dt
!
! Remove any particles lost during collisions
      if (flag_collision) then
        call remove_dead_bodies (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
        call output_codes (n,nbig,status,index,name,m)
        rcrit(1:n) = calc_rcrit  (dt,n,m,x,v)
        flag_accel = .true.
      end if
!
! Remove particles far from the central body, update close encounter limits
      if (abs(time - tfun) >= dtfun) then
        call ejections (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
        if (flag_ejections) then
          rcrit(1:n) = calc_rcrit (dt,n,m,x,v)
          flag_accel = .true.
        end if
        tfun  = tfun   +  sign(dtfun, dt)
      end if
!
! Write a progress report and update data dump files
      if (abs(time - tlog) >= dtdump) then
        call output_progress (n,nbig,m,x,v,s)
        call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
        tlog  = tlog   +  sign(dtdump, dt)
        tdump = tdump  +  sign(dtdump, dt)
      end if
!
! Stop the integration for any reason?
      if (flag_stop) return

! Go on to the next time step
    end do
!
    end subroutine driver_hybrid
!==============================================================================
! Does an integration using the mixed variable symplectic algorithm.
!
    subroutine driver_mvs (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    use constants;    use globals
    use interfaces, only: accel_mvs, calc_radius, calc_rce, &
      check_central, check_encounters, ejections, collide_central, &
      convert_helio_to_jacobi, convert_jacobi_to_helio, &
      corrector, inverse_corrector, drift, &
      output_coords, output_datadump, output_encounters, output_progress, &
      remove_dead_bodies
    implicit none
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::i,nclo,nhit
    real(R8)::dt,dtby2,gmcen,minside,msofar
    real(R8)::a(3,NMAX),rad(NMAX),x0(3,NMAX),v0(3,NMAX),gm(NMAX)
    type(encounter)::clo(NMAX),hit(NMAX)
!------------------------------------------------------------------------------
! Initialize variables
    gmcen = G * mcen;    flag_accel = .true.
    tdump = time;      tfun  = time;      tlog = time
    if (goal == 'start') dt = sign(dt0, tstart - time)
    if (goal == 'stop ') dt = sign(dt0, tstop  - time)
    dtby2 = dt * HALF
!
! Calculate physical radii and close encounter distances
    rad(1:n) = calc_radius (n,m,rho)
!
! Apply symplectic corrector if necessary
    if (algor == 9) call corrector (dt,n,nbig,m,x,v,ngf)
!
    do
      flag_stop = .false.
!
! See if we have reached the start time for the integration
      if (goal == 'start'.and.abs(time - tstart) <= HALF * dt0) then
        goal ='stop '
        dt = sign(dt0, tstop - time)
        dtby2 = dt * HALF
        tout = tstart  -  sign(dtout, dt)
        flag_accel = .true.
      end if
!
! Output the coords and velocities if it is time
      if (goal == 'stop '.and.abs(time - tout) >= dtout - HALF * dt0) then
        if (algor == 9) call inverse_corrector (dt,n,nbig,m,x,v,ngf)
        call output_coords   (n,nbig,m,x,v,s,rho,status,index,name)
        call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
        if (algor == 9) call corrector (dt,n,nbig,m,x,v,ngf)
        tout  = tout   +  sign(dtout,  dt)
        tdump = tdump  +  sign(dtdump, dt)
      end if
!
! See if the integration has finished
      if (goal == 'stop '.and.abs(time - tstop) <= HALF * dt0) then
        if (algor == 9) call inverse_corrector (dt,n,nbig,m,x,v,ngf)
        exit
      end if
!
!------------------------------------------------------------------------------
!  ADVANCE  ONE  TIMESTEP
!
! Save the current coordinates and velocities
      x0(:,1:n) = x(:,1:n)
      v0(:,1:n) = v(:,1:n)
!
      do
!
! If necessary, recalculate accelerations and central masses for drift step
        if (flag_accel) then
          minside = mcen
          do i = 1, nbig
            msofar = minside  +  m(i)
            gm(i) = G * mcen * msofar / minside
            minside = msofar
          end do
          gm(nbig+1:n) = G * mcen
!
          a(:,1:n) = accel_mvs (n,nbig,m,x,v,ngf)
          flag_accel = .false.
        end if
!
! Advance interaction Hamiltonian for half a timestep
        v(:,1:n) = v(:,1:n)  +  dtby2 * a(:,1:n)
!
! Advance Keplerian Hamiltonian (Jacobi/heliocentric coords for Big/Small bodies)
        call convert_helio_to_jacobi (n,nbig,m,x,v)
        do i = 1, n
          call drift (dt,gm(i),x(:,i),v(:,i))
        end do
        call convert_jacobi_to_helio (n,nbig,m,x,v)
!
! Advance interaction Hamiltonian for half a timestep
        a(:,1:n) = accel_mvs (n,nbig,m,x,v,ngf)
        v(:,1:n) = v(:,1:n)  +  dtby2 * a(:,1:n)
!
! If something hit the central body, restore saved coordinates
        call check_central (dt,n,nbig,m,x0,v0,x,v,nhit,hit)
        if (nhit == 0) exit
        x(:,1:n) = x0(:,1:n)
        v(:,1:n) = v0(:,1:n)
!
! Merge colliding particles with central body, and redo the Keplerian step
        call collide_central (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index, &
          name,nhit,hit)
        flag_accel = .true.
      end do
!
      time = time  +  dt
!
! Check for close-encounter minima
      call check_encounters(time,dt,n,nbig,m,x0,v0,x,v,rad, &
        rce_hill,name,nclo,clo,nhit,hit)
      if (nclo > 0) call output_encounters (n,nbig,rho,status,index,name,nclo,clo,m)
!
! Remove particles far from the central body, update close encounter limits
      if (abs(time - tfun) >= dtfun) then
        call ejections (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
        if (flag_ejections) flag_accel = .true.
        tfun  = tfun   +  sign(dtfun, dt)
      end if
!
! Write a progress report and update data dump files
      if (abs(time - tlog) >= dtdump) then
        if (algor == 9) call inverse_corrector (dt,n,nbig,m,x,v,ngf)
        call output_progress (n,nbig,m,x,v,s)
        call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
        if (algor == 9) call corrector (dt,n,nbig,m,x,v,ngf)
        tlog  = tlog   +  sign(dtdump, dt)
        tdump = tdump  +  sign(dtdump, dt)
      end if
!
! Stop the integration for any reason?
      if (flag_stop) return
!
! Go on to the next time step
    end do
!
    end subroutine driver_mvs
!==============================================================================
! Does an integration using one of the variable stepsize algorithms.
!
    subroutine driver_variable_step (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    use constants;    use globals
    use interfaces, only: calc_radius, output_coords, output_datadump, &
      bs1, bs2, check_central, collide_central, remove_dead_bodies, &
      check_encounters, output_encounters, collide_bodies, output_progress, &
      ejections, ra15, calc_rce, output_codes
    implicit none
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    type(encounter)::clo(NMAX),hit(NMAX)
    integer(I4)::i,j,k,nclo,nhit
    real(R8)::rad(NMAX),rcrit(NMAX),x0(3,NMAX),v0(3,NMAX)
    real(R8)::temp,dt_nextcall,tsmall,t,dt
!------------------------------------------------------------------------------
! Initialize variables
    dt = dt0;      tsmall = dt0 * 1.0e-8_R8
    flag_accel = .true.
    if (goal == 'start') dt = sign(dt0, tstart - time)
    if (goal == 'stop ') dt = sign(dt0, tstop  - time)
    tdump = time;      tfun  = time;      tlog = time
!
! Calculate physical radii and close encounter distances
    rad(1:n) = calc_radius (n,m,rho)
!
    do
      flag_stop = .false.
!
! See if we have reached the start time for the integration
      if (goal == 'start'.and.abs(time - tstart) <= tsmall) then
        goal = 'stop '
        dt = sign(dt, tstop - time)
        tout = tstart  -  sign(dtout, dt)
        flag_accel = .true.
      end if
!
! Output the coords and velocities if it is time
      if (goal == 'stop '.and.abs(time - tout) >= dtout - tsmall) then
        call output_coords   (n,nbig,m,x,v,s,rho,status,index,name)
        call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
        tout  = tout   +  sign(dtout,  dt)
        tdump = tdump  +  sign(dtdump, dt)
      end if
!
! See if the integration has finished
      if (goal == 'stop '.and.abs(time - tstop ) <= tsmall) exit
!
! Set the timestep
      if (goal == 'start') then
        temp = min(abs(dt), abs(time - tstart))
      else
        temp = min(abs(dt), abs(time - tstop), abs(time - tout - sign(dtout, dt)))
      end if
      dt = sign(temp, dt)
!
!------------------------------------------------------------------------------
!  ADVANCE  ONE  TIMESTEP
!
! Save the current coordinates and velocities
      x0(:,1:n) = x(:,1:n)
      v0(:,1:n) = v(:,1:n)
!
! Advance one timestep
      if (algor == 2) call bs1  (dt,dt_nextcall,n,nbig,m,x,v,ngf)
      if (algor == 3) call bs2  (dt,dt_nextcall,n,nbig,m,x,v,ngf,rcrit,0)
      if (algor == 4) call ra15 (dt,dt_nextcall,n,nbig,m,x,v,ngf)
      time = time  +  dt
!
! Resolve any collisions with the central body
      call check_central (dt,n,nbig,m,x0,v0,x,v,nhit,hit)
      if (nhit > 0) then
        call collide_central (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index, &
          name,nhit,hit)
        flag_accel = .true.
      end if
!
! Check for close-encounter minima and collisions
      call check_encounters (time,dt,n,nbig,m,x0,v0,x,v,rad,rce_hill, &
        name,nclo,clo,nhit,hit)
      if (nclo > 0) call output_encounters (n,nbig,rho,status,index,name,nclo,clo,m)
!
! Remove any particles lost during collisions
      if (nhit > 0.and.opt_collisions) then
        do k = 1, nhit
          i = hit(k) % i;      j = hit(k) % j;      t = hit(k) % t
          call collide_bodies (t,i,j,n,nbig,m,x0,v0,x,v,s,ngf,rho,rce_hill,rad, &
            rcrit,status,index,name)
        end do
        call remove_dead_bodies (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
        call output_codes (n,nbig,status,index,name,m)
        flag_accel = .true.
      end if
!
! Remove particles far from the central body, update close encounter limits
      if (abs(time - tfun) >= abs(dtfun)) then
        call ejections (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
        if (flag_ejections) flag_accel = .true.
        tfun = tfun  +  sign(dtfun, dt)
      end if
!
! Write a progress report and update data dump files
      if (abs(time - tlog) >= dtdump) then
        call output_progress (n,nbig,m,x,v,s)
        call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
        tlog  = tlog   +  sign(dtdump, dt)
        tdump = tdump  +  sign(dtdump, dt)
      end if
!
! Stop the integration for any reason?
      if (flag_stop) return
!
! Go on to the next time step
      dt = dt_nextcall
    end do
!
    end subroutine driver_variable_step
!==============================================================================
! Removes all particles that are more than RMAX from the central body
! and writes a message to the information file recording each ejection event.
! The coordinates need to be with respect to the central body.
!
    subroutine ejections (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
!
    use constants;    use globals
    use interfaces, only: calc_date, calc_energy, remove_dead_bodies
    implicit none
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:)
    real(R8),     intent(inout)::ngf(:,:),rho(:),rce_hill(:),rad(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    integer(I4)::i,year1,month1
    real(R8)::r2,rmax2,t1,energy1,energy2,angmom1,angmom2
!------------------------------------------------------------------------------
    flag_ejections = .false.
    rmax2 = rmax * rmax
!
! Flag each particle that needs to be ejected (make INDEX negative)
    do i = 1, n
      r2 = dot_product(x(:,i), x(:,i))
      if (r2 > rmax2) then
        status(i) = 'dead '
        flag_ejections = .true.
      end if
    end do
!
! Record ejection events in the information file
    if (flag_ejections) then
  20  open  (23,file=outfile(3),status='old',access='append',err=20)
!
      do i = 1, n
        if (status(i) == 'dead ') then
          if (opt_time == 1) then
            call calc_date (time / DAY,year1,month1,t1)
            write (23,'(1x,a8,a,i10,1x,i2,1x,f8.5)') name(i),' ejected at ', &
              year1,month1,t1
          else
            if (opt_time == 0) t1 = time / DAY
            if (opt_time == 2) t1 = (time  -  tstart) / DAY
            if (opt_time == 3) t1 = (time  -  tstart) / YEAR
            write (23,'(1x,a8,a,f20.7,a)') name(i),' ejected at ',t1,time_string
          end if
        end if
      end do
!
      close (23)
    end if
!
! Remove ejected objects from global arrays, update energy and angular momentum
    if (flag_ejections) then
      call calc_energy (n,nbig,m,x,v,s,energy1,angmom1)
!
      call remove_dead_bodies (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status,index,name)
!
      call calc_energy (n,nbig,m,x,v,s,energy2,angmom2)
      denergy = denergy  +  (energy1  -  energy2)
      dangmom = dangmom  +  (angmom1  -  angmom2)
    end if
!
    end subroutine ejections
!==============================================================================
! Tidy up at the end of an integration and do any final outputs.
!
    subroutine finish (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    use constants;    use globals
    use interfaces, only: calc_energy, output_datadump
    implicit none
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    real(R8)::energy,angmom
!------------------------------------------------------------------------------
! Do a data dump and record the overall change in energy and angular momentum
    call output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
    call calc_energy (n,nbig,m,x,v,s,energy,angmom)
!
 50 open  (23, file=outfile(3), status='old', access='append', err=50)
    write (23,'(/,a)') '   Integration complete.'
    write (*,'(a)')   '    Integration complete.'
    write (23,231) '   Fractional energy change due to integrator: ', &
                  abs((energy  +  denergy  -  energy0) / energy0)
    write (23,232) '   Fractional angular momentum change:         ', &
                  abs((angmom  +  dangmom  -  angmom0) / angmom0)
    write (23,231) '   Fractional energy change due to collisions/ejections: ', &
                  abs(denergy / energy0)
    write (23,232) '   Fractional angular momentum change:                   ', &
                  abs(dangmom / angmom0)
    close (23)
!
 231  format (/,a,es12.5)
 232  format (a,es12.5)
!
    end subroutine finish
!==============================================================================
! Applies an inverse symplectic corrector for a time DT to a set of particles.
!
    subroutine inverse_corrector (dt,n,nbig,m,x,v,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_mvs, drift, &
      convert_helio_to_jacobi, convert_jacobi_to_helio
    implicit none
    real*8,parameter::HA(2) = (/ ROOT10 * 3.0_R8 / 10.0_R8, ROOT10 / 5.0_R8/)
    real*8,parameter::HB(2) = (/-ROOT10 / 72.0_R8,  ROOT10 / 24.0_R8/)
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::dt,m(:),ngf(:,:)
    real(R8),    intent(inout)::x(:,:),v(:,:)
!
    integer(I4)::i,k
    real(R8)::a(3,NMAX),gm(NMAX),minside,msofar,temp
!------------------------------------------------------------------------------
! Calculate effective central masses for Kepler drifts
    minside = mcen
    do i = 1, nbig
      msofar = minside  +  m(i)
      gm(i) = G * mcen * msofar / minside
      minside = msofar
    end do
    gm(nbig+1:n) = G * mcen
!
! Two step corrector
    do k = 1, 2
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call convert_helio_to_jacobi (n,nbig,m,x,v)
      temp = dt * HA(k)
      do i = 1, n
        call drift (temp,gm(i),x(:,i),v(:,i))
      end do
      call convert_jacobi_to_helio (n,nbig,m,x,v)
!
! Advance Interaction Hamiltonian
      a(:,1:n) = accel_mvs (n,nbig,m,x,v,ngf)
      temp = -dt * HB(k)
      v(:,1:n) = v(:,1:n)  +  temp * a(:,1:n)
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call convert_helio_to_jacobi (n,nbig,m,x,v)
      temp = -2.0_R8 * dt * HA(k)
      do i = 1, n
        call drift (temp,gm(i),x(:,i),v(:,i))
      end do
      call convert_jacobi_to_helio (n,nbig,m,x,v)
!
! Advance Interaction Hamiltonian
      a(:,1:n) = accel_mvs (n,nbig,m,x,v,ngf)
      temp = dt * HB(k)
      v(:,1:n) = v(:,1:n)  +  temp * a(:,1:n)
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call convert_helio_to_jacobi (n,nbig,m,x,v)
      temp = dt * HA(k)
      do i = 1, n
        call drift (temp,gm(i),x(:,i),v(:,i))
      end do
      call convert_jacobi_to_helio (n,nbig,m,x,v)
!
    end do
!
    end subroutine inverse_corrector
!==============================================================================
! Applies the solar Hamiltonian of the Duncan et al. symplectic integrator for
! a time DT, causing all the coordinates to jump by an equal amount.
! The number of active particles is N.
!
    subroutine jump (dt,n,m,x,v)
    use constants;      use globals
    implicit none
!
    integer(I4), intent(in)::n
    real(R8),    intent(in)::dt,m(:),v(:,:)
    real(R8),    intent(inout)::x(:,:)
    real(R8)::mvsum(3)
!------------------------------------------------------------------------------
    mvsum(1) = sum ( m(1:n) * v(1,1:n) )
    mvsum(2) = sum ( m(1:n) * v(2,1:n) )
    mvsum(3) = sum ( m(1:n) * v(3,1:n) )
    mvsum = mvsum * dt / mcen
!
    x(1,1:n) = x(1,1:n)  +  mvsum(1)
    x(2,1:n) = x(2,1:n)  +  mvsum(2)
    x(3,1:n) = x(3,1:n)  +  mvsum(3)
!
    end subroutine jump
!==============================================================================
! Merges a target particle ITARG with a projectile particle IPROJ to form
! a single object conserving mass, momentum and angular momentum. 
!
    subroutine merge_bodies (itarg,iproj,m,x,v,s,status,name)
!
    use constants;      use globals
    use interfaces, only: calc_date, cross_product
    implicit none
    integer(I4),  intent(in)::itarg,iproj
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
!
    real(R8)::msum,mredu,dx(3),dv(3),temp
!------------------------------------------------------------------------------
    msum  = m(itarg)  +  m(iproj)
    mredu = m(itarg) * m(iproj) / msum
    dx(:) = x(:,itarg)  -  x(:,iproj)
    dv(:) = v(:,itarg)  -  v(:,iproj)
!
! Calculate energy loss due to the collision
    temp = HALF * mredu * dot_product(dv, dv) &
                 -  G * m(itarg) * m(iproj) / sqrt(dot_product(dx, dx))
    denergy = denergy  +  temp
!
! New coordinates, velocities, spin angular momentum and mass
    x(:,itarg) = (m(itarg) * x(:,itarg)  +  m(iproj) * x(:,iproj)) / msum
    v(:,itarg) = (m(itarg) * v(:,itarg)  +  m(iproj) * v(:,iproj)) / msum
    s(:,itarg) = s(:,itarg)  +  s(:,iproj)  +  mredu * cross_product(dx, dv)
    m(itarg) = msum
!
! Flag the lost body for removal, and move it away from the merged body
    m(iproj) = ZERO;             x(:,iproj) = -x(:,iproj)
    v(:,iproj) = -v(:,iproj);    s(:,iproj) = ZERO
    status(iproj) = 'dead '
!
    end subroutine merge_bodies
!==============================================================================
! Checks whether any new particles have been added since last call (they have
! INDEX = 0). The routine writes a list of all new particles with their names
! and indices to the main output file and the close encounter output file.
! Also outputs some global variables relating to the central body.
! All output quantities are expressed in CGS units.
!
    subroutine output_codes (n,nbig,status,index,name,m)
!
    use constants;    use globals
    use interfaces, only: calc_float_string, calc_real_string
    implicit none
    integer(I4),  intent(in)::n,nbig
    character(8), intent(in)::name(:)
    character(5), intent(in)::status(:)
    integer(I4),  intent(inout)::index(:)
!
    integer(I4)::i,nnew_big,nnew_sml
    logical::flag
    character(80)::header,c
!
    real(R8),     intent(in)::m(:)
!------------------------------------------------------------------------------
! Determine if any new bodies are present (they will have INDEX = 0)
    nnew_big = 0;      nnew_sml = 0;    flag = .false.
!
    do i = 1, n
      if (index(i) == 0) then
        if (status(i) == 'big  ') nnew_big = nnew_big  +  1
        if (status(i) == 'small') nnew_sml = nnew_sml  +  1
        flag = .true.
      end if
    end do
!
! If new bodies are present, update list of output codes in output files
    if (flag) then

 120 open (23, file=outfile(3), status='old', access='append', err=120)
    write (23,*)
    write (23,'(a,f14.3,a)') ' Updating output codes at ',time/YEAR,' years'

 10   open (21, file=outfile(1), status='old', access='append', err=10)
 20   open (22, file=outfile(2), status='old', access='append', err=20)
!
! Write a header line with time, number of bodies and relevant parameters
      header(1:8)   = calc_float_string (time/DAY) !DAY converts to years from who knows what
      header(9:16)  = calc_real_string (dble(nnew_big), ZERO, INDEX_MAX)
      header(12:19) = calc_real_string (dble(nnew_sml), ZERO, INDEX_MAX)
      header(15:22) = calc_float_string (mcen/MSUN) !MSUN converts to solar masses from grams
      header(23:30) = calc_float_string (j2 / rcen**2)
      header(31:38) = calc_float_string (j4 / rcen**4)
      header(39:46) = calc_float_string (j6 / rcen**6)
      header(47:54) = calc_float_string (rcen/AU) !AU converts to AU from cm, next line too
      header(55:62) = calc_float_string (rmax/AU)
      write (21,'(a1,a2,i2,a62,i1)') char(12),'7a',algor,header(1:62),3 !The 3 here on this line and the next are to specify precision for the element6 and close6 codes
      write (22,'(a1,a2,i2,a62,i1)') char(12),'7a',algor,header(1:62),3
!
! For each new body, write its index number, name and density 
      do i = 1, n
        if (index(i) == 0) then
          nindex = nindex  +  1
          index(i) = nindex

   write (23,*) '  adding: ',name(i),index(i)
          c(1:8)   = calc_real_string (dble(index(i)), ZERO, INDEX_MAX)
          c(4:11)  = name(i)
!          c(12:19) = calc_float_string (m(i)/MSUN) !MSUN converts to solar masses from grams
!          c(20:27) = calc_float_string (0.0_R8) !zeroes because the spins don't matter for me
!          c(28:35) = calc_float_string (0.0_R8)
!          c(36:43) = calc_float_string (0.0_R8)
!          c(44:51) = calc_float_string (0.0_R8)!rho(i)) Already know rho, don't actually need to print
          write (21,'(a11)') c(1:11)
          write (22,'(a11)') c(1:11)
!          write(21,'(a51)') c(1:51)
!          write(22,'(a51)') c(1:51)
        end if
      end do
!
      close (21);      close (22)
    close (23)
    end if
!
    end subroutine output_codes
!==============================================================================
! Writes compressed output for all active particles listing the current mass,
! position, velocity and spin angular momenta. All quantities are output in
! CGS units. The number of active particles is N.
!
    subroutine output_coords (n,nbig,m,x,v,s,rho,status,index,name)
!
    use constants;    use globals
    use interfaces, only: calc_float_string, calc_real_string, output_codes, &
      calc_spherical_polars
    implicit none
    integer(I4),  intent(in)::n,nbig
    real(R8),     intent(in)::m(:),x(:,:),v(:,:),s(:,:),rho(:)
    character(8), intent(in)::name(:)
    character(5), intent(in)::status(:)
    integer(I4),  intent(inout)::index(:)
!
    integer(I4)::i
    real(R8)::r1,v1,s1,theta,phi,vtheta,vphi,stheta,sphi
    character(85)::header,c
    real(R8)::fr,fv,rfac

    rfac = log10(rmax/rcen)
!------------------------------------------------------------------------------
! Update list of output codes if necessary
    call output_codes (n,nbig,status,index,name,m)
!
! Open the orbital elements output file
 10 open (21, file=outfile(1), status='old', access='append', err=10)
!
! Write a header line containing the time and number of objects
    header(1:8)   = calc_float_string (time/DAY) !Dividing by DAY puts in years, instead of don't know what
    header(9:16)  = calc_real_string (dble(nbig),     ZERO, INDEX_MAX)
    header(12:19) = calc_real_string (dble(n - nbig), ZERO, INDEX_MAX)
    write (21,'(a1,a2,a14)') char(12),'7b',header(1:14)
!
! For each body, write its index number, mass, coords, velocities, and spin momenta
    do i = 1, n
      call calc_spherical_polars (x(:,i), r1,  theta,  phi)
      call calc_spherical_polars (v(:,i), v1, vtheta, vphi)
      call calc_spherical_polars (s(:,i), s1, stheta, sphi)
      call mco_x2ov(m,x(:,i),v(:,i),fr,fv)
!
      c(1:8)   = calc_real_string  (dble(index(i)), ZERO, INDEX_MAX)
!      c(4:11)  = calc_float_string (m(i))
!      c(11:18) = calc_float_string (rho(i))
!
      c(4:11) = calc_real_string (fr,    ZERO,rfac)!r1/AU) !AU converts to AU from cm
      c(11:18) = calc_real_string  (theta,  ZERO, PI)
      c(18:25) = calc_real_string  (phi,    ZERO, TWOPI)
!
      c(25:32) = calc_real_string (fv,   ZERO,ONE)!v1/AU*DAY) !AU*DAY converts to AU/day from cm/s
      c(32:39) = calc_real_string  (vtheta, ZERO, PI)
      c(39:46) = calc_real_string  (vphi,   ZERO, TWOPI)
      c(46:53) = calc_float_string (m(i)/MSUN) !MSUN converts to solar masses from gram
!      c(18:25) = calc_float_string (r1)
!      c(25:32) = calc_real_string  (theta,  -PIBY2, PIBY2)
!      c(32:39) = calc_real_string  (phi,    ZERO, TWOPI)
!
!      c(39:46) = calc_float_string (v1)
!      c(46:53) = calc_real_string  (vtheta, -PIBY2, PIBY2)
!      c(53:60) = calc_real_string  (vphi,   ZERO, TWOPI)
!
!      c(64:71) = calc_float_string (s1)
!      c(72:79) = calc_real_string  (stheta, -PIBY2, PIBY2)
!      c(79:85) = calc_real_string  (sphi,   ZERO, TWOPI)
      write (21,'(a53)') c(1:53)
    end do
!
    close (21)
!
    end subroutine output_coords
!==============================================================================
! Writes the current state of the entire system to a set of dump files.
!
    subroutine output_datadump (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name)
    use constants;    use globals;    use frag_globals
    implicit none
    integer(I4),  intent(in)::n,nbig,index(:)
    real(R8),     intent(in)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
    character(8), intent(in)::name(:)
    integer(I4)::idp,i,j,k,len,jj1,jj2
    character(8)::c8
    character(70)::cline
    character(150)::c
!------------------------------------------------------------------------------
    cline = ')---------------------------------------------------------------------'
!
! Dump to temporary files (idp = 1) and real dump files (idp = 2)
    do idp = 1, 2
!
! Dump data for the Big (i=1) and Small (i=2) bodies
      do i = 1, 2
        if (idp == 1) then
          if (i == 1) c(1:12) = 'big.tmp     '
          if (i == 2) c(1:12) = 'small.tmp   '
  20      open (31, file=c(1:12), status='unknown', err=20)
        else
  25      open (31, file=dumpfile(i), status='old', err=25)
        end if
!
! Write header lines, data style (and epoch for Big bodies)
        if (i == 1) write (31,'(2a)') ')O+_06 Big-body initial data  ', &
          '(WARNING: Do not delete this line!!)'
        if (i == 2) write (31,'(2a)') ')O+_06 Small-body initial data  ', &
          '(WARNING: Do not delete this line!!)'
        write (31,'(a)') ') Lines beginning with `)'' are ignored.'
        write (31,'(a)') cline
        write (31,'(a)') ' style (Cartesian, Asteroidal, Cometary) = Cartesian'
        if (i == 1) write (31,*) ' epoch (in days) = ',time / DAY
        write (31,'(a)') cline
!
! For each body...
        if (i == 1) then
          jj1 = 1;               jj2 = nbig
        else
          jj1 = nbig + 1;        jj2 = n
        end if
!
        do j = jj1, jj2
          c(1:8) = name(j)
          write (c(9:37),'(a3,es11.5,a3,es11.5)') ' r=',rce_hill(j),' d=',rho(j)
          write (c(38:62),'(a3,es22.15)') ' m=',m(j) / MSUN
          len = 62
          do k = 1, 3
            if (ngf(k,j) /= 0) then
              write (c(len+1:len+16),'(a2,i1,a1,es12.5)') ' a',k,'=',ngf(k,j)
              len = len + 16
            end if
          end do
          write (31,'(a)') c(1:len)
          write (31,'(3(1x,es22.15))') x(:,j) / AU
          write (31,'(3(1x,es22.15))') v(:,j) * DAY / AU
          write (31,'(3(1x,es22.15))') s(:,j) * DAY / (MSUN * AU * AU)
        enddo
        close (31)
      end do
!
! Dump the integration parameters
  40  if (idp == 1) open (33,file='param.tmp',status='unknown',err=40)
  45  if (idp == 2) open (33, file=dumpfile(3), status='old', err=45)
!
! Important parameters
      write (33,'(2a)') ')O+_06 Integration parameters  ', &
        '(WARNING: Do not delete this line!!)'
      write (33,'(a)')  ') Lines beginning with `)'' are ignored.'
      write (33,'(a)') cline
      write (33,'(a)')  ') Important integration parameters:'
      write (33,'(a)') cline
      c8 = '0       '
      if (algor ==  1) c8 = 'MVS     '
      if (algor ==  2) c8 = 'BS      '
      if (algor ==  3) c8 = 'BS2     '
      if (algor ==  4) c8 = 'RADAU   '
      if (algor ==  9) c8 = 'COR     '
      if (algor == 10) c8 = 'HYBRID  '
      write (33,*) ' algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = ',c8
      write (33,*) ' start time (days) =      ',tstart / DAY
      write (33,*) ' stop time (days) =       ',tstop / DAY
      write (33,*) ' output interval (days) = ',dtout / DAY
      write (33,*) ' timestep (days) =        ',dt0 / DAY
      write (33,*) ' accuracy parameter =     ',tol
!
! Integration options
      write (33,'(a)') cline
      write (33,'(a)')  ') Integration options:'
      write (33,'(a)') cline
      if (.not.opt_no_encounters) write (33,'(2a)') ' stop after ', &
        'a close encounter = no '
      if (opt_no_encounters)      write (33,'(2a)') ' stop after ', &
        'a close encounter = yes'
!
      if (.not.opt_collisions) write (33,'(a)') ' allow collisions = no'
      if (opt_collisions)      write (33,'(a)') ' allow collisions = yes'
!
      if (.not.opt_fragment) write (33,'(a)') ' include fragmentation = no'
      if (opt_fragment)      write (33,'(a)') ' include fragmentation = yes'
!
      if (opt_time == 0.or.opt_time == 2) then
        write (33,'(2a)') ' express time in days or years =  days '
      else
        write (33,'(2a)') ' express time in days or years = years'
      end if
      if (opt_time == 2.or.opt_time == 3) then
        write (33,'(a)') ' express time relative to integration start time = yes'
      else
        write (33,'(a)') ' express time relative to integration start time = no '
      end if
!
      write (33,'(a)') ' < Not used at present > '
      write (33,'(a)') ' < Not used at present > '
      write (33,'(a)') ' < Not used at present > '
      if (opt_user_force) write (33,'(a)')      ' include user-defined force = yes'
      if (.not.opt_user_force) write (33,'(a)') ' include user-defined force = no '
!
! Infrequently-changed parameters
      write (33,'(a)') cline
      write (33,'(a)') ') These parameters do not need to be adjusted often:'
      write (33,'(a)') cline
      write (33,*) ' ejection distance (AU) =      ',rmax / AU
      write (33,*) ' radius of central body (AU) = ',rcen / AU
      write (33,*) ' central mass (solar masses) = ',mcen / MSUN
      write (33,*) ' central J2 = ',j2 / rcen**2
      write (33,*) ' central J4 = ',j4 / rcen**4
      write (33,*) ' central J6 = ',j6 / rcen**6
      write (33,*) ' < Not used at present > '
      write (33,*) ' < Not used at present > '
      write (33,*) ' Hybrid integrator changeover (Hill radii) =    ',cefac
      write (33,*) ' number of timesteps between data dumps =       ',ndump
      write (33,*) ' number of timesteps between periodic effects = ',nfun
      write (33,*) ' minimum fragment mass =                        ',mfrag_min / MSUN
      close (33)
!
! Create new version of the restart file
  60  if (idp == 1) open (35, file='restart.tmp', status='unknown',err=60)
  65  if (idp == 2) open (35, file=dumpfile(4), status='old', err=65)
      if (goal == 'start') write (35,'(1x,i2)') -1
      if (goal == 'stop ') write (35,'(1x,i2)') 0
      write (35,*) energy0 / CONVERT_EN
      write (35,*) angmom0 / CONVERT_ANG
      write (35,*) denergy / CONVERT_EN
      write (35,*) dangmom / CONVERT_ANG
      write (35,*) scen(1) / CONVERT_ANG
      write (35,*) scen(2) / CONVERT_ANG
      write (35,*) scen(3) / CONVERT_ANG
      write (35,*) nfrag
      write (35,*) tout / DAY 
      close (35)
    end do
!
    end subroutine output_datadump
!==============================================================================
! Writes details of NCLO close encounters to the close-encounter output file.
! All output quantities are expressed in CGS units.
!
    subroutine output_encounters (n,nbig,rho,status,index,name,nclo,clo,m)
!
    use constants;    use globals
    use interfaces, only: calc_real_string, calc_float_string, &
      calc_spherical_polars, output_codes
    implicit none
    integer(I4) , intent(in)::n,nbig,nclo
    real(R8),     intent(in)::rho(:)
    character(8), intent(in)::name(:)
    character(5), intent(in)::status(:)
    type(encounter), intent(in)::clo(:)
    integer(I4),  intent(inout)::index(:)
!
    real(R8),     intent(in)::m(:)
!
    integer(I4)::i,j,k
    real(R8)::r1,v1,theta,phi,vtheta,vphi
    character(120)::c
    real(R8)::fr,fv,rfac

    rfac = log10(rmax/rcen)
!------------------------------------------------------------------------------
! Update list of output codes if necessary
    call output_codes (n,nbig,status,index,name,m)
!
! Open the close encounter output file
 10 open (22, file=outfile(2), status='old', access='append', err=10)
!
! Output details of each close encounter
    do k = 1, nclo
!       write(*,*) clo(k) % im
!       write(*,*) clo(k) % jm
!       write(*,*) clo(k) % ix
!       write(*,*) clo(k) % ix(1)
!       write(*,*) clo(k) % ix(2)
!       write(*,*) clo(k) % iv
!       write(*,*) clo(k) % iv(1)
!       write(*,*) clo(k) % iv(2)
!       write(*,*) clo(k) % jx
!       write(*,*) clo(k) % jx(1)
!       write(*,*) clo(k) % jx(2)
!       write(*,*) clo(k) % jv
!       write(*,*) clo(k) % jv(1)
!       write(*,*) clo(k) % jv(2)
      i = index(clo(k) % i);      j = index(clo(k) % j)
      c(1:8)   = calc_float_string(clo(k) % t/DAY)
      c(9:16)  = calc_real_string (dble(i), ZERO, INDEX_MAX)
      c(12:19) = calc_real_string (dble(j), ZERO, INDEX_MAX)
      c(15:22) = calc_float_string( (clo(k) % d)/AU)
!      if (k.eq.1) write(*,*) (clo(k) % d)/AU
!
      call calc_spherical_polars (clo(k) % ix, r1,  theta,  phi)
      call calc_spherical_polars (clo(k) % iv, v1, vtheta, vphi)
      call mco_x2ov(clo(k) % im,clo(k) % ix(:),clo(k) % iv(:),fr,fv)
      c(23:30) = calc_real_string (fr,      ZERO, rfac)
      c(27:34) = calc_real_string  (theta,  ZERO, PI)
      c(31:38) = calc_real_string  (phi,    ZERO, TWOPI)
      c(35:42) = calc_real_string (fv,     ZERO,ONE)
      c(39:46) = calc_real_string  (vtheta, ZERO, PI)
      c(43:50) = calc_real_string  (vphi,   ZERO, TWOPI)
!      if (k.eq.1) write(*,*) fr
!      if (k.eq.1) write(*,*) theta
!      if (k.eq.1) write(*,*) phi
!      if (k.eq.1) write(*,*) fv
!
      call calc_spherical_polars (clo(k) % jx, r1,  theta,  phi)
      call calc_spherical_polars (clo(k) % jv, v1, vtheta, vphi)
      call mco_x2ov(clo(k) % jm,clo(k) % jx(:),clo(k) % jv(:),fr,fv)
      c(47:54)   = calc_real_string (fr,     ZERO,rfac)
      c(51:58)   = calc_real_string  (theta,  ZERO, PI)
      c(55:62)   = calc_real_string  (phi,    ZERO, TWOPI)
      c(59:66)   = calc_real_string (fv,     ZERO,ONE)
      c(63:74)  = calc_real_string  (vtheta, ZERO, PI)
      c(67:78) = calc_real_string  (vphi,   ZERO, TWOPI)
!      if (k.eq.1) write(*,*) fr
!      if (k.eq.1) write(*,*) theta
!      if (k.eq.1) write(*,*) phi
!      if (k.eq.1) write(*,*) fv
!
! Output the compressed close encounter details
      write (22,'(a1,a2,a70)') char(12),'7b',c(1:110)
    end do
    close (22)
!
    end subroutine output_encounters
!==============================================================================
! Writes a status report to the screen giving the current time, and energy and
! angular momentum errors.
!
    subroutine output_progress (n,nbig,m,x,v,s)
!
    use constants;    use globals
    use interfaces, only: calc_date, calc_energy
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),x(:,:),v(:,:),s(:,:)
!
    integer(I4)::year1,month1
    real(R8)::tmp0,tmp1,energy,angmom,t1
!------------------------------------------------------------------------------
! Get current energy and angular momentum of the system
    call calc_energy (n,nbig,m,x,v,s,energy,angmom)
!
    tmp0 = ZERO;      tmp1 = ZERO
    if (energy0 /= 0) tmp0 = (energy  +  denergy  -  energy0) / abs(energy0)
    if (angmom0 /= 0) tmp1 = (angmom  +  dangmom  -  angmom0) / abs(angmom0)
!
    if (opt_time == 1) then
      call calc_date (time/DAY,year1,month1,t1)
      write (*,'(1x,a,i10,1x,i2,1x,f4.1,2(a,es12.5))') 'Date: ',year1,month1,t1, &
        '   dE/E: ',tmp0,'   dL/L: ',tmp1
    else
      if (opt_time == 0) t1 = time / DAY
      if (opt_time == 2) t1 = (time - tstart) / DAY
      if (opt_time == 3) t1 = (time - tstart) / YEAR
      write (*,'(1x,a,f16.3,a,2(a,es12.5))') 'Time: ',t1,time_string, &
        '   dE/E: ',tmp0,'   dL/L: ',tmp1
    end if
!
    end subroutine output_progress
!==============================================================================
! Advances all particles for a time DT using Everhart's RADAU algorithm.
! On entry: DT is the desired timestep
! On exit:  DT is the actual timestep used
!           DT_NEXTCALL is the recommend timestep next time the routine is called
!
    subroutine ra15 (dt,dt_nextcall,n,nbig,m,x0,v0,ngf)
!
    use constants;    use globals
    use interfaces, only: accel_all
    implicit none
    integer(I4), intent(in)::n,nbig
    real(R8),    intent(in)::m(:),ngf(:,:)
    real(R8),    intent(inout)::dt,dt_nextcall,x0(:,:),v0(:,:)
!
    integer(I4)::j,k,ns,niter
    real(R8)::a0(3,NMAX),x(3,NMAX),v(3,NMAX),a(3,NMAX)
    real(R8)::gg(3,7,NMAX),b(3,7,NMAX),e(3,7,NMAX),xscal(NMAX)
    real(R8)::q,q2,q3,q4,q5,q6,q7,temp,errmax
    real(R8)::c(21),d(21),r(28),s(9),ss(3,7),tempv(3),gk(3)
    logical::firstflag = .true.
!
! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
! integrator. The sum of the H values should be 3.733333333333333
    real(R8), parameter::H(8) = (/         ZERO, 0.0562625605369221_R8, &
      0.1802406917368924_R8, 0.3526247171131696_R8, 0.5471536263305554_R8, &
      0.7342101772154105_R8, 0.8853209468390958_R8, 0.9775206135612875_R8/)
!
! Constant coefficients used in series expansions for X and V
    real(R8), parameter::XC(8) = (/HALF, SIXTH, ONE / 12.0_R8, &
      ONE / 20.0_R8, ONE / 30.0_R8, ONE / 42.0_R8, ONE / 56.0_R8, ONE / 72.0_R8/)
    real(R8), parameter::VC(7) = (/HALF, THIRD, ONE / FOUR, ONE / 5.0_R8, &
      SIXTH, ONE / 7.0_R8, ONE / 8.0_R8/)
    real(R8), parameter::NINTH = ONE / 9.0_R8
!
    save c,d,r,b,e,firstflag
!------------------------------------------------------------------------------
! If this is first call to the subroutine, set values of the constant arrays
! (R = R21, R31, R32, R41, R42, R43 in Everhart 1985.)
    if (firstflag) then
        ns = 0
        do j = 2, 8
          do k = 1, j - 1
            ns = ns + 1
            r(ns) = ONE / (H(j) - H(k))
          end do
        end do
!
! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
        c(1) = - H(2)
        d(1) =   H(2)
        ns = 1
        do j = 3, 7
          ns = ns + 1
          c(ns) = -H(j) * c(ns-j+2)
          d(ns) =  H(2) * d(ns-j+2)
          do k = 3, j - 1
            ns = ns + 1
            c(ns) = c(ns-j+1)  -  H(j) * c(ns-j+2)
            d(ns) = d(ns-j+1)  +  H(k) * d(ns-j+2)
          end do
          ns = ns + 1
          c(ns) = c(ns-j+1) - H(j)
          d(ns) = d(ns-j+1) + H(j)
        end do
!
        firstflag = .false.
      end if
!
  100 continue
!
! If this is first call to subroutine since number/masses of objects changed
! do 6 iterations and initialize B, E arrays, otherwise do 2 iterations.
      if (flag_accel) then
        niter = 6
        b(:,:,1:n) = ZERO;        e(:,:,1:n) = ZERO
      else
        niter = 2
      end if
!
! Calculate forces at the start of the sequence
      a0(:,1:n) = accel_all (n,nbig,m,x0,v0,ngf)
!
! Arrays used to scale the relative error (1 / absolute distance or velocity)
    do k = 1, n
      xscal(k) = ONE / dot_product (x0(:,k), x0(:,k))
    end do
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
      do k = 1, n
        gg(:,1,k) = b(:,7,k) * d(16)  +  b(:,6,k) * d(11)  +  b(:,5,k) * d(7) &
                  + b(:,4,k) * d(4)   +  b(:,3,k) * d(2)   +  b(:,2,k) * d(1) &
                  + b(:,1,k)
        gg(:,2,k) = b(:,7,k) * d(17)  +  b(:,6,k) * d(12)  +  b(:,5,k) * d(8) &
                  + b(:,4,k) * d(5)   +  b(:,3,k) * d(3)   +  b(:,2,k)
        gg(:,3,k) = b(:,7,k) * d(18)  +  b(:,6,k) * d(13)  +  b(:,5,k) * d(9) &
                  + b(:,4,k) * d(6)   +  b(:,3,k)
        gg(:,4,k) = b(:,7,k) * d(19)  +  b(:,6,k) * d(14)  +  b(:,5,k) * d(10) &
                  + b(:,4,k)
        gg(:,5,k) = b(:,7,k) * d(20)  +  b(:,6,k) * d(15)  +  b(:,5,k)
        gg(:,6,k) = b(:,7,k) * d(21)  +  b(:,6,k)
        gg(:,7,k) = b(:,7,k)
      end do
!
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for first call to subroutine, two otherwise)...
      do ns = 1, niter
!
! For each substep within a sequence...
        do j = 2, 8
!
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = dt * H(j)
          s(2) = s(1) * s(1) * HALF
          s(3) = s(2) * H(j) * THIRD
          s(4) = s(3) * H(j) * HALF
          s(5) = s(4) * H(j) * 0.6_R8
          s(6) = s(5) * H(j) * 0.6666666666666667_R8
          s(7) = s(6) * H(j) * 0.7142857142857143_R8
          s(8) = s(7) * H(j) * 0.75_R8
          s(9) = s(8) * H(j) * 0.7777777777777778_R8
!
          do k = 1, n
            x(:,k) = s(9) * b(:,7,k)  +  s(8) * b(:,6,k)  +  s(7) * b(:,5,k) &
                   + s(6) * b(:,4,k)  +  s(5) * b(:,3,k)  +  s(4) * b(:,2,k) &
                   + s(3) * b(:,1,k)  +  s(2) * a0(:,k)   +  s(1) * v0(:,k)  +  x0(:,k)
          end do
!
! If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
          if (flag_ngf) then
            s(1) = dt * H(j)
            s(2) = s(1) * H(j) * HALF
            s(3) = s(2) * H(j) * 0.6666666666666667_R8
            s(4) = s(3) * H(j) * 0.75_R8
            s(5) = s(4) * H(j) * 0.8_R8
            s(6) = s(5) * H(j) * 0.8333333333333333_R8
            s(7) = s(6) * H(j) * 0.8571428571428571_R8
            s(8) = s(7) * H(j) * 0.875_R8
!
            do k = 1, n
              v(:,k) = s(8) * b(:,7,k)  +  s(7) * b(:,6,k)  +  s(6) * b(:,5,k) &
                     + s(5) * b(:,4,k)  +  s(4) * b(:,3,k)  +  s(3) * b(:,2,k) &
                     + s(2) * b(:,1,k)  +  s(1) * a0(:,k)   +  v0(:,k)
            end do
          end if
!
! Calculate forces at the current substep
      a(:,1:n) = accel_all (n,nbig,m,x,v,ngf)
!
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
          if (j == 2) then
            do k = 1, n
              tempv(:)  = gg(:,1,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,1,k) = gk(:) * r(1)
              tempv(:) = gg(:,1,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:)
            end do
          else if (j == 3) then
            do k = 1, n
              tempv(:) = gg(:,2,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,2,k) = (gk(:)*r(2) - gg(:,1,k))*r(3)
              tempv(:) = gg(:,2,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:) * c(1)
              b(:,2,k) =  b(:,2,k)  +  tempv(:)
            end do
          else if (j == 4) then
            do k = 1, n
              tempv(:) = gg(:,3,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,3,k) = ((gk(:)*r(4) - gg(:,1,k))*r(5) - gg(:,2,k))*r(6)
              tempv(:) = gg(:,3,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:) * c(2)
              b(:,2,k) =  b(:,2,k)  +  tempv(:) * c(3)
              b(:,3,k) =  b(:,3,k)  +  tempv(:)
            end do
          else if (j == 5) then
            do k = 1, n
              tempv(:) = gg(:,4,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,4,k) = (((gk(:)*r(7) - gg(:,1,k))*r(8) - gg(:,2,k))*r(9) &
                     - gg(:,3,k))*r(10)
              tempv(:) = gg(:,4,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:) * c(4)
              b(:,2,k) =  b(:,2,k)  +  tempv(:) * c(5)
              b(:,3,k) =  b(:,3,k)  +  tempv(:) * c(6)
              b(:,4,k) =  b(:,4,k)  +  tempv(:)
            end do
          else if (j == 6) then
            do k = 1, n
              tempv(:) = gg(:,5,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,5,k) =  ((((gk(:)*r(11) - gg(:,1,k))*r(12) - gg(:,2,k))*r(13) &
                        - gg(:,3,k))*r(14) - gg(:,4,k))*r(15)
              tempv(:) = gg(:,5,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:) * c(7)
              b(:,2,k) =  b(:,2,k)  +  tempv(:) * c(8)
              b(:,3,k) =  b(:,3,k)  +  tempv(:) * c(9)
              b(:,4,k) =  b(:,4,k)  +  tempv(:) * c(10)
              b(:,5,k) =  b(:,5,k)  +  tempv(:)
            end do
          else if (j == 7) then
            do k = 1, n
              tempv(:) = gg(:,6,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,6,k) = (((((gk(:)*r(16) - gg(:,1,k))*r(17) - gg(:,2,k))*r(18) &
                        - gg(:,3,k))*r(19) - gg(:,4,k))*r(20) - gg(:,5,k))*r(21)
              tempv(:) = gg(:,6,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:) * c(11)
              b(:,2,k) =  b(:,2,k)  +  tempv(:) * c(12)
              b(:,3,k) =  b(:,3,k)  +  tempv(:) * c(13)
              b(:,4,k) =  b(:,4,k)  +  tempv(:) * c(14)
              b(:,5,k) =  b(:,5,k)  +  tempv(:) * c(15)
              b(:,6,k) =  b(:,6,k)  +  tempv(:)
            end do
          else if (j == 8) then
            do k = 1, n
              tempv(:) = gg(:,7,k)
              gk(:) = a(:,k) - a0(:,k)
              gg(:,7,k) = ((((((gk(:)*r(22) - gg(:,1,k))*r(23) - gg(:,2,k))*r(24) &
                        - gg(:,3,k))*r(25) - gg(:,4,k))*r(26) - gg(:,5,k))*r(27) &
                        - gg(:,6,k))*r(28)
              tempv(:) = gg(:,7,k)  -  tempv(:)
              b(:,1,k) =  b(:,1,k)  +  tempv(:) * c(16)
              b(:,2,k) =  b(:,2,k)  +  tempv(:) * c(17)
              b(:,3,k) =  b(:,3,k)  +  tempv(:) * c(18)
              b(:,4,k) =  b(:,4,k)  +  tempv(:) * c(19)
              b(:,5,k) =  b(:,5,k)  +  tempv(:) * c(20)
              b(:,6,k) =  b(:,6,k)  +  tempv(:) * c(21)
              b(:,7,k) =  b(:,7,k)  +  tempv(:)
            end do
          end if
        end do
      end do
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!
! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
      errmax = ZERO
      do k = 1, n
        temp = max( b(1,7,k)*b(1,7,k), b(2,7,k)*b(2,7,k), b(3,7,k)*b(3,7,k) )
        errmax = max(errmax, temp*xscal(k))
      end do
      errmax = sqrt(errmax)
!
      if (errmax == 0) then
        dt_nextcall = dt * 1.4_R8
      else
        dt_nextcall = sign( (tol * abs(dt)**7 * 72.0_R8 / errmax)**NINTH, dt)
      end if
!
! If sequence size for the first subroutine call is too big, go back and redo
! the sequence using a smaller size.
      if (flag_accel.and.abs(dt_nextcall/dt).lt.1) then
        dt = dt * 0.8_R8
        goto 100
      end if
!
! If new sequence size is much bigger than the current one, reduce it
      if (abs(dt_nextcall/dt).gt.1.4_R8) dt_nextcall = dt * 1.4_R8
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)
      temp = dt * dt
      do k = 1 , n
        x0(:,k) = (XC(8) * b(:,7,k)  +  XC(7) * b(:,6,k)  +  XC(6) * b(:,5,k) &
                +  XC(5) * b(:,4,k)  +  XC(4) * b(:,3,k)  +  XC(3) * b(:,2,k) &
                +  XC(2) * b(:,1,k)  +  XC(1) * a0(:,k))*temp  +  v0(:,k) * dt &
                +  x0(:,k)
!
        v0(:,k) = (VC(7) * b(:,7,k)  +  VC(6) * b(:,6,k)  +  VC(5) * b(:,5,k) &
                +  VC(4) * b(:,4,k)  +  VC(3) * b(:,3,k)  +  VC(2) * b(:,2,k) &
                +  VC(1) * b(:,1,k)  +  a0(:,k)) * dt  +  v0(:,k)
      end do
!
! Predict new B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
      q = dt_nextcall / dt
      q2 = q  * q;       q3 = q  * q2
      q4 = q2 * q2;      q5 = q2 * q3
      q6 = q3 * q3;      q7 = q3 * q4
!
      do k = 1, n
        ss(:,1:7) = b(:,1:7,k) - e(:,1:7,k)
!
! Estimate B values for the next sequence (Eqs. 13 of Everhart).
        e(:,1,k) = q* (b(:,7,k) *  7.0_R8  +  b(:,6,k) * 6.0_R8  +  b(:,5,k) *  5.0_R8 &
                 +     b(:,4,k) *  4.0_R8  +  b(:,3,k) * 3.0_R8  +  b(:,2,k) *  2.0_R8 &
                 +     b(:,1,k))
        e(:,2,k) = q2*(b(:,7,k) * 21.0_R8  +  b(:,6,k) *15.0_R8  +  b(:,5,k) * 10.0_R8 &
                 +     b(:,4,k) *  6.0_R8  +  b(:,3,k) * 3.0_R8  +  b(:,2,k))
        e(:,3,k) = q3*(b(:,7,k) * 35.0_R8  +  b(:,6,k) *20.0_R8  +  b(:,5,k) * 10.0_R8 &
                 +     b(:,4,k) *  4.0_R8  +  b(:,3,k))
        e(:,4,k) = q4*(b(:,7,k) * 35.0_R8  +  b(:,6,k) *15.0_R8  +  b(:,5,k) *  5.0_R8 &
                 +     b(:,4,k))
        e(:,5,k) = q5*(b(:,7,k) * 21.0_R8  +  b(:,6,k) * 6.0_R8   +  b(:,5,k))
        e(:,6,k) = q6*(b(:,7,k) *  7.0_R8  +  b(:,6,k))
        e(:,7,k) = q7* b(:,7,k)
!
        b(:,1:7,k) = e(:,1:7,k)  +  ss(:,1:7)
      end do
      flag_accel = .false.
!
      end subroutine ra15
!==============================================================================
! Eliminates particles flagged for removal from global particle arrays.
! Particles that are meant to be removed should have negative values of INDEX.
!
    subroutine remove_dead_bodies (n,nbig,m,x,v,s,ngf,rho,rce_hill,rad,status, &
      index,name)
    use constants;    use globals
    implicit none
!
    integer(I4),  intent(inout)::n,nbig,index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:)
    real(R8),     intent(inout)::rce_hill(:),rad(:)
    character(8), intent(inout)::name(:)
    character(5), intent(inout)::status(:)
    integer(I4)::i,j,k,nelim,nbigelim,elim(NMAX+1)
!------------------------------------------------------------------------------
! Find out how many objects are to be removed
    nelim = 0;      nbigelim = 0
!
    do j = 1, n
      if (status(j) == 'dead ') then
        nelim = nelim + 1
        elim(nelim) = j
        if (j <= nbig) nbigelim = nbigelim + 1
      end if
    end do
    elim(nelim+1) = n  +  1
!
! Eliminate unwanted objects
    do i = 1, nelim
      do j = elim(i) - i + 1, elim(i+1) - i - 1
        k = j + i
        m(j)        = m(k);           x(:,j)   = x(:,k)
        v(:,j)      = v(:,k);         s(:,j)   = s(:,k)
        ngf(:,j)    = ngf(:,k);       rho(j)   = rho(k)
        rce_hill(j) = rce_hill(k);    rad(j)   = rad(k)
        status(j)   = status(k);      index(j) = index(k)
        name(j)     = name(k)
      end do
    end do
!
! If no more big bodies remain, write a warning message
    if (nbig >= 1.and.nbig - nbigelim <= 0) then
 10   open (23,file=outfile(3),status='old',access='append',err=10)
      write (23,'(2a)') ' WARNING: No more Big bodies are present.'
      close (23)
    end if
!
! Change status of newly empty particles
    status((n+1-nelim):n) = 'empty'
!
! Update total number of bodies and big bodies
    n    = n     -  nelim
    nbig = nbig  -  nbigelim
!
    end subroutine remove_dead_bodies
!==============================================================================
! Reads the initial conditions and parameters from files and sets up other
! variables needed for an integration.
!
    subroutine setup (n,nbig,m,x,v,s,ngf,rho,rce_hill,status,index,name)
!
    use constants;    use globals;    use frag_globals
    use interfaces, only: calc_cartesian_coords, calc_energy, calc_substring, &
      synch_epochs
    implicit none
    integer(I4),  intent(out)::n,nbig,index(:)
    real(R8),     intent(out)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:),rho(:),rce_hill(:)
    character(8), intent(out)::name(:)
    character(5), intent(out)::status(:)
!
    integer::j,k,nsub,lim(2,10),lineno,informat,itemp,month1,year1
    real(R8)::q,a,e,incl,p,node,mean,temp,t1,epoch(NMAX)
    logical::test,oldflag,flag1,flag2
    character(1)::c1
    character(3)::c3
    character(50)::c50
    character(80)::filnam,infile(3),files(10),c80
    character(150)::string
    character(3), dimension(60)::alg = &
        (/'MVS','Mvs','mvs','mvs','mvs','BS ','Bs ','bs ','BS1','bs1', &
          'BS2','Bs2','bs2','bs2','bs2','RAD','Rad','rad','RA ','ra ', &
          'xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx', &
          'xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx', &
          'COR','Cor','cor','cor','cor','HYB','Hyb','hyb','HY ','hy ', &
          'CLO','Clo','clo','CB ','cb ','WID','Wid','wid','WB ','wb '/)
!------------------------------------------------------------------------------
    filnam(1:40)  = '                                        '
    filnam(41:80) = '                                        '
    infile(1:3) = filnam;    outfile(1:3) = filnam;    dumpfile(1:4) = filnam
!
! Read names of input and output file and check for duplicate names
    inquire (file='files.in', exist=test)
    if (.not.test) then
      write (*,'(/,a)') ' ERROR: This file is needed to continue:  files.in'
      stop
    end if
!
    open (15, file='files.in', status='old')
    do j = 1, 10
      read (15,'(a150)') string
      call calc_substring (string,nsub,lim)
      files(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
      do k = 1, j - 1
        if (files(j) == files(k)) then
          write (*,'(/a/a)') ' ERROR: This file name is duplicated in files.in: ', &
            files(j)
          stop
        end if
      end do
!
      if (j <= 3)           infile(j)     = files(j)
      if (j > 3.and.j <= 6) outfile(j-3)  = files(j)
      if (j > 6)            dumpfile(j-6) = files(j)
    end do
    close (15)
!
! Find out if this is an old integration (i.e. does the restart file exist)
    inquire (file=dumpfile(4), exist=oldflag)
!
! Make sure output and dump files exist for an old integration
    do j = 4, 10
      inquire (file=files(j), exist=test)
      if (oldflag) then
        if (.not.test) then
          write (*,'(/a/a)') ' ERROR: This file is needed to continue: ',files(j)
          stop
        end if
!
! Create output and dump files for a new integration
      else
        if (test) then
          write (*,'(/a/a)') ' ERROR: This file already exists: ',files(j)
          stop
        end if
        open (11, file = files(j), status='new')
        close (11)
      end if
    end do
!
!------------------------------------------------------------------------------
!  READ  IN  INTEGRATION  PARAMETERS
!
    if (oldflag)      open (13, file = dumpfile(3), status='old')
    if (.not.oldflag) open (13, file = infile(3),   status='old')
!
! Default options
    opt_no_encounters = .false.
    opt_collisions    = .true.
    opt_fragment      = .false.
    opt_user_force    = .false.
!
    lineno = 0
    do j = 1, 27
      do
        lineno = lineno  +  1
        read (13,'(a150)') string
        if (string(1:1) /= ')') exit
      end do
!      
      call calc_substring (string,nsub,lim)
      c80(1:3) = '   '
      c80 = string(lim(1,nsub):lim(2,nsub))
      c1 = c80(1:1)
!
      if (j == 1) then
        algor = 0
        do k = 1, 60
          if (c80(1:3) == alg(k)) algor = (k + 4) / 5
        end do
        if (algor == 0) goto 661
      end if
      if (j ==  2) read (c80,*,err=661) tstart
      if (j ==  3) read (c80,*,err=661) tstop
      if (j ==  4) read (c80,*,err=661) dtout
      if (j ==  5) read (c80,*,err=661) dt0
      if (j ==  6) read (c80,*,err=661) tol
      if (j ==  7.and.(c1 == 'y'.or.c1 == 'Y')) opt_no_encounters = .true.
      if (j ==  8.and.(c1 == 'n'.or.c1 == 'N')) opt_collisions    = .false.
      if (j ==  9.and.(c1 == 'y'.or.c1 == 'Y')) opt_fragment      = .true.
      if (j == 10.and.(c1 == 'y'.or.c1 == 'Y')) opt_time = 1
      if (j == 11.and.(c1 == 'y'.or.c1 == 'Y')) opt_time = opt_time + 2
      if (j == 15.and.(c1 == 'y'.or.c1 == 'Y')) opt_user_force    = .true.
      if (j == 16) read (c80,*,err=661) rmax
      if (j == 17) read (c80,*,err=661) rcen
      if (j == 18) read (c80,*,err=661) mcen
      if (j == 19) read (c80,*,err=661) j2
      if (j == 20) read (c80,*,err=661) j4
      if (j == 21) read (c80,*,err=661) j6
      if (j == 24) read (c80,*,err=661) cefac
      if (j == 25) read (c80,*,err=661) ndump
      if (j == 26) read (c80,*,err=661) nfun
!
! Read in minimum fragment mass
      if (j == 27) read (c80,*,err=661) mfrag_min
    end do
    close (13)
!
! Convert to CGS units etc
    tstart = tstart * DAY;     tstop = tstop * DAY
    dtout = abs(dtout) * DAY;  dt0 = dt0 * DAY
    rmax = rmax * AU;          rcen = rcen * AU;        mcen = mcen * MSUN
    dt0 = abs(dt0);            tol = abs(tol)
    rmax = abs(rmax);          rcen = abs(rcen);        cefac = abs(cefac)
    mcen = abs(mcen);          scen = ZERO
    j2 = j2 * rcen**2;         j4 = j4 * rcen**4;      j6 = j6 * rcen**6
    dtdump = dt0 * ndump;      dtfun = dt0 * nfun
    mfrag_min = mfrag_min * MSUN
!
!------------------------------------------------------------------------------
!  READ  IN  DATA  FOR  BIG  AND  SMALL  BODIES
!
    n = 0;      nbig = 0;      nindex = 0
    do j = 1, 2
      if (oldflag)      open(11, file = dumpfile(j), status = 'old')
      if (.not.oldflag) open(11, file = infile(j), status = 'old')
!
! Read data style
      do
        read (11,'(a150)') string
        if (string(1:1) /= ')') exit
      end do
      call calc_substring (string,nsub,lim)
      c3 = string(lim(1,nsub):(lim(1,nsub)+2))
      informat = 0
      if (c3 == 'Car'.or.c3 == 'car'.or.c3 == 'CAR') informat = 1
      if (c3 == 'Ast'.or.c3 == 'ast'.or.c3 == 'AST') informat = 2
      if (c3 == 'Com'.or.c3 == 'com'.or.c3 == 'COM') informat = 3
      if (informat == 0) then
        write (*,'(/,2a)') ' ERROR: Data style not recognized in ',filnam
        stop
      end if
!
! Read epoch of Big bodies
        if (j == 1) then 
          do
            read (11,'(a150)') string
            if (string(1:1) /= ')') exit
          end do
          call calc_substring (string,nsub,lim)
          read (string(lim(1,nsub):lim(2,nsub)),*,err=667) time
        end if
!
! Read information for each object
        do
          read (11,'(a)',end=140) string
          if (string(1:1) == ')') cycle
          call calc_substring (string,nsub,lim)
!
          n = n  +  1
          if (n > NMAX) then
            write (*,'(/,2a)') ' ERROR: Number of bodies exceeds NMAX. Modify ', &
              'and recompile source code.'
            stop
          end if
!
! Status (big or small)
          if (j == 1) status(n) = 'big  '
          if (j == 2) status(n) = 'small'
!
! Determine the name of the object
          if ((lim(2,1)-lim(1,1)) > 7) then
            write (*,'(/,3a)') ' WARNING: ', &
              'Truncating the name of this object to 8 characters: ', &
              string( lim(1,1):lim(2,1) )
          end if
          name(n) = string( lim(1,1):min(7+lim(1,1),lim(2,1)) )
!
! Check whether another object has the same name
          do k = 1, n - 1
            if (name(k) == name(n)) then
              write (*,'(/,2a)') ' ERROR:  Two bodies have this name: ',name(n)
              stop
            end if
          end do
!
! Default values of mass, close-encounter limit, density etc.
          m(n)     = ZERO;          rho(n) = ONE
          ngf(:,n) = ZERO;          rce_hill(n) = 0.1_R8
          index(n) = 0;             epoch(n) = time
!
! Read values of mass, close-encounter limit, density etc.
          do k = 3, nsub, 2
            c80 = string(lim(1,k-1):lim(2,k-1))
            read (string(lim(1,k):lim(2,k)),*,err=666) temp
!
            if (c80(1:1) == 'm'.or.c80(1:1) == 'M') m(n)   = temp
            if (c80(1:1) == 'r'.or.c80(1:1) == 'R') rce_hill(n) = abs(temp)
            if (c80(1:1) == 'd'.or.c80(1:1) == 'D') rho(n) = abs(temp)
            if (c80(1:2) == 'ep'.or.c80(1:2) == 'EP'.or.c80(1:2) == 'Ep') &
              epoch(n) = temp
            if (c80(1:2) == 'a1'.or.c80(1:2) == 'A1') ngf(1,n) = temp
            if (c80(1:2) == 'a2'.or.c80(1:2) == 'A2') ngf(2,n) = temp
            if (c80(1:2) == 'a3'.or.c80(1:2) == 'A3') ngf(3,n) = temp
          end do
!
! Convert to CGS units
          m(n)        = m(n) * MSUN
          epoch(n)    = epoch(n) * DAY
!
! Read Cartesian coordinates, velocities and spins of the bodies
          do
            read (11,'(a150)',end=666) string
            if (string(1:1) /= ')') exit
          end do
          backspace 11
          if (informat == 1) then
            read (11,*,err=666) x(:,n), v(:,n), s(:,n)
            x(:,n) = x(:,n) * AU
            v(:,n) = v(:,n) * AU / DAY
            s(:,n) = s(:,n) * MSUN * AU * AU / DAY
! Alternatively, read cometary or asteroidal elements
          else
            read (11,*,err=666) a,e,incl,p,node,mean,s(:,n)
            a = a * AU
            s(:,n) = s(:,n) * MSUN * AU * AU / DAY
            incl = incl * DR
            p = (p + node) * DR
            node = node * DR
            temp = G * (mcen  +  m(n))
!
            if (informat == 3) then
              q = a;        a = q / (ONE  -  e)
              mean = mod (sqrt(temp/(abs(a*a*a))) * (epoch(n) - mean), TWOPI)
            else
              q = a * (ONE  -  e);        mean = mean * DR
            end if
            call calc_cartesian_coords (temp,q,e,incl,p,node,mean,x(:,n),v(:,n))
          end if
        end do
 140    close (11)
        if (j == 1) nbig = n
    end do
    time = time * DAY
!
!------------------------------------------------------------------------------
!  CHECK  FOR  ATTEMPTS  TO  DO  INCOMPATIBLE  THINGS
!
! Check whether RCEN >= RMAX
    if (rcen >= rmax) then
      write (*,'(/,2a)') ' ERROR: Central-body radius must be less than ', &
        'the maximum radius.'
      stop
    end if
!
! Check if non-grav forces are being used with an incompatible algorithm
    if (flag_ngf.and.algor == 3) then
      write (*,'(/,3a)') ' ERROR: You cannot use non-gravitational forces ', &
        'with this algorithm.'
      stop
    end if
!
! Check whether MVS is being used to integrate massive Small bodies,
! or whether massive Small bodies have different epochs than Big bodies.
    flag1 = .false.;    flag2 = .false.
    do j = nbig  +  1, n
      if (m(j) /= 0) flag1 = .true.
      if (epoch(j) /= time) flag2 = .true.
    end do
!
    if (flag1) then
      if (algor == 1.or.algor == 9) then
        write (*,'(/,2a)') ' ERROR: You cannot integrate massive small bodies ', &
          'using this algorithm.'
        stop
      end if
      if (flag2) then
        write (*,'(/,2a)') ' ERROR: Massive small bodies must have the same ', &
          'epoch as the big bodies.'
        stop
      end if
    end if
!
!------------------------------------------------------------------------------
!  SET  FLAGS
!
! Set cometary non-gravitational-forces flag
    flag_ngf = .false.
    do j = 1, n
      if (ngf(1,j) /= 0.or.ngf(2,j) /= 0.or.ngf(3,j) /= 0) flag_ngf = .true.
    end do
!
! Set oblateness flag
    flag_oblate = .false.
    if (j2 /= ZERO.or.j4 /= ZERO.or.j6 /= ZERO) flag_oblate = .true.
!
! Set up a string containing the appropriate time unit for output messages
    time_string = '      '
    if (opt_time == 3) time_string = ' years'
    if (opt_time == 0.or.opt_time == 2) time_string = ' days '
!
!------------------------------------------------------------------------------
!  IF  CONTINUING  AN  OLD  INTEGRATION
!
    if (oldflag) then
      open(23, file = outfile(3), status = 'old', access = 'append')
      if (opt_time == 1) then
        call calc_date (time/DAY,year1,month1,t1)
        write (23,'(/,a,i10,i2,f8.5,/)') &
          '   Continuing integration from dump files at ',year1,month1,t1
      else
        if (opt_time == 0) t1 = time / DAY
        if (opt_time == 2) t1 = (time  -  tstart) / DAY
        if (opt_time == 3) t1 = (time  -  tstart) / YEAR
        write (23,'(/,a,f20.7,a,/)') &
          '   Continuing integration from dump files at ',t1,time_string
      end if
!
! Read in energy and angular momentum variables, and convert to internal units
 330  open (35, file=dumpfile(4), status='old', err=330)
        read (35,*) itemp
        if (itemp <  0) goal = 'start'
        if (itemp >= 0) goal = 'stop '
        read (35,*) energy0,angmom0,denergy,dangmom
        read (35,*) scen(:)
        read (35,*) nfrag
        read (35,*) tout
        energy0 = energy0 * CONVERT_EN;    angmom0 = angmom0 * CONVERT_ANG
        denergy = denergy * CONVERT_EN;    dangmom = dangmom * CONVERT_ANG
        scen = scen * CONVERT_ANG
        tout = tout * DAY
      close (35)
!
!------------------------------------------------------------------------------
!  IF  STARTING  A  NEW  INTEGRATION
!
    else
      open(23, file = outfile(3), status = 'old')
      goal = 'start'
      tout = tstart
!
! Write integration parameters to information file
      write (23,'(/,a)') '                MERCURY 7_0'
      write (23,'(/,a)') '           Integration parameters'
      write (23,'(a)')   '           ----------------------'
      if (algor ==  1) c50 = 'Second-order mixed-variable symplectic (MVS)      '
      if (algor ==  2) c50 = 'Bulirsch-Stoer (general)                          '
      if (algor ==  3) c50 = 'Bulirsch-Stoer (conservative systems)             '
      if (algor ==  4) c50 = '15th-order RADAU                                  '
      if (algor ==  9) c50 = 'Second-order MVS with symplectic corrector        '
      if (algor == 10) c50 = 'Hybrid symplectic integrator (mixed coordinates)  '
      write (23,'(/,2a)') '   Algorithm: ',c50
!
      if (abs(tstart / DAY) >= 1.d10) then
        write (23,'(/,a,es19.13,a)') '   Integration start epoch:     ', &
        tstart / DAY,' days'
      else
        write (23,'(/,a,f19.7,a)')   '   Integration start epoch:     ', &
        tstart / DAY,' days'
      end if
      if (abs(tstop / DAY) > 1.d10) then
        write (23,'(a,es19.13)') '   Integration stop  epoch:     ',tstop / DAY
      else
        write (23,'(a,f19.7)')   '   Integration stop  epoch:     ',tstop / DAY
      end if
      write (23,'(a,f15.3)') '   Output interval:             ',dtout / DAY
!
      write (23,'(/,a,f10.4,a)')  '   Initial timestep:      ',dt0 / DAY,' days '
      write (23,'(a,es10.4)')     '   Accuracy parameter:       ',tol
      write (23,'(a,es10.4,a)')   '   Central mass:             ',mcen / MSUN, &
       ' solar masses'
      write (23,'(a,es10.4)')     '   J2:                       ',j2 / rcen**2
      write (23,'(a,es10.4)')     '   J4:                       ',j4 / rcen**4
      write (23,'(a,es10.4)')     '   J6:                       ',j6 / rcen**6
      write (23,'(a,es10.4,a)')   '   Ejection distance:        ',rmax / AU,' AU'
      write (23,'(a,es10.4,a)')   '   Radius of central body:   ',rcen / AU,' AU'
      write (23,*)
!
      if (.not.opt_no_encounters) write (23,'(2a)') '   Stop after ', &
        'a close encounter:        no '
      if (opt_no_encounters)      write (23,'(2a)') '   Stop after ', &
        'a close encounter:        yes'
!
      if (.not.opt_collisions) write (23,'(2a)') '   Allow collisions: ', &
       '                   no'
      if (opt_collisions)      write (23,'(2a)') '   Allow collisions: ', &
       '                   yes'
!
      if (.not.opt_fragment) write (23,'(2a)') '   Include fragmentation: ', &
       '              no'
      if (opt_fragment)      write (23,'(2a)') '   Include fragmentation: ', &
       '              yes'
!
      if (opt_user_force)      write (23,'(2a)') '   Include user-defined ', &
        'force routine:  yes'
      if (.not.opt_user_force) write (23,'(2a)') '   Include user-defined ', &
        'force routine:  no '
!
      write (23,'(/,a,i4)') '   Number of Big bodies:     ', nbig
      write (23,'(a,i4)')   '   Number of Small bodies:   ', n - nbig
!
! Set up the initial number of fragments and the minimum fragment mass
      nfrag = 0
      write (23,'(a,es11.4,a)') '   Minimum fragment mass:      ', &
        mfrag_min/MSUN,' solar masses'
!
! Calculate initial energy and angular momentum and write to information file
      denergy = ZERO;      dangmom = ZERO
      call calc_energy (n,nbig,m,x,v,s,energy0,angmom0)
      write (23,'(//,a)') '           Integration details'
      write (23,'(a)')    '           -------------------'
      write (23,'(/,a,es12.5,a)') '   Initial energy:           ', &
        energy0 / CONVERT_EN,  ' solar masses AU^2 day^-2'
      write (23,'(a,es12.5,a)')   '   Initial angular momentum: ', &
        angmom0 / CONVERT_ANG, ' solar masses AU^2 day^-1'
!
! Write warning messages if necessary
      if (tstop < tstart) write (23,'(/,2a)') ' WARNING: ', &
       'Main integration is backwards.'
      if (nbig <= 0) write (23,'(/,2a)') ' WARNING: No Big bodies are present.'
      if (nbig == n) write (23,'(/,2a)') ' WARNING: No Small bodies are present.'
!
! Integrate all the bodies to a common epoch
      write (23,'(/,2a)') '   Integrating massive bodies and particles to the ', &
        'same epoch.'
      call synch_epochs (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name,epoch)
      write (23,'(/,a,/)') '   Beginning the main integration.'
    end if
!
    close (23)
    return
!

! Error reading from the input file containing integration parameters
 661  write (*,'(/,a,i4,2a)') ' ERROR: A problem occurred reading from line ', &
        lineno,' of ',filnam
      stop
!
! Error reading from the input file for Big or Small bodies
 666  write (*,'(/,2a)') ' ERROR: A problem occurred reading data for this body: ', &
        name(n)
      stop
!
! Error reading epoch of Big bodies
 667  write (*,'(/,a)') ' ERROR: A problem occurred reading the epoch for the big bodies.'
      stop
!
    end subroutine setup
!==============================================================================
! Sorts an array X of size N in increasing order.
!
    subroutine sort_array (n,x,order)
!
    use constants
    implicit none
    integer(I4), intent(in)::n
    real(R8),    intent(inout)::x(n)
    integer(I4), intent(out)::order(n)
!
    integer(I4)::i,j,k,l,m,inc,iy
    integer(I4)::incarr(9) = (/1,4,13,40,121,364,1093,3280,9841/)
    real(R8)::y
!------------------------------------------------------------------------------
    do i = 1, n
      order(i) = i
    end do
!
    m = 0
 10 m = m + 1
    if (incarr(m) < n) goto 10
    m = m - 1
!
    do i = m, 1, -1
      inc = incarr(i)
      do j = 1, inc
        do k = inc, n - j, inc
          y = x(j+k)
          iy = order(j+k)
          do l = j + k - inc, j, -inc
            if (x(l) <= y) goto 20
            x(l+inc) = x(l)
            order(l+inc) = order(l)
          end do
 20       x(l+inc) = y
          order(l+inc) = iy
        end do
      end do
    end do
!
    end subroutine sort_array
!==============================================================================
! Integrates all the particles to a common epoch.
!
    subroutine synch_epochs (n,nbig,m,x,v,s,ngf,rho,rce_hill,index,name,epoch)
!
    use constants;    use globals
    use interfaces, only: bs1, sort_array
    implicit none
    integer(I4),  intent(in)::n,nbig
    integer(I4),  intent(inout)::index(:)
    real(R8),     intent(inout)::m(:),x(:,:),v(:,:),s(:,:),ngf(:,:)
    real(R8),     intent(inout)::rho(:),rce_hill(:),epoch(:)
    character(8), intent(inout)::name(:)
!
    integer(I4)::i,j,nsml,nsofar,order(NMAX),indextemp(NMAX)
    real(R8)::temp,epsml(NMAX),dt,dt_nextcall,tsmall
    real(R8)::mtemp(NMAX),xtemp(3,NMAX),vtemp(3,NMAX),stemp(3,NMAX)
    real(R8)::ngftemp(3,NMAX),rhotemp(NMAX),rcetemp(NMAX),epochtemp(NMAX)
    character(8)::nametemp(NMAX)
!------------------------------------------------------------------------------
! Sort the small bodies by epoch
    nsml = n  -  nbig
    do i = nbig  +  1, n
      epsml(i-nbig) = epoch(i)
    end do
!
    call sort_array (nsml,epsml,order)
!
    do i = nbig + 1, n
      j = nbig  +  order(i - nbig)
      mtemp(i)     = m(j)
      xtemp(:,i)   = x(:,j)
      vtemp(:,i)   = v(:,j)
      stemp(:,i)   = s(:,j)
      ngftemp(:,i) = ngf(:,j)
      rhotemp(i)   = rho(j)
      rcetemp(i)   = rce_hill(j)
      epochtemp(i) = epoch(j)
      indextemp(i) = index(j)
      nametemp(i)  = name(j)
    end do
    m(nbig+1:n)     = mtemp(nbig+1:n)
    x(:,nbig+1:n)   = xtemp(:,nbig+1:n)
    v(:,nbig+1:n)   = vtemp(:,nbig+1:n)
    s(:,nbig+1:n)   = stemp(:,nbig+1:n)
    ngf(:,nbig+1:n) = ngftemp(:,nbig+1:n)
    rho(nbig+1:n)   = rhotemp(nbig+1:n)
    rce_hill(nbig+1:n) = rcetemp(nbig+1:n)
    epoch(nbig+1:n)    = epochtemp(nbig+1:n)
    index(nbig+1:n)    = indextemp(nbig+1:n)
    name(nbig+1:n)     = nametemp(nbig+1:n)
!
! Integrate all the small bodies up to the same epoch
    dt = dt0;      tsmall = dt0 * SMALL_NUMBER
!
    do nsofar = nbig, n - 1
      do while (abs(time  -  epoch(nsofar+1)) > tsmall)
        temp = epoch(nsofar+1)  -  time
        dt = max(min(abs(temp),abs(dt)),tsmall)
        dt = sign(dt, temp)
        call bs1 (dt,dt_nextcall,nsofar,nbig,m,x,v,ngf)
        time = time  +  dt
        dt   = dt_nextcall
      end do
    end do
!
    end subroutine synch_epochs
!==============================================================================


!*************************************************************************
!                        DRIFT_DAN.F
!*************************************************************************
! This subroutine does the Danby and decides which vbles to use
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 x0,y0,z0         ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 vx0,vy0,vz0      ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 dt0            ==>  time step
!             Output:
!                 x0,y0,z0         ==>  final position in jacobi coord 
!                                       (real scalars)
!                 vx0,vy0,vz0      ==>  final position in jacobi coord 
!                                       (real scalars)
!                 iflg             ==>  integer flag (zero if satisfactory)
!					      (non-zero if nonconvergence)
!
! Authors:  Hal Levison & Martin Duncan  
! Date:    2/10/93
! Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged

      subroutine drift_dan(mu,x0,y0,z0,vx0,vy0,vz0,dt0,iflg)
      use swift
      implicit none
!
      real(R8),    intent(in)::mu,dt0
      real(R8),    intent(inout)::x0,y0,z0,vx0,vy0,vz0
      integer(I4), intent(out)::iflg
      real(R8)::x,y,z,vx,vy,vz,dt
      real(R8)::f,g,fdot,c1,c2
      real(R8)::c3,gdot
      real(R8)::u,alpha,fp,r0,v0s
      real(R8)::a,asq,en
      real(R8)::dm,ec,es,esq,xkep
      real(R8)::fchk,s,c
!------------------------------------------------------------------------------
!...  Set dt = dt0 to be sure timestep is not altered while solving
!...  for new coords.
	dt = dt0
	iflg = 0
        r0 = sqrt(x0*x0 + y0*y0 + z0*z0)
        v0s = vx0*vx0 + vy0*vy0 + vz0*vz0
        u = x0*vx0 + y0*vy0 + z0*vz0
        alpha = 2.0*mu/r0 - v0s
        
	if (alpha > 0.d0) then
           a = mu/alpha
           asq = a*a
           en = sqrt(mu/(a*asq))
           ec = 1.0d0 - r0/a
           es = u/(en*asq)
	   esq = ec*ec + es*es
	   dm = dt*en - int(dt*en/TWOPI)*TWOPI
	   dt = dm/en
	   if((dm*dm  >  0.16d0) .or. (esq > 0.36d0)) goto 100

	   if(esq*dm*dm  <  0.0016) then

               call drift_kepmd(dm,es,ec,xkep,s,c)
	       fchk = (xkep - ec*s +es*(1.-c) - dm)

	       if(fchk*fchk  >  DANBYB) then
		  iflg = 1
		  return
	       endif

               fp = 1. - ec*c + es*s
               f = (a/r0) * (c-1.) + 1.
               g = dt + (s-xkep)/en
               fdot = - (a/(r0*fp))*en*s
               gdot = (c-1.)/fp + 1.

               x = x0*f + vx0*g
               y = y0*f + vy0*g
               z = z0*f + vz0*g
               vx = x0*fdot + vx0*gdot
               vy = y0*fdot + vy0*gdot
               vz = z0*fdot + vz0*gdot

               x0 = x
               y0 = y
               z0 = z
               vx0 = vx
               vy0 = vy
               vz0 = vz

	       iflg = 0
	       return

	   endif

         endif
             
100      call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)

         if(iflg  == 0) then
           f = 1.0 - (mu/r0)*c2
           g = dt - mu*c3
           fdot = -(mu/(fp*r0))*c1
           gdot = 1. - (mu/fp)*c2

           x = x0*f + vx0*g
           y = y0*f + vy0*g
           z = z0*f + vz0*g
           vx = x0*fdot + vx0*gdot
           vy = y0*fdot + vy0*gdot
           vz = z0*fdot + vz0*gdot

           x0 = x
           y0 = y
           z0 = z
           vx0 = vx
           vy0 = vy
           vz0 = vz
	endif

        return
        end subroutine drift_dan
!
!********************************************************************#
!                  DRIFT_KEPMD
!********************************************************************#
!  Subroutine for solving kepler's equation in difference form for an
!  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
!  for the criteria.
!  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
!  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
!  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
!
!	Input:
!	    dm		==> increment in mean anomaly M (real*8 scalar)
!	    es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars)
!
!       Output:
!            x          ==> solution to Kepler's difference eqn (real*8 scalar)
!            s,c        ==> sin and cosine of x (real*8 scalars)
!

        subroutine drift_kepmd(dm,es,ec,x,s,c)

	implicit none

!...    Inputs
	real*8 dm,es,ec
	
!...	Outputs
	real*8 x,s,c

!...    Internals
	real*8 A0, A1, A2, A3, A4
        parameter(A0 = 39916800.d0, A1 = 6652800.d0, A2 = 332640.d0)
	parameter(A3 = 7920.d0, A4 = 110.d0)
	real*8 dx
	real*8 fac1,fac2,q,y
	real*8 f,fp,fpp,fppp


!...    calc initial guess for root
	fac1 = 1.d0/(1.d0 - ec)
	q = fac1*dm
	fac2 = es*es*fac1 - ec/3.d0
	x = q*(1.d0 -0.5d0*fac1*q*(es -q*fac2))

!...  excellent approx. to sin and cos of x for small x.
	y = x*x
	s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
        c = sqrt(1.d0 - s*s)

!...    Compute better value for the root using quartic Newton method
        f = x - ec*s + es*(1.-c) - dm
        fp = 1. - ec*c + es*s
        fpp = ec*s + es*c
        fppp = ec*c - es*s
        dx = -f/fp
        dx = -f/(fp + 0.5*dx*fpp)
        dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666*dx*dx*fppp)
        x = x + dx
     
!...  excellent approx. to sin and cos of x for small x.
	y = x*x
	s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
        c = sqrt(1.d0 - s*s)

	return
	end subroutine drift_kepmd
!-----------------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU.F
!*************************************************************************
! subroutine for solving kepler's equation using universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflg          ==>  =0 if converged; !=0 if not
!
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 2/3/93

      subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
      use swift

!...  Inputs: 
      real*8 dt,r0,mu,alpha,u

!...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflg

!...  Internals:
      real*8 s,st,fo,fn

!----
!...  Executable code 

        call drift_kepu_guess(dt,r0,mu,alpha,u,s)
         
        st = s
!..     store initial guess for possible use later in
!..     laguerre's method, in case newton's method fails.

        call drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
        if(iflg /= 0) then
           call drift_kepu_fchk(dt,r0,mu,alpha,u,st,fo)
           call drift_kepu_fchk(dt,r0,mu,alpha,u,s,fn)
           if(abs(fo) < abs(fn)) then
               s = st 
           endif
           call drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
        endif

        return
        end subroutine drift_kepu
!----------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU_FCHK.F
!*************************************************************************
! Returns the value of the function f of which we are trying to find the root
! in universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalar)
!                 r0            ==>  Distance between `Sun' and particle
!                                     (real scalar)
!                 mu            ==>  Reduced mass of system (real scalar)
!                 alpha         ==>  Twice the binding energy (real scalar)
!                 u             ==>  Vel. dot radial vector (real scalar)
!                 s             ==>  Approx. root of f 
!             Output:
!                 f             ==>  function value ( = 0 if O.K.) (integer)
!
! Author:  Martin Duncan  
! Date:    March 12/93
! Last revision: March 12/93

      subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)

!...  Inputs: 
      real*8 dt,r0,mu,alpha,u,s

!...  Outputs:
      real*8 f

!...  Internals:
      real*8  x,c0,c1,c2,c3

!----
!...  Executable code 

        x=s*s*alpha
        call drift_kepu_stumpff(x,c0,c1,c2,c3)
        c1=c1*s
        c2 = c2*s*s
        c3 = c3*s*s*s
        f = r0*c1 + u*c2 + mu*c3 - dt

        return
        end subroutine drift_kepu_fchk
!-------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU_GUESS.F
!*************************************************************************
! Initial guess for solving kepler's equation using universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  initial guess for the value of 
!                                    universal variable
!
! Author:  Hal Levison & Martin Duncan 
! Date:    3/12/93
! Last revision: April 6/93
! Modified by JEC: 8/6/98

      subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)
      use swift
      use interfaces, only: calc_sin_cos

!...  Inputs: 
      real*8 dt,r0,mu,alpha,u

!...  Inputs and Outputs:
      real*8 s

!...  Internals:
      integer iflg
      real*8 y,sy,cy,sigma,es
      real*8 x,a
      real*8 en,ec,e

!----
!...  Executable code 

        if (alpha > 0.0) then 
!...       find initial guess for elliptic motion

            if( dt/r0  <=  0.4)  then
              s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0)
	      return
            else
              a = mu/alpha
              en = sqrt(mu/(a*a*a))
              ec = 1.0 - r0/a
              es = u/(en*a*a)
              e = sqrt(ec*ec + es*es)
              y = en*dt - es
!
              call calc_sin_cos (y,sy,cy)
!
              sigma = dsign(1.d0,(es*cy + ec*sy))
              x = y + sigma*.85*e
              s = x/sqrt(alpha)
	    endif

        else
!...       find initial guess for hyperbolic motion.
	   call drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
	   if(iflg /= 0) then
	      s = dt/r0
	   endif
        endif

        return
        end subroutine drift_kepu_guess
!-------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU_LAG.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using LAGUERRE'S METHOD
!
!             Input:
!                 s             ==>  inital value of universal variable
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  final value of universal variable
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflgn          ==>  =0 if converged; !=0 if not
!
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 4/21/93

      subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
      use swift

!...  Inputs: 
      real*8 s,dt,r0,mu,alpha,u

!...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflg

!...  Internals:
      integer nc,ncmax
      real*8 ln
      real*8 x,fpp,ds,c0,f
      real*8 fdt

      integer NTMP
      parameter(NTMP=NLAG2+1)

!----
!...  Executable code 

!...    To get close approch needed to take lots of iterations if alpha<0
        if(alpha < 0.0) then
           ncmax = NLAG2
        else
           ncmax = NLAG2
        endif

        ln = 5.0
!...    start laguere's method
        do nc =0,ncmax
           x = s*s*alpha
           call drift_kepu_stumpff(x,c0,c1,c2,c3)
           c1 = c1*s 
           c2 = c2*s*s 
           c3 = c3*s*s*s
           f = r0*c1 + u*c2 + mu*c3 - dt
           fp = r0*c0 + u*c1 + mu*c2
           fpp = (-40.0*alpha + mu)*c1 + u*c0
           ds = - ln*f/(fp + dsign(1.d0,fp)*sqrt(abs((ln - 1.0)* &
             (ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
           s = s + ds

           fdt = f/dt

!..        quartic convergence
           if( fdt*fdt < DANBYB*DANBYB) then 
             iflg = 0
             return
           endif
!...      Laguerre's method succeeded
        enddo

        iflg = 2

        return

        end subroutine drift_kepu_lag
!-----------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU_NEW.F
!*************************************************************************
! subroutine for solving kepler's equation in universal variables.
! using NEWTON'S METHOD
!
!             Input:
!                 s             ==>  inital value of universal variable
!                 dt            ==>  time step (real scalor)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalor)
!                 mu            ==>  Reduced mass of system (real scalor)
!                 alpha         ==>  energy (real scalor)
!                 u             ==>  angular momentun  (real scalor)
!             Output:
!                 s             ==>  final value of universal variable
!                 fp            ==>  f' from p170  
!                                       (real scalors)
!                 c1,c2,c3      ==>  c's from p171-172
!                                       (real scalors)
!                 iflgn          ==>  =0 if converged; !=0 if not
!
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 4/21/93
! Modified by JEC: 31/3/98

      subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflgn)
      use swift

!...  Inputs: 
      real*8 s,dt,r0,mu,alpha,u

!...  Outputs:
      real*8 fp,c1,c2,c3
      integer iflgn

!...  Internals:
      integer nc
      real*8 x,c0,ds,s2
      real*8 f,fpp,fppp,fdt

!----
!...  Executable code 

      do nc=0,6
         s2 = s * s
         x = s2*alpha
         call drift_kepu_stumpff(x,c0,c1,c2,c3)
         c1 = c1*s 
         c2 = c2*s2 
         c3 = c3*s*s2
         f = r0*c1 + u*c2 + mu*c3 - dt
         fp = r0*c0 + u*c1 + mu*c2
         fpp = (mu - r0*alpha)*c1 + u*c0
         fppp = (mu - r0*alpha)*c0 - u*alpha*c1
         ds = - f/fp
         ds = - f/(fp + .5d0*ds*fpp)
         ds = -f/(fp + .5d0*ds*fpp + ds*ds*fppp*.1666666666666667d0)
         s = s + ds
         fdt = f/dt

!..      quartic convergence
         if( fdt*fdt < DANBYB*DANBYB) then 
             iflgn = 0
             return
         endif
!...     newton's method succeeded

        enddo

!..     newton's method failed
        iflgn = 1
        return

        end subroutine drift_kepu_new
!----------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU_P3SOLVE.F
!*************************************************************************
! Returns the real root of cubic often found in solving kepler
! problem in universal variables.
!
!             Input:
!                 dt            ==>  time step (real scalar)
!                 r0            ==>  Distance between `Sun' and paritcle
!                                     (real scalar)
!                 mu            ==>  Reduced mass of system (real scalar)
!                 alpha         ==>  Twice the binding energy (real scalar)
!                 u             ==>  Vel. dot radial vector (real scalar)
!             Output:
!                 s             ==>  solution of cubic eqn for the  
!                                    universal variable
!                 iflg          ==>  success flag ( = 0 if O.K.) (integer)
!
! Author:  Martin Duncan  
! Date:    March 12/93
! Last revision: March 12/93

      subroutine drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)

!...  Inputs: 
      real*8 dt,r0,mu,alpha,u

!...  Outputs:
      integer iflg
      real*8 s

!...  Internals:
      real*8 denom,a0,a1,a2,q,r,sq2,sq,p1,p2

!----
!...  Executable code 

	denom = (mu - alpha*r0)/6.d0
	a2 = 0.5*u/denom
	a1 = r0/denom
	a0 =-dt/denom

	q = (a1 - a2*a2/3.d0)/3.d0
	r = (a1*a2 -3.d0*a0)/6.d0 - (a2**3)/27.d0
	sq2 = q**3 + r**2

	if( sq2  >=  0.d0) then
	   sq = sqrt(sq2)

	   if ((r+sq)  <=  0.d0) then
	      p1 =  -(-(r + sq))**(1.d0/3.d0)
	   else
	      p1 = (r + sq)**(1.d0/3.d0)
	   endif
	   if ((r-sq)  <=  0.d0) then
	      p2 =  -(-(r - sq))**(1.d0/3.d0)
	   else
	      p2 = (r - sq)**(1.d0/3.d0)
	   endif

	   iflg = 0
	   s = p1 + p2 - a2/3.d0

	else
	   iflg = 1
	   s = 0
	endif

        return
        end subroutine drift_kepu_p3solve
!-------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_KEPU_STUMPFF.F
!*************************************************************************
! subroutine for the calculation of stumpff functions
! see Danby p.172  equations 6.9.15
!
!             Input:
!                 x             ==>  argument
!             Output:
!                 c0,c1,c2,c3   ==>  c's from p171-172
!                                       (real scalors)
! Author:  Hal Levison  
! Date:    2/3/93
! Last revision: 2/3/93
! Modified by JEC: 31/3/98
!
      subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)
      use swift

!...  Inputs: 
      real*8 x

!...  Outputs:
      real*8 c0,c1,c2,c3

!...  Internals:
      integer n,i
      real*8 xm,x2,x3,x4,x5,x6

!----
!...  Executable code 

      n = 0
      xm = 0.1
      do while(abs(x) >= xm)
         n = n + 1
         x = x * .25d0
      enddo
!
      x2 = x  * x
      x3 = x  * x2
      x4 = x2 * x2
      x5 = x2 * x3
      x6 = x3 * x3
!
      c2 = 1.147074559772972d-11*x6 - 2.087675698786810d-9*x5 &
         + 2.755731922398589d-7*x4  - 2.480158730158730d-5*x3 &
         + 1.388888888888889d-3*x2  - 4.166666666666667d-2*x + .5d0
!
      c3 = 7.647163731819816d-13*x6 - 1.605904383682161d-10*x5 &
         + 2.505210838544172d-8*x4  - 2.755731922398589d-6*x3 &
         + 1.984126984126984d-4*x2  - 8.333333333333333d-3*x &
         + 1.666666666666667d-1
!
      c1 = 1. - x*c3
      c0 = 1. - x*c2
!
      if(n /= 0) then
         do i=n,1,-1
            c3 = (c2 + c0*c3)*.25d0
            c2 = c1*c1*.5d0
            c1 = c0*c1
            c0 = 2.*c0*c0 - 1.
            x = x * 4.
          enddo
       endif

       return
       end subroutine drift_kepu_stumpff
!------------------------------------------------------------------
!
!*************************************************************************
!                        DRIFT_ONE.F
!*************************************************************************
! This subroutine does the danby-type drift for one particle, using 
! appropriate vbles and redoing a drift if the accuracy is too poor 
! (as flagged by the integer iflg).
!
!             Input:
!                 nbod          ==>  number of massive bodies (int scalar)
!                 mass          ==>  mass of bodies (real array)
!                 x,y,z         ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 vx,vy,vz      ==>  initial position in jacobi coord 
!                                    (real scalar)
!                 dt            ==>  time step
!             Output:
!                 x,y,z         ==>  final position in jacobi coord 
!                                       (real scalars)
!                 vx,vy,vz      ==>  final position in jacobi coord 
!                                       (real scalars)
!                 iflg          ==>  integer (zero for successful step)
!
! Authors:  Hal Levison & Martin Duncan 
! Date:    2/10/93
! Last revision: 2/10/93
!

      subroutine drift_one(mu,x,y,z,vx,vy,vz,dt,iflg)
      use swift

!...  Inputs Only: 
      real*8 mu,dt

!...  Inputs and Outputs:
      real*8 x,y,z
      real*8 vx,vy,vz

!...  Output
	integer iflg
	
!...  Internals:
	integer i
	real*8 dttmp

!----
!...  Executable code 

           call drift_dan(mu,x,y,z,vx,vy,vz,dt,iflg)

	   if(iflg  /=  0) then
	    
	     do i = 1,10
	       dttmp = dt/10.d0
               call drift_dan(mu,x,y,z,vx,vy,vz,dttmp,iflg)
	       if(iflg  /=  0) return
	     enddo

	   endif

        return
        end subroutine drift_one
!-------------------------------------------------------------------
!
!**********************************************************************
!                    ORBEL_FGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_fget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!           Cel. Mech. ".  Quartic convergence from Danby's book.
!     REMARKS: 
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: 2/26/93 hfl
!     Modified by JEC
!***********************************************************************

      real*8 function orbel_fget(e,capn)
      use swift
      use interfaces, only: calc_sinh_cosh

!...  Inputs Only: 
	real*8 e,capn

!...  Internals:
	integer i,IMAX
	real*8 tmp,x,shx,chx
	real*8 esh,ech,f,fp,fpp,fppp,dx
	PARAMETER (IMAX = 10)

!----
!...  Executable code 

! Function to solve "Kepler's eqn" for F (here called
! x) for given e and CAPN. 

!  begin with a guess proposed by Danby	
	if( capn  <  0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0
	   x = -log(tmp)
	else
	   tmp = +2.d0*capn/e + 1.8d0
	   x = log( tmp)
	endif

	orbel_fget = x

	do i = 1,IMAX
          call calc_sinh_cosh (x,shx,chx)
	  esh = e*shx
	  ech = e*chx
	  f = esh - x - capn
!	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_fget = x + dx
!   If we have converged here there's no point in going on
	  if(abs(dx)  <=  SMALL_NUMBER) RETURN
	  x = orbel_fget
	enddo	

	write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	return
	end function orbel_fget
!------------------------------------------------------------------
!
!**********************************************************************
!                    ORBEL_FHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           n ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
!	         For larger N, uses FGET
!     REMARKS: 
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 26,1992.
!     REVISIONS: 
!     REVISIONS: 2/26/93 hfl
!**********************************************************************

	real*8 function orbel_fhybrid(e,n)
        use swift

!...  Inputs Only: 
	real*8 e,n

!...  Internals:
	real*8 abn

!----
!...  Executable code 

	abn = n
	if(n < 0.d0) abn = -abn

	if(abn  <  0.636d0*e -0.6d0) then
	  orbel_fhybrid = orbel_flon(e,n)
	else 
	  orbel_fhybrid = orbel_fget(e,n)
	endif   

	return
	end function orbel_fhybrid
!-------------------------------------------------------------------
!
!**********************************************************************
!                    ORBEL_FLON.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_flon ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Uses power series for N in terms of F and Newton,s method
!     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 26, 1992.
!     REVISIONS: 
!**********************************************************************

	real*8 function orbel_flon(e,capn)
        use swift

!...  Inputs Only: 
	real*8 e,capn

!...  Internals:
	integer iflag,i,IMAX
	real*8 a,b,sq,biga,bigb
	real*8 x,x2
	real*8 f,fp,dx
	real*8 diff
	real*8 a0,a1,a3,a5,a7,a9,a11
	real*8 b1,b3,b5,b7,b9,b11
	PARAMETER (IMAX = 10)
	PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
	PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
	PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
	PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

!----
!...  Executable code 


! Function to solve "Kepler's eqn" for F (here called
! x) for given e and CAPN. Only good for smallish CAPN 

	iflag = 0
	if( capn  <  0.d0) then
	   iflag = 1
	   capn = -capn
	endif

	a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
	a0 = -6227020800.d0*capn/e
	b1 = a1

!  Set iflag nonzero if capn < 0., in which case solve for -capn
!  and change the sign of the final answer for F.
!  Begin with a reasonable guess based on solving the cubic for small F	


	a = 6.d0*(e-1.d0)/e
	b = -6.d0*capn/e
	sq = sqrt(0.25*b*b +a*a*a/27.d0)
	biga = (-0.5*b + sq)**0.3333333333333333d0
	bigb = -(+0.5*b + sq)**0.3333333333333333d0
	x = biga + bigb
!	write(6,*) 'cubic = ',x**3 +a*x +b
	orbel_flon = x
! If capn is tiny (or zero) no need to go further than cubic even for
! e =1.
	if( capn  <  TINY) go to 100

	do i = 1,IMAX
	  x2 = x*x
	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
	  dx = -f/fp
!	  write(6,*) 'i,dx,x,f : '
!	  write(6,432) i,dx,x,f
432	  format(1x,i3,3(2x,1p,1e22.15))
	  orbel_flon = x + dx
!   If we have converged here there's no point in going on
	  if(abs(dx)  <=  TINY) go to 100
	  x = orbel_flon
	enddo	

! Abnormal return here - we've gone thru the loop 
! IMAX times without convergence
	if(iflag  ==  1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif
	write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	  diff = e*sinh(orbel_flon) - orbel_flon - capn
	  write(6,*) 'N, F, ecc*sinh(F) - F - N : '
	  write(6,*) capn,orbel_flon,diff
	return

!  Normal return here, but check if capn was originally negative
100	if(iflag  ==  1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif

	return
	end function orbel_flon
!
!**********************************************************************
!                    ORBEL_ZGET.F
!**********************************************************************
!     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
!          given Q (Fitz. notation.)
!
!             Input:
!                           q ==>  parabola mean anomaly. (real scalar)
!             Returns:
!                  orbel_zget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
!     REMARKS: For a parabola we can solve analytically.
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: May 27 - corrected it for negative Q and use power
!	      series for small Q.
!**********************************************************************

	real*8 function orbel_zget(q)
        use swift

!...  Inputs Only: 
	real*8 q

!...  Internals:
	integer iflag
	real*8 x,tmp

!----
!...  Executable code 

	iflag = 0
	if(q < 0.d0) then
	  iflag = 1
	  q = -q
	endif

	if (q < 1.d-3) then
	   orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
	else
	   x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
	   tmp = x**(1.d0/3.d0)
	   orbel_zget = tmp - 1.d0/tmp
	endif

	if(iflag  == 1) then
           orbel_zget = -orbel_zget
	   q = -q
	endif
	
	return
	end function orbel_zget




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!
!!!!   Added in by Joshua Wallace
!!!!
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!==========================================================================
! The following is borrowed from Mercury 6.2, function mco_x2ov
! This allows the proper fr and fv values to be calculated 
! which can be read by element6 and produce the correct output.

      subroutine mco_x2ov (m,x,v,fr,fv)
      use constants; use globals 
      implicit none
! Input/Output
      real(R8),  intent(in)::m
      real(R8),  intent(out)::fr,fv
      real(R8),  intent(in)::x(3),v(3)
!
! Local
      real(R8) r,v2,be,ke,temp
!
!------------------------------------------------------------------------------
!
        r = sqrt(x(1) * x(1)  +  x(2) * x(2) + x(3) * x(3))
        v2 =     v(1) * v(1)  +  v(2) * v(2) +  v(3) * v(3)
        be = G*(mcen + m) / r
        ke = HALF * v2
!
        fr = log10 (min(max(r, rcen), rmax) / rcen)
        temp = ke / be
        fv = ONE / (ONE + TWO*temp*temp)
!
!------------------------------------------------------------------------------
!
      end subroutine mco_x2ov


