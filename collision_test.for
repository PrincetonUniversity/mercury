c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 
c   Created by Joshua Wallace
c     on 11 Dec 2015
c 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c This is to test the collision outcomes before
c implementing it


c First the main loop

      program yayhooray
      implicit none
      include 'mercury.inc'

      integer j,algor,nbod,nbig,opt(8),stat(NMAX),lmem(NMESS),k
      integer opflag,ngflag,ndump,nfun
      integer coltype_num, num_tests, num_success
      real*8 m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
      real*8 rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
      real*8 cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
      real*8 rphys(NMAX), elost
      character*8 id(NMAX)
      character*80 outfile(3), dumpfile(4), mem(NMESS)
      character*80 outfilename,testname
      external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy
      external mco_dh2h,mco_h2dh
      external mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden
c
      data opt/1,2,3,2,0,1,0,0/

c     mem, lmem initialization
      mem(67) = " collided with the central body at "
      lmem(67) = 35
      mem(69) = " was hit by "
      lmem(69) = 12
      mem(71) = " at "
      lmem(71) = 4
      mem(2) = " years "
      mem(1) = " days  "

c     central object
      m(1) = 0.1
      vh(1,1) = 0.
      vh(2,1) = 0.
      vh(3,1) = 0.
      xh(1,1) = 0.
      xh(2,1) = 0.
      xh(3,1) = 0.
      rcen = 0.001

      nbod = 20
      nbig = nbod
      do k=2,nbig
         m(k) = 1.0E-5
      end do

      do k=2,nbig
         xh(1,k) = 0.1 + k*0.001
         xh(2,k) = 0.1 + k*0.001
         xh(3,k) = 0.1 + k*0.001
         vh(1,k) = 1.0E-19
         vh(2,k) = 0.
         vh(3,k) = 0.
         id(k) = char(k-1)
      end do


c Lets' test a basic merger.  Two objects nearly on top of each other, no relative vel

      xh(1,3) = xh(1,2)
      xh(2,3) = xh(2,2)
      xh(3,3) = xh(3,2) + 1.0E-8

      vh(1,3) = vh(1,2)
      vh(2,3) = vh(2,2)
      vh(3,3) = vh(3,2) 

      call phys_radii(nbod,rphys,m,rho_forall*( AU * AU * AU  / MSUN))

      outfilename = "testing_output.candel.txt"
      
      write(*,"(A20)")  "Beginning testing..."

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,1,3,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "collision with central object"
      call int_checker(coltype_num,-1,num_tests,
     %     num_success,testname)

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,2,3,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "basicmerger"
      call int_checker(coltype_num,1,num_tests,
     %     num_success,testname)



      m(4) = 5.682E-7
      m(5) = m(4) * 0.2
      vh(1,4) = 0.00577548327 ! AU/day
      vh(1,5) = 0.0
      xh(1,4) = xh(1,5) - 8.0E-6
      xh(2,4) = xh(2,5) - 7.2202E-6
      xh(3,4) = xh(3,5)

      call phys_radii(nbod,rphys,m,rho_forall*( AU * AU * AU  / MSUN))

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "hit and run"
      call int_checker(coltype_num,2,num_tests,
     %     num_success,testname)

      xh(1,4) = xh(1,5) - 1.0E-5
      xh(2,4) = xh(2,5) - 3.0E-7
      xh(3,4) = xh(3,5)

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "partial accretion"
      call int_checker(coltype_num,5,num_tests,
     %     num_success,testname)
      

      m(4) = 5.682E-7
      m(5) = m(4) * 0.2
      vh(1,4) = 0.0577548327 ! AU/day
      vh(1,5) = 0.0
      xh(1,4) = xh(1,5) - 8.0E-6
      xh(2,4) = xh(2,5) - 7.2202E-6
      xh(3,4) = xh(3,5)

      call phys_radii(nbod,rphys,m,rho_forall*( AU * AU * AU  / MSUN))

c      write(*,"(A20,E20.7, E15.7)") "  phys radii: ", rphys(4), m(4)

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "supercatastrophic"
      call int_checker(coltype_num,3,num_tests,
     %     num_success,testname)

      vh(1,4) = 0.01477548327 ! AU/day
      vh(1,5) = 0.0
      xh(1,4) = xh(1,5) - 1.35E-5
      xh(2,4) = xh(2,5) - 1.0E-5
      xh(3,4) = xh(3,5)

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "erosive disruption"
      call int_checker(coltype_num,4,num_tests,
     %     num_success,testname)

      vh(1,4) = 0.01377548327 ! AU/day

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "partial accretion"
      call int_checker(coltype_num,5,num_tests,
     %     num_success,testname)

      vh(1,4) = 0.01507548327 ! AU/day
      vh(1,5) = 0.0
      xh(1,4) = xh(1,5) - 1.245E-5

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "partial accretion, triple point"
      call int_checker(coltype_num,5,num_tests,
     %     num_success,testname)

      xh(1,4) = xh(1,5) - 1.210E-5

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "hit and run, triple point"
      call int_checker(coltype_num,2,num_tests,
     %     num_success,testname)

      vh(1,4) = 0.01597548327 ! AU/day

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "erosive disruption, triple point"
      call int_checker(coltype_num,4,num_tests,
     %     num_success,testname)

      vh(1,4) = 0.0371 ! AU/day

      write(*,"(A30)") "---------------------------"
      call mce_coll_frag(0.,0.,elost,0.,4,5,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfilename,coltype_num)
      testname = "supercat, close to boundary"
      call int_checker(coltype_num,3,num_tests,
     %     num_success,testname)




      write(*,"(A3)")  "   "
      write(*,"(A39)")  "************************************** "
      write(*,"(A39)")    "*********   Test results   *********** "

      if (num_success.eq.num_tests) then
         write(*,"(A27)") "   all passed, no problems!"
      else
         write(*,"(A15)") "   Warning! ###"
         write(*,"(A8,I3,A10,I3,A7)") " Failed ",  
     %   num_tests - num_success, " tests of ", num_tests, " total."
         write(*,"(A3)") "   "
         write(*,"(A45)") "Reminder: -1 cent. collision, 1 perf. merger"
         write(*,"(A44)") "2 hit & run, 3 supercatastrophic disruption"
         write(*,"(A45)") "4 is erosive disruption, 5 partial accretion"
      end if

      stop
      end


c Now the function

c available variables
c time,tstart,elost,jcen,i,j,nbod,nbig,m,xh,
c     %  vh,s,rphys,stat,id,opt,mem,lmem,outfile

c time - time of impact relative to end of timestep (from mce_cent)
c tstart - start time of integration
c elost  - energy lost    en(3)
c jcen   - non-sphericity of central object
c i      - first particle involved, index    ihit(k)
c j      - second particle involved, index   jhit(k)
c nbod   - total number of bodies
c nbig   - number of big bodies
c m      - mass array
c xh     - x array, relative to central body
c vh     - v array, relative to central body
c s      - spin angular momentum of the bodies
c rphys  - physical radii of objects       local array
c stat   - tells code whether to keep particle or not 
c                      (0 is keep, other is not keep)
c id     - the id array of the particles
c opt    - array storing various flags set when code begins
c                 opt(2) 0=none 1=merge  2=merge+fragment
c mem    - messages
c lmem   - messages
c outfile-  the outfile, outfile(3)

c if not the central body colliding, all the arrays that are passed
c  are passed as their "bs" versions.  No, I am not making that up.
c  that is literally what they are called. xbs, sbs, statbs...
c  probably stands for Bulrisch-Stoer

c The "bs" versions are abbreviated versions of the oriiginal versions of
c each array, chosen if their ce(j) value is not equal to 0.
c except for ibs and jbs, which are integer values referring to positions
c in the bs arrays, not the original arrays

c ############ In the above case, then collision won't have access
c  #### to the m, x, v, etc. arrays themselves for to adding new particles
c  ### I may have to create some special feature to allow for adding
c   ### new particles.
c   ### Or some new way of checking for collisions so that the full
c   ### arrays can be passed

c this function is called in a loop from do k=1, nhit


      subroutine mce_coll_frag (time,tstart,elost,jcen,i,j,nbod,nbig,
     %  m,xh,vh,s,rphys,stat,id,opt,mem,lmem,outfile,coltype_num)

      implicit none
      include 'mercury.inc'

c Input/Output
      integer i,j,nbod,nbig,stat(nbod),opt(8),lmem(NMESS),coltype_num
      real*8 time,tstart,elost,jcen(3)
      real*8 m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),rphys(nbod)
      character*80 outfile,mem(NMESS)
      character*8 id(nbod)

c Local
      integer year,month,itmp
      real*8 t1
      character*38 flost,fcol
      character*6 tstring
c
      real*8 gamma, mtot, l_, xrel(3), vrel(3), costheta_squared
      real*8 vrel_magnitude_squared, rhoforall_mercunits
      real*8 b_, alpha, M_, R_, vesc_squared, rhocgs, bcrit
      real*8 rc1,qpd,vpd,mu,muint,qrdstar,vstar,qrdstarprime
      real*8 vstarprime,qrer,ver_squred,qsupercat,vsupercat_squred
      real*8 rho1
      integer collision_type, graze, filestatus
c     -1 is central collision, 1 is perfect merger, 2 is hit & run
c     3 is supercatastrophic disruption, 4 is erosive disruption, 5 is partial accretion

      real*8 m_target, m_proj,r_target,r_proj

c     calculate rho in proper units, mercury units
      rhocgs = AU * AU * AU  / MSUN
      rhoforall_mercunits = rho_forall * rhocgs

c     rho1 also needs to be put into correct unit
      rho1 = rho1_cgs * rhocgs
      


c If two bodies collided, check that the less massive one is given 
c the larger index (unless the more massive one is a Small body)
      if (i.ne.0) then
        if (m(j).gt.m(i).and.j.le.nbig) then
          itmp = i
          i = j
          j = itmp
        end if
      end if
      

c  Write message to info file (I=0 implies collision with the central body)
c    Hey, I could add my own output to ce.out here!!!  And remove
c    the other output
 10   open (23, file=outfile, status='old', access='append', 
     %  iostat = filestatus, err=10)
c 11   write(*,*) "cant open file, error code ", filestatus
c
      if (opt(3).eq.1) then

        call mio_jd2y (time,year,month,t1)
        if (i.eq.0) then ! lost to the outer reaches of space
          flost = '(1x,a8,a,i10,1x,i2,1x,f8.5)'
c  mem(67) = 'collided with the central body at'
          write (23,flost) id(j),mem(67)(1:lmem(67)),year,month,t1
        else
          fcol  = '(1x,a8,a,a8,a,i10,1x,i2,1x,f4.1)'
c  mem(69) = 'was hit by', mem(71) = 'at'
          write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j),
     %      mem(71)(1:lmem(71)),year,month,t1
        end if
      else
        if (opt(3).eq.3) then
          t1 = (time - tstart) / 365.25d0
          tstring = mem(2)
          flost = '(1x,a8,a,f18.7,a)'
          fcol  = '(1x,a8,a,a8,a,1x,f14.3,a)'
        else
          if (opt(3).eq.0) t1 = time
          if (opt(3).eq.2) t1 = time - tstart
          tstring = mem(1)
          flost = '(1x,a8,a,f18.5,a)'
          fcol  = '(1x,a8,a,a8,a,1x,f14.1,a)'
        end if
        if (i.eq.0.or.i.eq.1) then ! Why 0 and 1?  When would 0 come up?
          write (23,flost) id(j),mem(67)(1:lmem(67)),t1,tstring
        else
          write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j),
     %      mem(71)(1:lmem(71)),t1,tstring
        end if
      end if
      close (23)

c Now we start to perform the analysis of Leinhardt & Stewart, and Stewart & Leinhardt

      if (i.eq.1) then  ! If central object
         call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
         collision_type = -1
c         goto 654 !Skip all the fragmentation nonsense that comes next
c      endif
      else

        gamma = m(j)/m(i)       ! mass particle / mass target greater
        mtot = (m(j) + m(i)) !* MSUN ! total mass
c        m_target = m(i) * MSUN
c        m_proj = m(j) * MSUN
c        r_target = rphys(i) * AU
c        r_proj = rphys(j) * AU

c        write(*,"(ES14.7E2)") m_target
c        write(*,"(ES14.7E2)") m_proj
c        write(*,"(ES14.7E2)") r_target
c        write(*,"(ES14.7E2)") r_proj


        xrel(1) = (xh(1,j) - xh(1,i))!*AU
        xrel(2) = (xh(2,j) - xh(2,i))!*AU
        xrel(3) = (xh(3,j) - xh(3,i))!*AU
        vrel(1) = (vh(1,j) - vh(1,i))!*AU/(24.0*3600.0)
        vrel(2) = (vh(2,j) - vh(2,i))!*AU/(24.0*3600.0)
        vrel(3) = (vh(3,j) - vh(3,i))!*AU/(24.0*3600.0)
        vrel_magnitude_squared = vrel(1)*vrel(1) + vrel(2)*vrel(2) + 
     %       vrel(3)*vrel(3)

        if (vrel_magnitude_squared.eq.(0.0)) then
           write(*,"(A41)") "Somehow, the relative velocity is zero!"
           stop 987
        end if

        costheta_squared = (xrel(1)*vrel(1) + xrel(2)*vrel(2) + xrel(3)*
     %   vrel(3))**2 / ( (xrel(1)*xrel(1) + xrel(2)*xrel(2) + xrel(3)*
     %   xrel(3) ) * vrel_magnitude_squared ) ! Used for next calculation

        b_ = sqrt(1.0 - costheta_squared) ! Impact parameter, sin(theta)
c        write(*,"(A6,F8.6)") "  b_: ", b_
        l_ = (rphys(i) + rphys(j)) * (1.0 - b_) ! length (absolute in CGS) of projecticle that overlaps the target
        alpha = (3.0*rphys(j)*(l_*l_) - (l_*l_*l_))/(4.0*(rphys(j)*
     %       rphys(j)*rphys(j) )) ! Intersecting mass fraction

        if (rphys(i).gt.(b_*(rphys(i)+rphys(j) ) +rphys(j)) ) alpha=1.0 !Maximum alpha
        M_ = m(i) + alpha*m(j)  !mass of target + interacting mass of projectile
        R_ = ( (3.*M_)/(4.*PI*rhoforall_mercunits) )**(1./3.) !radius that target + interacting mass would have
        vesc_squared = 2.*K2*M_/R_ !escape velocity of target + interacting mass eq. 53
        write(*,"(A6,F12.1)") "Vesc: ",
     %   sqrt(vesc_squared)*AU/(24.0*3600.0)
        if (vrel_magnitude_squared.lt.vesc_squared) then
           call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
           collision_type = 1

        else
           bcrit = rphys(i)/(rphys(i) + rphys(j))

           if (b_ < bcrit) then
              graze = 0
           else
              graze = 1
           endif
           write(*,"(I3)") graze
           
           rc1 = (3.0*mtot/(4.*PI*rho1))**(1./3.) !radius of all mass if rho = 1
           qpd = cstar*4.0/5.0*PI*rho1*K2*rc1*rc1 !eq. 28, specific
c           ! impact energy when mt=mp
           vpd = sqrt(32.*PI*cstar*rho1*K2/5.)*rc1 !eq. 30, vel version of above
           mu  = (m(i) * m(j) / mtot) ! reduced mass
           muint = alpha*m(i)*m(j)/(m(i) + alpha*m(j))

           qrdstar = qpd*(   (((gamma+1.0)*(gamma+1.0))/(4.0*
     %          gamma))**(2.0/(3.0*mubar)-1.0)   )
c           eq. 23, specific impact energy-catastrophic disruption threshold
           vstar = vpd*(    (((gamma+1.0)*(gamma+1.0))/(4.0*
     %          gamma))**(1.0/(3.0*mubar))    )
c     eq. 22, velocity at catastrophic disruption threshold

           qrdstarprime = qrdstar* (mu/muint)**(2.-(1.5*mubar))
c eq. 15, specific impact energy at catastrophic disruption threshold for oblique (b>0) impacts
           vstarprime = sqrt(2.*qrdstarprime*mtot/mu)
c     eq. 16, impact velocity at catastrophic disruption threshold for oblique (b>0) impacts

           qrer = qrdstarprime*((-2.*m(i)/mtot)+2.0)
c  eq. 5(rearranged), specific impact energy at erosion threshold
           ver_squred = 2.0*qrer*mtot/mu
c     eq. 1(rearranged), velocity at erosion threshold

           if (graze.eq.1.and.vrel_magnitude_squared.lt.ver_squred) then
c     In this case, hit and run regime. Target intact but projectile may be disrupted
              collision_type = 2

           else  ! If not hit and run regime
              qsupercat = 1.8*qrdstarprime ! super-cat specific impact energy
              vsupercat_squred = 2.0*qsupercat*mtot/mu ! super-cat impact velocity
              write(*,"(A19,E20.7)") "V_supercat: ", 
     %         sqrt(vsupercat_squred)*AU/(24.0*3600.0)
              if (vrel_magnitude_squared.gt.ver_squred) then
                 if (vrel_magnitude_squared.gt.vsupercat_squred) then
                    collision_type = 3
                 else
                    collision_type = 4
                 end if         ! this if is checking erosion case, supercat or not
                 
              else
                 collision_type = 5 ! partial accretion

              end if ! this if is checking if erosion regime or not
              
           end if ! this if is checking hit and run
           
        end if ! this if is checking if perfect merger
      end if ! This if is checking if central object
 654  continue
      coltype_num = collision_type
      return
      end










c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------
c ----------------------------------------------------------------









c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCE_MERG.FOR    (ErikSoft   2 October 2000)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c
c Author: John E. Chambers
c
c Merges objects I and J inelastically to produce a single new body by 
c conserving mass and linear momentum.
c   If J <= NBIG, then J is a Big body
c   If J >  NBIG, then J is a Small body
c   If I = 0, then I is the central body
c
c N.B. All coordinates and velocities must be with respect to central body.
c ===
c
c------------------------------------------------------------------------------
c
      subroutine mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer i, j, nbod, nbig, stat(nbod)
      real*8 jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),elost
c
c Local
      integer k
      real*8 tmp1, tmp2, dx, dy, dz, du, dv, dw, msum, mredu, msum_1
      real*8 e0, e1, l2
c
c------------------------------------------------------------------------------
c
c If a body hits the central body
      if (i.le.1) then
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e0,l2)
c
c If a body hit the central body...
        msum   = m(1) + m(j)
        msum_1 = 1.d0 / msum
        mredu  = m(1) * m(j) * msum_1
        dx = xh(1,j)
        dy = xh(2,j)
        dz = xh(3,j)
        du = vh(1,j)
        dv = vh(2,j)
        dw = vh(3,j)
c
c Calculate new spin angular momentum of the central body
        s(1,1) = s(1,1)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
        s(2,1) = s(2,1)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
        s(3,1) = s(3,1)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
c
c Calculate shift in barycentric coords and velocities of central body
        tmp2 = m(j) * msum_1
        xh(1,1) = tmp2 * xh(1,j)
        xh(2,1) = tmp2 * xh(2,j)
        xh(3,1) = tmp2 * xh(3,j)
        vh(1,1) = tmp2 * vh(1,j)
        vh(2,1) = tmp2 * vh(2,j)
        vh(3,1) = tmp2 * vh(3,j)
        m(1) = msum
        m(j) = 0.d0
        s(1,j) = 0.d0
        s(2,j) = 0.d0
        s(3,j) = 0.d0
c
c Shift the heliocentric coordinates and velocities of all bodies
        do k = 2, nbod
          xh(1,k) = xh(1,k) - xh(1,1)
          xh(2,k) = xh(2,k) - xh(2,1)
          xh(3,k) = xh(3,k) - xh(3,1)
          vh(1,k) = vh(1,k) - vh(1,1)
          vh(2,k) = vh(2,k) - vh(2,1)
          vh(3,k) = vh(3,k) - vh(3,1)
        end do
c
c Calculate energy loss due to the collision
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e1,l2)
        elost = elost + (e0 - e1)
      else
c
c If two bodies collided...
        msum   = m(i) + m(j)
        msum_1 = 1.d0 / msum
        mredu  = m(i) * m(j) * msum_1
        dx = xh(1,i) - xh(1,j)
        dy = xh(2,i) - xh(2,j)
        dz = xh(3,i) - xh(3,j)
        du = vh(1,i) - vh(1,j)
        dv = vh(2,i) - vh(2,j)
        dw = vh(3,i) - vh(3,j)
c
c Calculate energy loss due to the collision
        elost = elost  +  .5d0 * mredu * (du*du + dv*dv + dw*dw)
     %        -  m(i) * m(j) / sqrt(dx*dx + dy*dy + dz*dz)
c
c Calculate spin angular momentum of the new body
        s(1,i) = s(1,i)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
        s(2,i) = s(2,i)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
        s(3,i) = s(3,i)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
c
c Calculate new coords and velocities by conserving centre of mass & momentum
        tmp1 = m(i) * msum_1
        tmp2 = m(j) * msum_1
        xh(1,i) = xh(1,i) * tmp1  +  xh(1,j) * tmp2
        xh(2,i) = xh(2,i) * tmp1  +  xh(2,j) * tmp2
        xh(3,i) = xh(3,i) * tmp1  +  xh(3,j) * tmp2
        vh(1,i) = vh(1,i) * tmp1  +  vh(1,j) * tmp2
        vh(2,i) = vh(2,i) * tmp1  +  vh(2,j) * tmp2
        vh(3,i) = vh(3,i) * tmp1  +  vh(3,j) * tmp2
        m(i) = msum
      end if
c
c Flag the lost body for removal, and move it away from the new body
      stat(j) = -2
      xh(1,j) = -xh(1,j)
      xh(2,j) = -xh(2,j)
      xh(3,j) = -xh(3,j)
      vh(1,j) = -vh(1,j)
      vh(2,j) = -vh(2,j)
      vh(3,j) = -vh(3,j)
      m(j)   = 0.d0
      s(1,j) = 0.d0
      s(2,j) = 0.d0
      s(3,j) = 0.d0
c
c------------------------------------------------------------------------------
c
      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     MXX_EN.FOR    (ErikSoft   21 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates the total energy and angular-momentum for a system of objects
c with masses M, coordinates X, velocities V and spin angular momenta S.
c
c N.B. All coordinates and velocities must be with respect to the central
c ===  body.
c
c------------------------------------------------------------------------------
c
      subroutine mxx_en  (jcen,nbod,nbig,m,xh,vh,s,e,l2)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      integer nbod,nbig
      real*8 jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),e,l2
c
c Local
      integer j,k,iflag,itmp(8)
      real*8 x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r2,tmp,ke,pe,l(3)
      real*8 r_1,r_2,r_4,r_6,u2,u4,u6,tmp2(4,NMAX)
c
c------------------------------------------------------------------------------
c
      ke = 0.d0
      pe = 0.d0
      l(1) = 0.d0
      l(2) = 0.d0
      l(3) = 0.d0
c
c Convert to barycentric coordinates and velocities
      call mco_h2b(temp,jcen,nbod,nbig,temp,m,xh,vh,x,v,tmp2,iflag,itmp)
c
c Do the spin angular momenta first (probably the smallest terms)
      do j = 1, nbod
        l(1) = l(1) + s(1,j)
        l(2) = l(2) + s(2,j)
        l(3) = l(3) + s(3,j)
      end do
c
c Orbital angular momentum and kinetic energy terms
      do j = 1, nbod
        l(1) = l(1)  +  m(j)*(x(2,j) * v(3,j)  -  x(3,j) * v(2,j))
        l(2) = l(2)  +  m(j)*(x(3,j) * v(1,j)  -  x(1,j) * v(3,j))
        l(3) = l(3)  +  m(j)*(x(1,j) * v(2,j)  -  x(2,j) * v(1,j))
        ke = ke + m(j)*(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j))
      end do
c
c Potential energy terms due to pairs of bodies
      do j = 2, nbod
        tmp = 0.d0
        do k = j + 1, nbod
          dx = x(1,k) - x(1,j)
          dy = x(2,k) - x(2,j)
          dz = x(3,k) - x(3,j)
          r2 = dx*dx + dy*dy + dz*dz
          if (r2.ne.0) tmp = tmp + m(k) / sqrt(r2)
        end do
        pe = pe  -  tmp * m(j)
      end do
c
c Potential energy terms involving the central body
      do j = 2, nbod
        dx = x(1,j) - x(1,1)
        dy = x(2,j) - x(2,1)
        dz = x(3,j) - x(3,1)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2.ne.0) pe = pe  -  m(1) * m(j) / sqrt(r2)
      end do
c
c Corrections for oblateness
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        do j = 2, nbod
          r2 = xh(1,j)*xh(1,j) + xh(2,j)*xh(2,j) + xh(3,j)*xh(3,j)
          r_1 = 1.d0 / sqrt(r2)
          r_2 = r_1 * r_1
          r_4 = r_2 * r_2
          r_6 = r_4 * r_2
          u2 = xh(3,j) * xh(3,j) * r_2
          u4 = u2 * u2
          u6 = u4 * u2
          pe = pe + m(1) * m(j) * r_1
     %       * (jcen(1) * r_2 * (1.5d0*u2 - 0.5d0)
     %       +  jcen(2) * r_4 * (4.375d0*u4 - 3.75d0*u2 + .375d0)
     %       +  jcen(3) * r_6
     %       *(14.4375d0*u6 - 19.6875d0*u4 + 6.5625d0*u2 - .3125d0))
        end do
      end if
c
      e = .5d0 * ke  +  pe
      l2 = sqrt(l(1)*l(1) + l(2)*l(2) + l(3)*l(3))
c
c------------------------------------------------------------------------------
c
      return	
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_H2B.FOR    (ErikSoft   2 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts coordinates with respect to the central body to barycentric
c coordinates.
c
c------------------------------------------------------------------------------
c
      subroutine mco_h2b (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag,
     %  opt)
c
      implicit none
c
c Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real*8 time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
      real*8 v(3,nbod),ngf(4,nbod)
c
c Local
      integer j
      real*8 mtot,temp
c
c------------------------------------------------------------------------------
c
      mtot = 0.d0
      x(1,1) = 0.d0
      x(2,1) = 0.d0
      x(3,1) = 0.d0
      v(1,1) = 0.d0
      v(2,1) = 0.d0
      v(3,1) = 0.d0
c
c Calculate coordinates and velocities of the central body
      do j = 2, nbod
        mtot = mtot  +  m(j)
        x(1,1) = x(1,1)  +  m(j) * xh(1,j)
        x(2,1) = x(2,1)  +  m(j) * xh(2,j)
        x(3,1) = x(3,1)  +  m(j) * xh(3,j)
        v(1,1) = v(1,1)  +  m(j) * vh(1,j)
        v(2,1) = v(2,1)  +  m(j) * vh(2,j)
        v(3,1) = v(3,1)  +  m(j) * vh(3,j)
      enddo
c
      temp = -1.d0 / (mtot + m(1))
      x(1,1) = temp * x(1,1)
      x(2,1) = temp * x(2,1)
      x(3,1) = temp * x(3,1)
      v(1,1) = temp * v(1,1)
      v(2,1) = temp * v(2,1)
      v(3,1) = temp * v(3,1)
c
c Calculate the barycentric coordinates and velocities
      do j = 2, nbod
        x(1,j) = xh(1,j) + x(1,1)
        x(2,j) = xh(2,j) + x(2,1)
        x(3,j) = xh(3,j) + x(3,1)
        v(1,j) = vh(1,j) + v(1,1)
        v(2,j) = vh(2,j) + v(2,1)
        v(3,j) = vh(3,j) + v(3,1)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_JD2Y.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts from Julian day number to Julian/Gregorian Calendar dates, assuming
c the dates are those used by the English calendar.
c
c Algorithm taken from `Practical Astronomy with your calculator' (1988)
c by Peter Duffett-Smith, 3rd edition, C.U.P.
c
c Algorithm for negative Julian day numbers (Julian calendar assumed) by
c J. E. Chambers.
c
c N.B. The output date is with respect to the Julian Calendar on or before
c ===  4th October 1582, and with respect to the Gregorian Calendar on or 
c      after 15th October 1582.
c
c
c------------------------------------------------------------------------------
c
      subroutine mio_jd2y (jd0,year,month,day)
c
      implicit none
c
c Input/Output
      integer year,month
      real*8 jd0,day
c
c Local
      integer i,a,b,c,d,e,g
      real*8 jd,f,temp,x,y,z
c
c------------------------------------------------------------------------------
c
      if (jd0.le.0) goto 50
c
      jd = jd0 + 0.5d0
      i = sign( dint(dabs(jd)), jd )
      f = jd - 1.d0*i
c
c If on or after 15th October 1582
      if (i.gt.2299160) then
        temp = (1.d0*i - 1867216.25d0) / 36524.25d0
        a = sign( dint(dabs(temp)), temp )
        temp = .25d0 * a
        b = i + 1 + a - sign( dint(dabs(temp)), temp )
      else
        b = i
      end if
c
      c = b + 1524
      temp = (1.d0*c - 122.1d0) / 365.25d0
      d = sign( dint(dabs(temp)), temp )
      temp = 365.25d0 * d
      e = sign( dint(dabs(temp)), temp )
      temp = (c-e) / 30.6001d0
      g = sign( dint(dabs(temp)), temp )
c
      temp = 30.6001d0 * g
      day = 1.d0*(c-e) + f - 1.d0*sign( dint(dabs(temp)), temp )
c
      if (g.le.13) month = g - 1
      if (g.gt.13) month = g - 13
c
      if (month.gt.2) year = d - 4716
      if (month.le.2) year = d - 4715
c
      if (day.gt.32) then
        day = day - 32
        month = month + 1
      end if
c
      if (month.gt.12) then
        month = month - 12
        year = year + 1
      end if
      return
c
  50  continue
c
c Algorithm for negative Julian day numbers (Duffett-Smith doesn't work)
      x = jd0 - 2232101.5
      f = x - dint(x)
      if (f.lt.0) f = f + 1.d0
      y = dint(mod(x,1461.d0) + 1461.d0)
      z = dint(mod(y,365.25d0))
      month = int((z + 0.5d0) / 30.61d0)
      day = dint(z + 1.5d0 - 30.61d0*dble(month)) + f
      month = mod(month + 2, 12) + 1
c
      year = 1399 + int (x / 365.25d0)
      if (x.lt.0) year = year - 1
      if (month.lt.3) year = year + 1
c
c------------------------------------------------------------------------------
c
      return
      end


c
c
c    calculate the physical radii of all the objects
      subroutine phys_radii(nbod,rphys,m,rho)

      implicit none
      include 'mercury.inc'

      integer nbod
      real*8 rphys(nbod),m(nbod),rho
      integer k ! local

      do k=2, nbod
         rphys(k) = (3.*m(k)/(4.*PI*rho))**(1./3.)
      end do

c      write(*,"(A20,E20.5)") "   rho is: ", rho
      return
      end

c
c
c
c
c
c       For the checker test, the function that does the checking
c


      subroutine int_checker(inttocheck,inttocheckagainst,numberoftests,
     % numberofsuccess,testname)

      implicit none

      integer inttocheck,inttocheckagainst,numberoftests,numberofsuccess
      character*100 testname

      if (inttocheck.eq.inttocheckagainst) then
         numberoftests = numberoftests + 1
         numberofsuccess = numberofsuccess + 1
      else
         numberoftests = numberoftests + 1
         write(*,"(A26,A20)" ) "*****  Failed test! Name: ", testname
c         write(*,"(A8)") "--------"
         write(*,"(A4,I1,A8,I1)") "    ",inttocheck," is not ", 
     %   inttocheckagainst
      end if


      return
      end
