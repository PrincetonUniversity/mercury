

c first call in main loop of mal_hcon,
c under "Convert to heliocentric coordinates and output data for all bodies"
c which is inside an if statement, asking "Is it time for output ?"
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,
     %    0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,
     %    outfile,nstored,0)


c second call in main loop of mal_hcon,
c under "If encounter minima occurred, output details and decide whethertostop"
c I decided to include whole loop it was in
c this loop is ran immediately after onestep is called
      if (nclo.gt.0.and.opflag.ge.-1) then
        itmp = 1
        if (colflag.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo,
     %    iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,
     %    outfile,nstored,itmp)
        if (stopflag.eq.1) return
      end if

c  third call in main loop of mal_hcon, again
c  under "DATA  DUMP  AND  PROGRESS  REPORT
c
c Convert to heliocentric coords and do the data dump"
      if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
c some irrelevant code skipped
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,
     %    0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem,
     %    outfile,nstored,0)


c  Ejections and periodic effects calculations follow this,
c  followed by goto 100
