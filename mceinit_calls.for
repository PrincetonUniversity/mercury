
c First call right in mal_hcon
c Right after initializing variables
c Calculate close-encounter limits and physical radii
      call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %  m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)


c after onestep

c Second call under "Collisions"
c
c If collisions occurred, output details and remove lost objects
      if (colflag.ne.0) then

c .... skip some code

        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %    m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
c ....
      end if



c then check for collisions with the central body


c and then the code says:
c "Remove lost objects, reset flags and recompute Hill and physical radii"

c   mxx_elim is called first, presumably to remove lost objects

        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig,
     %    m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
