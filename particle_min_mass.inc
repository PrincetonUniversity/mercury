c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      PARTICLE_MIN_MASS.INC    (Joshua Wallace   30 Dec 2015)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 
c
c  This defines a minimum mass for particles to be made
c  after a fragmentation occurs.
c  Similar to Chambers 2013
c   units of solar mass


      real*8 MINMASS

      parameter (MINMASS = 1.4d-8)