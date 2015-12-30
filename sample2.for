      program test
      implicit none
      integer i
      character*8 name

      i = 1
      write(name,"(I0)") i
      write(*,"(A8,A1)") name,"H"

      end

