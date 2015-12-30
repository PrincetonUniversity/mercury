      program test
      implicit none
      include 'sample.inc'

      integer num2,i
c      num = 20
      real*8 mass(num)
      character name(80)

      num2=10
      
      do i = 1,num
         mass(i) = 1.
      end do
      
      call foo(num2,mass)

      write(*,"(A30,E15.9)") "19th one, in main function:",mass(19)

      call state_saver(i)
      write(*,"(A30,I8)") "statesaver:",i
      call state_saver(i)
      write(*,"(A30,I8)") "statesaver:",i
      call state_saver(i)
      write(*,"(A30,I8)") "statesaver:",i
      call state_saver(i)
      write(*,"(A30,I8)") "statesaver:",i


      i = 1
      write(name,"(I50)") i
      write(*,"(A32)") name

      end


      subroutine foo(num2,mass)
      
      integer num2
      real*8 mass(num2)
      write(*,"(A20,E15.9)") "first one:",mass(1)
      write(*,"(A20,E15.9)") "tenth one:",mass(10)
      write(*,"(A20,E15.9)") "twentieth one:",mass(20)
c      write(*,"(A20,E15.9)") "twenty-second one:",mass(22)

      mass(19) = 2.0


      continue

      end

      subroutine state_saver(in)

      integer in,state
      data state /250/
      save state 
      in = state
      state = state+1

      return
      end
