      program test
      implicit none

      integer num,num2,i
      num = 20
      num2=10
      real*8 mass(num)
      
      do i = 1,num
         mass(i) = 1.
      end do
      
      call foo(num2,mass)

      end


      subroutine foo(num2,mass)
      
      integer num2
      real*8 mass(num2)

      continue

      end
