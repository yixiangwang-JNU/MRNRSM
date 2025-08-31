*****************************************************************
*     Listplot the 3D dispersion of the nodal-ring semimetals
*****************************************************************
      module value
      implicit double precision(a-h,o-z)
	real(8), parameter :: Pi=3.1416
      real(8), parameter :: xM=0.05         ! Dirac mass
      real(8), parameter :: epsilonr=0.094 
      real(8), parameter :: xMbar=0.53      ! xM/epsilonr 
      real(8), parameter :: gamma=0.7 
cc
      end module value
cc
**********************
*     main program
**********************
      program main
	use value
      implicit double precision (a-h,o-z)
      open(10,file='10.dat')
cc
      aa=sqrt(xMbar**2+1)  
      xn=0 
      do 100 xmu=0,3,0.01
         if(xmu>xMbar .and. xmu<aa) then
            xn=xmu**2-xMbar**2
         else if(xmu>=aa) then
            xn1=xmu**2-xMbar**2
            xn2=sqrt(xmu**2-aa**2)/Pi
     &          -(xmu**2-xMbar**2)*atan(sqrt(xmu**2-aa**2))/Pi
            xn=xn1+xn2  
         end if 
cc
         xn=xn/(8*Pi*gamma)
         write(10,'(100f16.8)') xmu, xn 
100   continue  
cc
      end program main
cc