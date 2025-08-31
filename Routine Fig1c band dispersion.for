*****************************************************************
*     Listplot the 3D dispersion of the nodal-ring semimetals
*****************************************************************
      module value
      implicit double precision(a-h,o-z)
	real(8), parameter :: Pi=3.14159
      real(8), parameter :: xM=0.05           ! Dirac mass
      real(8), parameter :: epsilonr=0.094 
      real(8), parameter :: xMbar=xM/epsilonr 
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
      open(1,file='11.dat')
      open(2,file='12.dat')
cc
      do 100 xk=-3,3,0.01
         temp=(xk*xk-1)**2+xMbar**2
         en1=sqrt(temp)
         en2=-en1         
         write(1,'(100f16.8)') xk,en1,en2 
100   continue      
cc
      end program main
cc
cc