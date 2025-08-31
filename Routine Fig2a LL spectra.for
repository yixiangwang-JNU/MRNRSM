*****************************************************
*     Calculate the LLs of a nodal-ring semimetal 
*****************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xM=0.05         ! Dirac mass
      real(8), parameter :: epsilonr=0.094 
      real(8), parameter :: xMbar=xM/epsilonr 
      real(8), parameter :: gamma=0.7 
      integer, parameter :: Ncut=15
      real(8) en(Ncut+1)
cc      
      end module value
cc      
**********************
*     main program
**********************
      program main
      use value
      implicit double precision (a-h,o-z)
      character(100) filename 
      write(filename,'(a100)') 'LLs.dat'
      open(10,file=filename)        
cc
      do 100 xB=0,1.0001,0.001 
      do 200 n=0,Ncut    
         aa=2*n*xB+xB-1
         en(n+1)=sqrt(aa*aa+xMbar*xMbar)          
200   continue      
cc 
      write(10,'(100f16.8)') xB,(en(n),n=1,Ncut+1)
100   continue          
cc      
      end program main
cc