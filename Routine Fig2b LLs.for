*****************************************************
*     Calculate the LLs of a nodal-ring semimetal 
*****************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xMbar=0.53        ! Dirac mass
      real(8), parameter :: gamma=0.7         ! dimensionless Fermi velocity 
      real(8), parameter :: xB=1./40          ! magnetic field
      integer, parameter :: Ncut=50
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
      write(filename,'(a2,f6.3,a4)') 'B=',xB,'.dat'
      open(10,file=filename)        
cc
      do 100 zk=-2,2,0.01
      do 200 n=0,Ncut    
         aa=2*n*xB+xB-1
         en(n+1)=sqrt(aa*aa+xMbar**2+gamma*gamma*zk*zk)          
200   continue      
cc 
      write(10,'(1000f16.8)') zk,(en(n),n=1,Ncut+1)
100   continue          
cc      
      end program main
cc