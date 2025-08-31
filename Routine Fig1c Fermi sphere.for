*****************************************************************
*     Listplot the 3D dispersion of the nodal-ring semimetals
*****************************************************************
      module value
      implicit double precision(a-h,o-z)
	real(8), parameter :: Pi=3.14159
      real(8), parameter :: xMbar=0.53  
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
      open(10,file='en=0.8a.dat')
      open(11,file='en=0.8b.dat')
      open(20,file='en=1.4a.dat')
      open(21,file='en=1.4b.dat')      
cc
      en=0.8
      zk=0 
      do 100 xk=0,2,0.01
         aa=en**2-xMbar**2-gamma*gamma*zk*zk
         if(aa>=0) then
            bb=sqrt(aa)+1-xk**2 
            if(bb>=0) then 
               yk=sqrt(bb)
               write(10,'(3f16.8)') xk,yk,zk  
cc               write(10,'(3f16.8)') xk,-yk,zk                 
            end if    
         end if   
cc
100   continue      
cc
      end program main
cc