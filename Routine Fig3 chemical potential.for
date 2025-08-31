********************************************************************
*     Calculate the chemical potential in a nodal-ring semimetal 
********************************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xMbar=0.53        ! Dirac mass 
      real(8), parameter :: gamma=0.7 
cc    real(8), parameter :: xmu=0.8           ! chemical potential 
      real(8), parameter :: density0=894.6    ! the chosen value 204.2, 894.6 
      integer, parameter :: Ncut=3000         ! cutoff
      real(8) en(Ncut+1)                      ! store energy 
      integer enindex(Ncut+1)                 ! store index  
cc      
      end module value
cc      
**********************
*     main program
**********************
      program main
      use value
      implicit double precision (a-h,o-z)
      character(100) filename1,filename2,filename3 
      write(filename1,'(a100)') 'LLs1.dat'
      write(filename2,'(a100)') 'LLs2.dat'
      write(filename3,'(a100)') 'mu=1.4b.dat'
      open(10,file=filename1)
      open(20,file=filename2)
      open(30,file=filename3)
cc
      do 101 xB=0.001,0.011,0.001   
      write(*,*) xB 
      do i=0,Ncut    
         aa=2*i*xB+xB-1
         en(i+1)=sqrt(aa*aa+xMbar*xMbar)          
         enindex(i+1)=i 
      end do 
cc
      do i=0,Ncut-1                         ! order the LL energy from low to high  
      do j=i+1,Ncut
         if(en(j+1)<en(i+1)) then
            temp1=en(j+1) 
            en(j+1)=en(i+1)
            en(i+1)=temp1
cc
cc          temp2=enindex(j+1)
cc          enindex(j+1)=enindex(i+1) 
cc          enindex(i+1)=temp2
         end if
      end do
      end do
cc    
      do 100 xmu=xMbar+0.001,5.1,0.0001
      density=0
      do 200 i=0,Ncut
         if(en(i+1)<xmu) then
            bb=sqrt(xmu**2-en(i+1)**2)
            density=density+1.105*xB*655.36*bb
         end if
cc         
       if(density>density0) exit         
       if(en(i+1)>xmu) exit
200   continue
cc
      if(density>density0) then 
         write(30,'(2f16.8,I4)') xB,xmu,i
         exit 
      end if 
100   continue 
cc
101   continue
cc       
      end program main
cc