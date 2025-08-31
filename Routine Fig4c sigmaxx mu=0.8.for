*****************************************************
*     Calculate sigmaxx of a nodal-ring semimetal 
*****************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xMbar=0.53          ! Dirac mass
      real(8), parameter :: gamma=0.7 
      real(8), parameter :: xlower=-2           ! lower limit of integration
      real(8), parameter :: xupper=2            ! upper limit of integration      
      complex(8), parameter :: Im=(0,1.d0)     
cc     
      real(8) sigmaxx(600),etaArray(8) 
      complex(8) cni1,cnj1 
      common cni1,cni2,cnj1,cnj2,
     &       Binv,xB,xmu,eta, 
     &       eni,enj,  
     &       ni,nj,is,js,Ncut                    ! index for the initial and final LLs  
      end module value
cc      
**********************
*     main program
**********************
      program main
      use value
      implicit double precision (a-h,o-z)
      integer, parameter :: Num=19971 
      character(100) filename 
      real(8) Array(Num,3)               ! store inverseB and chemical potential  
      open(10,file='mu=0.8b.dat')
cc
      etaArray(1)=0.01
      etaArray(2)=0.4
      etaArray(3)=0.6
      etaArray(4)=0.7 
      etaArray(5)=0.8
      etaArray(6)=0.9
      etaArray(7)=1       
cc             
      do i=1,Num
         read(10,*) Array(i,1),Array(i,2),Array(i,3)     ! read chemical potential and Ncut
      end do      
cc
      do 101 j=1,1
      eta=etaArray(j) 
      write(filename,'(a12,f6.3,a4)') 'sigmaxx,eta=',eta,'.dat' 
      open(20,file=filename)       
cc
      do 100 i=1,Num                         ! number of magnetic field 
      Binv=Array(i,1)
      xB=1./Binv 
      xmu=Array(i,2)
      Ncut=Array(i,3)
cc 
cc    transition: n,s->n+1,s'  
      do ni=0,Ncut+200                       ! from the zeroth LL  
         nj=ni+1
         call solve_fn(Res)
         sigmaxx(ni+1)=Res 
      end do    
cc
      sum=0
      do ni=0,Ncut+200
         sum=sum+sigmaxx(ni+1) 
      end do
      write(20,'(100f16.8)') Binv,sum
cc 
100   continue     
101   continue 
cc 
      end program main
cc
*****************************
*     subroutine solve_fn
*****************************
      subroutine solve_fn(Res)
      include 'link_fnl_shared.h'
      use QDAG_INT      
	use value
      implicit double precision (a-h,o-z)
      external fn 
cc   
      Res=0
      do 100 is=-1,1,2
      do 200 js=-1,1,2
         call QDAG(fn,xlower,xupper,Resn)          
         Resn=Resn*2*xB*xB*eta*eta/(pi*pi) 
         Res=Res+Resn
200   continue
100   continue
cc      
      end subroutine solve_fn
cc      
**********************
*     function: fn
**********************
      real function fn(zk)
      use value
      implicit double precision (a-h,o-z)  
      complex(8) temp1
cc
      call solve_ni(zk)
      call solve_nj(zk) 
cc            
      temp1=cni2*cnj1+conjg(cni1)*cnj2
      temp2=(ni+1)*abs(temp1)**2        
      temp3=(xmu-eni)**2+eta**2
      temp4=(xmu-enj)**2+eta**2
      fn=temp2/(temp3*temp4)
cc
      end function fn      
cc 
*********************************************
*     subroutine: solve the initial ni LL
********************************************* 
      subroutine solve_ni(zk)
      use value
      implicit double precision (a-h,o-z)      
cc
      aa=(2*ni*xB+xB-1)**2+xMbar**2+(gamma*zk)**2
      en=sqrt(aa)
      eni=is*en 
cc
      cni1=2*ni*xB+xB-1-Im*gamma*zk
      cni2=-(xMbar+is*en)
      bb=sqrt(abs(cni1)**2+cni2**2) 
      cni1=cni1/bb
      cni2=cni2/bb       
cc      
      end subroutine solve_ni      
cc       
********************************************
*     subroutine: solve the final nj LL
******************************************** 
      subroutine solve_nj(zk)
      use value
      implicit double precision (a-h,o-z)      
cc
      aa=(2*nj*xB+xB-1)**2+xMbar**2+(gamma*zk)**2
      en=sqrt(aa)
      enj=js*en      
cc
      cnj1=2*nj*xB+xB-1-Im*gamma*zk
      cnj2=-(xMbar+js*en)
      bb=sqrt(abs(cnj1)**2+cnj2**2) 
      cnj1=cnj1/bb
      cnj2=cnj2/bb       
cc      
      end subroutine solve_nj           
cc