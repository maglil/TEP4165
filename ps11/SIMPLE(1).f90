module declarations
  implicit none
  private
  double precision,public, allocatable, dimension(:) :: x,x_u,y,y_v !grid coordinates?
  double precision,public :: xl,yl
  integer,public :: ipref,jpref
  double precision,public, allocatable, dimension(:,:) :: u,v,pc,p,T,rho,mu,gamma,cp !pc - pressure correction
  double precision,public, parameter :: pi=3.1415927
  double precision,public :: Tamb,radius,perim,h,Tpar ! not relevant here?
  double precision,public, allocatable, dimension(:,:) :: f_u,f_v ! convective fluxes?
  double precision,public, allocatable, dimension(:,:) :: d_u,d_v ! Pressure correction d's
  double precision,public :: m_in,m_out ! mass flow in out
  real ,public :: relax(6) ! relaxation parameters
  integer,public :: iter,last,npi,npj
end module declarations

module sub
  contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine init()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!**** purpose: to set the number of grid points and initialize all parameters.
!
! Number of grid points in each direciton:
  npi=52
  npj=22 ! must be even, see pressure reference

  allocate(u(npi,npj),v(npi,npj),p(npi,npj),pc(npi,npj),T(npi,npj),rho(npi,npj))
  allocate(mu(npi,npj),gamma(npi,npj),cp(npi,npj),d_u(npi,npj),d_v(npi,npj)) ! why not f??

!
!---- reference for zero pressure
!
  ipref=npi-1
  jpref=npj/2

!**** purpose: to initialise all parameters.

!
!---- initilising all variables 
!
  m_in=1.
  m_out=1.
!
  u = 0.05
  v = 0.00
  p = 0.
  pc=0.  !pressure correction (equivalet to p´ in ref. 1).
  T=423 !temperature, inital guess somewhere in the middle 
  d_u=0. !variable d(i,j) to calculate pc defined in 6.23
  d_v=0. !variable d(i,j) to calculate pc defined in 6.23

!
!---- set constant values
!
  rho=800  !density
  mu=2.168e-3   !dynamic viscosity 
  gamma=0.145!thermal conductivety
  cp=2010   !specific heat - assumed constant for this problem
  
!      
!---- setting the relaxation parameters
!
  relax(1)=0.8 !relax u, see eq. 6.36
  relax(2)=0.8 !relax v, see eq. 6.37
  relax(4)=0.25 !relax p, see eq. 6.33
  relax(5)=1.0 !relaxation factor for temperature
  relax(6)=0.1 !relaxation factor for density

end subroutine init

subroutine grid()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  use declarations
  implicit none
  integer :: i,j
  double precision :: dx,dy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: defining the grid 
!**** see fig. 6.2-6.4 in ref. 1 
!
!
  allocate(x(npi),x_u(npi))
  allocate(y(npj),y_v(npj))
!
!**** purpose: definding the geometrical variables 
!**** see fig. 6.2-6.4 in ref. 1. 
!
!---- length of the area in the x- and y direction 
!
  xl=100e-6
  yl=20e-6
!
!---- length of volume element
!
  dx=xl/real(npi-2) ! include also edge points dx/2 
  dy=yl/real(npj-2)
!
!
!---- length variable for the scalar points in the x direction
!
  x(1)=0.
  x(2)=0.5*dx
  do i=3,npi-1
    x(i)=x(i-1)+dx
  end do
  x(npi)=x(npi-1)+0.5*dx
!
!---- length variable for the scalar points fi(i,j) in the y direction
!
  y(1)=0.
  y(2)=0.5*dy
  do j=3,npj-1
    y(j)=y(j-1)+dy
  end do
  y(npj)=y(npj-1)+0.5*dy
!
!---- length variable for the velocity components u(i,j) in the x direction
!
  x_u(1)=0.
  x_u(2)=0.
  do i=3,npi
    x_u(i)=x_u(i-1)+dx
  end do
!
!---- length variable for the velocity components v(i,j) in the y direction
!
  y_v(1)=0.
  y_v(2)=0.
  do j=3,npj
    y_v(j)=y_v(j-1)+dy
  end do   
!
end subroutine grid
      

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fixedbound()
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: specify the boundary values that remain untouched during 
!**** the iteration process     
!
!---- fixed temperature at the upper and lower wall
!
  T(:,1)=300   !Temperature in kelvin 
  T(:,npj)=573  
!
!---- fixed temperature at the left wall and the incomming fluid
!
  T(1,:)=273
!
!---- setting the velocity at inlet    
!
  u(1,2:npj-1)=0.05
  u(2,2:npj-1)=0.05
  v(1,:)=0.0
!
!---- setting the velocity at the upper and lower wall
!
  v(:,1)=0. 
  v(:,2)=0. 
  v(:,npj)=0.
  u(:,1)=0. 
  u(:,npj)=0.

end subroutine fixedbound


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine bound()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  integer :: j
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: specify boundary conditions for a calculation
!
  call globcont()
!
!---- velocity and temperature gradient at outlet = zero:
!
  do j=2,npj-1
    u(npi,j)=u(npi-1,j)*m_in/m_out
    v(npi,j)=v(npi-1,j)
    T(npi,j)=T(npi-1,j)
  end do

end subroutine bound


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine globcont()
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  integer :: j
  real :: areaw  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: calculate mass in and out of the calculation domain to
!**** correct for the continuety at outlet.
!
  call conv()
!      
  m_in=0.
  m_out=0.
!
  do j=2,npj-1
    areaw=y_v(j+1)-y_v(j) !see fig. 6.3  
    m_in=m_in+f_u(2,j)*areaw
    m_out=m_out+f_u(npi-1,j)*areaw        
  end do
!
end subroutine globcont


subroutine output()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: creating result table
!
  if(iter.eq.last) call print()
!
end subroutine output


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine velcorr()
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  integer :: i,j
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!****  purpose: to correct the pressure and the velocities by eq. 6.24, 6.25
!****  and a modified version of eq. 6.33. 
!
  do i=2,npi-1
    do j=2,npj-1
!
!---- eq. 6.33 modified for dynamical pressure:
!---- this formulation makes the pressure maintain the zero value 
!---- at i=ipref and j=jpref, so the pressure anywhere else is calculated
!---- in reference to the pressure here.
!
      p(i,j)=p(i,j)+relax(4)*(pc(i,j)-pc(ipref,jpref)) 
!
!---- velocity correction
!
      if(i.ne.2) u(i,j)=u(i,j)+d_u(i,j)*(pc(i-1,j)-pc(i,j)) !eq. 6.24
      if(j.ne.2) v(i,j)=v(i,j)+d_v(i,j)* (pc(i,j-1)-pc(i,j)) !eq. 6.25
    end do
  end do   
  p(1,:) = p(2,:)  
  p(:,1) = p(:,2) 
  p(npi,:) = p(npi-1,:)
  p(:,npj) = p(:,npj-1)
!      
end subroutine velcorr


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine print()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
implicit none
integer :: i,j
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c**** purpose: print out result table to file
!c
      open(100,file='u.dat',status='unknown')
      open(101,file='v.dat',status='unknown')

      open(110,file='x.dat',status='unknown')
      open(111,file='y.dat',status='unknown')

      open(12,file='t.dat',status='unknown')
      open(13,file='p.dat',status='unknown')

      do i = 1,npi
      	write(110,*) x(i)
      end do

      do j = 1,npj
        write(111,*) y(j)
      end do 
 

      write(100, '(*(F10.3 : ", "))') (u(2,j), j=1,npj)
      do i=2,npi-1!
	 write(100, '(*(F10.3 : ", "))') ((u(i,j)+u(i+1,j))/2, j=1,npj)
      end do
      write(100, '(*(F10.3 : ", "))') (u(npi,2), j=1,npj)


      write(101, '(*(F10.3 : ", "))') (v(i,2), i=1,npi)
      do j=2,npj-1
	 write(101, '(*(F10.3 : ", "))') ((v(i,j)+v(i,j+1))/2, i=1,npi)
      end do
      write(101, '(*(F10.3 : ", "))') (v(i,npj), i=1,npi)



      do j = 1, npj
	 write(12, '(*(F10.3 : ", "))') (t(i,j), i=1,npi)
      end do

      do j=1,npj
	write(13, '(*(F10.3 : ", "))') (p(i,j), i=1,npi)
      end do

      close(100)
      close(101)
      close(110)
      close(111)
      close(12)
      close(13)

!
 100  format(' ',3e15.6e3)
 101  format(' ',4e15.6e3)
!
end subroutine print

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine solve(fi,b,ae,aw,an,as,ap,istart,iend,jstart,jend)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  double precision, dimension(:,:), intent(in out) :: fi
  double precision, dimension(:,:), intent(in) :: b
  double precision, dimension(:,:), intent(in) :: ae,aw,an,as,ap
  integer, intent(in) :: istart,iend,jstart,jend
  integer :: i,j
  double precision, allocatable, dimension(:) :: ath, cmth
  double precision :: cth
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: to solve the algebraic equation 7.7.
!

  allocate(ath(iend),cmth(iend))


!---- tdma along a horisontal row from west to east. equation to solve:
!----
!---- - aw*fiw + ap*fip - ae*fie = as*fis + an*fin + b
!----
!---- equivalences with variables in eq. 7.1-7.6:
!---- beta=aw(i,j)   def. in eq. 7.2
!---- d=ap(i,j)      def. in eq. 7.2 
!---- alfa=ae(i,j)   def. in eq. 7.2
!---- a=ath(i)        def. in eq. 7.6b
!---- c=cth           the right side assumed temporarily known (see eq. 7.8)
!---- c´=cmth(i) def. in eq. 7.6c
!---- b=b(i,j)       def. in eq. 7.7
!
!
!---- solving the (e-w) lines from the south
!
  do j=jstart+1,jend-1 
!---- at the inlet boundary:
    ath(istart)=0.   
    cmth(istart)=fi(istart,j) 
!
    do i=istart+1,iend-1 !forward substitution
      ath(i)=ae(i,j)/(ap(i,j)-aw(i,j)*ath(i-1)) !eq. 7.6b
      cth=an(i,j)*fi(i,j+1)+as(i,j)*fi(i,j-1)+b(i,j)   
      cmth(i)=(aw(i,j)*cmth(i-1)+cth)/(ap(i,j)-aw(i,j)*ath(i-1)) !eq. 7.6c
    end do   
!
    do i=iend-1,istart+1,-1 !back substitution  
      fi(i,j)=ath(i)*fi(i+1,j)+cmth(i) !eq. 7.6a
    end do
  end do
!
!  
!---- solving the (e-w) lines from the north
!
  do j=jend-2,jstart+1,-1 
!---- at the inlet boundary:
    ath(istart)=0. 
    cmth(istart)=fi(istart,j)
!
    do i=istart+1,iend-1 !forward substitution
      ath(i)=ae(i,j)/(ap(i,j)-aw(i,j)*ath(i-1)) !eq. 7.6b
      cth=an(i,j)*fi(i,j+1)+as(i,j)*fi(i,j-1)+b(i,j)  
      cmth(i)=(aw(i,j)*cmth(i-1)+cth)/(ap(i,j)-aw(i,j)*ath(i-1))  !eq. 7.6c
    end do   
!
    do i=iend-1,istart+1,-1 !back substitution  
      fi(i,j)=ath(i)*fi(i+1,j)+cmth(i) !eq. 7.6a
    end do
  end do      
!

!---- tdma along a vertical column from south to north. equation to solve:
!----
!---- - as*fiw + ap*fip - an*fie = aw*fis + ae*fin + b !eq. 7.8
!----
!---- equivalences with variables in eq. 7.1-7.6:
!---- beta=as(i,j)   def. in eq. 7.2
!---- d=ap(i,j)      def. in eq. 7.2 
!---- alfa=an(i,j)   def. in eq. 7.2
!---- a=ath(i)        def. in eq. 7.6b
!---- c=cth           the right side assumed temporarily known (see eq. 7.8)
!---- c´=cmth(i) def. in eq. 7.6c
!---- b=b(i,j)       def. in eq. 7.7
!
!
  deallocate(ath,cmth)
!
!
!---- solving (n-s) lines from the west
!
  allocate(ath(jend),cmth(jend))
!
  do i=istart+1,iend-1
!---- at the bottom boundary:
    ath(jstart)=0.   
    cmth(jstart)=fi(i,jstart) 
!
    do j=jstart+1,jend-1
      ath(j)=an(i,j)/(ap(i,j)-as(i,j)*ath(j-1)) !eq. 7.6b
      cth=ae(i,j)*fi(i+1,j)+aw(i,j)*fi(i-1,j)+b(i,j)  
      cmth(j)=(as(i,j)*cmth(j-1)+cth)/(ap(i,j)-as(i,j)*ath(j-1)) !eq. 7.6c
    end do   
!
    do j=jend-1,jstart+1,-1 !back substitution
      fi(i,j)=ath(j)*fi(i,j+1)+cmth(j) !eq. 7.6a 
    end do   
  end do
!
!
!---- solving (n-s) lines from the east

  do i=iend-2,istart+1,-1
!---- at the bottom boundary:
    ath(jstart)=0. 
    cmth(jstart)=fi(i,jstart)
!         
    do j=jstart+1,jend-1 !foreward substitution
      ath(j)=an(i,j)/(ap(i,j)-as(i,j)*ath(j-1)) !eq. 7.6b
      cth=ae(i,j)*fi(i+1,j)+aw(i,j)*fi(i-1,j)+b(i,j)  
      cmth(j)=(as(i,j)*cmth(j-1)+cth)/(ap(i,j)-as(i,j)*ath(j-1)) !eq. 7.6c
    end do   
!
    do j=jend-1,jstart+1,-1 !back substitution
      fi(i,j)=ath(j)*fi(i,j+1)+cmth(j) !eq. 7.6a 
    end do   
  end do   
!
  deallocate(ath,cmth)
!
end subroutine solve


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ucoeff(ae,aw,an,as,ap,b)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
  integer :: i,j
  real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: to calculate the coefficients for the u-equation
!
!
  call conv()
!
  do i=3,npi-1
    do j=2,npj-1
!
!---- geometrical parameters
!---- areas of the cell faces
!
      areaw=y_v(j+1)-y_v(j) !see fig. 6.3  
      areae=areaw
      areas=x(i)-x(i-1)
      arean=areas
!
!---- eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a 
!---- note:  f=rho*u but fw=(rho*u)w=rho*u*areaw per definition.
!
      fw=((f_u(i,j)+f_u(i-1,j))/2)*areaw
      fe=((f_u(i+1,j)+f_u(i,j))/2)*areae
      fs=((f_v(i,j)+f_v(i-1,j))/2)*areas
      fn=((f_v(i,j+1)+f_v(i-1,j+1))/2)*arean
!
!---- eq. 6.9e-6.9h - the diffusion conductance defined in eq. 5.8b 
!---  note: d=mu/dx but dw=(mu/dx)*areaw per definition
!
      dw=(mu(i-1,j)/(x_u(i)-x_u(i-1)))*areaw
      de=(mu(i,j)/(x_u(i+1)-x_u(i)))*areae
      ds=((mu(i-1,j)+mu(i,j)+mu(i-1,j-1)+mu(i,j-1))/(4*(y(j)-y(j-1))))*areas
      dn=((mu(i-1,j+1)+mu(i,j+1)+mu(i-1,j)+mu(i,j))/(4*(y(j+1)-y(j))))*arean    
!
!---- the source terms
!
      sp=0.   
      su=0.
!
!---- the coefficients (upwind differencing sheme)----------------------
!
      aw(i,j)=dw+max(fw,0.)
      ae(i,j)=de+max(0.,-fe)
      as(i,j)=ds+max(fs,0.)
      an(i,j)=dn+max(0.,-fn)
!
!---- eq. 8.31 without time dependent terms (see also eq. 5.14):
! 
      ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp 
!
!-------------------------------------------------------------------------
!
!---- calculation of d(i,j)=d_u(i,j) defined in eq. 6.23 for use in the 
!---- equation for pression correction (eq. 6.32). see subroutine pccoeff.
!                      
      d_u(i,j)=areaw*relax(1)/ap(i,j)
!
!---- putting the integrated pressure gradient into the source term b(i,j)
!---- the reason is to get an equation on the generelised form 
!---- (eq. 7.7 ) to be solved by the thomas algorithm. 
!---- note: in reality b=a0p*fip+su=0. 
!
      b(i,j)=(p(i-1,j)-p(i,j))*areaw+su
!            
!---- introducing relaxation by eq. 6.36 . and putting also the last 
!---- term on the right side into the source term b(i,j)
!            
      ap(i,j)=ap(i,j)/relax(1)
      b(i,j)=b(i,j)+(1-relax(1))*ap(i,j)*u(i,j)
!
!---- now we have implemented eq. 6.36  on the form of eq. 7.7
!---- and the thomas algorithm can be called to solve it. this is done 
!---- in the next step of the main program.            
!
    end do
  end do   
!
end subroutine ucoeff


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine vcoeff(ae,aw,an,as,ap,b)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
  integer :: i,j
  real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: to calculate the coefficients for the v-equation
!
  call conv()
!
  do i=2,npi-1
    do j=3,npj-1
!
!---- geometrical parameters
!---- areas of the cell faces
!
      areaw=y(j)-y(j-1) !see fig. 6.4
      areae=areaw
      areas=x_u(i+1)-x_u(i)
      arean=areas
!
!---- eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a 
!---- note:  f=rho*u but fw=(rho*u)w=rho*u*areaw per definition.
!
      fw=((f_u(i,j)+f_u(i,j-1))/2)*areaw
      fe=((f_u(i+1,j)+f_u(i+1,j-1))/2)*areae
      fs=((f_v(i,j)+f_v(i,j-1))/2)*areas
      fn=((f_v(i,j)+f_v(i,j+1))/2)*arean
!
!---- eq. 6.11e-6.11h - the diffusion conductance defined in eq. 5.8b 
!---  note: d=mu/dx but dw=(mu/dx)*areaw per definition
!
      dw=((mu(i-1,j-1)+mu(i,j-1)+mu(i-1,j)+mu(i,j))/(4*(x(i)-x(i-1))))*areaw
      de=((mu(i,j-1)+mu(i+1,j-1)+mu(i,j)+mu(i+1,j))/(4*(x(i+1)-x(i))))*areae    
      ds=(mu(i,j-1)/(y_v(j)-y_v(j-1)))*areas
      dn=(mu(i,j)/(y_v(j+1)-y_v(j)))*arean
!
!---- the source terms
!
      sp=0.   
      su=0.
!
!---- the coefficients (upwind differencing sheme)----------------------
!
      aw(i,j)=dw+max(fw,0.)
      ae(i,j)=de+max(0.,-fe)
      as(i,j)=ds+max(fs,0.)
      an(i,j)=dn+max(0.,-fn)
!
!---- eq. 8.31 without time dependent terms (see also eq. 5.14):
! 
      ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp
!
!-------------------------------------------------------------------------
!
!---- calculation of d(i,j)=d_v(i,j) defined in eq. 6.23 for use in the 
!---- equation for pression correction (eq. 6.32) (see subroutine pccoeff).
!                      
      d_v(i,j)=areas*relax(2)/ap(i,j)
!
!---- putting the integrated pressure gradient into the source term b(i,j)
!---- the reason is to get an equation on the generelised form 
!---- (eq. 7.7 ) to be solved by the thomas algorithm. 
!---- note: in reality b=a0p*fip+su=0. 
!
      b(i,j)=(p(i,j-1)-p(i,j))*areas+su
!            
!---- introducing relaxation by eq. 6.37 . and putting also the last 
!---- term on the right side into the source term b(i,j)
!            
      ap(i,j)=ap(i,j)/relax(2)
      b(i,j)=b(i,j)+(1-relax(2))*ap(i,j)*v(i,j)
!            
!---- now we have implemented eq. 6.36  on the form of eq. 7.7
!---- and the thomas algorithm can be called to solve it. this is done 
!---- in the next step of the main program.            
!
    end do
  end do   
!
end subroutine vcoeff


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine tcoeff(ae,aw,an,as,ap,b)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
  integer :: i,j
  real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: to calculate the coefficients for the t equation.
!
!
  call conv()
!
  do i=2,npi-1
    do j=2,npj-1
!---- geometrical parameters
!---- areas of the cell faces
!
      areaw=y_v(j+1)-y_v(j) !=a(i,j) see fig. 6.2 or fig. 6.5
      areae=areaw
      areas=x_u(i+1)-x_u(i)!=a(i,j)
      arean=areas
!
!---- the convective mass flux defined in eq. 5.8a 
!---- note:  f=rho*u but fw=(rho*u)w=rho*u*areaw per definition.
!
      fw=f_u(i,j)*cp(i,j)*areaw
      fe=f_u(i+1,j)*cp(i,j)*areae
      fs=f_v(i,j)*cp(i,j)*areas
      fn=f_v(i,j+1)*cp(i,j)*arean
!
!---- the diffusion conductance defined in eq. 5.8b 
!---  note: d=mu/dx but dw=(mu/dx)*areaw per definition
!
!---- the conductivity, gamma,
!---- at the interface is calculated by use of harmonic mean.  
!
      dw=((gamma(i-1,j)*gamma(i,j))/(gamma(i-1,j)*    &
            (x(i)-x_u(i))+gamma(i,j)*(x_u(i)-x(i-1))))*areaw
      de=((gamma(i,j)*gamma(i+1,j))/(gamma(i,j)*       &
            (x(i+1)-x_u(i+1))+gamma(i+1,j)*(x_u(i+1)-x(i))))*areae
      ds=((gamma(i,j-1)*gamma(i,j))/(gamma(i,j-1)*      &
            (y(j)-y_v(j))+gamma(i,j)*(y_v(j)-y(j-1))))*areas
      dn=((gamma(i,j)*gamma(i,j+1))/(gamma(i,j)*         &
            (y(j+1)-y_v(j+1))+gamma(i,j+1)*(y_v(j+1)-y(j))))*arean
!
!---- the source terms
!
      sp=0.
      su=0.
!
!
!---- the coefficients (upwind differencing sheme)----------------------
!
      aw(i,j)=dw+max(fw,0.)
      ae(i,j)=de+max(0.,-fe)
      as(i,j)=ds+max(fs,0.)
      an(i,j)=dn+max(0.,-fn)
!
!---- eq. 8.31 without time dependent terms (see also eq. 5.14):
! 
      ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp 
!
!-------------------------------------------------------------------------
!
!---- setting the source term equal to b
!
      b(i,j)=su
!
!---- now the thomas algorithm can be called to solve the equation. 
!---- this is done in the next step of the main program.            
!
    end do
  end do   
!
end subroutine tcoeff


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine conv()
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
use declarations
  implicit none
  integer :: i,j
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: to calculate the convective mass flux component pr. unit
!**** area defined in eq. 5.7
!
  do i=2,npi
    do j=2,npj
      f_u(i,j)= (rho(i-1,j)*(x(i)-x_u(i))+           &
            rho(i,j)*(x_u(i)-x(i-1)))*u(i,j)/(x(i)-x(i-1)) ! =f(i,j)
      f_v(i,j)= (rho(i,j-1)*(y(j)-y_v(j))+           &
            rho(i,j)*(y_v(j)-y(j-1)))*v(i,j)/(y(j)-y(j-1)) ! =f(i,j)
    end do
 end do   
!      
end subroutine conv


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine pccoeff(ae,aw,an,as,ap,b)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
  implicit none
  double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
  integer :: i,j
  real ::  areaw,areae,areas,arean,sp
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**** purpose: to calculate the coefficients for the pressure correction 
!**** equation, eq. 6.32
!
  call conv()
!
  do i=2,npi-1
    do j=2,npj-1
!
!---- geometrical parameters
!---- areas of the cell faces
!
      areaw=y_v(j+1)-y_v(j) !=a(i,j) see fig. 6.2 or fig. 6.5
      areae=areaw
      areas=x_u(i+1)-x_u(i)!=a(i,j)
      arean=areas
!
!---- the constant b´ in eq. 6.32
!
      b(i,j)=f_u(i,j)*areaw-f_u(i+1,j)*areae+   &
                f_v(i,j)*areas-f_v(i,j+1)*arean
      sp=0.   
!
!---- the coefficients ---------------------------------------------------
!               
      ae(i,j)=(rho(i,j)+rho(i+1,j))*d_u(i+1,j)*areae/2.
      aw(i,j)=(rho(i-1,j)+rho(i,j))*d_u(i,j)*areaw/2.
      an(i,j)=(rho(i,j)+rho(i,j+1))*d_v(i,j+1)*arean/2.
      as(i,j)=(rho(i,j-1)+rho(i,j))*d_v(i,j)*areas/2.
!
      ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)-sp
!
      pc(i,j)=0. !important!! old pressure corrections must be set
                 !to zero to avoid influence on the new ones

!-------------------------------------------------------------------------
!---- note: at the points nearest the boundaries, some coefficients are
!---- necesserily zero. for instance at i=2 and j=2, the coefficients
!---- as and aw are zero since they are on the outside of the calculation
!---- domain. this is automatically satisfied through the initialisation
!---- where d_u(i,j) and d_v(i,j) are set to zero at these points.   
!
    end do
  end do
!
end subroutine pccoeff

end module sub


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 program SIMPLE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!**** solves: non transient, compressible  convection-diffusion problems.  
!**** version: central differences.
!**** documentation: educational (overdocumented). all equations cited are 
!**** from ref. 1
!**** description:
!**** this program solves non transient convection-diffusion problems      
!**** using the SIMPLE algorithm described in ch. 6.4 in "computational 
!**** fluid dynamics" by h.k. versteeg and w. malalasekera. symbols and 
!**** variables follows exactly the notations in this reference, and all 
!**** equations cited are from this reference unless somthing else is said.
!
!**** references: 1. computational fluid dynamics, h.k. versteeg and w. 
!****                malalasekera, longman group ltd, 1995
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use declarations
use sub
implicit none
double precision, allocatable, dimension(:,:)   :: ae,aw,an,as,ap,b
integer :: k,it,jt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  call init()
  call grid()    
  call fixedbound()

!--- allocate remaining arrays: 
  allocate(ae(npi,npj),aw(npi,npj),an(npi,npj),as(npi,npj),ap(npi,npj),b(npi,npj))
  allocate(f_u(npi,npj),f_v(npi,npj))

  it=npi/2
  jt=npj/2
  write(*,*) 'Node for convergence history:',it,jt

! Set the maximum number of iterations
  last=1000
!
! SIMPLE-loop
  do iter=1,last
! 
! ENTER YOUR CODE HERE
! 
! Print convergence history
    if (mod(iter,10) == 0) then
      write (*,'(i6,5g15.5)')  iter, &
              & u(it,jt),v(it,jt),p(it,jt),pc(it,jt),T(it,jt) 
    end if
  end do
!       
end program SIMPLE
