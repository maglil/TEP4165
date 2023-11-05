! Program to compute temperature distribution in a 2D system
! undergoing Pueseill flow, including diffusion and prescribed
! temperatures at boundaries.

module declarations

	private
	! Physical variables
	double precision,public :: rho, thermal_const, heat_cap, umax

	! Geometric variables
	double precision,public, allocatable, dimension(:) :: x, x_face, y, y_face ! x and y coordinates of nodes and cell faces. Left/south face same index.
	double precision,public :: dx,dy  ! Cell size
	double precision,public :: xl, yl ! Extent of simulation domain
	integer,public :: npi,npj ! Number of nodes (cells + 2)

	! System variables
	double precision,public, allocatable, dimension(:,:) :: T,Told,diffT, thermal	! Temperature, thermal conductivity
	double precision,public, allocatable, dimension(:,:) :: Sp,Su
	double precision,public, allocatable, dimension(:,:) :: aW,aE,aP,aN,aS ! Matrix coefficients

	! Algoritmic parameters
	integer,public :: iter,last
	double precision, public :: eps ! convergence criterion

end module

module procedures
! Procedures
implicit none

contains 
	subroutine gauss_seidel(phi, b, istart, iend, jstart, jend) !do we need to pass in phi, if so shouldn aX's also be passed?
	
	use declarations
	implicit none
	
	double precision,dimension(:,:), intent(in out) :: phi
	double precision,dimension(:,:), intent(in) :: b
	integer :: i, j
	integer, intent(in) :: istart,iend, jstart, jend
	
	do i=istart+1,iend-1 !i start, iend refers to node grid. computer array also includes boundary conditions which are excluded. 
		do j = jstart+1, jend-1
			! gauss_seidel
			phi(i,j) = (aE(i,j)*phi(i+1,j) + aW(i,j)*phi(i-1,j) + aS(i,j)*phi(i,j-1) + aN(i,j)*phi(i,j+1) + b(i,j))/aP(i,j)
		end do
	end do
	
	end subroutine
end module

program convdiff2d

	use declarations
	use procedures
	implicit none

!---set physical variables
	call system_params()

!---initalize variables
	call init()
	
!---define grid
	call grid()

!---Set Boundary conditions
	call bound()

!---Set linear system coefficients
	 ! Place inside loop if linearized nonlinear propblem
	
!---Solve system iteratively
	do iter=1,last		
		call tcoeff()
		call gauss_seidel(T,Su,1,npi,1,npj)
		
		if(mod(iter,200)==0)then			
			write (*,*) 'iteration no:',iter,'  Temperatur:',(T(int((npi)/2),int((npj)/2))+T(int((npi)/2+1),int((npj)/2+1)))/2			
		end if
	end do

	call printout()
write(*,*) 'Completed program'

end program convdiff2d

subroutine system_params()

	use declarations
	implicit none
	
!---Set physical quantities
    thermal_const = 0.145
	heat_cap = 2010
	rho = 800

!---System properties
	umax = 0.05 

!---Define geometry
    xl = 100E-6
	yl = 20e-6

end subroutine system_params

subroutine init()
	
	use declarations
	implicit none
	
	integer :: i,j
	
!---Ask user for number of nodes---
	write(*,*) 'Number of nodal points in the grid in x-direction:'
    read(*,*) npi
	write(*,*) 'Number of nodal points in the grid in y-direction:'
    read(*,*) npj
    write(*,*) 'Here is the numbers:',npi, ' ', npj
	
!---Allocate dynamic arrays
	allocate(T(npi,npj),Told(npi,npj), diffT(npi, npj),thermal(npi,npj))
	allocate(Sp(npi,npj),Su(npi,npj))
	allocate(aE(npi,npj),aW(npi,npj),aN(npi,npj), aS(npi,npj), aP(npi,npj))
	
!---Initalize arrays (except coefficients)
	do i=1,npi
		do j=1,npj
		  thermal(i,j)=thermal_const !thermal conductivity		  
		  T(i,j)=0. ! guess
		  Sp(i,j)=0.
		  Su(i,j)=0.
	  end do
    end do
	
!---Set number of iterations---
	write(*,*) 'Set number of iterations:'
    read(*,*) last
	
!---Set convergence criterion
	eps = 1.0e-6

end subroutine init

subroutine grid()

	use declarations
    implicit none
	
    integer :: i,j
    	
	allocate(x(npi),x_face(npi),y(npj), y_face(npj))

	!---Size of cell
	dx=xl/real(npi-2.)
	dy=yl/real(npj-2.)
	
!---Set node coordinates
	x(1)=0.
    x(2)=0.5*dx
    do i=3,npi-1
      x(i)=x(i-1)+dx
    end do
    x(npi)=x(npi-1)+0.5*dx
	
	y(1)= -yl/2 ! y=0 at midpoint of y-range
    y(2)= y(1) + 0.5*dy
    do j=3,npj-1
      y(j)=y(j-1)+dy
    end do
	
    y(npj)=y(npj-1)+0.5*dy
	
!---Set face coordinates
	x_face(1)=0.
    x_face(2)=0.
    do i=3,npi
      x_face(i)=x_face(i-1)+dx
    end do
	
	y_face(1)=0.
    y_face(2)=0.
    do i=3,npj
      y_face(i)=y_face(i-1)+dy
    end do

end subroutine grid

subroutine bound()
	
	use declarations
	implicit none
	
	integer :: i,j
		
	do j = 1,npj
		T(1,j)=900. ! left boundary
	end do
	
	do i = 1,npi
		T(i,1)=250. ! Lower bondary
		T(i,npj) = 450. ! Upper boundary
	end do

end subroutine bound
	
subroutine tcoeff()
	use declarations
    implicit none
	
    integer :: i, j
    double precision :: areaw,areae,arean,areas, Dw,De,Ds,Dn, Fe,Fw,Fs,Fn ,u , h
	
	! Cell face areas
	! Assume length of slab is 1 m.
	! Move inside loop if dy,dx variable
	areaw = dy*1.
	areae = areaw
	areas = dx*1.
	arean = areas
	
	do i = 3,npi-1
		do j = 2, npj-1
		! Diffusion coefficients k/c_p * A/dx
		Dw = (thermal(i-1,j) + thermal(i,j))/(2*(x(i) - x(i-1))) * areaw 
		De = (thermal(i+1,j) + thermal(i,j))/(2*(x(i+1) - x(i))) * areae 
		Ds = (thermal(i,j-1) + thermal(i,j))/(2*(y(j) - y(j-1))) * areas 
		Dn = (thermal(i,j+1) + thermal(i,j))/(2*(y(j+1) - y(j))) * arean 			
		
		! Convection coefficient rho u A
		h = yl/2.
		u = umax*(1-y(j)**2/h**2)
		Fe = rho * u * areae * heat_cap
		Fw = rho * u * areaw * heat_cap
		!write(*,*) Fe, De
		
		! Set coefficients
		aW(i,j) = Dw + Fw/2		
		aE(i,j) = De - Fe/2
		aS(i,j) = Ds
		aN(i,j) = Dn		
		aP(i,j) = aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)			
		end do
		!write(*,*) Dw, De, Dn, Ds, Fw, Fe
		!write(*,*) Dw/De
	end do
	
	
	! Left boundary condition, i = 2
	do j=2,npj-1
		i = 2
		! Diffusion coefficients k/c_p * A/dx
		Dw = (thermal(i-1,j) + thermal(i,j))/(2*(x(i) - x(i-1))) * areaw 
		De = (thermal(i+1,j) + thermal(i,j))/(2*(x(i+1) - x(i))) * areae 
		Ds = (thermal(i,j-1) + thermal(i,j))/(2*(y(j) - y(j-1))) * areas 
		Dn = (thermal(i,j+1) + thermal(i,j))/(2*(y(j+1) - y(j))) * arean 		
		
		! Convection coefficient rho u A
		h = yl/2.
		u = umax*(1-y(j)**2/h**2)
		
		Fe = rho * u * areae * heat_cap
		Fw = rho * u * areaw * heat_cap
		!write(*,*) Fe, De
		
		Su(i,j) = (Dw + Fw)*T(1,j)
		!aW(i,j) = Dw
		aE(i,j) = De - Fe/2
		aS(i,j) = Ds
		aN(i,j) = Dn		
		!aP(i,j) = aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)
		aP(i,j) = (Dw + Fw)+ aE(i,j)+aS(i,j)+aN(i,j)
	end do
	
	! Right boundary
	do j=2,npj-1
		i = npi
		T(i,j) = T(i-1,j)
	end do
	
	
end subroutine tcoeff