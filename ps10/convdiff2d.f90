! Program to compute temperature distribution in a 2D system
! undergoing Pueseill flow, including diffusion and prescribed
! temperatures at boundaries.

module declarations

	private
	! Physical variables
	double precision,public :: rho, thermal_const, heat_cap

	! Geometric variables
	double precision,public, allocatable, dimension(:) :: x, x_face, y, y_face ! x and y coordinates of nodes and cell faces. Left/south face same index.
	double precision,public :: dx,dy  ! Cell size
	double precision,public :: xl, yl ! Extent of simulation domain
	integer,public :: npi,npj ! Number of nodes (cells + 2)

	! Solution variables
	double precision,public, allocatable, dimension(:,:) :: T,Told,diffT, thermal	! Temperature, thermal conductivity
	double precision,public, allocatable, dimension(:,:) :: aW,aE,aP,aN,aS ! Matrix coefficients

	! Algoritmic parameters
	integer,public :: iter,last
	double precision, public :: eps ! convergence criterion

end module

program convdiff2d

	use declarations
	implicit none

!---set physical variables
	call system_params()

!---initalize variables
	call init()
	
!---define grid
	call grid()

! Boundary conditions
	call bound()


! Define equations system coefficients

! Solve linear system
write(*,*) 'Completed program'

end program convdiff2d

subroutine system_params()

	use declarations
	implicit none
	
!---Set physical quantities
    thermal_const=50.
	
!---Define geometry
    xl = 0.5
	yl = 0.2

end subroutine system_params

subroutine init()
	
	use declarations
	implicit none
	
	integer :: i,j
	
!---Set number of nodal points and allocate arrays---
	write(*,*) 'Number of nodal points in the grid in x-direction:'
    read(*,*) npi
	write(*,*) 'Number of nodal points in the grid in y-direction:'
    read(*,*) npj
    write(*,*) 'Here is the numbers:',npi, ' ', npj
	
!---Allocate dynamic arrays
	allocate(T(npi,npj),Told(npi,npj), diffT(npi, npj),thermal(npi,npj))
	allocate(aE(npi,npj),aW(npi,npj),aN(npi,npj), aS(npi,npj), aP(npi,npj))
	
!---Initalize arrays (except coefficients)
	do i=1,npi
		do j=1,npj
		  !thermal(i,j)=thermal_const !thermal conductivity		  
		  T(i,j)=0. ! guess
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
    x(npi)=x(npi-1)+0.5*dy
	
	y(1)=0.
    y(2)=0.5*dy
    do i=3,npi-1
      y(i)=y(i-1)+dy
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
		T(i,npj) = 250. ! Upper boundary
	end do

end subroutine bound