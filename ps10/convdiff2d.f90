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
! set physical variables

! initalize variables

! define grid

! Boundary conditions


! Define equations system coefficients

! Solve linear system
write(*,*) 'Completed program'

end program convdiff2d

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

end subroutine init