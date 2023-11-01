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

subroutine init

end subroutine init