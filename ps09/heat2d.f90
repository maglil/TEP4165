! To solve 2D heat equation

module declarations
! Decleare the necessary variables.
	implicit none

	private
	double precision,public, parameter :: pi=3.141592653589793
	double precision,public, allocatable, dimension(:) :: x, x_face, y, y_face ! x and y coordinates of nodes and cell faces. Left/south face same index as node. (should these be 2D?)
	double precision,public :: xl, yl ! extent of simulation domain
	
	double precision,public, allocatable, dimension(:,:) :: aW,aE,aP,aN,aS
	double precision,public, allocatable, dimension(:,:) :: T,thermal	! Temperature, thermal conductivity
	double precision,public, allocatable, dimension(:,:) :: Sp,Su	
	
	double precision,public :: rho, Tamb, radius, perim, thermal_const, h, q 
	
	double precision,public :: relax(1)
	integer,public :: iter,last,npi,npj
  
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
			phi(i,j) = (aE(i,j)*phi(i-1,j) + aW(i,j)*phi(i+1,j) + aS(i,j)*phi(i,j-1) + aN(i,j)*phi(i,j+1) + b(i,j))/aP(i,j)
		end do
	end do
	
	end subroutine
end module


program heat2d
! main program

	use declarations
	use procedures
	implicit none
	
	call init()
	call grid()
	
	call bound() 
	call Tcoeff()
	

end program

! external subroutines
! subroutines *use* module declarations to access variables.

! Should have one subroutine to set system physical parameters.

! Initalize variables
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
	
	allocate(T(npi,npj),thermal(npi,npj),Su(npi,npj),Sp(npi,npj))
	
!---Set number of iterations---
	write(*,*) 'Set number of iterations:'
    read(*,*) last
	
!---Set physical quantities
	Tamb=290
    thermal_const=50.
    radius=0.01
    perim=pi*2*radius
	h=80. ! heat loss coefficient
	q = 100E3
	
!---Initalize arrays (except coefficients)
	do i=1,npi
		do j=1,npj
		  thermal(i,j)=thermal_const !thermal conductivity
		  SP(i,j)=0.
		  Su(i,j)=0.
		  T(i,j)=0.
	  end do
    end do
	
!---Set relaxation parameter
end subroutine

! Initalize grid
subroutine grid()

	use declarations
    implicit none
	
    integer :: i,j
    double precision :: dx,dy
	
	allocate(x(npi),x_face(npi),y(npj), y_face(npj))
	
!---Set physical extent
	xl = 1
	yl = 1
	
!---Size of cell
	dx=xl/real(npi-2.)
	dy=yl/real(npi-2.)
	
!---Set node coordinates
	x(1)=0.
    x(2)=0.5*dx
    do i=3,npi-1
      x(i)=x(i-1)+dx
    end do
    y(npj)=x(npj-1)+0.5*dy
	y(1)=0.
    y(2)=0.5*dy
    do i=3,npi-1
      y(i)=y(i-1)+dy
    end do
    y(npi)=y(npi-1)+0.5*dy
	
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
	
end subroutine

! Boundary conditions
subroutine bound()

	use declarations
	implicit none
	
	integer :: i,j
!---- Constant temp at each end of the bar
	do j = 1,npj
		T(1,j)=300.
		T(npi,j) = 600.
	end do
	do i = 1,npi
		T(i,1)=300.
		T(i,npj) = 600.
	end do

end subroutine

! Coefficients
subroutine Tcoeff()
	
	use declarations
    implicit none
    integer :: i
    double precision :: AREAw,AREAe,Dw,De

end subroutine
