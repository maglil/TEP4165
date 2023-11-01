! To solve 2D heat equation

module declarations
! Decleare the necessary variables.
	implicit none

	private
	double precision,public, parameter :: pi=3.141592653589793
	double precision,public, allocatable, dimension(:) :: x, x_face, y, y_face ! x and y coordinates of nodes and cell faces. Left/south face same index as node. (should these be 2D?)
	double precision,public :: dx,dy
	double precision,public :: xl, yl ! extent of simulation domain
	
	double precision,public, allocatable, dimension(:,:) :: aW,aE,aP,aN,aS
	double precision,public, allocatable, dimension(:,:) :: T,Told,diffT, thermal	! Temperature, thermal conductivity
	double precision,public, allocatable, dimension(:,:) :: Sp,Su	
	
	double precision,public :: rho, Tamb, perim, thermal_const, h, q 
	double precision,public :: width, height
	
	double precision,public :: relax(1)
	integer,public :: iter,last,npi,npj
	double precision, public :: norm0,norm1,dnorm0, n ! norms
	double precision, public :: eps ! convergence criterion
	
	integer, public :: sys !Shift between different system parameters.
	
	
	double precision, public :: norm
  
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
	
	sys = 2 ! 1 or 2 to switch between systems
	
	call init()
	call grid()
	
	allocate(aE(npi,npj),aW(npi,npj),aN(npi,npj), aS(npi,npj), aP(npi,npj))
	
	call bound() 
	
	call calc_norm()
	norm0 = norm
	write(*,*) norm0
	
	Told = T
	do iter=1,last
		
		call Tcoeff()
		call gauss_seidel(T,Su,1,npi,1,npj)
		
		! Alternating gauss seidel
		! For each j
		! b = 
		!tdma(T,b
		
		if (iter==1) then
			call calc_norm()
			norm0 = norm
			write(*,*) norm0, 'first'
		end if
		
		if(mod(iter,200)==0)then
			call calc_norm()		
			write(*,*) norm/norm0
			write (*,*) 'iteration no:',iter,'  Temperatur:',(T(int((npi)/2),int((npj)/2))+T(int((npi)/2+1),int((npj)/2+1)))/2
			if (norm/norm0 < eps) then
				write(*,*) 'exiting'
				exit
			end if
		end if
		Told = T
	end do
	
	call printout()
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
	
	allocate(T(npi,npj),Told(npi,npj), diffT(npi, npj),thermal(npi,npj),Su(npi,npj),Sp(npi,npj))
	
!---Set number of iterations---
	write(*,*) 'Set number of iterations:'
    read(*,*) last
	
!---Set physical quantities
	Tamb=290
    thermal_const=50.
    width = 0.5
	height = 0.2
    !perim=pi*2*radius
	h=80. ! heat loss coefficient
	q = 42E3
	
!---Initalize arrays (except coefficients)
	do i=1,npi
		do j=1,npj
		  thermal(i,j)=thermal_const !thermal conductivity
		  SP(i,j)=0.
		  Su(i,j)=0.
		  Told(i,j) = 0.
		  diffT(i,j) = 0.
		  T(i,j)=0. ! guess
	  end do
    end do
	
!---Set convergence criterion
	eps = 1.0e-6
	
!---Set relaxation parameter
end subroutine

! Initalize grid
subroutine grid()

	use declarations
    implicit none
	
    integer :: i,j
    
	
	allocate(x(npi),x_face(npi),y(npj), y_face(npj))
	
!---Set physical extent
	xl = 0.5
	yl = 0.2
	
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
	
end subroutine

! Boundary conditions
subroutine bound()

	use declarations
	implicit none
	
	integer :: i,j
!---- Problem a) Constant temp at faces.
	if (sys==1) then
		do j = 1,npj
			T(1,j)=300.
			T(npi,j) = 600.
		end do
		do i = 1,npi
			T(i,1)=300.
			T(i,npj) = 600.
		end do
	end if

!---Problem d) southern, western faces insulated, northern fixed temperature, easter fixed flux
	if (sys==2) then
		do i = 1,npi
			T(i,npj) = 400.
		end do
	end if

end subroutine

! Coefficients
subroutine Tcoeff()
	
	use declarations
    implicit none
    integer :: i, j
    double precision :: areaw,areae,arean,areas,Dw,De,Ds,Dn
	
	! Cell face areas
	! Assume length of slab is 1 m.
	! Move inside loop if dy,dx variable
	areaw = dy*1.
	areae = areaw
	areas = dx*1.
	arean = areas
	
	! D = k*A/dx Diffusion conductance
	
	do i = 2,npi-1
		do j = 2, npj-1
		Dw = (thermal(i-1,j) + thermal(i,j))/2*(x(i) - x(i-1)) * areaw
		De = (thermal(i+1,j) + thermal(i,j))/2*(x(i+1) - x(i)) * areae
		Ds = (thermal(i,j-1) + thermal(i,j))/2*(y(j) - y(j-1)) * areas
		Dn = (thermal(i,j+1) + thermal(i,j))/2*(y(j+1) - y(j)) * arean
		
		! Set source if present
		
		! Set coefficients
		aW(i,j) = Dw		
		aE(i,j) = De
		aS(i,j) = Ds
		aN(i,j) = Dn		
		aP(i,j) = aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j) + Sp(i,j)				
		end do
	end do
	
	if (sys==2) then
		do j = 2, npj-1
			aW(2,j) = 0. ! Western insulated
			Su(2,j) = 0.
		end do
		do i = 2, npi-1
			aS(i,2) = 0. ! Southern insulated
			Su(i,2) = 0.
		end do
		do j = 2, npj-1
			aE(npi-1,j) = 0. !Eastern fixed flux
			Su(npi-1,j) = q*dy
		end do
	end if
	
end subroutine

subroutine calc_norm()

	use declarations
	implicit none
	
	!double precision, intent(in) :: T(npi,npj)
	!double precision, intent(in) :: dx,dy
	!integer, intent(in) :: npi, npj
	integer :: i,j
	
	diffT = T-Told
	norm = 0;
	
	do i=2,npi-1
		do j=2,npi-1
			norm = norm + (diffT(i,j)**2)			
		end do
	end do	
	norm = sqrt(dx*dy*norm)
	return
	
end subroutine calc_norm
