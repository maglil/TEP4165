!**** FORTRAN90 file skeleton.f90
!**** Consists two modules (declarations and subroutines) and a main program (simple).
!
!**** Results written to file 'res.dat'. Can be imported to, e.g., MATLAB for plotting.
!
!	Overview of the Skeleton file:
!		1. A module 'declarations': collects a number of useful variables
!		2. A module 'LES_solver': here you put the TDMA algorithm
!		   as a solver for linear equation systems (LES)
!		3. Main program 'simple'
!		4. Subroutines needed by 'simple'
!**************************************************************************************


module declarations
  implicit none
  private
  double precision,public, allocatable, dimension(:) :: x,x_face
  double precision,public :: xl
  double precision,public, allocatable, dimension(:) :: aW,aE,aP
  double precision,public, allocatable, dimension(:) :: T,thermal
  double precision,public, parameter :: pi=3.141592653589793
  double precision,public :: rho,Tamb,radius,perim,thermal_const,h
  double precision,public, allocatable, dimension(:) :: Sp,Su
  double precision,public :: relax(1)
  integer,public :: iter,last,npi
end module declarations

module procedure
	implicit none
	
	contains
		subroutine solve(phi,b,istart,iend)
!
!**** Purpose: To solve the algebraic equation 7.7 in ref. 1.
!
   		 use declarations
   		 implicit none
   		 double precision,dimension(:), intent(in out) :: phi
    	double precision,dimension(:), intent(in) :: b
   		 integer :: i
   		 integer, intent(in) :: istart,iend

!---- Solving from left to right
   		 do i=istart+1,iend-1
    	  phi(i)=(aE(i)*phi(i+1)+aW(i)*phi(i-1)+b(i))/aP(i)
   		 end do
!---- Solving from right to left
  
		end subroutine solve


end module procedure





!module LES_solver
!
!**** This is a module where you put the TDMA routine (a solver for linear equation systems (LES)).
!****
!
!  contains

!
!
!end module LES_solver






program simple
!
!**** Solves: Steady state heat conduction.
!****
!**** References: 
!**** [1] H.K. Versteeg & W. Malalasekera: 
!****     An Introduction to Computational Fluid Dynamics: The Finite Volume Method, 
!****     2nd edition, Pearson / Prentice Hall, 2007.
!**** [2] S.V. Patankar: Numerical Heat Transfer and Fluid Flow, McGraw-Hill, 1980.
!
  use declarations
  use procedure
  !use LES_solver
  implicit none
  call init()
  call grid()

  !---Allocating remaining arrays:
  allocate(aE(npi),aW(npi),aP(npi))

  do iter=1,last
    call bound()
    call Tcoeff()
    call solve(T,Su,1,npi)
	!call tdma(T,0, aW, aE, aP, 1, npi)

    if(mod(iter,200)==0)then
      write (*,*) 'iteration no:',iter,'  Temperatur:',T(npi/2)
    end if

    call output()
  end do
       
end program simple




subroutine init()
!
!**** Purpose: To initialize all parameters.
!
    use declarations
    implicit none
    integer :: i

!---Number of nodal points in the grid:
    write(*,*) 'Number of nodal points in the grid:'
    read(*,*) npi
    write(*,*) 'Here is the number:',npi


  allocate(T(npi),thermal(npi),Su(npi),Sp(npi))

!---Number of iterations
    write(*,*) 'Set number of iterations:'
    read(*,*) last

!---Setting the constant thermal conductivity, ambient temperature and the radius and 
!---perimeter of the bar (SI units):
    !Tamb=
    thermal_const=50.
    radius=0.01
    perim=pi*2*radius

!---- Initilising all other variables 
    do i=1,npi
      thermal(i)=thermal_const !thermal conductivity
      SP(i)=0.
      Su(i)=0.
      T(i)=0.
    end do

!---Relaxation parameters
    !relax(1)=1.0
end subroutine init



subroutine grid()
!
!**** Purpose: Defining the grid 
!**** See fig. 6.2-6.4 in ref. 1 
!
    use declarations
    implicit none
    integer :: i
    double precision :: dx

    allocate(x(npi),x_face(npi))

!---Length of the calculation domain in the x-direction
    xl=1.
!----Length of volume element
    dx=xl/real(npi-2.)
!----Length variable for the scalar points in the x direction
    x(1)=0.
    x(2)=0.5*Dx
    do i=3,npi-1
      x(i)=x(i-1)+dx
    end do
    x(npi)=x(npi-1)+0.5*Dx

!---Length variable for the interfaces between the scalar points in the x direction
    x_face(1)=0.
    x_face(2)=0.
    do i=3,npi
      x_face(i)=x_face(i-1)+dx
    end do

end subroutine grid



subroutine bound()
!
!**** Purpose: Specify boundary conditions for a calculation
!
  use declarations
  implicit none

!---- Constant temp at each end of the bar
    T(1)=350.
    T(npi)=1100.

end subroutine bound



subroutine Tcoeff()
!
!**** Purpose: To calculate the coefficients for the T equation.
!
    use declarations
    implicit none
    integer :: i
    double precision :: AREAw,AREAe,Dw,De

    do i=2,npi-1

!---- Geometrical parameters
!---- Areas of the cell faces
      AREAw=pi*radius**2
      AREAe=AREAw

!---- The diffusion conductance D=(thermal/Dx)*AREA defined eq. 5.8b in ref. 1
!---- for one dimentional cases. (Note: AREA=1 for one dimentional cases). 
      Dw=((thermal(i-1)+thermal(i))/(2*(x(i)-x(i-1))))*AREAw
      De=((thermal(i)+thermal(i+1))/(2*(x(i+1)-x(i))))*AREAe    

!---- The source terms
      !SP(i)=-h*perim*(x_face(i+1)-x_face(i))
      !Su(i)=h*perim*(x_face(i+1)-x_face(i))*Tamb

!---- The coefficients (central differencing sheme)----------------------
      aW(i)=Dw
      aE(i)=De
      aP(i)=aW(i)+aE(i)-SP(i)

!---- Introducing relaxation by eq. 4.55b in ref. 2 and putting also the 
!---- last term on the right side into the source term b(i,J)
      !aP(i)=aP(i)/relax(1)
      !b(i)=b(i)+(1-relax(1))*aP(i)*T(i)

    end do
end subroutine Tcoeff







subroutine printout()
!
!**** Purpose: To print the numerical result to files
!
    use declarations
    implicit none
    integer :: i

    open(10,file='res.dat',status='unknown')

    do i=1,npi
      write(10,100) x(i),T(i)
    end do  
    100  format(' ',2E17.8)

    close(10)

end subroutine printout


subroutine output()
!
!**** Purpose: To call the subroutine "printout" at the end of the 
!***           calculation
!
    use declarations
    implicit none

    if(iter.eq.last) call printout()

end subroutine output




