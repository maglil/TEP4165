program heat2d
  implicit none
  integer, parameter :: nx = 100 ! number of grid points in x direction
  integer, parameter :: ny = 100 ! number of grid points in y direction
  integer, parameter :: nt = 1000 ! number of time steps
  real, parameter :: lx = 2.0 ! length of the domain in x direction
  real, parameter :: ly = 2.0 ! length of the domain in y direction
  real, parameter :: alpha = 0.1 ! thermal diffusivity constant
  real, parameter :: u0 = 0.0 ! initial and boundary temperature
  real, parameter :: u1 = 100.0 ! temperature at the top boundary
  real :: dx, dy, dt ! grid spacing and time step size
  real :: u(nx,ny) ! temperature field
  integer :: i, j, m ! loop indices
  character(len=20) :: filename ! file name for output
  open(unit=10, file='heat2d.dat', status='replace', action='write') ! open a file for output
  dx = lx / (nx - 1) ! compute the grid spacing in x direction
  dy = ly / (ny - 1) ! compute the grid spacing in y direction
  dt = 0.25 * (dx * dy)**2 / (alpha * (dx**2 + dy**2)) ! compute the time step size using the stability condition
  u = u0 ! initialize the temperature field to u0
  u(:,ny) = u1 ! set the temperature at the top boundary to u1
  do m = 1, nt ! loop over the time steps
    do j = 2, ny - 1 ! loop over the y direction
      do i = 2, nx - 1 ! loop over the x direction
        ! update the temperature using the finite difference scheme
        u(i,j) = u(i,j) + dt * alpha * ((u(i-1,j) - 2.0 * u(i,j) + u(i+1,j)) / dx**2 &
               + (u(i,j-1) - 2.0 * u(i,j) + u(i,j+1)) / dy**2)
      end do
    end do
    write(filename, '(A,I4.4,A)') 'heat2d_', m, '.txt' ! generate a file name for output
    open(unit=11, file=filename, status='replace', action='write') ! open a file for output
    do j = 1, ny ! loop over the y direction
      write(11, *) u(:,j) ! write the temperature field to the file
    end do
    close(unit=11) ! close the file
  end do
  close(unit=10) ! close the file
end program heat2d