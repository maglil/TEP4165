subroutine printout()
!
!**** Purpose: To print the numerical result to files
!
    use declarations
    implicit none
    integer :: i, j


    open(10,file='temperature.dat',status='unknown')
    do i=1,npi
    	write(10, '(*(F10.3 : ", "))') (T(i,j), j=1,npj)
    end do  
    close(10)
    
    open(11,file='x.dat',status='unknown')
    do i=1,npi
      write(11,*) x(i)
    end do  
    close(11)

    open(12,file='y.dat',status='unknown')
    do i=1,npj
      write(12,*) y(i)
    end do  
    close(12)

end subroutine printout
