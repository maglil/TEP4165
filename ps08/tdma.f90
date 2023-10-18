subroutine TDMA(fi,b,aE,aW,aP,Istart,Iend)
!
  implicit none
  double precision, dimension(:), intent(in out) :: fi
  double precision, dimension(:), intent(in) :: b
  double precision, dimension(:), intent(in) :: aE,aW,aP
  integer, intent(in) :: Istart,Iend
  integer :: I
  double precision, dimension(Iend) :: Ath, Cmth
  double precision :: Cth
!
!---- TDMA from west to east. equation to solve:
!----
!---- - aW*fiW + aP*fiP - aE*fiE = b
!----
!---- equivalences with variables in eq. 7.1-7.6:
!---- BETA=aW(I,J)   Def. in eq. 7.2
!---- D=aP(I,J)      Def. in eq. 7.2 
!---- ALFA=aE(I,J)   Def. in eq. 7.2
!---- A=Ath(I)       Def. in eq. 7.6b
!---- C=Cth          Same as the constant b=Su def. in eq 4.11
!---- C´=Cmth(I)     Def. in eq. 7.6c
!
!
!---- Solving from east to west
!
!---- At the west boundary:
  Ath(Istart)=0.   
  Cmth(Istart)=fi(Istart) 

!---- Forward substitution
  do I=Istart+1,Iend-1 
    Ath(I)=aE(I)/(aP(I)-aW(I)*Ath(I-1)) !eq. 7.6b
    Cth=b(I)   
    Cmth(I)=(aW(I)*Cmth(I-1)+Cth)/(aP(I)-aW(I)*Ath(I-1))
  end do   

!---- Back substitution  
  do I=Iend-1,Istart+1,-1 
    fi(I)=Ath(I)*fi(I+1)+Cmth(I)
  end do

end subroutine TDMA
