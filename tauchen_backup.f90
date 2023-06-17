MODULE TAUCHEN

!USE imsl_libraries

!USE RNNOF_INT
USE ANORDF_INT
    
IMPLICIT NONE

Contains

Subroutine tauchen_hans(STDE,RHO,ZZ,SZ,PZEE,pi)
!=============================================================
!DISCRETIZE A NORMAL DISTRIBUTION BY USING THE TAUCHEN(1986)
!==============================================================
!	Input:  STDW(standard deviation)
!		    Log(z_t+1)=RHO*Log(z_t)+e_t+1; e~N(0,STDE^2).
!			M is the number of st. dev. within the distribution is approximated
!			ZZ is the number of shocks
!	Ouput:  SZ is the set of shock values,
!	        PZEE is the transition matrix of z
!		    pi is a stationary distribution of z

Implicit NONE

Real(8), Intent(In) :: STDE,RHO
Integer, Intent(In) :: ZZ
Real(8), Intent(Out) :: SZ(ZZ),PZEE(ZZ,ZZ),pi(ZZ)
INTEGER i,j,k
Integer, parameter :: draws_number=1000000
Real(8)  STDW, step, M, temp1,dum
Real(8), Dimension(draws_number) :: normal_draws

M=1.5d0

!==================================
!USE PZW TO SET THE VALUE OF SHOCKS
!==================================

STDW=STDE/sqrt(1-RHO*RHO)
SZ(ZZ)=M*STDW

SZ(1)=-SZ(ZZ)
step=2*STDW*M/FLOAT(ZZ-1)
DO i=2,ZZ-1
   SZ(i)=SZ(i-1)+step
Enddo

dum=(5D-1)*step
    
do i=1,ZZ
    PZEE(i,1)=D_ANORDF((SZ(1)-RHO*SZ(i)+dum)/STDE)
end do

do i=2,ZZ-1
    do j=1,ZZ
       PZEE(j,i)=D_ANORDF((SZ(i)-RHO*SZ(j)+dum)/STDE)-D_ANORDF((SZ(i)-RHO*SZ(j)-dum)/STDE)
    end do
end do

do i=1,ZZ
    PZEE(i,ZZ)=1D0-D_ANORDF((SZ(ZZ)-RHO*SZ(i)-dum)/STDE)
end do

!===========================================================
!FIND TRANSITION MATRIX BY TAUCHEN(1986):
!Use Monte Carlo Integration
!===========================================================
	
!Draw a vector normal_draws from the normal distribution.

!DO  i = 1, draws_number
!	normal_draws(i) = D_RNNOF()
!END DO
!
!
!normal_draws=normal_draws*stde

!Compute a share of draws in each interval.
!PZEE=0
!Do j=1,ZZ
!	temp1=rho*SZ(j)
!	Do i=1,draws_number
!		If (normal_draws(i)<=SZ(1)+step/2-temp1) PZEE(j,1)=PZEE(j,1)+1
!		If (normal_draws(i)> SZ(ZZ)-step/2-temp1) PZEE(j,ZZ)=PZEE(j,ZZ)+1
!	Enddo
!Enddo
!
!Do j=1,ZZ
!	temp1=rho*SZ(j)
!	Do i=1,draws_number
!		Do k=2, ZZ-1
!			If ((normal_draws(i)> SZ(k)-step/2-temp1) .and. (normal_draws(i)<=SZ(k)+step/2-temp1)) &
!				 PZEE(j,k)=PZEE(j,k)+1
!		Enddo
!	Enddo
!Enddo
!
!
!PZEE=PZEE/draws_number
!
!DO i=1,ZZ
!   SZ(i)=EXP(SZ(i))
!Enddo



!=============
!Find the stationary distribution
!=============
pi=1d0/ZZ
Do i=1,100
	pi=Matmul(pi,PZEE)
Enddo

RETURN

END Subroutine tauchen_hans

END MODULE TAUCHEN