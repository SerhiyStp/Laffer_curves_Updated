subroutine Statistics_to_file(file_id)

!This subroutine computes aggregate statistics from the simulation

use Model_Parameters
use PolicyFunctions
use Utilities
USE RLSE_INT
use CORVC_int
implicit none
integer, intent(in) :: file_id
integer :: i,country,ia,ia2,iu,ix,um,it2,ik,ifc, NVAR=2,it3,it4,ICOPT=2,ik2
real(8) :: dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,dum12,dum13,dum14,dum15,dum16,dum17,SST,SSE,COV(2,2)
real(8) :: dum18,dum19,dum20,dum21,dum22,dum23,dum24,dum25,r_ret,ss_tax,ss_expense
real(8), dimension (:,:), allocatable :: XVARS
real(8), dimension (:), allocatable :: YVAR, BREG
real(8), allocatable :: spousewage(:,:), spousewage2(:,:)

allocate(XVARS(nsim2*nsim*T,1))
allocate(YVAR(nsim2*nsim*T))
allocate(BREG(2))

!write(file_id, *)'Gamma_redistr',Gamma_redistr/2d0

!Computing the weight of each generation

WeightActive(1)=1d0
do i=2,T
    WeightActive(i)=WeightActive(i-1)*OmegaActive(i-1)
end do

WeightRet(1)=WeightActive(T)*OmegaActive(T)
do i=2,Tret
    WeightRet(i)=WeightRet(i-1)*OmegaRet(i-1)
end do

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

!Labor Supply before 65

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
    dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Labor supply below 65 is',dum2



dum2=0d0
dum3=0d0

!Male Labor Supply

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
    dum3=dum3+1d0*WeightActive(i)
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Male labor supply is',dum2

!Single Male Labor Supply

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)<0.5) then
        dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Single Male labor supply is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.260d0)/0.260d0)**2d0
dum10=((dum2-0.260d0)/0.260d0)**2d0

!Single male labor force participation

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)<0.5) then
        if(Sim1m(it2,it,i,4)>1d-3) then
            dum2=dum2+(1d0)*WeightActive(i)
        end if
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

dum10=dum10+((dum2-0.727d0)/0.727d0)**2d0

write(file_id, *)'Single male labor force participation is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.727d0)/0.727d0)**2d0

YVAR=sqrt(-1.0)
XVARS=sqrt(-1.0)

do i=2,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)<0.5) then
        if(Sim1m(it2,it,i-1,4)>1d-3) then
            XVARs((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        else
            XVARs((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        end if    
        if(Sim1m(it2,it,i,4)>1d-3) then
            YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        else
            YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        end if      
    end if
end do
end do

end do

CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)

dum2=1d0-SSE/SST

dum10=dum10+((dum2-0.408)/0.408)**2

write(file_id, *)'Persistence of single male LFP is',BREG(2)

write(file_id, *)'Single male LFP R2 is',dum2


!Married Male Labor Supply

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        dum2=dum2+Sim1m(it2,it,i,4)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Married Male labor supply is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.349d0)/0.349d0)**2d0
dum10=dum10+((dum2-0.349d0)/0.349d0)**2d0


!Married male labor force participation

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        if(Sim1m(it2,it,i,4)>1d-3) then
            dum2=dum2+(1d0)*WeightActive(i)
        end if
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

dum10=dum10+((dum2-0.871d0)/0.871d0)**2d0

write(file_id, *)'Married male labor force participation is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.871d0)/0.871d0)**2d0
YVAR=sqrt(-1.0)
XVARS=sqrt(-1.0)

do i=2,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        if(Sim1m(it2,it,i-1,4)>1d-3) then
            XVARs((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        else
            XVARs((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        end if    
        if(Sim1m(it2,it,i,4)>1d-3) then
            YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        else
            YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        end if      
    end if
end do
end do

end do

CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)

dum2=1d0-SSE/SST

dum10=dum10+((dum2-0.457)/0.457)**2

write(file_id, *)'Persistence of married male LFP is',BREG(2)

write(file_id, *)'Married male LFP R2 is',dum2


!Female labor supply

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
    dum3=dum3+1d0*WeightActive(i)
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Female labor supply is',dum2

!Single female labor supply

dum2=0d0
dum3=0d0
dum4=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)<0.5) then
        dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Single female labor supply is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.236d0)/0.236d0)**2d0
dum10=dum10+((dum2-0.236d0)/0.236d0)**2d0

!Variance of single female hours

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)<0.5) then
        dum4=dum4+((Sim1f(it2,it,i,4)-dum2)**2)*WeightActive(i)
    end if
end do
end do

end do

dum4=dum4/dum3

write(file_id, *)'Stdev single female labor supply is',SQRT(dum4)


!Variance of single female hours
! >>>>>
!do i=10,10
!
!do it2=1,nsim2
!do it=1,nsim
!    if(Sim1f(it2,it,i,10)<0.5) then
!        dum4=dum4+((Sim1f(it2,it,i,4)-dum2)**2)*WeightActive(i)
!    end if
!end do
!end do
!
!end do
!
!dum4=dum4/dum3
!
!write(file_id, *)'Stdev single female labor supply at age 30 is',SQRT(dum4)
! >>>>>

!Married female labor supply

dum2=0d0
dum3=0d0
dum4=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)>0.5) then
        dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Married female labor supply is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.231d0)/0.231d0)**2d0
dum10=dum10+((dum2-0.231d0)/0.231d0)**2d0

!Variance of married female hours

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)>0.5) then
        dum4=dum4+((Sim1f(it2,it,i,4)-dum2)**2)*WeightActive(i)
    end if
end do
end do

end do

dum4=dum4/dum3

write(file_id, *)'Stdev married female labor supply is',SQRT(dum4)

!Female labor force participation

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,4)>1d-3) then
        dum2=dum2+(1d0)*WeightActive(i)
    end if
    dum3=dum3+1d0*WeightActive(i)
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Female labor force participation is',dum2

!Single female labor force participation

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)<0.5) then
        if(Sim1f(it2,it,i,4)>1d-3) then
            dum2=dum2+(1d0)*WeightActive(i)
        end if
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Single female labor force participation is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.694d0)/0.694d0)**2d0
dum10=dum10+((dum2-0.694d0)/0.694d0)**2d0

YVAR=sqrt(-1.0)
XVARS=sqrt(-1.0)

do i=2,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)<0.5) then
        if(Sim1f(it2,it,i-1,4)>1d-3) then
            XVARs((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        else
            XVARs((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        end if    
        if(Sim1f(it2,it,i,4)>1d-3) then
            YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        else
            YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        end if      
    end if
end do
end do

end do

CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)

dum2=1d0-SSE/SST

dum10=dum10+((dum2-0.463)/0.463)**2

write(file_id, *)'Persistence of single female LFP is',BREG(2)

write(file_id, *)'Single female LFP R2 is',dum2

!Married female labor force participation

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)>0.5) then
        if(Sim1f(it2,it,i,4)>1d-3) then
            dum2=dum2+(1d0)*WeightActive(i)
        end if
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Married female labor force participation is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.668d0)/0.668d0)**2d0
dum10=dum10+((dum2-0.668d0)/0.668d0)**2d0

YVAR=sqrt(-1.0)
XVARS=sqrt(-1.0)

do i=2,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)>0.5) then
        if(Sim1f(it2,it,i-1,4)>1d-3) then
            XVARS((i-1)*it*it2+(it2-1)*it+it,1)=1d0
        else
            XVARS((i-1)*it*it2+(it2-1)*it+it,1)=0d0
        end if      
        if(Sim1f(it2,it,i,4)>1d-3) then
            YVAR((i-1)*it*it2+(it2-1)*it+it)=1d0
        else
            YVAR((i-1)*it*it2+(it2-1)*it+it)=0d0
        end if      
    end if
end do
end do

end do

CALL RLSE (YVAR, XVARS, BREG, SST=SST, SSE=SSE)

dum2=1d0-SSE/SST

dum10=dum10+((dum2-0.553)/0.553)**2

write(file_id, *)'Persistence of married female LFP is',BREG(2)

write(file_id, *)'Married female LFP R2 is',dum2

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1f(it2,it,i,4)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Female intensive margin is',dum2

!Variance of log male earnings

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,5)>0d0) then
    dum2=dum2+log(Sim1m(it2,it,i,5))*WeightActive(i)
    dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

do i=1,T

do it2=1,nsim2    
do it=1,nsim
    if(Sim1m(it2,it,i,5)>0d0) then
    dum4=dum4+((log(Sim1m(it2,it,i,5))-dum2)**2)*WeightActive(i)
    end if
end do
end do

end do

dum4=dum4/dum3

write(file_id, *)'Stdev of log male earnings is',SQRT(dum4)

dum2=0d0
dum3=0d0
dum4=0d0

!Variance of log female earnings

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,5)>1d-3) then
        dum2=dum2+log(Sim1f(it2,it,i,5))*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,5)>1d-3) then
        dum4=dum4+((log(Sim1f(it2,it,i,5))-dum2)**2)*WeightActive(i)
    end if
end do
end do

end do

dum4=dum4/dum3

write(file_id, *)'Stdev of log female earnings is',SQRT(dum4)

!Variance of log male wage

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,5)>0d0) then
    dum2=dum2+log(Sim1m(it2,it,i,3))*WeightActive(i)
    dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

do i=1,T

do it2=1,nsim2    
do it=1,nsim
    if(Sim1m(it2,it,i,5)>0d0) then
    dum4=dum4+((log(Sim1m(it2,it,i,3))-dum2)**2)*WeightActive(i)
    end if
end do
end do

end do

dum4=dum4/dum3
dum4=SQRT(dum4)

write(file_id, *)'Stdev of log male wage is',dum4

!dum10=dum10+((dum4-0.776d0)/0.776d0)**2d0

dum2=0d0
dum3=0d0
dum4=0d0

!Variance of log female wage

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,5)>1d-3) then
        dum2=dum2+log(Sim1f(it2,it,i,3))*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,5)>1d-3) then
        dum4=dum4+((log(Sim1f(it2,it,i,3))-dum2)**2)*WeightActive(i)
    end if
end do
end do

end do

dum4=dum4/dum3

dum4=SQRT(dum4)

write(file_id, *)'Stdev of log female wage is',dum4

!dum10=dum10+((dum4-0.724d0)/0.724d0)**2d0

!Male experience

dum2=0d0
dum3=0d0
dum4=0d0

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+1d0*(exp2m(it2,it,T+1,1)-1)
    dum3=dum3+1d0
end do
end do

dum2=dum2/dum3

write(file_id, *)'Male experience at age 65 is',dum2

!Female experience

dum2=0d0
dum3=0d0
dum4=0d0

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+1d0*(exp2f(it2,it,T+1,1)-1)
    dum3=dum3+1d0
end do
end do

dum2=dum2/dum3

write(file_id, *)'Female experience at age 65 is',dum2

!Male wage

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1m(it2,it,i,3)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Average male wage is',dum2

dum7=dum2

!Female wage

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1f(it2,it,i,3)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Average female wage is',dum2

write(file_id, *)'Male_wage/Female_wage is',dum7/dum2

write(file_id, *)'Average male wage is',dum2

!Wages at 0.5AE and 2AE

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    
    if((Sim1m(it2,it,i,5)>0.4d0).AND.(Sim1m(it2,it,i,5)<0.6d0)) then
        dum2=dum2+Sim1m(it2,it,i,3)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
    
    if((Sim1f(it2,it,i,5)>0.4d0).AND.(Sim1f(it2,it,i,5)<0.6d0)) then
        dum2=dum2+Sim1f(it2,it,i,3)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
    
    if((Sim1m(it2,it,i,5)>1.9d0).AND.(Sim1m(it2,it,i,5)<2.1d0)) then
        dum4=dum4+Sim1m(it2,it,i,3)*WeightActive(i)
        dum5=dum5+1d0*WeightActive(i)
    end if
    
    if((Sim1f(it2,it,i,5)>1.9d0).AND.(Sim1f(it2,it,i,5)<2.1d0)) then
        dum4=dum4+Sim1f(it2,it,i,3)*WeightActive(i)
        dum5=dum5+1d0*WeightActive(i)
    end if
    
end do
end do

end do


write(file_id, *)'Wage rate of people making 0.5AE is',dum2/dum3
write(file_id, *)'Wage rate of people making 2AE is',dum4/dum5

!Male earnings

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Average male earnings is',dum2

dum7=dum2

!Female earnings

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1f(it2,it,i,5)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

dum7=dum7/dum2

write(file_id, *)'Average female earnings is',dum2

write(file_id, *)'Male_earnings/Female_earnings is',dum7
!write(file_id, *)'contribution to FCN is is',((dum7-1.432d0)/1.432d0)**2d0
dum10=dum10+((dum7-1.432d0)/1.432d0)**2d0

!Male earnings Age 35-45

dum2=0d0
dum3=0d0
! >>>>>
!do i=16,26
!
!do it2=1,nsim2
!do it=1,nsim
!    if(Sim1m(it2,it,i,4)>1d-3) then
!        dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
!        dum3=dum3+1d0*WeightActive(i)
!    end if
!end do
!end do
!
!end do
!
!dum2=dum2/dum3
!
!write(file_id, *)'Male earnings 35-45 is',dum2
! >>>>>
!write(file_id, *)'contribution to FCN is is',((dum2-1.258d0)/1.258d0)**2d0
!dum10=dum10+((dum2-1.258d0)/1.258d0)**2d0

!Male earnings Age 55-64

dum2=0d0
dum3=0d0

do i=Tret,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Male earnings 55-64 is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-1.289d0)/1.289d0)**2d0
!dum10=dum10+((dum2-1.289d0)/1.289d0)**2d0


!Female earnings age 35-45

dum2=0d0
dum3=0d0

! >>>>>
!do i=16,26
!
!do it2=1,nsim2
!do it=1,nsim
!    if(Sim1f(it2,it,i,4)>1d-3) then
!        dum2=dum2+Sim1f(it2,it,i,5)*WeightActive(i)
!        dum3=dum3+1d0*WeightActive(i)
!    end if
!end do
!end do
!
!end do
!
!dum2=dum2/dum3
!
!write(file_id, *)'Female earnings 35-45 is',dum2
! >>>>>

!write(file_id, *)'contribution to FCN is is',((dum2-0.87d0)/0.87d0)**2d0
!dum10=dum10+((dum2-0.87d0)/0.87d0)**2d0

!Female earnings age 55-64

dum2=0d0
dum3=0d0

do i=Tret,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,4)>1d-3) then
        dum2=dum2+Sim1f(it2,it,i,5)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Female earnings 55-64 is',dum2
!write(file_id, *)'contribution to FCN is is',((dum2-0.868d0)/0.868d0)**2d0
!dum10=dum10+((dum2-0.868d0)/0.868d0)**2d0

!Correlation in spousal ability
it4=0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        !if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it3,it,i,4)>0.001d0)) then
            it4=it4+1
        !end if
    end if
end do
end do

end do

allocate(spousewage(it4,2))

it4=0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        
        it3=exp1m(it2,it,i,4)
        !if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it3,it,i,4)>0.001d0)) then
            it4=it4+1
            Spousewage(it4,1)=A(1,exp1m(it2,it,i,2))
            Spousewage(it4,2)=A(2,exp1f(it2,it3,i,2))
        !end if
    end if
end do
end do

end do

CALL D_CORVC(NVAR, Spousewage, COV, ICOPT=ICOPT)

write(file_id, *)'Correlation of spousal ability is',COV(1,2)

!Correlation in spousal education

it4=0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        
        it3=exp1m(it2,it,i,4)
        !if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it3,it,i,4)>0.001d0)) then
            it4=it4+1
            Spousewage(it4,1)=exp1m(it2,it,i,2)*1d0
            Spousewage(it4,2)=exp1f(it2,it3,i,2)*1d0
        !end if
    end if
end do
end do

end do

CALL D_CORVC(NVAR, Spousewage, COV, ICOPT=ICOPT)

write(file_id, *)'Correlation of spousal education is',COV(1,2)
!write(file_id, *)'contribution to FCN is is',((COV(1,2)-0.646)/0.646)**2
!dum10=dum10+((COV(1,2)-0.646)/0.646)**2

it4=0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        
        it3=exp1m(it2,it,i,4)
        if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
            it4=it4+1
        end if
    end if
end do
end do

end do

allocate(spousewage2(it4,2))

it4=0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        
        it3=exp1m(it2,it,i,4)
        if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
            it4=it4+1
            Spousewage2(it4,1)=Sim1m(it2,it,i,3)
            Spousewage2(it4,2)=Sim1f(it2,it3,i,3)
        end if
    end if
end do
end do

end do

CALL D_CORVC(NVAR, Spousewage2, COV, ICOPT=ICOPT)
!
write(file_id, *)'Correlation of spousal wages is',COV(1,2)


dum10=dum10+((COV(1,2)-0.287)/0.287)**2

!it4=0
!
!do i=1,T
!
!do it2=1,nsim2
!do it=1,nsim
!    if(Sim1m(it2,it,i,10)>0.5) then
!        
!        it3=exp1m(it2,it,i,4)
!        if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
!            it4=it4+1
!            Spousewage2(it4,1)=Sim1m(it2,it,i,5)
!            Spousewage2(it4,2)=Sim1f(it2,it3,i,5)
!        end if
!    end if
!end do
!end do
!
!end do
!
!CALL D_CORVC(NVAR, Spousewage2, COV, ICOPT=ICOPT)
!
!write(file_id, *)'Correlation of spousal earnings is',COV(1,2)

!Fractions of 2- and 1- earner households

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum7=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5) then
        
        it3=exp1m(it2,it,i,4)
        dum3=dum3+1d0
        
        if((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
            dum2=dum2+1d0
        elseif((Sim1m(it2,it,i,4)>0.001).AND.(Sim1f(it2,it3,i,4)<0.001d0)) then
            dum4=dum4+1d0
        elseif((Sim1m(it2,it,i,4)<0.001).AND.(Sim1f(it2,it3,i,4)>0.001d0)) then
            dum5=dum5+1d0
        else
            dum7=dum7+1d0
        end if
    end if
    
end do
end do

end do

write(file_id, *)'Fraction of married households with 2 earners',dum2/dum3
write(file_id, *)'Fraction of married households with only male earner',dum4/dum3
write(file_id, *)'Fraction of married households with only female earner',dum5/dum3
write(file_id, *)'Fraction of married households with no earners',dum7/dum3

!Total labor income taxes and tax revenue

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum7=0d0
dum15=0d0
dum16=0d0
dum17=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum5=dum5+2d0*WeightActive(i)
    dum2=dum2+Sim1m(it2,it,i,6)*(1d0+t_employer)*WeightActive(i)
    dum7=dum7+Sim1m(it2,it,i,7)*WeightActive(i)
    dum15=dum15+Sim1m(it2,it,i,9)*WeightActive(i)
    dum4=dum4+Sim1m(it2,it,i,9)*WeightActive(i)+Sim1m(it2,it,i,7)*WeightActive(i)
    dum3=dum3+Sim1m(it2,it,i,9)*WeightActive(i)+Sim1m(it2,it,i,7)*WeightActive(i)+Sim1m(it2,it,i,8)*WeightActive(i)+(Sim1m(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1f(it2,it,i,6)*(1d0+t_employer)*WeightActive(i)
        dum7=dum7+Sim1f(it2,it,i,7)*WeightActive(i)
        dum15=dum15+Sim1f(it2,it,i,9)*WeightActive(i)
        dum4=dum4+Sim1f(it2,it,i,9)*WeightActive(i)+Sim1f(it2,it,i,7)*WeightActive(i)
        dum3=dum3+Sim1f(it2,it,i,9)*WeightActive(i)+Sim1f(it2,it,i,7)*WeightActive(i)+Sim1f(it2,it,i,8)*WeightActive(i)+(Sim1f(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
    end if
end do
end do

end do

do i=1,Tret
    
r_ret=r

do it2=1,nsim2
do it=1,nsim
    dum5=dum5+2d0*WeightRet(i)
    dum3=dum3+SimR1m(it2,it,i,3)*WeightRet(i)+(SimR1m(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum3=dum3+SimR1f(it2,it,i,3)*WeightRet(i)+(SimR1f(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
    end if
end do
end do

end do

write(file_id, *)'Labor income tax rate including TSS is',dum4/dum2

write(file_id, *)'Tax revenue per capita including TSS is',dum3/dum5

!write(file_id, *)'Labor Income tax per capita is',dum7/dum5

!write(file_id, *)'Social security tax per capita is',dum15/dum5


!Social security

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum7=0d0
dum15=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum3=dum3+2d0*WeightActive(i)
    dum4=dum4+Sim1m(it2,it,i,9)*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum4=dum4+Sim1f(it2,it,i,9)*WeightActive(i)
    end if
    
end do
end do

end do
    

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    if(SimR1m(it2,it,i,5)<1d-3) then
         dum2=dum2+1d0*WeightRet(i)
         dum5=dum5+SimR1m(it2,it,i,14)*WeightRet(i)
    end if
    
    if(SimR1f(it2,it,i,5)<1d-3) then
         dum2=dum2+1d0*WeightRet(i)
         dum5=dum5+SimR1f(it2,it,i,14)*WeightRet(i)
    end if
    
end do
end do

end do

dum4=dum4/dum2
write(file_id, *)'SS tax per retiree is',dum4
write(file_id, *)'Average pension is',dum5/dum2

dum5=dum5/dum3

ss_expense=dum5

write(file_id, *)'Social Security expenses per capita is',dum5
!write(file_id, *)'Pension',Psi_pension/2d0

epsilon=Psi0-dum4

!epsilon=0d0
Psi0=Psi0-0.1d0*(Psi0-dum4)

!Savings

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,1)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1f(it2,it,i,1)*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+SimR1m(it2,it,i,1)*WeightRet(i)
    dum3=dum3+2d0*WeightRet(i)
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+SimR1f(it2,it,i,1)*WeightRet(i)
    end if
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Savings per capita is',dum2
dum6=dum2

! Assets for redistribution

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i+1,1)*WeightActive(i)*(1d0-OmegaActive(i))
    dum3=dum3+1d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1f(it2,it,i+1,1)*WeightActive(i)*(1d0-OmegaActive(i))
    end if
    
end do
end do

end do

do i=1,36

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+SimR1m(it2,it,i+1,1)*WeightRet(i)*(1d0-OmegaRet(i))
    dum3=dum3+1d0*WeightRet(i)
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+SimR1f(it2,it,i+1,1)*WeightRet(i)*(1d0-OmegaRet(i))
    end if
    
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Assets per capita to be redistributed',dum2/2d0
write(file_id, *)'Gamma_redistr',Gamma_redistr/2d0

epsilon2=Gamma_redistr-dum2
!epsilon2=0d0
Gamma_redistr=Gamma_redistr-0.3d0*(Gamma_redistr-dum2)

! Capital tax

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2 
do it=1,nsim
    dum2=dum2+(Sim1m(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
    dum3=dum3+2d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+(Sim1f(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
    end if
    
end do
end do

end do

do i=1,Tret

r_ret=r    
    
do it2=1,nsim2
do it=1,nsim
    dum2=dum2+(SimR1m(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
    dum3=dum3+2d0*WeightRet(i)
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+(SimR1f(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
    end if
    
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Capital tax per capita is',dum2

dum7=dum2

!Consumption

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,2)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1f(it2,it,i,2)*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+SimR1m(it2,it,i,2)*WeightRet(i)
    dum3=dum3+2d0*WeightRet(i)
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+SimR1f(it2,it,i,2)*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Consumption per capita is',dum2

dum2=0d0
dum3=0d0

!Consumption tax

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,8)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1f(it2,it,i,8)*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+SimR1m(it2,it,i,3)*WeightRet(i)
    dum3=dum3+2d0*WeightRet(i)
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+SimR1f(it2,it,i,3)*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Consumption tax per capita is',dum2

!Labor income tax

dum4=0d0
dum5=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum4=dum4+Sim1m(it2,it,i,7)*WeightActive(i)
    dum5=dum5+2d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum4=dum4+Sim1f(it2,it,i,7)*WeightActive(i)
    end if
    
end do
end do

end do
 

dum4=dum4/dum5

write(file_id, *)'Labor income tax per capita is',dum4

write(file_id, *)'Tax revenue per capita is',dum4+dum2+dum7

dum5=dum4+dum2

dum2=0d0
dum3=0d0
dum15=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum3=dum3+2d0*WeightActive(i)
    if(Sim1f(it2,it,i,4)<0.001) then
        dum15=dum15+Unemp_benefit*WeightActive(i)
    end if
end do
end do

end do


dum15=dum15/dum3

!GDP per capita

dum9=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum9=dum9+(Sim1m(it2,it,i,6)*(1d0+t_employer)/w)*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum9=dum9+(Sim1f(it2,it,i,6)*(1d0+t_employer)/w)*WeightActive(i)
    end if
    
end do
end do

end do

    
    
!write(file_id, *)'Ltot is',dum9

dum3=((ratio*dum9)**alpha)*(dum9**(1-alpha))/dum3

write(file_id, *)'GDP per capita is',dum3

write(file_id, *)'Lumpsum is',lumpsum/2d0


!Government Budget

lumpsumdum=(dum5+dum7)+mu*debttoGDP*dum3-(dum15+r*debttoGDP*dum3+2d0*milspendtoGDP*dum3)

write(file_id, *)'Net revenue is',lumpsumdum

lumpsumdum=lumpsumdum*2d0

epsilon3=lumpsum-lumpsumdum

!epsilon3=0d0
lumpsum=lumpsum-0.1d0*(lumpsum-lumpsumdum)



!Labor income tax level

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,6)*WeightActive(i)
    dum4=dum4+Sim1m(it2,it,i,7)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1f(it2,it,i,6)*WeightActive(i)
        dum4=dum4+Sim1f(it2,it,i,7)*WeightActive(i)
    end if
    
end do
end do

end do

    
write(file_id, *)'Average labor income tax rate is',dum4/dum2

write(file_id, *)'Average individual earnings is',dum2/dum3

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,4)>0.001) then
        dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
        dum3=dum3+WeightActive(i)
    end if
    if(Sim1f(it2,it,i,4)>0.001) then
        dum2=dum2+Sim1f(it2,it,i,5)*WeightActive(i)
        dum3=dum3+WeightActive(i)
    end if
end do
end do

end do

    
write(file_id, *)'Average individual earnings for working people is',dum2/dum3

write(file_id, *)'AE is',AE
!write(file_id, *)'contribution to FCN is is',((dum2/dum3-1d0)/1d0)**2
dum10=dum10+((dum2/dum3-1d0)/1d0)**2

epsilon5=AE-dum2/dum3

AE=AE-0.1d0*(AE-dum2/dum3)

Unemp_benefit=0.201795*AE

!Filling the average earnings matrix to be used for pensions


do ia=1,na

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

do i=11,T

do it2=1,nsim2
do it=1,nsim

if(exp1m(it2,it,i,2)==ia) then
    if(Sim1m(it2,it,i,10)>0.5) then    
        if(Sim1m(it2,it,i,4)>0.001) then
            dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
            dum3=dum3+WeightActive(i)
        end if
    end if
end if
    
if(exp1f(it2,it,i,2)==ia) then
    if(Sim1f(it2,it,i,10)>0.5) then
        if(Sim1f(it2,it,i,4)>0.001) then
            dum4=dum4+Sim1f(it2,it,i,5)*WeightActive(i)
            dum5=dum5+WeightActive(i)
        end if
    end if
end if
    
    
end do
end do

end do

!av_earnings(1,1,ia)=dum2/dum3
!av_earnings(2,1,ia)=dum4/dum5

end do


do ia=1,na

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

do i=11,T

do it2=1,nsim2
do it=1,nsim

if(exp1m(it2,it,i,2)==ia) then
    if(Sim1m(it2,it,i,10)<0.5) then    
        if(Sim1m(it2,it,i,4)>0.001) then
            dum2=dum2+Sim1m(it2,it,i,5)*WeightActive(i)
            dum3=dum3+WeightActive(i)
        end if
    end if
end if
    
if(exp1f(it2,it,i,2)==ia) then
    if(Sim1f(it2,it,i,10)<0.5) then
        if(Sim1f(it2,it,i,4)>0.001) then
            dum4=dum4+Sim1f(it2,it,i,5)*WeightActive(i)
            dum5=dum5+WeightActive(i)
        end if
    end if
end if
    
    
end do
end do

end do

!av_earnings(1,2,ia)=dum2/dum3
!av_earnings(2,2,ia)=dum4/dum5

end do


!Prices

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+(Sim1m(it2,it,i,1)+Gamma_redistr)*WeightActive(i)
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum2=dum2+(Sim1f(it2,it,i,1)+Gamma_redistr)*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+(SimR1m(it2,it,i,1)+Gamma_redistr)*WeightRet(i)
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+(SimR1f(it2,it,i,1)+Gamma_redistr)*WeightRet(i)
    end if
    
end do
end do

end do

dum6=dum2-debttoGDP*((ratio*dum9)**(alpha))*(dum9**(1d0-alpha))

!write(file_id, *)'Ktot is',dum6

ratiodum=dum6/dum9

write(file_id, *)'Ratio between capital and labor is',ratiodum

write(file_id, *)'Implied wage is',(1d0-alpha)*ratiodum**alpha
write(file_id, *)'Implied interest is',alpha*ratiodum**(alpha-1d0)-delta

write(file_id, *)'Wage is',w
write(file_id, *)'Interest is',r

write(file_id, *)'K/Y is',ratiodum**(1d0-alpha)
!write(file_id, *)'contribution to FCN is is',(((ratiodum**(1d0-alpha))-2.683d0)/2.683d0)**2d0
dum10=dum10+(((ratiodum**(1d0-alpha))-2.683d0)/2.683d0)**2d0

write(file_id, *)'FCN is',dum10

! beta=1.00251
! Fw=0.0125
! Chi=26.1d0
! gamma0=0.3288d0
! FCN=0.00000559

!write(file_id, *)'I/Y is',delta*ratiodum**(1d0-alpha)

dum9=max(maxval(Sim1m(:,:,:,1)),maxval(Sim1f(:,:,:,1)))

dum9=max(dum9,maxval(SimR1m(:,:,:,1)),maxval(SimR1f(:,:,:,1)))

write(file_id, *)'Max savings is',dum9

dum2=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,i,10)>0.5d0) then
        dum2=dum2+1d0*WeightActive(i)
    end if
    
end do
end do

end do

write(file_id, *)'Fraction of Married females is',dum2/(1d0*T*nsim2*nsim)

dum2=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    if(Sim1m(it2,it,i,10)>0.5d0) then
        dum2=dum2+1d0*WeightActive(i)
    end if
    
end do
end do

end do

write(file_id, *)'Fraction of Married males is',dum2/(1d0*T*nsim2*nsim)

!dum2=0d0
!dum3=0d0
!dum4=0d0
!dum5=0d0
!
!do i=1,T
!
!do it2=1,nsim2
!do it=1,nsim
!    if((Sim1f(it2,it,i,10)>0.5d0).AND.(exp1f(it2,it,i,6)==1)) then
!        dum2=dum2+1d0*WeightActive(i)
!        dum5=dum5+1d0*WeightActive(i)
!    elseif((Sim1f(it2,it,i,10)>0.5d0).AND.(exp1f(it2,it,i,6)==2)) then
!        dum3=dum3+1d0*WeightActive(i)
!        dum5=dum5+1d0*WeightActive(i)
!    elseif((Sim1f(it2,it,i,10)>0.5d0).AND.(exp1f(it2,it,i,6)==3)) then
!        dum4=dum4+1d0*WeightActive(i)
!        dum5=dum5+1d0*WeightActive(i)
!    end if
!end do
!end do
!
!end do
!
!write(file_id, *)'Fraction of married fixed cost 1 is',dum2/dum5
!write(file_id, *)'Fraction of married fixed cost 2 is',dum3/dum5
!write(file_id, *)'Fraction of married fixed cost 3 is',dum4/dum5
!
!dum2=0d0
!dum3=0d0
!dum4=0d0
!dum5=0d0
!
!do i=1,T
!
!do it2=1,nsim2
!do it=1,nsim
!    if((Sim1f(it2,it,i,10)<0.5d0).AND.(exp1f(it2,it,i,6)==1)) then
!        dum2=dum2+1d0*WeightActive(i)
!        dum5=dum5+1d0*WeightActive(i)
!    elseif((Sim1f(it2,it,i,10)<0.5d0).AND.(exp1f(it2,it,i,6)==2)) then
!        dum3=dum3+1d0*WeightActive(i)
!        dum5=dum5+1d0*WeightActive(i)
!    elseif((Sim1f(it2,it,i,10)<0.5d0).AND.(exp1f(it2,it,i,6)==3)) then
!        dum4=dum4+1d0*WeightActive(i)
!        dum5=dum5+1d0*WeightActive(i)
!    end if
!end do
!end do
!
!end do
!
!write(file_id, *)'Fraction of single fixed cost 1 is',dum2/dum5
!write(file_id, *)'Fraction of single fixed cost 2 is',dum3/dum5
!write(file_id, *)'Fraction of single fixed cost 3 is',dum4/dum5


mpartnerdum=0d0
fpartnerdum=0d0

do ia2=1,na

do i=1,T

dum3=0d0    
    
do it2=1,nsim2
do it=1,nsim
    
if((Sim1m(it2,it,i,10)<0.5).AND.(ia2==exp1m(it2,it,i,2))) then
    dum2=Sim1m(it2,it,i,1)
    ia=ia2
    iu=exp1m(it2,it,i,3)
    ifc=exp1m(it2,it,i,6)
    dum5=exp2m(it2,it,i,1)
    
    if(i==1) then
        ix=1
    end if
    
    if(dum5<exp_grid(2,i)/2d0) then
        ix=1
    end if
    
    do ik=2,nexp-1  
    if(((exp_grid(ik,i)-(exp_grid(ik,i)-exp_grid(ik-1,i))/2d0)<dum5).AND.(dum5<(exp_grid(ik,i)+(exp_grid(ik+1,i)-exp_grid(ik,i))/2d0))) then
        ix=ik
    end if
    end do

    if(dum5>(exp_grid(nexp,i)-(exp_grid(nexp,i)-exp_grid(nexp-1,i))/2d0)) then
        ix=nexp
    end if
    
    if(dum2<k_grid(2)/2d0) then
        mpartnerdum(1,ix,ia,iu,i,ifc)=mpartnerdum(1,ix,ia,iu,i,ifc)+1d0
    end if
    
    do ik=2,nk-1  
    if(((k_grid(ik)-(k_grid(ik)-k_grid(ik-1))/2d0)<dum2).AND.(dum2<(k_grid(ik)+(k_grid(ik+1)-k_grid(ik))/2d0))) then
        mpartnerdum(ik,ix,ia,iu,i,ifc)=mpartnerdum(ik,ix,ia,iu,i,ifc)+1d0
    end if
    end do

    if(dum2>(k_grid(nk)-(k_grid(nk)-k_grid(nk-1))/2d0)) then
        mpartnerdum(nk,ix,ia,iu,i,ifc)=mpartnerdum(nk,ix,ia,iu,i,ifc)+1d0
    end if
    
    dum3=dum3+1d0
end if
    
end do
end do

mpartnerdum(:,:,ia2,:,i,:)=mpartnerdum(:,:,ia2,:,i,:)/dum3

end do

end do

do ia2=1,na
    
do i=1,T

dum3=0d0     
        
do it2=1,nsim2
do it=1,nsim
    
if((Sim1f(it2,it,i,10)<0.5).AND.(ia2==exp1f(it2,it,i,2))) then
    dum2=Sim1f(it2,it,i,1)
    ia=ia2
    iu=exp1f(it2,it,i,3)
    ifc=exp1f(it2,it,i,6)
    dum5=exp2f(it2,it,i,1)
    dum4=0d0
    
    if(i==1) then
        ix=1
    end if
    
    if(dum5<exp_grid(2,i)/2d0) then
        ix=1
    end if
    
    do ik=2,nexp-1  
    if(((exp_grid(ik,i)-(exp_grid(ik,i)-exp_grid(ik-1,i))/2d0)<dum5).AND.(dum5<(exp_grid(ik,i)+(exp_grid(ik+1,i)-exp_grid(ik,i))/2d0))) then
        ix=ik
    end if
    end do

    if(dum5>(exp_grid(nexp,i)-(exp_grid(nexp,i)-exp_grid(nexp-1,i))/2d0)) then
        ix=nexp
    end if
    
    if(dum2<k_grid(2)/2d0) then
        fpartnerdum(1,ix,ia,iu,i,ifc)=fpartnerdum(1,ix,ia,iu,i,ifc)+1d0
        dum4=1d0
    end if
    
    do ik=2,nk-1
    if(((k_grid(ik)-(k_grid(ik)-k_grid(ik-1))/2d0)<dum2).AND.(dum2<(k_grid(ik)+(k_grid(ik+1)-k_grid(ik))/2d0))) then
        fpartnerdum(ik,ix,ia,iu,i,ifc)=fpartnerdum(ik,ix,ia,iu,i,ifc)+1d0
        dum4=1d0
    end if
    end do

    if(dum2>(k_grid(nk)-(k_grid(nk)-k_grid(nk-1))/2d0)) then
        fpartnerdum(nk,ix,ia,iu,i,ifc)=fpartnerdum(nk,ix,ia,iu,i,ifc)+1d0
        dum4=1d0
    end if
    
    dum3=dum3+1d0
end if

end do
end do
    
     fpartnerdum(:,:,ia2,:,i,:)=fpartnerdum(:,:,ia2,:,i,:)/dum3

end do

end do


fpartnerdum2=abs(fpartnerdum-fpartner)

mpartnerdum2=abs(mpartnerdum-mpartner)

epsilon6=max(maxval(fpartnerdum2),maxval(mpartnerdum2))

write(file_id, *)'max distance between distribution of singles is',epsilon6

!if(iter>5) then
    epsilon6=0d0
!end if

fpartner=fpartnerdum

mpartner=mpartnerdum

dum2=0d0
! >>>>>
!do ifc=1,nfc
!do ik=1,nk
!    do ia=1,na
!        do ix=1,nexp
!            do iu=1,nu
!                dum2=dum2+fpartner(ik,ix,ia,iu,30,ifc)
!            end do
!        end do
!    end do
!end do
!end do
!
!write(file_id, *)'sum fpartner is',dum2
! >>>>>

!dum2=0d0
!
!do ifc=1,nfcm
!do ik=1,nk
!    do ia=1,na
!        do ix=1,nexp
!            do iu=1,nu
!                dum2=dum2+mpartner(ik,ix,ia,iu,30,ifc)
!            end do
!        end do
!    end do
!end do
!end do
!
!write(file_id, *)'sum mpartner is',dum2

!Here we are computing the matrix of  probabilities for marrying someone of ability a' if you have probability a.

ability_prob=0d0

do ia=1,na
    
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    
if((Sim1m(it2,it,i,10)>0.5).AND.(ia==exp1m(it2,it,i,2))) then
    
    it4=exp1m(it2,it,i,4)
    ia2=exp1f(it2,it4,i,2)
    
    ability_prob(ia,ia2)=ability_prob(ia,ia2)+1d0
    dum3=dum3+1d0
    
end if
    
end do
end do
    
end do

ability_prob(ia,:)=ability_prob(ia,:)/dum3

end do

write(file_id, *)'ability_prob is',ability_prob(2,:)

!dum2=0d0
!dum3=0d0
!dum4=0d0
!dum5=0d0
!dum6=0d0
!dum7=0d0
!
!do i=1,T
!
!do it=1,nsim
!    dum7=dum7+1d0*WeightActive(i)
!    if(exp1(it,i,2)==1) then
!        dum2=dum2+1d0*WeightActive(i)
!    elseif(exp1(it,i,2)==2) then
!        dum3=dum3+1d0*WeightActive(i)
!    elseif(exp1(it,i,2)==3) then
!        dum4=dum4+1d0*WeightActive(i)
!    elseif(exp1(it,i,2)==4) then
!        dum5=dum5+1d0*WeightActive(i)
!    else
!        dum6=dum6+1d0*WeightActive(i)
!    end if
!end do
!
!end do
!
!dum2=dum2/dum7
!dum3=dum3/dum7
!dum4=dum4/dum7
!dum5=dum5/dum7
!dum6=dum6/dum7
!
!write(file_id, *)'Distribution of male ability',dum2,dum3,dum4,dum5,dum6
!
!dum2=0d0
!dum3=0d0
!dum4=0d0
!dum5=0d0
!dum6=0d0
!dum7=0d0
!
!do i=1,T
!
!do it=1,nsim
!    dum7=dum7+1d0*WeightActive(i)
!    if(exp1(it,i,3)==1) then
!        dum2=dum2+1d0*WeightActive(i)
!    elseif(exp1(it,i,3)==2) then
!        dum3=dum3+1d0*WeightActive(i)
!    elseif(exp1(it,i,3)==3) then
!        dum4=dum4+1d0*WeightActive(i)
!    elseif(exp1(it,i,3)==4) then
!        dum5=dum5+1d0*WeightActive(i)
!    else
!        dum6=dum6+1d0*WeightActive(i)
!    end if
!end do
!
!end do
!
!dum2=dum2/dum7
!dum3=dum3/dum7
!dum4=dum4/dum7
!dum5=dum5/dum7
!dum6=dum6/dum7
!
!write(file_id, *)'Distribution of male idiosyncraic shock',dum2,dum3,dum4,dum5,dum6
!
!dum2=0d0
!dum3=0d0
!dum4=0d0
!dum5=0d0
!dum6=0d0
!dum7=0d0
!
!do i=1,T
!
!do it=1,nsim
!    dum7=dum7+1d0*WeightActive(i)
!    if(exp1(it,i,4)==1) then
!        dum2=dum2+1d0*WeightActive(i)
!    elseif(exp1(it,i,4)==2) then
!        dum3=dum3+1d0*WeightActive(i)
!    elseif(exp1(it,i,4)==3) then
!        dum4=dum4+1d0*WeightActive(i)
!    elseif(exp1(it,i,4)==4) then
!        dum5=dum5+1d0*WeightActive(i)
!    else
!        dum6=dum6+1d0*WeightActive(i)
!    end if
!end do
!
!end do
!
!dum2=dum2/dum7
!dum3=dum3/dum7
!dum4=dum4/dum7
!dum5=dum5/dum7
!dum6=dum6/dum7
!
!write(file_id, *)'Distribution of female ability',dum2,dum3,dum4,dum5,dum6
!
!dum2=0d0
!dum3=0d0
!dum4=0d0
!dum5=0d0
!dum6=0d0
!dum7=0d0
!
!do i=1,T
!
!do it=1,nsim
!    dum7=dum7+1d0*WeightActive(i)
!    if(exp1(it,i,5)==1) then
!        dum2=dum2+1d0*WeightActive(i)
!    elseif(exp1(it,i,5)==2) then
!        dum3=dum3+1d0*WeightActive(i)
!    elseif(exp1(it,i,5)==3) then
!        dum4=dum4+1d0*WeightActive(i)
!    elseif(exp1(it,i,5)==4) then
!        dum5=dum5+1d0*WeightActive(i)
!    else
!        dum6=dum6+1d0*WeightActive(i)
!    end if
!end do
!
!end do
!
!dum2=dum2/dum7
!dum3=dum3/dum7
!dum4=dum4/dum7
!dum5=dum5/dum7
!dum6=dum6/dum7
!
!write(file_id, *)'Distribution of female idiosyncraic shock',dum2,dum3,dum4,dum5,dum6

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum6=0d0
dum7=0d0
dum8=0d0
dum9=0d0
dum10=0d0
dum11=0d0
dum12=0d0
dum13=0d0
dum14=0d0
dum15=0d0
dum16=0d0
dum17=0d0
dum18=0d0
dum19=0d0
dum20=0d0
dum21=0d0
dum22=0d0
dum23=0d0
dum24=0d0
dum25=0d0


do i=1,T

do it2=1,nsim2
do it=1,nsim

if(Sim1m(it2,it,i,10)<0.5) then
    dum7=dum7+1d0*WeightActive(i)
    if(exp1m(it2,it,i,6)==1) then
        dum2=dum2+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==2) then
        dum3=dum3+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==3) then
        dum4=dum4+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==4) then
        dum5=dum5+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==5) then
        dum6=dum6+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==6) then
        dum8=dum8+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==7) then
        dum9=dum9+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==8) then
        dum10=dum10+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==9) then
        dum11=dum11+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==10) then
        dum12=dum12+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==11) then
        dum13=dum13+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==12) then
        dum14=dum14+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==13) then
        dum15=dum15+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==14) then
        dum16=dum16+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==15) then
        dum17=dum17+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==16) then
        dum18=dum18+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==17) then
        dum19=dum19+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==18) then
        dum20=dum20+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==19) then
        dum21=dum21+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==20) then
        dum22=dum22+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==21) then
        dum23=dum23+1d0*WeightActive(i)
    elseif(exp1m(it2,it,i,6)==22) then
        dum24=dum24+1d0*WeightActive(i)
    else
        dum25=dum25+1d0*WeightActive(i)
    end if

end if

end do
end do

end do

dum2=dum2/dum7
dum3=dum3/dum7
dum4=dum4/dum7
dum5=dum5/dum7
dum6=dum6/dum7
dum8=dum8/dum7
dum9=dum9/dum7
dum10=dum10/dum7
dum11=dum11/dum7
dum12=dum12/dum7
dum13=dum13/dum7
dum14=dum14/dum7
dum15=dum15/dum7
dum16=dum16/dum7
dum17=dum17/dum7
dum18=dum18/dum7
dum19=dum19/dum7
dum20=dum20/dum7
dum21=dum21/dum7
dum22=dum22/dum7
dum23=dum23/dum7
dum24=dum24/dum7
dum25=dum25/dum7



write(file_id, *)'Distribution of single male fixed cost shock is',dum2,dum3,dum4,dum5,dum6,dum8,dum9,dum10,dum11,dum12,dum13,dum14,dum15,dum16,dum17,dum18,dum19,dum20,dum21,dum22,dum23,dum24,dum25

! Taxes by demographic group

!Singles

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum6=0d0
dum7=0d0


do i=1,T

do it2=1,nsim2 
do it=1,nsim
    
    if(Sim1m(it2,it,i,10)<0.5d0) then
        dum2=dum2+(Sim1m(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
        dum3=dum3+1d0*WeightActive(i)
    end if
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum4=dum4+(Sim1f(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
        dum5=dum5+1d0*WeightActive(i)
    end if
    
end do
end do

end do



do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    r_ret=r
    
    if(Sim1m(it2,it,T,10)<0.5d0) then
        dum2=dum2+(SimR1m(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
        dum3=dum3+1d0*WeightRet(i)
    end if
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum4=dum4+(SimR1f(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
        dum5=dum5+1d0*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3
dum4=dum4/dum5

write(file_id, *)'Capital tax for singles is',(dum2+dum4)/2d0
write(file_id, *)'Capital tax for single men is',dum2
write(file_id, *)'Capital tax for single women is',dum4

dum6=dum2
dum7=dum4


!Consumption tax singles

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1m(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1m(it2,it,i,8)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum4=dum4+Sim1f(it2,it,i,8)*WeightActive(i)
        dum5=dum5+1d0*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum2=dum2+SimR1m(it2,it,i,3)*WeightRet(i)
        dum3=dum3+1d0*WeightRet(i)
    end if
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum4=dum4+SimR1f(it2,it,i,3)*WeightRet(i)
        dum5=dum5+1d0*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3
dum4=dum4/dum5

write(file_id, *)'Consumption tax for singles is',(dum2+dum4)/2d0
write(file_id, *)'Consumption tax for single men is',dum2
write(file_id, *)'Consumption tax for single women is',dum4

dum6=dum6+dum2
dum7=dum7+dum4

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum10=0d0
dum11=0d0
dum12=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1m(it2,it,i,10)<0.5d0) then
        dum2=dum2+Sim1m(it2,it,i,7)*WeightActive(i)
        dum10=dum10+Sim1m(it2,it,i,6)*WeightActive(i)
        dum3=dum3+1d0*WeightActive(i)
    end if
    
    if(Sim1f(it2,it,i,10)<0.5d0) then
        dum4=dum4+Sim1f(it2,it,i,7)*WeightActive(i)
        dum11=dum11+Sim1f(it2,it,i,6)*WeightActive(i)
        dum5=dum5+1d0*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum3=dum3+1d0*WeightRet(i)
    end if
    
    if(Sim1f(it2,it,T,10)<0.5d0) then
        dum5=dum5+1d0*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3
dum4=dum4/dum5

dum10=dum10/dum3
dum11=dum11/dum5

write(file_id, *)'Labor income tax for singles is',(dum2+dum4)/2d0
write(file_id, *)'Labor income tax for single men is',dum2
write(file_id, *)'Labor income tax for single women is',dum4

write(file_id, *)'Labor income tax rate for singles is',(dum2+dum4)/(dum10+dum11)
write(file_id, *)'Labor income tax rate for single men is',dum2/dum10
write(file_id, *)'Labor income tax rate for single women is',dum4/dum11

write(file_id, *)'Total tax revenue for singles is',(dum2+dum4+dum6+dum7)/2d0
write(file_id, *)'Total tax revenueor single men is',dum2+dum6
write(file_id, *)'Total tax revenue for single women is',dum4+dum7

!Married

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0
dum6=0d0
dum7=0d0


do i=1,T

do it2=1,nsim2 
do it=1,nsim
    
    if(Sim1m(it2,it,i,10)>0.5d0) then
        dum2=dum2+(Sim1m(it2,it,i,1)+Gamma_redistr)*WeightActive(i)*r*tk
        dum3=dum3+2d0*WeightActive(i)
    end if
    
end do
end do

end do



do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    r_ret=r
    
    if(Sim1m(it2,it,T,10)>0.5d0) then
        dum2=dum2+(SimR1m(it2,it,i,1)+Gamma_redistr)*WeightRet(i)*r_ret*tk
        dum3=dum3+2d0*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Capital tax for married is',dum2

dum6=dum2


!Consumption tax married

dum2=0d0
dum3=0d0
dum4=0d0
dum5=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1m(it2,it,i,10)>0.5d0) then
        dum2=dum2+Sim1m(it2,it,i,8)*WeightActive(i)
        dum3=dum3+2d0*WeightActive(i)
    end if
    
end do
end do

end do


do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1f(it2,it,T,10)>0.5d0) then
        dum2=dum2+SimR1m(it2,it,i,3)*WeightRet(i)
        dum3=dum3+2d0*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Consumption tax for married is',dum2+dum4

dum6=dum6+dum2

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1m(it2,it,i,10)>0.5d0) then
        dum2=dum2+Sim1m(it2,it,i,7)*WeightActive(i)
        dum12=dum12+Sim1m(it2,it,i,6)*WeightActive(i)
        dum3=dum3+2d0*WeightActive(i)
    end if
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    
    if(Sim1f(it2,it,T,10)>0.5d0) then
        dum3=dum3+2d0*WeightRet(i)
    end if
    
end do
end do

end do

dum2=dum2/dum3
dum12=dum12/dum3

write(file_id, *)'Labor income tax for married is',dum2

write(file_id, *)'Labor income tax rate for married is',dum2/dum12

write(file_id, *)'Total tax revenue for married is',dum2+dum6

!Social Welfare

dum2=0d0
dum3=0d0

do i=1,T

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,11)*WeightActive(i)
    dum2=dum2+Sim1f(it2,it,i,11)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
    
end do
end do

end do

do i=1,Tret

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+SimR1m(it2,it,i,11)*WeightRet(i)
    dum2=dum2+SimR1f(it2,it,i,11)*WeightRet(i)
    dum3=dum3+2d0*WeightRet(i)
end do
end do

end do

dum2=dum2/dum3

write(file_id, *)'Social welfare of everyone is',dum2

dum2=0d0
dum3=0d0

i=1

do it2=1,nsim2
do it=1,nsim
    dum2=dum2+Sim1m(it2,it,i,11)*WeightActive(i)
    dum2=dum2+Sim1f(it2,it,i,11)*WeightActive(i)
    dum3=dum3+2d0*WeightActive(i)
    
end do
end do

dum2=dum2/dum3

write(file_id, *)'Social welfare of 20-year olds is',dum2

!STOP

end subroutine Statistics_to_file