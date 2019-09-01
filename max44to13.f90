subroutine max44to13(it)
!This subroutine computes optimal policies for people aged 20-63
use defineparameters
use glob0
USE CSVAL_INT
USE CSINT_INT
USE BS2IN_INT
USE BS3IN_INT
USE BSNAK_INT
USE BS2VL_INT
USE BS3VL_INT
USE QDVAL_INT
USE QD2VL_INT
USE QD3VL_INT
USE LCONF_INT

implicit none
integer, INTENT(IN) :: it
integer :: it2,it3,it4,it5,it6,it7,it8,it9,expcount,expcount2,temp,count,count2,count3,count4,exper,experf
integer :: d,d3,d4,d5,d6,d1,MAXFCN,NEQ=0
real(8) :: hours,dum2,dum3,L1,L2,L3,dum4,dum5,diff,ip,dum6, dum7, dum8, dum9,labor1,labor2,dum1,ACC=0.001,OBJ
real(8) :: X1,X2,X3,X4,V1,V2,V3,V4,Y1,Y2,Y3,Y4,p1,p2,p3,p4,B0,B1
real(8), dimension (1:40) :: KNOT1, KNOT2
real(8), dimension (1:640) :: KNOT3
real(8) :: SOL(1), XGUESS(1), XLB(1), XUB(1), GRAD(1,1), RHS(1)
real(8) :: SOL2(2), XGUESS2(2), XLB2(2), XUB2(2), GRAD2(1,2)
real(8), dimension (:,:,:,:,:), allocatable :: XSDATA,VSDATA,CSDATA,GKSDATA,GNSDATA
real(8), dimension (:,:,:,:,:,:,:), allocatable :: XMDATA,YMDATA,VMDATA,CMDATA,GKMDATA,GNMDATA,GNFMDATA
EXTERNAL labors, laborm

!j is gender, it2 is experience, it3 is education

d6=int((i-1)/4.0)+1
if((d6-1)*4+1==i) then
    d6=d6-1
end if
j=1
    do expcount=1,d6
        do it3=1,q
            do it5=1,l
                        it2=(expcount-1)*4+1
                        exper=it2-1
                        if(it3==1) then
                            wage=exp(gam0+gam1*exper+gam2*(exper**2)+gam3*(exper**3)+U(1,it3,it5))
                        elseif(it3==2) then
                            wage=exp(gamsc0+gamsc1*exper+gamsc2*(exper**2)+gamsc3*(exper**3)+U(1,it3,it5))
                        end if
                        if(it3==1) then
                            it4=2
                            dum8=(1.0-(1.0-marr)/0.9346)
                        else
                            it4=1
                            dum8=(1.0-(1.0-marr)/1.0654)
                        end if
                        B1=wage-wage*(t0+t1*(((wage)/AE)**0.2)+t2*(((wage)/AE)**0.4)+t3*(((wage)/AE)**0.6)+t4*(((wage)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/wage)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/(1.0+tc))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fsm
                            do it8=1,l
                                V2=V2+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2+1,it3,it5),CSCOEFS(j,i+1,:,it2+1,it3,it5,:))
                            end do
                            dum9=0.0
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it3,it8,it9,:))
                                            V2=V2+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it4,it8,it9,:))
                                            dum9=dum9+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/(1.0+tc))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fsm
                            do it8=1,l
                                V3=V3+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2+1,it3,it5),CSCOEFS(j,i+1,:,it2+1,it3,it5,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it3,it8,it9,:))
                                            V3=V3+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it4,it8,it9,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=SOL(1)
                        end if
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc))
                            do it8=1,l
                                V2=V2+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2,it3,it5),CSCOEFS(j,i+1,:,it2,it3,it5,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it3,it8,it9,:))
                                            V2=V2+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it4,it8,it9,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P3)/(1.0+tc))
                            do it8=1,l
                                V3=V3+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2,it3,it5),CSCOEFS(j,i+1,:,it2,it3,it5,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it3,it8,it9,:))
                                            V3=V3+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it4,it8,it9,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=0.0
                        end if
            end do
        end do
    end do

j=2

    do expcount=1,d6
        do it3=1,q
            do it5=1,l
                        it2=(expcount-1)*4+1
                        exper=it2-1
                        if(it3==1) then
                            wage=exp(gamf0+gamf1*exper+gamf2*(exper**2)+gamf3*(exper**3)+U(2,it3,it5))
                        elseif(it3==2) then
                            wage=exp(gamscf0+gamscf1*exper+gamscf2*(exper**2)+gamscf3*(exper**3)+U(2,it3,it5))
                        end if
                        if(it3==1) then
                            it4=2
                            dum8=(1.0-(1.0-marr)/0.9346)
                        else
                            it4=1
                            dum8=(1.0-(1.0-marr)/1.0654)
                        end if
                        B1=wage-wage*(t0+t1*(((wage)/AE)**0.2)+t2*(((wage)/AE)**0.4)+t3*(((wage)/AE)**0.6)+t4*(((wage)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/wage)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/(1.0+tc))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fsf
                            do it8=1,l
                                V2=V2+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2+1,it3,it8),CSCOEFS(j,i+1,:,it2+1,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it3,it3,it9,it8,:))
                                            V2=V2+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
                            CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/(1.0+tc))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fsf
                            do it8=1,l
                                V3=V3+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2+1,it3,it8),CSCOEFS(j,i+1,:,it2+1,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it3,it3,it9,it8,:))
                                            V3=V3+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=SOL(1)
                        end if
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc))
                            do it8=1,l
                                V2=V2+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2,it3,it8),CSCOEFS(j,i+1,:,it2,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it3,it3,it9,it8,:))
                                            V2=V2+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P3)/(1.0+tc))
                            do it8=1,l
                                V3=V3+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2,it3,it8),CSCOEFS(j,i+1,:,it2,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it3,it3,it9,it8,:))
                                            V3=V3+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=0.0
                        end if
            end do
        end do
    end do


    do expcount=1,d6
        do expcount2=1,d6
            do it4=1,q
                do it5=1,q
                    do it6=1,l
                    do it7=1,l
                        it2=(expcount-1)*4+1
                        it3=(expcount2-1)*4+1
                        exper=it2-1
                        experf=it3-1
                        if(it4==1) then
                            wage=exp(gam0+gam1*exper+gam2*(exper**2)+gam3*(exper**3)+U(1,it4,it6))
                        elseif(it4==2) then
                            wage=exp(gamsc0+gamsc1*exper+gamsc2*(exper**2)+gamsc3*(exper**3)+U(1,it4,it6))
                        end if
                        if(it5==1) then
                            wagef=exp(gamf0+gamf1*experf+gamf2*(experf**2)+gamf3*(experf**3)+U(2,it5,it7))
                        elseif(it5==2) then
                            wagef=exp(gamscf0+gamscf1*experf+gamscf2*(experf**2)+gamscf3*(experf**3)+U(2,it5,it7))
                        end if
                        dum3=wage+wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB2(1)=0.04*(AE/wage)
                        XUB2(1)=1.0
                        XLB2(2)=0.04*(AE/wagef)
                        XUB2(2)=1.0
                        GRAD2=-1.0
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(2)
                        end if
                        
                        dum3=wage
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                        
                        dum3=wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=2
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                        end if
                        
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                    end do
                    end do
                end do           
            end do
        end do
    end do
    
    j=1
        do it3=1,q
            do it5=1,l
                        it2=i
                        exper=it2-1
                        if(it3==1) then
                            wage=exp(gam0+gam1*exper+gam2*(exper**2)+gam3*(exper**3)+U(1,it3,it5))
                        elseif(it3==2) then
                            wage=exp(gamsc0+gamsc1*exper+gamsc2*(exper**2)+gamsc3*(exper**3)+U(1,it3,it5))
                        end if
                        if(it3==1) then
                            it4=2
                            dum8=(1.0-(1.0-marr)/0.9346)
                        else
                            it4=1
                            dum8=(1.0-(1.0-marr)/1.0654)
                        end if
                        B1=wage-wage*(t0+t1*(((wage)/AE)**0.2)+t2*(((wage)/AE)**0.4)+t3*(((wage)/AE)**0.6)+t4*(((wage)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/wage)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/(1.0+tc))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fsm
                            do it8=1,l
                                V2=V2+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2+1,it3,it5),CSCOEFS(j,i+1,:,it2+1,it3,it5,:))
                            end do
                            dum9=0.0
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it3,it8,it9,:))
                                            V2=V2+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it4,it8,it9,:))
                                            dum9=dum9+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/(1.0+tc))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fsm
                            do it8=1,l
                                V3=V3+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2+1,it3,it5),CSCOEFS(j,i+1,:,it2+1,it3,it5,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it3,it8,it9,:))
                                            V3=V3+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2+1,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2+1,it6,it3,it4,it8,it9,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=SOL(1)
                        end if
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc))
                            do it8=1,l
                                V2=V2+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2,it3,it5),CSCOEFS(j,i+1,:,it2,it3,it5,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it3,it8,it9,:))
                                            V2=V2+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it4,it8,it9,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P3)/(1.0+tc))
                            do it8=1,l
                                V3=V3+transu(1,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2,it3,it5),CSCOEFS(j,i+1,:,it2,it3,it5,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3,K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(1,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it3,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it3,it8,it9,:))
                                            V3=V3+transu(1,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,2)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it2,it6,it3,it4,it8,it9),CSCOEFM(i+1,:,it2,it6,it3,it4,it8,it9,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=0.0
                        end if
            end do
        end do

j=2

        do it3=1,q
            do it5=1,l
                        it2=i
                        exper=it2-1
                        if(it3==1) then
                            wage=exp(gamf0+gamf1*exper+gamf2*(exper**2)+gamf3*(exper**3)+U(2,it3,it5))
                        elseif(it3==2) then
                            wage=exp(gamscf0+gamscf1*exper+gamscf2*(exper**2)+gamscf3*(exper**3)+U(2,it3,it5))
                        end if
                        if(it3==1) then
                            it4=2
                            dum8=(1.0-(1.0-marr)/0.9346)
                        else
                            it4=1
                            dum8=(1.0-(1.0-marr)/1.0654)
                        end if
                        B1=wage-wage*(t0+t1*(((wage)/AE)**0.2)+t2*(((wage)/AE)**0.4)+t3*(((wage)/AE)**0.6)+t4*(((wage)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/wage)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/(1.0+tc))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fsf
                            do it8=1,l
                                V2=V2+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2+1,it3,it8),CSCOEFS(j,i+1,:,it2+1,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it3,it3,it9,it8,:))
                                            V2=V2+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
                            CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/(1.0+tc))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fsf
                            do it8=1,l
                                V3=V3+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2+1,it3,it8),CSCOEFS(j,i+1,:,it2+1,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it3,it3,it9,it8,:))
                                            V3=V3+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2+1,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2+1,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=SOL(1)
                        end if
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc))
                            do it8=1,l
                                V2=V2+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P2,BREAKS(j,i+1,:,it2,it3,it8),CSCOEFS(j,i+1,:,it2,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P2, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P2,K(n))
                                            end if
                                            V2=V2+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it3,it3,it9,it8,:))
                                            V2=V2+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+unempbenefit+lumpsum-P3)/(1.0+tc))
                            do it8=1,l
                                V3=V3+transu(2,it3,it5,it8)*(1.0-Probm(i))*beta*D_CSVAL(P3,BREAKS(j,i+1,:,it2,it3,it8),CSCOEFS(j,i+1,:,it2,it3,it8,:))
                            end do
                            do it8=1,l
                                do it6=1,i+1
                                    do it7=1,8
                                        do it9=1,l
                                            dum2=min(K(2*it7-1)+P3, K(n))
                                            if(it7==1) then
                                                dum2=min(K(2)*0.5+P3,K(n))
                                            end if
                                            V3=V3+transu(2,it3,it5,it8)*dum8*pprob(i,it6,it7,it3,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it3,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it3,it3,it9,it8,:))
                                            V3=V3+transu(2,it3,it5,it8)*(1.0-dum8)*pprob(i,it6,it7,it4,it9,1)*Probm(i)*beta*D_CSVAL(dum2,BREAKM(i+1,:,it6,it2,it4,it3,it9,it8),CSCOEFM(i+1,:,it6,it2,it4,it3,it9,it8,:))
                                        end do
                                    end do
                                end do
                            end do
                            if(K(it)*(1.0+r)+unempbenefit+lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3S(j,i,it,it2,it3,it5)) then
                            temp3S(j,i,it,it2,it3,it5)=V2
                            gkS(j,i,it,it2,it3,it5)=P2
                            cS(j,i,it,it2,it3,it5)=(K(it)*(1.0+r)+unempbenefit+lumpsum-P2)/(1.0+tc)
                            gnS(j,i,it,it2,it3,it5)=0.0
                        end if
            end do
        end do


        do expcount2=1,d6
            do it4=1,q
                do it5=1,q
                    do it6=1,l
                    do it7=1,l
                        it2=i
                        it3=(expcount2-1)*4+1
                        exper=it2-1
                        experf=it3-1
                        if(it4==1) then
                            wage=exp(gam0+gam1*exper+gam2*(exper**2)+gam3*(exper**3)+U(1,it4,it6))
                        elseif(it4==2) then
                            wage=exp(gamsc0+gamsc1*exper+gamsc2*(exper**2)+gamsc3*(exper**3)+U(1,it4,it6))
                        end if
                        if(it5==1) then
                            wagef=exp(gamf0+gamf1*experf+gamf2*(experf**2)+gamf3*(experf**3)+U(2,it5,it7))
                        elseif(it5==2) then
                            wagef=exp(gamscf0+gamscf1*experf+gamscf2*(experf**2)+gamscf3*(experf**3)+U(2,it5,it7))
                        end if
                        dum3=wage+wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB2(1)=0.04*(AE/wage)
                        XUB2(1)=1.0
                        XLB2(2)=0.04*(AE/wagef)
                        XUB2(2)=1.0
                        GRAD2=-1.0
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(2)
                        end if
                        
                        dum3=wage
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                        
                        dum3=wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=2
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                        end if
                        
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                    end do
                    end do
                end do           
            end do
        end do
    
    do expcount=1,d6
            do it4=1,q
                do it5=1,q
                    do it6=1,l
                    do it7=1,l
                        it2=(expcount-1)*4+1
                        it3=i
                        exper=it2-1
                        experf=it3-1
                        if(it4==1) then
                            wage=exp(gam0+gam1*exper+gam2*(exper**2)+gam3*(exper**3)+U(1,it4,it6))
                        elseif(it4==2) then
                            wage=exp(gamsc0+gamsc1*exper+gamsc2*(exper**2)+gamsc3*(exper**3)+U(1,it4,it6))
                        end if
                        if(it5==1) then
                            wagef=exp(gamf0+gamf1*experf+gamf2*(experf**2)+gamf3*(experf**3)+U(2,it5,it7))
                        elseif(it5==2) then
                            wagef=exp(gamscf0+gamscf1*experf+gamscf2*(experf**2)+gamscf3*(experf**3)+U(2,it5,it7))
                        end if
                        dum3=wage+wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB2(1)=0.04*(AE/wage)
                        XUB2(1)=1.0
                        XLB2(2)=0.04*(AE/wagef)
                        XUB2(2)=1.0
                        GRAD2=-1.0
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(2)
                        end if
                        
                        dum3=wage
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                        
                        dum3=wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=2
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                        end if
                        
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                    end do
                    end do
                end do           
            end do
    end do
    
            do it4=1,q
                do it5=1,q
                    do it6=1,l
                    do it7=1,l
                        it2=i
                        it3=i
                        exper=it2-1
                        experf=it3-1
                        if(it4==1) then
                            wage=exp(gam0+gam1*exper+gam2*(exper**2)+gam3*(exper**3)+U(1,it4,it6))
                        elseif(it4==2) then
                            wage=exp(gamsc0+gamsc1*exper+gamsc2*(exper**2)+gamsc3*(exper**3)+U(1,it4,it6))
                        end if
                        if(it5==1) then
                            wagef=exp(gamf0+gamf1*experf+gamf2*(experf**2)+gamf3*(experf**3)+U(2,it5,it7))
                        elseif(it5==2) then
                            wagef=exp(gamscf0+gamscf1*experf+gamscf2*(experf**2)+gamscf3*(experf**3)+U(2,it5,it7))
                        end if
                        dum3=wage+wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB2(1)=0.04*(AE/wage)
                        XUB2(1)=1.0
                        XLB2(2)=0.04*(AE/wagef)
                        XUB2(2)=1.0
                        GRAD2=-1.0
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc-P2)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/wage)-0.04*(AE/wagef)
                            XGUESS2=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(laborm, NEQ, GRAD2, RHS, XLB2, XUB2, SOL2, XGUESS=XGUESS2,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc-P3)/((1.0+tc)*1.7))-chim*(SOL2(1)**(1.0+etam))/(1.0+etam)-chif*(SOL2(2)**(1.0+etaf))/(1.0+etaf)-Fmm-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL2(2)
                        end if
                        
                        dum3=wage
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=1
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chim*(SOL(1)**(1.0+etam))/(1.0+etam)-Fmm
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2+1,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2+1,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2+1,it4,it8),CSCOEFS(1,i+1,:,it2+1,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                        
                        dum3=wagef
                        B1=dum3-dum3*(t0M+t1M*(((dum3)/AE)**0.2)+t2M*(((dum3)/AE)**0.4)+t3M*(((dum3)/AE)**0.6)+t4M*(((dum3)/AE)**0.8))
                        P1=K(1)
                        P4=min(K(n),K(it)*(1.0+r)+B1-0.05)
                        itg=it
                        XLB(1)=0.04*(AE/dum3)
                        XUB(1)=1.0
                        GRAD=-1.0
                        ind=2
                        j=2
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P2
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V2=log((K(it)*(1.0+r)+netinc+unempbenefit-P2)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P2<0.01) then
                                V2=-999999999.0
                            end if
                            RHS(1)=-0.04*(AE/dum3)
                            XGUESS(1)=0.9
		                    MAXFCN=1000000
		                    dum=P3
		                    CALL D_LCONF(labors, NEQ, GRAD, RHS, XLB, XUB, SOL, XGUESS=XGUESS,MAXFCN=MAXFCN, ACC=ACC, OBJ=OBJ)
                            V3=log((K(it)*(1.0+r)+netinc+unempbenefit-P3)/((1.0+tc)*1.7))-chif*(SOL(1)**(1.0+etaf))/(1.0+etaf)-Fmf
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3+1,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3+1,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3+1,it5,it9),CSCOEFS(2,i+1,:,it3+1,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+netinc+unempbenefit-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)-P2+netinc+unempbenefit)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=SOL(1)
                        end if
                        
                        P1=K(1)
                        P4=K(n)
                        do
                            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                            V2=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V2=V2+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P2,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P2*0.5
                            do it8=1,l
                                V2=V2+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V2=V2+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2<0.01) then
                                V2=-999999999.0
                            end if
                            V3=log((K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3)/((1.0+tc)*1.7))
                            do it8=1,l
                                do it9=1,l
                                    V3=V3+transu(1,it4,it6,it8)*transu(2,it5,it7,it9)*(1.0-Probd(i))*beta*D_CSVAL(P3,BREAKM(i+1,:,it2,it3,it4,it5,it8,it9),CSCOEFM(i+1,:,it2,it3,it4,it5,it8,it9,:))
                                end do
                            end do
                            dum2=P3*0.5
                            do it8=1,l
                                V3=V3+transu(1,it4,it6,it8)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(1,i+1,:,it2,it4,it8),CSCOEFS(1,i+1,:,it2,it4,it8,:))
                            end do
                            do it9=1,l
                                V3=V3+transu(2,it5,it7,it9)*Probd(i)*beta*D_CSVAL(dum2,BREAKS(2,i+1,:,it3,it5,it9),CSCOEFS(2,i+1,:,it3,it5,it9,:))
                            end do
                            if(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P3<0.01) then
                                V3=-999999999.0
                            end if        
                            if (V2 < V3) then
                                P1=P2
                            else
                                P4=P3
                            end if
                            if((P4-P1)<1d-3) exit
                        end do
                        if(V2>temp3M(i,it,it2,it3,it4,it5,it6,it7)) then
                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=V2
                            gkM(i,it,it2,it3,it4,it5,it6,it7)=P2
                            cM(i,it,it2,it3,it4,it5,it6,it7)=(K(it)*(1.0+r)+2.0*unempbenefit+2.0*lumpsum-P2)/(1.0+tc)
                            gnM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=0.0
                        end if
                    end do
                    end do
                end do           
            end do
            
            
            expcount=int((i-1)/4.0)+2
            if(((expcount-2)*4+1)==i) then
                expcount=expcount-1
            end if
            allocate(XSDATA(1:2,1:n,1:expcount,1:q,1:l))
            allocate(VSDATA(1:2,1:n,1:expcount,1:q,1:l))
            allocate(CSDATA(1:2,1:n,1:expcount,1:q,1:l))
            allocate(GKSDATA(1:2,1:n,1:expcount,1:q,1:l))
            allocate(GNSDATA(1:2,1:n,1:expcount,1:q,1:l))
            do j=1,2
                    do it2=1,expcount
                        do it3=1,q
                            do it4=1,l
                                if (it2<expcount) then
                                    XSDATA(j,it,it2,it3,it4)=((it2-1)*4+1)*1D0
                                    VSDATA(j,it,it2,it3,it4)=temp3S(j,i,it,(it2-1)*4+1,it3,it4)
                                    CSDATA(j,it,it2,it3,it4)=cS(j,i,it,(it2-1)*4+1,it3,it4)
                                    GKSDATA(j,it,it2,it3,it4)=gkS(j,i,it,(it2-1)*4+1,it3,it4)
                                    GNSDATA(j,it,it2,it3,it4)=gnS(j,i,it,(it2-1)*4+1,it3,it4)
                                else
                                    XSDATA(j,it,it2,it3,it4)=i*1D0
                                    VSDATA(j,it,it2,it3,it4)=temp3S(j,i,it,i,it3,it4)
                                    CSDATA(j,it,it2,it3,it4)=cS(j,i,it,i,it3,it4)
                                    GKSDATA(j,it,it2,it3,it4)=gkS(j,i,it,i,it3,it4)
                                    GNSDATA(j,it,it2,it3,it4)=gnS(j,i,it,i,it3,it4)
                                end if
                            end do
                        end do
                    end do
            end do
            
            !do j=1,2
            !        do it2=1,i
            !            do it3=1,q
            !                do it4=1,l
            !                    dum2=it2*1D0
            !                    temp3S(j,i,it,it2,it3,it4)=D_QDVAL(dum2,XSDATA(j,it,:,it3,it4),VSDATA(j,it,:,it3,it4))
            !                    cS(j,i,it,it2,it3,it4)=D_QDVAL(dum2,XSDATA(j,it,:,it3,it4),CSDATA(j,it,:,it3,it4))
            !                    gkS(j,i,it,it2,it3,it4)=D_QDVAL(dum2,XSDATA(j,it,:,it3,it4),GKSDATA(j,it,:,it3,it4))
            !                    gnS(j,i,it,it2,it3,it4)=D_QDVAL(dum2,XSDATA(j,it,:,it3,it4),GNSDATA(j,it,:,it3,it4))
            !                end do
            !            end do
            !        end do
            !end do
            
            do j=1,2
                    do it2=1,i
                        do it3=1,q
                            do it4=1,l
                                dum2=it2*1D0
                                d=0
                                do
                                    d=d+1
                                    if(XSDATA(j,it,d,it3,it4)>dum2) exit
                                    if(d==expcount) exit
                                end do
                                d=d-1
                                d=max(d,1)
                                p1=(dum2-XSDATA(j,it,d,it3,it4))/(XSDATA(j,it,d+1,it3,it4)-XSDATA(j,it,d,it3,it4))
                                temp3S(j,i,it,it2,it3,it4)=VSDATA(j,it,d,it3,it4)+p1*(VSDATA(j,it,d+1,it3,it4)-VSDATA(j,it,d,it3,it4))
                                cS(j,i,it,it2,it3,it4)=CSDATA(j,it,d,it3,it4)+p1*(CSDATA(j,it,d+1,it3,it4)-CSDATA(j,it,d,it3,it4))
                                gkS(j,i,it,it2,it3,it4)=GKSDATA(j,it,d,it3,it4)+p1*(GKSDATA(j,it,d+1,it3,it4)-GKSDATA(j,it,d,it3,it4))
                                gnS(j,i,it,it2,it3,it4)=GNSDATA(j,it,d,it3,it4)+p1*(GNSDATA(j,it,d+1,it3,it4)-GNSDATA(j,it,d,it3,it4))
                            end do
                        end do
                    end do
            end do
            
            allocate(XMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            allocate(YMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            allocate(VMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            allocate(CMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            allocate(GKMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            allocate(GNMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            allocate(GNFMDATA(1:n,1:expcount,1:expcount,1:q,1:q,1:l,1:l))
            
                do it2=1,expcount
                    do it3=1,expcount
                        do it4=1,q
                            do it5=1,q
                                do it6=1,l
                                    do it7=1,l
                                        if ((it2<expcount).AND.(it3<expcount)) then
                                            XMDATA(it,it2,it3,it4,it5,it6,it7)=((it2-1)*4+1)*1D0
                                            YMDATA(it,it2,it3,it4,it5,it6,it7)=((it3-1)*4+1)*1D0
                                            VMDATA(it,it2,it3,it4,it5,it6,it7)=temp3M(i,it,(it2-1)*4+1,(it3-1)*4+1,it4,it5,it6,it7)
                                            CMDATA(it,it2,it3,it4,it5,it6,it7)=cM(i,it,(it2-1)*4+1,(it3-1)*4+1,it4,it5,it6,it7)
                                            GKMDATA(it,it2,it3,it4,it5,it6,it7)=gkM(i,it,(it2-1)*4+1,(it3-1)*4+1,it4,it5,it6,it7)
                                            GNMDATA(it,it2,it3,it4,it5,it6,it7)=gnM(i,it,(it2-1)*4+1,(it3-1)*4+1,it4,it5,it6,it7)
                                            GNFMDATA(it,it2,it3,it4,it5,it6,it7)=gnfM(i,it,(it2-1)*4+1,(it3-1)*4+1,it4,it5,it6,it7)
                                        elseif((it2==expcount).AND.(it3<expcount)) then
                                            XMDATA(it,it2,it3,it4,it5,it6,it7)=i*1D0
                                            YMDATA(it,it2,it3,it4,it5,it6,it7)=((it3-1)*4+1)*1D0
                                            VMDATA(it,it2,it3,it4,it5,it6,it7)=temp3M(i,it,i,(it3-1)*4+1,it4,it5,it6,it7)
                                            CMDATA(it,it2,it3,it4,it5,it6,it7)=cM(i,it,i,(it3-1)*4+1,it4,it5,it6,it7)
                                            GKMDATA(it,it2,it3,it4,it5,it6,it7)=gkM(i,it,i,(it3-1)*4+1,it4,it5,it6,it7)
                                            GNMDATA(it,it2,it3,it4,it5,it6,it7)=gnM(i,it,i,(it3-1)*4+1,it4,it5,it6,it7)
                                            GNFMDATA(it,it2,it3,it4,it5,it6,it7)=gnfM(i,it,i,(it3-1)*4+1,it4,it5,it6,it7)
                                        elseif ((it2<expcount).AND.(it3==expcount)) then
                                            XMDATA(it,it2,it3,it4,it5,it6,it7)=((it2-1)*4+1)*1D0
                                            YMDATA(it,it2,it3,it4,it5,it6,it7)=i*1D0
                                            VMDATA(it,it2,it3,it4,it5,it6,it7)=temp3M(i,it,(it2-1)*4+1,i,it4,it5,it6,it7)
                                            CMDATA(it,it2,it3,it4,it5,it6,it7)=cM(i,it,(it2-1)*4+1,i,it4,it5,it6,it7)
                                            GKMDATA(it,it2,it3,it4,it5,it6,it7)=gkM(i,it,(it2-1)*4+1,i,it4,it5,it6,it7)
                                            GNMDATA(it,it2,it3,it4,it5,it6,it7)=gnM(i,it,(it2-1)*4+1,i,it4,it5,it6,it7)
                                            GNFMDATA(it,it2,it3,it4,it5,it6,it7)=gnfM(i,it,(it2-1)*4+1,i,it4,it5,it6,it7)
                                        elseif ((it2==expcount).AND.(it3==expcount)) then
                                            XMDATA(it,it2,it3,it4,it5,it6,it7)=i*1D0
                                            YMDATA(it,it2,it3,it4,it5,it6,it7)=i*1D0
                                            VMDATA(it,it2,it3,it4,it5,it6,it7)=temp3M(i,it,i,i,it4,it5,it6,it7)
                                            CMDATA(it,it2,it3,it4,it5,it6,it7)=cM(i,it,i,i,it4,it5,it6,it7)
                                            GKMDATA(it,it2,it3,it4,it5,it6,it7)=gkM(i,it,i,i,it4,it5,it6,it7)
                                            GNMDATA(it,it2,it3,it4,it5,it6,it7)=gnM(i,it,i,i,it4,it5,it6,it7)
                                            GNFMDATA(it,it2,it3,it4,it5,it6,it7)=gnfM(i,it,i,i,it4,it5,it6,it7)
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            
                !do it2=1,i
                !    do it3=1,i
                !        do it4=1,q
                !            do it5=1,q
                !                do it6=1,l
                !                    do it7=1,l
                !                        dum2=it2*1D0
                !                        dum3=it3*1D0
                !                        temp3M(i,it,it2,it3,it4,it5,it6,it7)=D_QD2VL(dum2,dum3,XMDATA(it,:,1,it4,it5,it6,it7),YMDATA(it,1,:,it4,it5,it6,it7),VMDATA(it,:,:,it4,it5,it6,it7))
                !                        cM(i,it,it2,it3,it4,it5,it6,it7)=D_QD2VL(dum2,dum3,XMDATA(it,:,1,it4,it5,it6,it7),YMDATA(it,1,:,it4,it5,it6,it7),CMDATA(it,:,:,it4,it5,it6,it7))
                !                        gkM(i,it,it2,it3,it4,it5,it6,it7)=D_QD2VL(dum2,dum3,XMDATA(it,:,1,it4,it5,it6,it7),YMDATA(it,1,:,it4,it5,it6,it7),GKMDATA(it,:,:,it4,it5,it6,it7))
                !                        gnM(i,it,it2,it3,it4,it5,it6,it7)=D_QD2VL(dum2,dum3,XMDATA(it,:,1,it4,it5,it6,it7),YMDATA(it,1,:,it4,it5,it6,it7),GNMDATA(it,:,:,it4,it5,it6,it7))
                !                        gnfM(i,it,it2,it3,it4,it5,it6,it7)=D_QD2VL(dum2,dum3,XMDATA(it,:,1,it4,it5,it6,it7),YMDATA(it,1,:,it4,it5,it6,it7),GNFMDATA(it,:,:,it4,it5,it6,it7))
                !                    end do
                !                end do
                !            end do
                !        end do
                !    end do
                !end do
                
                do it2=1,i
                    do it3=1,i
                        do it4=1,q
                            do it5=1,q
                                do it6=1,l
                                    do it7=1,l
                                        dum2=it2*1D0
                                        dum3=it3*1D0
                                        
                                        d=0
                                        do
                                            d=d+1
                                            if(XMDATA(it,d,1,it4,it5,it6,it7)>dum2) exit
                                            if(d==expcount) exit
                                        end do
                                        d=d-1
                                        d=max(d,1)
                                        p1=(dum2-XMDATA(it,d,1,it4,it5,it6,it7))/(XMDATA(it,d+1,1,it4,it5,it6,it7)-XMDATA(it,d,1,it4,it5,it6,it7))
                                        
                                        d2=0
                                        do
                                            d2=d2+1
                                            if(YMDATA(it,1,d2,it4,it5,it6,it7)>dum3) exit
                                            if(d2==expcount) exit
                                        end do
                                        d2=d2-1
                                        d2=max(d2,1)
                                        p2=(dum3-YMDATA(it,1,d2,it4,it5,it6,it7))/(YMDATA(it,1,d2+1,it4,it5,it6,it7)-YMDATA(it,1,d2,it4,it5,it6,it7))
                                        
                                        if(p1+p2<1) then
                                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=VMDATA(it,d,d2,it4,it5,it6,it7)*(1.0-(p1+p2))
                                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=temp3M(i,it,it2,it3,it4,it5,it6,it7)+VMDATA(it,d+1,d2,it4,it5,it6,it7)*p1+VMDATA(it,d,d2+1,it4,it5,it6,it7)*p2
                                            
                                            cM(i,it,it2,it3,it4,it5,it6,it7)=CMDATA(it,d,d2,it4,it5,it6,it7)*(1.0-(p1+p2))
                                            cM(i,it,it2,it3,it4,it5,it6,it7)=cM(i,it,it2,it3,it4,it5,it6,it7)+CMDATA(it,d+1,d2,it4,it5,it6,it7)*p1+CMDATA(it,d,d2+1,it4,it5,it6,it7)*p2
                                            
                                            gkM(i,it,it2,it3,it4,it5,it6,it7)=GKMDATA(it,d,d2,it4,it5,it6,it7)*(1.0-(p1+p2))
                                            gkM(i,it,it2,it3,it4,it5,it6,it7)=gkM(i,it,it2,it3,it4,it5,it6,it7)+GKMDATA(it,d+1,d2,it4,it5,it6,it7)*p1+GKMDATA(it,d,d2+1,it4,it5,it6,it7)*p2
                                            
                                            gnM(i,it,it2,it3,it4,it5,it6,it7)=GNMDATA(it,d,d2,it4,it5,it6,it7)*(1.0-(p1+p2))
                                            gnM(i,it,it2,it3,it4,it5,it6,it7)=gnM(i,it,it2,it3,it4,it5,it6,it7)+GNMDATA(it,d+1,d2,it4,it5,it6,it7)*p1+GNMDATA(it,d,d2+1,it4,it5,it6,it7)*p2
                                            
                                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=GNFMDATA(it,d,d2,it4,it5,it6,it7)*(1.0-(p1+p2))
                                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=gnfM(i,it,it2,it3,it4,it5,it6,it7)+GNFMDATA(it,d+1,d2,it4,it5,it6,it7)*p1+GNFMDATA(it,d,d2+1,it4,it5,it6,it7)*p2
                                        else
                                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=VMDATA(it,d+1,d2+1,it4,it5,it6,it7)*(p1+p2-1.0)
                                            temp3M(i,it,it2,it3,it4,it5,it6,it7)=temp3M(i,it,it2,it3,it4,it5,it6,it7)+VMDATA(it,d+1,d2,it4,it5,it6,it7)*(1.0-p2)+VMDATA(it,d,d2+1,it4,it5,it6,it7)*(1.0-p1)
                                            
                                            cM(i,it,it2,it3,it4,it5,it6,it7)=CMDATA(it,d+1,d2+1,it4,it5,it6,it7)*(p1+p2-1.0)
                                            cM(i,it,it2,it3,it4,it5,it6,it7)=cM(i,it,it2,it3,it4,it5,it6,it7)+CMDATA(it,d+1,d2,it4,it5,it6,it7)*(1.0-p2)+CMDATA(it,d,d2+1,it4,it5,it6,it7)*(1.0-p1)
                                            
                                            gkM(i,it,it2,it3,it4,it5,it6,it7)=GKMDATA(it,d+1,d2+1,it4,it5,it6,it7)*(p1+p2-1.0)
                                            gkM(i,it,it2,it3,it4,it5,it6,it7)=gkM(i,it,it2,it3,it4,it5,it6,it7)+GKMDATA(it,d+1,d2,it4,it5,it6,it7)*(1.0-p2)+GKMDATA(it,d,d2+1,it4,it5,it6,it7)*(1.0-p1)
                                            
                                            gnM(i,it,it2,it3,it4,it5,it6,it7)=GNMDATA(it,d+1,d2+1,it4,it5,it6,it7)*(p1+p2-1.0)
                                            gnM(i,it,it2,it3,it4,it5,it6,it7)=gnM(i,it,it2,it3,it4,it5,it6,it7)+GNMDATA(it,d+1,d2,it4,it5,it6,it7)*(1.0-p2)+GNMDATA(it,d,d2+1,it4,it5,it6,it7)*(1.0-p1)
                                            
                                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=GNFMDATA(it,d+1,d2+1,it4,it5,it6,it7)*(p1+p2-1.0)
                                            gnfM(i,it,it2,it3,it4,it5,it6,it7)=gnfM(i,it,it2,it3,it4,it5,it6,it7)+GNFMDATA(it,d+1,d2,it4,it5,it6,it7)*(1.0-p2)+GNFMDATA(it,d,d2+1,it4,it5,it6,it7)*(1.0-p1)
                                        end if
                                        
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do

end subroutine max44to13