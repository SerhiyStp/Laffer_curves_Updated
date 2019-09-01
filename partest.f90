subroutine partest(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    USE QDVAL_INT
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

!Assigning the grid points
dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
ik=int(dum3)+1
dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
ium=int(dum3)+1
iam=counter-(ik-1)*na*nu-(ium-1)*nu  
    
    j=2
        do ifc=1,nfc
            do ix = 1, nexp
                    do iu2=1,nu
                        do ik2=1,nk
                            do iaf=1,na
                                do iuf=1,nu
                                        dum=k_grid(ik)+k_grid(ik2)
                                        if(dum<k_grid(nk)-0.001d0) then
                                            !evm(j,ik,ix,iam,ium,T-it)=evm(j,ik,ix,iam,ium,T-it)+trans_u(2,ium,iu2)*mpartner(ik2,iaf,iuf,T-it)*Probfam(j2,T-it)*D_CSVAL(dum,BREAK,v_spln_coefs(:,:,ix,iaf,iuf,iam,iu2,T-it,j2))
                                            vnext = D_QDVAL(dum, k_grid, v(:,ix,iaf,iuf,iam,iu2,T-it,ifc))
                                            evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,iaf,iuf,T-it)*vnext
                                        else
                                            evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,iaf,iuf,T-it)*LinInterp(dum,k_grid,v(:,ix,iaf,iuf,iam,iu2,T-it,ifc),nk)
                                        end if
                                end do
                            end do
                        end do
                    end do
            end do
        end do

end subroutine partest
    
subroutine partest2(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    USE QDVAL_INT
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

!Assigning the grid points
dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
ik=int(dum3)+1
dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
ium=int(dum3)+1
iam=counter-(ik-1)*na*nu-(ium-1)*nu

    
    j=1
            do iu2=1,nu
                do ik2=1,nk
                    do ix=1,nexp
                        do iaf=1,na
                            do iuf=1,nu
                                    do ifc=1,nfc
                                        dum=k_grid(ik)+k_grid(ik2)
                                        if(dum<k_grid(nk)-0.001d0) then
                                            !evm(j,ik,:,iam,ium,T-it)=evm(j,ik,:,iam,ium,T-it)+trans_u(1,ium,iu2)*fpartner(ik2,ix,iaf,iuf,T-it)*Probfam(j2,T-it)*D_CSVAL(dum,BREAK,v_spln_coefs(:,:,ix,iam,iu2,iaf,iuf,T-it,j2))
                                            vnext = D_QDVAL(dum, k_grid, v(:,ix,iam,iu2,iaf,iuf,T-it,ifc))
                                            evm(j,ik,:,iam,ium,T-it,ifc)=evm(j,ik,:,iam,ium,T-it,ifc)+trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext
                                        else
                                            evm(j,ik,:,iam,ium,T-it,ifc)=evm(j,ik,:,iam,ium,T-it,ifc)+trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*LinInterp(dum,k_grid,v(:,ix,iam,iu2,iaf,iuf,T-it,ifc),nk)
                                        end if
                                    end do
                            end do
                        end do
                    end do
                end do
            end do


end subroutine partest2
    
subroutine partest3(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    USE BS2IN_INT
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE CSINT_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ium
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    
!Assigning the grid points
dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
ium=int(dum3)+1
dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
ifc=int(dum3)+1
iam=counter-(ium-1)*nfc*na-(ifc-1)*na

            
        do iaf=1,na
            do iuf=1,nu
                !call d_csint(k_grid, ev(:,ix,iam,ium,iaf,iuf,T-it,j2), BREAK, ev_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T-it,j2))
                CALL D_BS2IN(k_grid, exp_grid(:,T-it), ev(:,:,iam,ium,iaf,iuf,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), ev_spln_coefs(:,:,iam,ium,iaf,iuf,T-it,ifc) , nk)
                CALL D_BS2IN(k_grid, exp_grid(:,T-it), v(:,:,iam,ium,iaf,iuf,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), v_spln_coefs(:,:,iam,ium,iaf,iuf,T-it,ifc) , nk)
            end do
        end do
            

end subroutine partest3

subroutine partest4(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE CSINT_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ium
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    
!Assigning the grid points
dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
ium=int(dum3)+1
dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
ifc=int(dum3)+1
iam=counter-(ium-1)*nfc*na-(ifc-1)*na

    
            do ix = 1, nexp
                    do iaf=1,na
                        do iuf=1,nu
                            !call d_csint(k_grid, v(:,ix,iam,ium,iaf,iuf,T-it,j2), BREAK, v_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T-it,j2))
                            !call ppp_csint(k_grid, v(:,ix,iam,ium,iaf,iuf,T-it,ifc), v_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T-it,ifc), nk)
                        end do
                    end do
            end do
            

end subroutine partest4
    
subroutine partest5(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    USE BS2IN_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ium
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    
!Assigning the grid points
dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
ium=int(dum3)+1
dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
ifc=int(dum3)+1
iam=counter-(ium-1)*nfc*na-(ifc-1)*na
    
    do j=1,2
        !call d_csint(k_grid, evs(j,:,ix,iam,ium,T-it), BREAK, evs_spln_coefs(j,:,:,ix,iam,ium,T-it))
        CALL D_BS2IN(k_grid, exp_grid(:,T-it), evs(j,:,:,iam,ium,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), evs_spln_coefs(j,:,:,iam,ium,T-it,ifc), nk)
    end do
            
    do j=1,2
        !call d_csint(k_grid, vs(j,:,ix,iam,ium,T-it), BREAK, vs_spln_coefs(j,:,:,ix,iam,ium,T-it))
        CALL D_BS2IN(k_grid, exp_grid(:,T-it), vs(j,:,:,iam,ium,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), vs_spln_coefs(j,:,:,iam,ium,T-it,ifc), nk)
    end do

end subroutine partest5   
    
subroutine partest6(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    USE BS2IN_INT
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE CSINT_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ium
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2

!Assigning the grid points
dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
ium=int(dum3)+1
dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
ifc=int(dum3)+1
iam=counter-(ium-1)*nfc*na-(ifc-1)*na

    do j=1,2
        !call d_csint(k_grid, evm(j,:,ix,iam,ium,T-it), BREAK, evm_spln_coefs(j,:,:,ix,iam,ium,T-it))
        CALL D_BS2IN(k_grid, exp_grid(:,T-it), evm(j,:,:,iam,ium,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), evm_spln_coefs(j,:,:,iam,ium,T-it,ifc), nk)
    end do

end subroutine partest6
    
subroutine partest7(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE CSINT_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT
    USE QD2VL_INT
    USE QDVAL_INT
    USE BS2VL_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

!Assigning the grid points
dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
ik=int(dum3)+1
dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
ium=int(dum3)+1
iam=counter-(ik-1)*na*nu-(ium-1)*nu
    
    if(T-it>1) then

        ev(ik,:,iam,ium,:,:,T-it,:)=0d0
        dum3=K_grid(ik)/2d0
        
            do ifc=1,nfc
                    do ix = 1, nexp
                            do iuf = 1, nu
                                do iaf = 1, na

                                do iu2=1,nu
                                    do iu3=1,nu
                                            ev(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ev(ik,ix,iam,ium,iaf,iuf,T-it,ifc)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V(ik,ix,iam,iu2,iaf,iu3,T-it,ifc)
                                    end do
                                end do

                                do iu3=1,nu
                                    !ev(ik,ix,iam,ium,iaf,iuf,T-it,j)=ev(ik,ix,iam,ium,iaf,iuf,T-it,j)+Probd(T-it-1)*trans_u(2,iuf,iu3)*0.5d0*D_CSVAL(dum3,BREAK,vs_spln_coefs(2,:,:,ix,iaf,iu3,T-it))
                                    !vnext = D_QDVAL(dum3,k_grid,vs(2,:,ix,iaf,iu3,T-it,ifc))
                                    vnext = D_BS2VL(dum3,exp_grid(ix,T-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it), nk, nexp, vs_spln_coefs(2,:,:,iaf,iu3,T-it,ifc))
                                    ev(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ev(ik,ix,iam,ium,iaf,iuf,T-it,ifc)+Probd(T-it-1)*trans_u(2,iuf,iu3)*0.5d0*vnext
                                end do

                                do iu2=1,nu
                                    !ev(ik,ix,iam,ium,iaf,iuf,T-it,j)=ev(ik,ix,iam,ium,iaf,iuf,T-it,j)+Probd(T-it-1)*trans_u(1,ium,iu2)*0.5d0*D_CSVAL(dum3,BREAK,vs_spln_coefs(1,:,:,1,iam,iu2,T-it))
                                    !vnext = D_QDVAL(dum3,k_grid,vs(1,:,1,iam,iu2,T-it,1))
                                    vnext = D_BS2VL(dum3,exp_grid(ix,T-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it), nk, nexp, vs_spln_coefs(1,:,:,iam,iu2,T-it,1))
                                    ev(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ev(ik,ix,iam,ium,iaf,iuf,T-it,ifc)+Probd(T-it-1)*trans_u(1,ium,iu2)*0.5d0*vnext
                                end do

                                end do
                            end do
                    end do
            end do

    end if
            
end subroutine partest7
   
    
subroutine partest8(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE CSINT_INT
    !USE NEQNF_INT
    !USE NEQBF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

!Assigning the grid points
dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
ik=int(dum3)+1
dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
ium=int(dum3)+1
iam=counter-(ik-1)*na*nu-(ium-1)*nu
    
evs(:,ik,:,iam,ium,T-it,:)=0d0
            
    j=2
    
        do ifc=1,nfc
            do ix = 1, nexp
                    do iu2=1,nu
                        evs(j,ik,ix,iam,ium,T-it,ifc)=evs(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*Vs(j,ik,ix,iam,iu2,T-it,ifc)
                    end do
            end do
        end do
        
            
    j=1
            do iu2=1,nu    
                evs(j,ik,:,iam,ium,T-it,:)=evs(j,ik,:,iam,ium,T-it,:)+trans_u(1,ium,iu2)*Vs(j,ik,1,iam,iu2,T-it,1)
            end do
            
            
            
end subroutine partest8