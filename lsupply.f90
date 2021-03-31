subroutine lsupply(ik)
    !This subroutine computes optimal policies durin active worklife.
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE NEQNF_INT
    !USE ERSET_INT
    use hybrd_wrapper, only: getHybrSoln

    implicit none

    integer, INTENT(IN) :: ik
    integer :: ix,ixd,ixm,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2

    real(8) :: sol_test(2)
    real(8) :: xguess_loc(2), sol_loc(2)
    real(8) :: xguess_loc_1(1), sol_loc_1(1)


    EXTERNAL labor2
    EXTERNAL labor1
    EXTERNAL labor3
    EXTERNAL labors
    EXTERNAL labor2_hybrd
    EXTERNAL labor3_hybrd
    EXTERNAL labor1_hybrd
    EXTERNAL labors_hybrd

    !CALL ERSET(IERSVR, IPACT, ISACT)

    dum=c_grid(ik)

    ! Married couples:
        do ium=1,nw
            do iuf=1,nw
                wagem=wage_grid(ium)
                wagef=wage_grid(iuf)
                XGUESS_loc(1)=0.99d0
                XGUESS_loc(2)=0.99d0
                ITMAX=1000000
                !For each level of capital we call D_LCONF to find optimal labor
                !CALL D_NEQNF(labor2, SOL_loc, ERREL,ITMAX=ITMAX,XGUESS=XGUESS_loc,FNORM=FNORM)
                call getHybrSoln(xguess_loc, 2, labor2_hybrd, sol_loc, fnorm)

                dum4=sol_loc(1)
                dum5=sol_loc(2)
                dum6=FNORM
                if(dum6>10d0*(ERREL**2)) then
                    XGUESS2(1)=0.01d0
                    XGUESS2(2)=0.01d0
                    ITMAX=1000000
                    !CALL D_NEQNF(labor2, sol_loc, ERREL,ITMAX=ITMAX,XGUESS=XGUESS2,FNORM=FNORM)
                    call getHybrSoln(xguess_loc, 2, labor2_hybrd, sol_loc, fnorm)
                    if(FNORM<dum6) then
                        dum4=sol_loc(1)
                        dum5=sol_loc(2)
                        dum6=FNORM
                    end if
                end if
                if(dum6>10d0*(ERREL**2)) then
                    XGUESS2(1)=0.3d0
                    XGUESS2(2)=0.3d0
                    ITMAX=1000000
                    !CALL D_NEQNF(labor2, sol_loc, ERREL,ITMAX=ITMAX,XGUESS=XGUESS2,FNORM=FNORM)
                    call getHybrSoln(xguess_loc, 2, labor2_hybrd, sol_loc, fnorm)
                    if(FNORM<dum6) then
                        dum4=sol_loc(1)
                        dum5=sol_loc(2)
                        dum6=FNORM
                    end if
                end if

                if((dum4>1d0).AND.(dum5>1d0)) then
                    dum4=1d0
                    dum5=1d0
                elseif((dum4>1d0).AND.(dum5<1d0)) then
                    XGUESS_loc_1(1)=0.99d0
                    ITMAX=1000000
                    ind=1
                    !CALL D_NEQNF(labor3, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                    call getHybrSoln(xguess_loc_1, 1, labor3_hybrd, sol_loc_1, fnorm)
                    dum4=1d0
                    dum5=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                    if(dum6>10d0*(ERREL**2)) then
                        xguess_loc_1(1)=0.01d0
                        ITMAX=1000000
                        !CALL D_NEQNF(labor3, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                        call getHybrSoln(xguess_loc_1, 1, labor3_hybrd, sol_loc_1, fnorm)
                        if(FNORM<dum6) then
                            dum5=min(sol_loc_1(1),1d0)
                            dum6=FNORM
                        end if
                    end if
                    if(dum6>10d0*(ERREL**2)) then
                        xguess_loc_1(1)=0.3d0
                        ITMAX=1000000
                        !CALL D_NEQNF(labor3, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                        call getHybrSoln(xguess_loc_1, 1, labor3_hybrd, sol_loc_1, fnorm)
                        if(FNORM<dum6) then
                            dum5=min(sol_loc_1(1),1d0)
                            dum6=FNORM
                        end if
                        !if(dum6>10d0*(ERREL**2)) then
                        !    Print *,'ik is',ik
                        !    Print *,'iam is',iam
                        !    Print *,'ium is',ium
                        !    Print *,'iaf is',iaf
                        !    Print *,'iuf is',iuf
                        !    Print *,'dum6 is',dum6
                        !    Print *,'ERREL is',ERREL
                        !end if
                    end if
                elseif((dum4<1d0).AND.(dum5>1d0)) then
                    xguess_loc_1(1)=0.99d0
                    ITMAX=1000000
                    ind=2
                    !CALL D_NEQNF(labor3, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                    call getHybrSoln(xguess_loc_1, 1, labor3_hybrd, sol_loc_1, fnorm)
                    dum5=1d0
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                    if(dum6>10d0*(ERREL**2)) then
                        xguess_loc_1(1)=0.01d0
                        ITMAX=1000000
                        !CALL D_NEQNF(labor3, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                        call getHybrSoln(xguess_loc_1, 1, labor3_hybrd, sol_loc_1, fnorm)
                        if(FNORM<dum6) then
                            dum4=min(sol_loc_1(1),1d0)
                            dum6=FNORM
                        end if
                    end if
                    if(dum6>10d0*(ERREL**2)) then
                        xguess_loc_1(1)=0.3d0
                        ITMAX=1000000
                        !CALL D_NEQNF(labor3, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                        call getHybrSoln(xguess_loc_1, 1, labor3_hybrd, sol_loc_1, fnorm)
                        if(FNORM<dum6) then
                            dum4=min(sol_loc_1(1),1d0)
                            dum6=FNORM
                        end if
                        !if(dum6>10d0*(ERREL**2)) then
                        !    Print *,'ik is',ik
                        !    Print *,'iam is',iam
                        !    Print *,'ium is',ium
                        !    Print *,'iaf is',iaf
                        !    Print *,'iuf is',iuf
                        !    Print *,'dum6 is',dum6
                        !    Print *,'ERREL is',ERREL
                        !end if
                    end if
                end if

                laborm(ik,ium,iuf)=dum4
                laborf(ik,ium,iuf)=dum5
            end do
        end do

        !Male works
        
        ind2=1
        
        do ium=1,nw
            wagem=wage_grid(ium)
            xguess_loc_1(1)=0.99d0
            ITMAX=1000000
            !CALL D_NEQNF(labor1, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS)
            call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
            dum6=FNORM
            dum4=min(sol_loc_1(1),1d0)
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.01d0
                ITMAX=1000000
                !CALL D_NEQNF(labor1, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.3d0
                ITMAX=1000000
                !CALL D_NEQNF(labor1, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
                !if(dum6>10d0*(ERREL**2)) then
                !    Print *,'ik is',ik
                !    Print *,'iam is',iam
                !    Print *,'ium is',ium
                !    Print *,'iaf is',iaf
                !    Print *,'iuf is',iuf
                !    Print *,'dum6 is',dum6
                !    Print *,'ERREL is',ERREL
                !end if
            end if
            labormwork(ik,ium)=dum4
        end do

        
     !Female works
     
     ind2=2
        
        do ium=1,nw
            wagem=wage_grid(ium)
            xguess_loc_1(1)=0.99d0
            ITMAX=1000000
            !CALL D_NEQNF(labor1, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS)
            call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
            dum6=FNORM
            dum4=min(sol_loc_1(1),1d0)
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.01d0
                ITMAX=1000000
                !CALL D_NEQNF(labor1, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.3d0
                ITMAX=1000000
                !CALL D_NEQNF(labor1, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
                !if(dum6>10d0*(ERREL**2)) then
                !    Print *,'ik is',ik
                !    Print *,'iam is',iam
                !    Print *,'ium is',ium
                !    Print *,'iaf is',iaf
                !    Print *,'iuf is',iuf
                !    Print *,'dum6 is',dum6
                !    Print *,'ERREL is',ERREL
                !end if
            end if
            laborfwork(ik,ium)=dum4
        end do   
        
    ! Singles:

    ind2=1

    do ium=1,nw
        wagem=wage_grid(ium)
        xguess_loc_1(1)=0.99d0
        ITMAX=1000000
        !For each level of capital we call D_LCONF to find optimal labor
        !CALL D_NEQNF(labors, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
        call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
        dum4=min(sol_loc_1(1),1d0)
        dum6=FNORM
        if(dum6>10d0*(ERREL**2)) then
            xguess_loc_1(1)=0.01d0
            ITMAX=1000000
            !CALL D_NEQNF(labors, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
            call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
            if(FNORM<dum6) then
                dum4=min(sol_loc_1(1),1d0)
                dum6=FNORM
            end if
        end if
        if(dum6>10d0*(ERREL**2)) then
            xguess_loc_1(1)=0.3d0
            ITMAX=1000000
            !CALL D_NEQNF(labors, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
            call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
            if(FNORM<dum6) then
                dum4=min(sol_loc_1(1),1d0)
                dum6=FNORM
            end if
            !if(dum6>10d0*(ERREL**2)) then
            !    Print *,'ik is',ik
            !    Print *,'iam is',iam
            !    Print *,'ium is',ium
            !    Print *,'iaf is',iaf
            !    Print *,'iuf is',iuf
            !    Print *,'dum6 is',dum6
            !    Print *,'ERREL is',ERREL
            !end if
        end if

        laborsinglem(ik,ium)=dum4
    end do

    ind2=2

    do ium=1,nw
        wagem=wage_grid(ium)
        xguess_loc_1(1)=0.99d0
        ITMAX=1000000
        !For each level of capital we call D_LCONF to find optimal labor
        !CALL D_NEQNF(labors, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
        call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
        dum4=min(sol_loc_1(1),1d0)
        dum6=FNORM
        if(dum6>10d0*(ERREL**2)) then
            xguess_loc_1(1)=0.01d0
            ITMAX=1000000
            !CALL D_NEQNF(labors, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
            call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
            if(FNORM<dum6) then
                dum4=min(sol_loc_1(1),1d0)
                dum6=FNORM
            end if
        end if
        if(dum6>10d0*(ERREL**2)) then
            xguess_loc_1(1)=0.3d0
            ITMAX=1000000
            !CALL D_NEQNF(labors, SOL, ERREL,ITMAX=ITMAX,XGUESS=XGUESS,FNORM=FNORM)
            call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
            if(FNORM<dum6) then
                dum4=min(sol_loc_1(1),1d0)
                dum6=FNORM
            end if
            !if(dum6>10d0*(ERREL**2)) then
            !    Print *,'ik is',ik
            !    Print *,'iam is',iam
            !    Print *,'ium is',ium
            !    Print *,'iaf is',iaf
            !    Print *,'iuf is',iuf
            !    Print *,'dum6 is',dum6
            !    Print *,'ERREL is',ERREL
            !end if
        end if

        laborsinglef(ik,ium)=dum4
    end do

end subroutine lsupply
