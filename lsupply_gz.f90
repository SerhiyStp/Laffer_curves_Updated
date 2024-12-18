subroutine lsupply(ik)
    !This subroutine computes optimal policies durin active worklife.
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
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

    real(8) :: sol_test(2), f2(2)
    real(8) :: xguess_loc(2), sol_loc(2)
    real(8) :: xguess_loc_1(1), sol_loc_1(1)
    integer :: i_reopt
    logical :: verbose
    
    real(8) :: sol_loc_gz(2), hm_gz, hf_gz, fnorm_gz
    
    integer :: ii
    
    integer, parameter :: nh_test = 50
    integer :: i
    real(8) :: hgrid_test(nh_test), f_test(nh_test), foc_test, h_test
    real(8), parameter :: dh = 1d0/(nh_test-1)    

    EXTERNAL labor2_hybrd
    EXTERNAL labor2_gz_hybrd
    EXTERNAL labor3_hybrd
    EXTERNAL labor1_hybrd
    EXTERNAL labors_hybrd

    !CALL ERSET(IERSVR, IPACT, ISACT)
    hgrid_test = [(dh*i, i = 0, nh_test-1)]
    
    verbose = .true.  

    dum=c_grid(ik)

    ! Married couples:
    do ium=1,nw
        do iuf=1,nw
            wagem=wage_grid(ium)
            wagef=wage_grid(iuf)
            
            if (iuf > 1) then
                xguess_loc = [laborm(ik,ium,iuf-1), laborf(ik,ium,iuf-1)]
            else if (ium > 1) then
                xguess_loc = [laborm(ik,ium-1,iuf), laborf(ik,ium-1,iuf)]
            else 
                XGUESS_loc(1)=0.99d0
                XGUESS_loc(2)=0.99d0
            end if
            ITMAX=1000000
            xguess_loc(1) = xguess_loc(1) - 1.0d0
            xguess_loc(2) = xguess_loc(2) - 1.0d0
            call getHybrSoln(xguess_loc, 2, labor2_gz_hybrd, sol_loc, fnorm)
            dum4 = 1.0d0 - max(0.0d0,-sol_loc(1))
            dum5 = 1.0d0 - max(0.0d0,-sol_loc(2))
            dum6=FNORM
            
            if(dum6>10d0*(ERREL**2)) then
                XGUESS2(1)=0.01d0
                XGUESS2(2)=0.01d0
                ITMAX=1000000
                xguess2(1) = xguess2(1) - 1.0d0
                xguess2(2) = xguess2(2) - 1.0d0                
                call getHybrSoln(XGUESS2, 2, labor2_gz_hybrd, sol_loc, fnorm)
                if(FNORM<dum6) then
                    dum4=1.0d0 - max(0.0d0,-sol_loc(1))
                    dum5=1.0d0 - max(0.0d0,-sol_loc(2))
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                XGUESS2(1)=0.3d0
                XGUESS2(2)=0.3d0
                ITMAX=1000000
                xguess2(1) = xguess2(1) - 1.0d0
                xguess2(2) = xguess2(2) - 1.0d0                
                call getHybrSoln(XGUESS2, 2, labor2_gz_hybrd, sol_loc, fnorm)
                if(FNORM<dum6) then
                    dum4=1.0d0 - max(0.0d0,-sol_loc(1))
                    dum5=1.0d0 - max(0.0d0,-sol_loc(2))
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                XGUESS2(1)=0.1d0
                XGUESS2(2)=0.6d0
                ITMAX=1000000
                xguess2(1) = xguess2(1) - 1.0d0
                xguess2(2) = xguess2(2) - 1.0d0                
                call getHybrSoln(XGUESS2, 2, labor2_gz_hybrd, sol_loc, fnorm)
                if(FNORM<dum6) then
                    dum4=1.0d0 - max(0.0d0,-sol_loc(1))
                    dum5=1.0d0 - max(0.0d0,-sol_loc(2))
                    dum6=FNORM
                end if
            end if   
            if(dum6>10d0*(ERREL**2)) then
                XGUESS2(1)=0.6d0
                XGUESS2(2)=0.1d0
                ITMAX=1000000
                xguess2(1) = xguess2(1) - 1.0d0
                xguess2(2) = xguess2(2) - 1.0d0                
                call getHybrSoln(XGUESS2, 2, labor2_gz_hybrd, sol_loc, fnorm)
                if(FNORM<dum6) then
                    dum4=1.0d0 - max(0.0d0,-sol_loc(1))
                    dum5=1.0d0 - max(0.0d0,-sol_loc(2))
                    dum6=FNORM
                end if
            end if    
            
            if(dum6>1d-4) then
                i_reopt = 1
                do while (dum6>1d-4 .and. i_reopt < 100)
                    call RANDOM_NUMBER(XGUESS2(1))
                    call RANDOM_NUMBER(XGUESS2(2))
                    ITMAX=1000000
                    xguess2(1) = xguess2(1) - 1.0d0
                    xguess2(2) = xguess2(2) - 1.0d0                
                    call getHybrSoln(XGUESS2, 2, labor2_gz_hybrd, sol_loc, fnorm)
                    if(FNORM<dum6) then
                        dum4=1.0d0 - max(0.0d0,-sol_loc(1))
                        dum5=1.0d0 - max(0.0d0,-sol_loc(2))
                        dum6=FNORM
                    end if 
                    i_reopt = i_reopt + 1
                end do
                if (verbose) then
                    if (dum6>1d-4) then
                        print *, 'reoptimisation failed (married)'
                    else
                        print '(a, i4, a, ES9.2)', '(married) successfully reoptimized in ', i_reopt-1, ' tries, fnorm = ', dum6
                    end if
                end if
            end if
            
            
            if (dum6>1d-4) then
                print *, 'WARNING: hours choice for married couple with 2 spouses working not solved'
            end if

            laborm(ik,ium,iuf)=dum4
            laborf(ik,ium,iuf)=dum5
        end do
    end do

    !Male works
    ind2=1
    do ium=1,nw
        wagem=wage_grid(ium)
        h_test = 1d0
        call labor1(h_test, foc_test, 1)
        
        if (foc_test > 0d0) then
            labormwork(ik,ium)=h_test
        else                
        
            if (ium > 1) then
                xguess_loc_1(1)=labormwork(ik,ium-1)
            else
                xguess_loc_1(1)=0.99d0
            end if
            ITMAX=1000000
            call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
            dum6=FNORM
            dum4=min(sol_loc_1(1),1d0)
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.01d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.3d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
        
            if(dum6>1d-4) then
                i_reopt = 1
                do while (dum6>1d-4 .and. i_reopt < 15)
                    call RANDOM_NUMBER(xguess_loc_1(1))
                    ITMAX=1000000
                    call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                    if(FNORM<dum6) then
                        dum4=min(sol_loc_1(1),1d0)
                        dum6=FNORM
                    end if   
                    i_reopt = i_reopt + 1
                end do
            end if                    
            
            if (dum6>1d-4) then
                open(11, file='htest.txt')
                do i = 1, nh_test
                    call labor1 (hgrid_test(i), f_test(i), 1)    
                    write(11, '(2f12.6)') hgrid_test(i), f_test(i)
                end do
                close(11)              
            
            
                print *, 'WARNING: hours choice for married couple with only man working not solved'
            end if        
            labormwork(ik,ium)=dum4
        end if
    end do

        
    !Female works
    ind2=2
    do ium=1,nw
        wagem=wage_grid(ium)
        
        h_test = 1d0
        call labor1(h_test, foc_test, 1)
        
        if (foc_test > 0d0) then
            laborfwork(ik,ium)=h_test
        else        
        
            if (ium > 1) then
                xguess_loc_1(1)=labormwork(ik,ium-1)
            else        
                xguess_loc_1(1)=0.99d0
            end if
            ITMAX=1000000
            call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
            dum6=FNORM
            dum4=min(sol_loc_1(1),1d0)
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.01d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.3d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
        
            if(dum6>1d-4) then
                i_reopt = 1
                do while (dum6>1d-4 .and. i_reopt < 15)
                    call RANDOM_NUMBER(xguess_loc_1(1))
                    ITMAX=1000000
                    call getHybrSoln(xguess_loc_1, 1, labor1_hybrd, sol_loc_1, fnorm)
                    if(FNORM<dum6) then
                        dum4=min(sol_loc_1(1),1d0)
                        dum6=FNORM
                    end if   
                    i_reopt = i_reopt + 1
                end do
            end if    
        
            if (dum6>1d-4) then
                open(11, file='htest.txt')
                do i = 1, nh_test
                    call labor1 (hgrid_test(i), f_test(i), 1)    
                    write(11, '(2f12.6)') hgrid_test(i), f_test(i)
                end do
                close(11)               
            
                print *, 'WARNING: hours choice for married couple with only woman working not solved'
            end if        
            laborfwork(ik,ium)=dum4
        end if
    end do   
        
    ! Singles:
    ind2=1
    do ium=1,nw
        wagem=wage_grid(ium)
        
        h_test = 1d0
        call labors(h_test, foc_test, 1)
        
        if (foc_test > 0d0) then
            laborsinglem(ik,ium)=h_test
        else
        
            if (ium > 1) then
                xguess_loc_1(1)=laborsinglem(ik,ium-1)
            else
                xguess_loc_1(1)=0.99d0
            end if
            ITMAX=1000000
            call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
            dum4=min(sol_loc_1(1),1d0)
            dum6=FNORM
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.01d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.3d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
        
            if(dum6>1d-4) then
                i_reopt = 1
                do while (dum6>1d-4 .and. i_reopt < 15)
                    call RANDOM_NUMBER(xguess_loc_1(1))
                    ITMAX=1000000
                    call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
                    if(FNORM<dum6) then
                        dum4=min(sol_loc_1(1),1d0)
                        dum6=FNORM
                    end if   
                    i_reopt = i_reopt + 1
                end do
                if (verbose) then
                    if (dum6>1d-4) then
                    
                        open(11, file='htest.txt')
                        do i = 1, nh_test
                            call labors (hgrid_test(i), f_test(i), 1)    
                            write(11, '(2f12.6)') hgrid_test(i), f_test(i)
                        end do
                        close(11)                     
                    
                        print *, 'reoptimisation failed (single M)'
                    
                    else
                        print '(a, i4, a, ES9.2)', '(single M) successfully reoptimized in ', i_reopt-1, ' tries, fnorm = ', dum6
                    end if 
                end if            
            end if    
        
            if (dum6>1d-4) then
                print *, 'WARNING: hours choice for single men not solved'
            end if         

            laborsinglem(ik,ium)=dum4
        end if
    end do

    ind2=2
    do ium=1,nw
        wagem=wage_grid(ium)
        
        h_test = 1d0
        call labors(h_test, foc_test, 1)
        
        if (foc_test > 0d0) then
            laborsinglef(ik,ium)=h_test
        else        
        
            if (ium > 1) then
                xguess_loc_1(1)=laborsinglef(ik,ium-1)
            else
                xguess_loc_1(1)=0.99d0
            end if
            ITMAX=1000000
            call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
            dum4=min(sol_loc_1(1),1d0)
            dum6=FNORM
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.01d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
            if(dum6>10d0*(ERREL**2)) then
                xguess_loc_1(1)=0.3d0
                ITMAX=1000000
                call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
                if(FNORM<dum6) then
                    dum4=min(sol_loc_1(1),1d0)
                    dum6=FNORM
                end if
            end if
        
            if(dum6>1d-4) then
                i_reopt = 1
                do while (dum6>1d-4 .and. i_reopt < 15)
                    call RANDOM_NUMBER(xguess_loc_1(1))
                    ITMAX=1000000
                    call getHybrSoln(xguess_loc_1, 1, labors_hybrd, sol_loc_1, fnorm)
                    if(FNORM<dum6) then
                        dum4=min(sol_loc_1(1),1d0)
                        dum6=FNORM
                    end if   
                    i_reopt = i_reopt + 1
                end do
                if (verbose) then
                    if (dum6>1d-4) then
                    
                        open(11, file='htest.txt')
                        do i = 1, nh_test
                            call labors (hgrid_test(i), f_test(i), 1)    
                            write(11, '(2f12.6)') hgrid_test(i), f_test(i)
                        end do
                        close(11)                         
                    
                        print *, 'reoptimisation failed (single F)'
                    else
                        print '(a, i4, a, ES9.2)', '(single F) successfully reoptimized in ', i_reopt-1, ' tries, fnorm = ', dum6
                    end if 
                end if             
            end if    
        
            if (dum6>1d-4) then
                print *, 'WARNING: hours choice for single women not solved'
            end if         

            laborsinglef(ik,ium)=dum4
        end if
    end do
    

    
end subroutine lsupply
