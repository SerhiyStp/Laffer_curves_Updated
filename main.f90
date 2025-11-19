program Laffer

    use Utilities
    use Model_Parameters
    use PolicyFunctions
    use Tauchen
    use hybrd_wrapper, only: setHybrParams
    use, intrinsic :: iso_fortran_env, only: output_unit
    use partests
    use SolveActive

    implicit none
    integer :: ik,tprint,it2,it3,it4,it6,it7,it8,ium,iam,iuf,iaf,ix,j,iu2,ik2,ifc,counter,iter_ratio
    real(8) :: dum,dum2,dum3,dum4,dum5,dum6,epsilon_ratio=1d0,epsilon_ratio_old=1d0,step_ratio=0.05
    integer :: i_k, i, ii
    EXTERNAL labor2
    EXTERNAL labor1
    EXTERNAL labor3
    EXTERNAL labors

    !call OMP_SET_NUM_THREADS(106)

    write(output_unit, *) "Program started"
    call Initialize
    call setHybrParams(2)

    print *, tSS_employee(0.1d0)
    print *, tSS_employee(0.0d0)


    iter_ratio=1

    open(61, file="Laffer_Results.txt")

    !tax_level_scale = 1.1d0
    tax_level_scale = 1.0d0

    do while (tax_level_scale > 0.2d0)

        write(61, *) "========================="
        write(61, "(a, f10.6)") "tax_level_scale = ", tax_level_scale
        write(61, *) "========================="

        theta(1) = 0.975d0*tax_level_scale
        thetas(1) = 0.895d0*tax_level_scale
        epsilon_ratio=1d0

        do while(abs(epsilon_ratio)>0.003d0)

            epsilon=1d0
            epsilon2=1d0
            epsilon3=1d0
            epsilon5=1d0

            iter=0

            do while((abs(epsilon3)>0.001d0).OR.(abs(epsilon)>0.001d0).OR.(abs(epsilon2)>0.001d0).OR.(abs(epsilon5)>0.001d0))

                iter=iter+1

                dum2= 1.5d0*wage(1,a(1,na),dble(T),u(1,nu))/(1d0+t_employer)
                call MakeGrid(nw,wage_grid,0.01d0,dum2,2d0)
                !!$OMP PARALLEL PRIVATE(ik)
                !!$OMP DO SCHEDULE(DYNAMIC)
                do ik=1,nc
                    call lsupply(ik)
                end do
                !!$OMP END DO    
                !!$OMP END PARALLEL

                !open(22, file='laborm.txt')
                !open(23, file='laborf.txt')
                !ium = 40
                !iuf = 40
                !do ii = 1, nk
                !    write(22, '(2f12.6)') c_grid(ii), laborm(ii,10,10) !, laborm(ii,10,80), labor
                !    write(23, '(2f12.6)') c_grid(ii), laborf(ii,ium,iuf)
                !end do
                !close(22)
                !close(23)

                print *, 'hours worked problem solved'            

                do it=1,Tret

                    Print *,'t is',T+Tret+1-it

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*nexp*na
                        call SolveInRetirement(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL


                    ev_ret(:,:,:,:,:,Tret-it+1)=v_ret(:,:,:,:,:,Tret-it+1)
                    edc_ret(:,:,:,:,:,Tret-it+1)=Uprime_ret(:,:,:,:,:,Tret-it+1)
                    evs_ret(:,:,:,:,Tret-it+1)=vs_ret(:,:,:,:,Tret-it+1)
                    edcs_ret(:,:,:,:,Tret-it+1)=Uprimes_ret(:,:,:,:,Tret-it+1)

                    ! Compute spline coefficients:

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter = 1, na
                        call partest10(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL               

                end do

                !Compute optimal policies at age 64

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call Solvelastactive(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                it=0
                Print *,'t is',T-it        

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call partest8(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nu*na*nfc
                    call partest5(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call partest7(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                ! Compute spline coefficients:

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nu*na*nfc
                    call partest3(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                ! Compute spline coefficients:


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call partest(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*nexp*na*nu
                    call partest2(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nu*na*nfc
                    call partest6(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                call update_ev_aux(T)
                call update_lfp_policies(T)


                !Compute optimal policies for age 2-63

                do it=1,T-2

                    Print *,'t is',T-it
                    !Print *,'Gamma_redistr',Gamma_redistr/2d0

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*nexp*na*nu
                        call SolveActiveLife(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL


                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*na*nu
                        call partest8(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL


                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter = 1, nu*na*nfc
                        call partest5(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter = 1, nk*nexp*na*nu
                        call partest7(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL   

                    ! Compute spline coefficients:

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter = 1, nfc*na*nu
                        call partest3(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL


                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*nexp*na*nu
                        call partest(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL


                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*nexp*na*nu
                        call partest2(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)
                    do counter = 1, nfc*na*nu
                        call partest6(counter)
                    end do
                    !$OMP END DO    
                    !$OMP END PARALLEL

                    call update_ev_aux(T-it)
                    call update_lfp_policies(T-it)

                end do


                it=T-1               

                !Print *,'t is',T-it

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call Solvefirstactive(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL  

                !$OMP PARALLEL PRIVATE(ik)
                !$OMP DO SCHEDULE(DYNAMIC)
                do ik=1,nsim2
                    call Simulation(ik)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

                call Statistics(output_unit)

                !open(1, file='singledist.txt')
                !write (1, *) fpartner,mpartner
                !close(1)
                !
                !open(1, file='abilityprob.txt')
                !write (1, *) ability_prob
                !close(1)

                Print *,'epsilon is',epsilon
                Print *,'epsilon2 is',epsilon2
                Print *,'epsilon3 is',epsilon3
                Print *,'epsilon5 is',epsilon5

                open(41, file='last_results.txt')
                write(41, *) AE, lumpsum, psi0, gamma_redistr 
                close(41)
                !Variables are age, gender, ID number, weight, marital status, asset holdings, household labor income, Household_Labor_Income_Tax_Paid,  Household_consumption_Tax_Paid, ability, hours, earnings
                !open(1, file='Simulation_output.txt')
                !do it2=1,T
                !    do it3=1,nsim2
                !        do it4=1,nsim
                !            write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') it2*1d0, 1d0, nsim*(it3-1)*1d0+it4*1d0, 1d0, Sim1m(it3,it4,it2,10), Sim1m(it3,it4,it2,1), Sim1m(it3,it4,it2,6), Sim1m(it3,it4,it2,7), Sim1m(it3,it4,it2,8), exp1m(it3,it4,it2,2)*1d0, Sim1m(it3,it4,it2,4), Sim1m(it3,it4,it2,5)  
                !            write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') it2*1d0, 2d0, nsim*(it3-1)*1d0+it4*1d0, 1d0, Sim1f(it3,it4,it2,10), Sim1f(it3,it4,it2,1), Sim1f(it3,it4,it2,6), Sim1f(it3,it4,it2,7), Sim1f(it3,it4,it2,8), exp1f(it3,it4,it2,2)*1d0, Sim1f(it3,it4,it2,4), Sim1f(it3,it4,it2,5)
                !        end do
                !    end do
                !end do
                !close(1)

                write(output_unit, *) 'Program finished'
                stop

            end do

            !epsilon_ratio=ratiodum-ratio
            epsilon_ratio=0d0
            Print *,'epsilon_ratio is',epsilon_ratio

        end do

        call Statistics(61)
        !tax_level_scale = tax_level_scale - 0.01d0

        open(1, file='Simulation_output.txt')
        do it2=1,T
            do it3=1,nsim2
                do it4=1,nsim
                    write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') it2*1d0, 1d0, nsim*(it3-1)*1d0+it4*1d0, 1d0, Sim1m(it3,it4,it2,10), Sim1m(it3,it4,it2,1), Sim1m(it3,it4,it2,6), Sim1m(it3,it4,it2,7), Sim1m(it3,it4,it2,8), exp1m(it3,it4,it2,2)*1d0, Sim1m(it3,it4,it2,4), Sim1m(it3,it4,it2,5)  
                    write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') it2*1d0, 2d0, nsim*(it3-1)*1d0+it4*1d0, 1d0, Sim1f(it3,it4,it2,10), Sim1f(it3,it4,it2,1), Sim1f(it3,it4,it2,6), Sim1f(it3,it4,it2,7), Sim1f(it3,it4,it2,8), exp1f(it3,it4,it2,2)*1d0, Sim1f(it3,it4,it2,4), Sim1f(it3,it4,it2,5)
                end do
            end do
        end do
        close(1)

        write(output_unit, *) 'Program finished'
        stop

    end do

    close(61)

    contains

    subroutine Initialize()

        !USE ANORDF_INT
        implicit none

        print *, Tret
        print *, T

        allocate(v(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(ev(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(c(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(k(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(nm(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(nf(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(v_aux(nk,nexp,nexp,na,nu,na,nu,nfc,nfcm,2))
        allocate(ev_aux(nk,nexp,nexp,na,nu,na,nu,nfc,nfcm,2))

        allocate(lfpm(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(lfpf(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm))
        allocate(c_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))
        allocate(v_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2,2))
        allocate(k_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))  
        allocate(nm_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))
        allocate(nf_lfp(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfcm,2,2))           

        allocate(vs(2,nk,nexp,na,nu,T,nfc))
        allocate(evs(2,nk,nexp,na,nu,T,nfc))
        allocate(evm(2,nk,nexp,na,nu,T,nfc))
        allocate(cs(2,nk,nexp,na,nu,T,nfc))
        allocate(edcs(2,nk,nexp,na,nu,T,nfc))
        allocate(Uprimes(2,nk,nexp,na,nu,T,nfc))
        allocate(ks(2,nk,nexp,na,nu,T,nfc))
        allocate(ns(2,nk,nexp,na,nu,T,nfc))

        allocate(lfps(2,nk,nexp,na,nu,T,nfc))
        allocate(cs_lfp(2,nk,nexp,na,nu,T,nfc,2))
        allocate(ns_lfp(2,nk,nexp,na,nu,T,nfc,2))
        allocate(ks_lfp(2,nk,nexp,na,nu,T,nfc,2))  
        allocate(vs_lfp(2,nk,nexp,na,nu,T,nfc,2))          

        allocate(Sim1m(nsim2,nsim,T+1,11))
        allocate(Sim1f(nsim2,nsim,T+1,11))
        allocate(exp1m(nsim2,nsim,T+1,6))
        allocate(exp1f(nsim2,nsim,T+1,6))
        allocate(exp2m(nsim2,nsim,T+1,6))
        allocate(exp2f(nsim2,nsim,T+1,6))
        allocate(SimR1m(nsim2,nsim,Tret+1,14))
        allocate(SimR1f(nsim2,nsim,Tret+1,14))
        allocate(expR1m(nsim2,nsim,Tret+1,4))
        allocate(expR1f(nsim2,nsim,Tret+1,4))
        allocate(Random3m(nsim2,nsim,T+Tret))
        allocate(Random3f(nsim2,nsim,T+Tret))
        allocate(marstatm(nsim2,nsim,T))
        allocate(marstatf(nsim2,nsim,T))
        allocate(marstatm_init(nsim2,nsim))
        allocate(marstatf_init(nsim2,nsim))
        allocate(partshock(nsim2,nsim,1))
        allocate(partshock2(nsim2,nsim,1))
        allocate(Random1m(nsim2,nsim))
        allocate(Random1f(nsim2,nsim))
        allocate(Random2m(nsim2,nsim))
        allocate(Random2f(nsim2,nsim))

        allocate(a(2,na))
        allocate(Prob_a(2,na))
        allocate(fc(2,nfc))
        allocate(fcage(2,2))
        allocate(Prob_fc(2,nfc))
        allocate(fcm(2,nfcm))
        allocate(Prob_fcm(2,nfcm))
        allocate(u(2,nu))
        allocate(Prob_u(2,nu))
        allocate(trans_u(2,nu,nu))
        allocate(trans_a(2,na,na))
        allocate(trans_fc(2,nfc,nfc))
        allocate(trans_fcm(2,nfcm,nfcm))
        allocate(OmegaRet(Tret))
        allocate(OmegaRet2(Tret))
        allocate(OmegaActive(T))
        allocate(Probm(T))
        allocate(Probd(T))
        allocate(WeightRet(Tret))
        allocate(WeightActive(T))

        allocate(fpartner(nk,nexp,na,nu,T,nfc))
        allocate(mpartner(nk,nexp,na,nu,T,nfcm))
        allocate(fpartnerdum(nk,nexp,na,nu,T,nfc))
        allocate(mpartnerdum(nk,nexp,na,nu,T,nfcm))
        allocate(fpartnerdum2(nk,nexp,na,nu,T,nfc))
        allocate(mpartnerdum2(nk,nexp,na,nu,T,nfcm))
        allocate(ability_prob(na,na))
        allocate(av_earnings(2,2,na))
        allocate(laborm(nc,nw,nw))
        allocate(laborf(nc,nw,nw))
        allocate(labormwork(nc,nw))
        allocate(laborfwork(nc,nw))
        allocate(laborsinglem(nc,nw))
        allocate(laborsinglef(nc,nw))

        allocate(c_grid(nc))
        allocate(wage_grid(nw))
        allocate(k_grid(nk))
        allocate(exp_grid(nexp,T+Tret))
        !allocate(K_KNOT(nk+KORDER))
        !allocate(EXP_KNOT(nexp+EXPORDER,T+Tret))

        allocate(c_ret(nk,nexp,nexp,na,na,Tret))
        allocate(edc_ret(nk,nexp,nexp,na,na,Tret))
        allocate(Uprime_ret(nk,nexp,nexp,na,na,Tret))
        allocate(v_ret(nk,nexp,nexp,na,na,Tret))
        allocate(ev_ret(nk,nexp,nexp,na,na,Tret))
        allocate(k_ret(nk,nexp,nexp,na,na,Tret))
        allocate(vs_ret(2,nk,nexp,na,Tret))
        allocate(evs_ret(2,nk,nexp,na,Tret))
        allocate(cs_ret(2,nk,nexp,na,Tret))
        allocate(Eulers_ret(2,nk,nexp,na,Tret))
        allocate(edcs_ret(2,nk,nexp,na,Tret))
        allocate(Uprimes_ret(2,nk,nexp,na,Tret))
        allocate(ks_ret(2,nk,nexp,na,Tret))
        allocate(break(nk))

        allocate(pol_v_mar_lfp(na, nu, na, nu, T, nfc, nfcm, 2, 2, 2))


        call init_Tax_ss(test=.true.)

        !print *, t_employee1
        !print *, t_employee2


        open(1, file='singledist.txt')

        read (1, *) fpartner,mpartner

        close(1)
        !
        open(1, file='abilityprob.txt')

        read (1, *) ability_prob

        close(1)

        !open(1, file='av_earnings.txt')
        !
        !read (1, *) av_earnings
        !
        !close(1)


        !Print *,ability_prob(5,:)
        !STOP
        !ability_prob=1d0/7d0
        !
        !fpartner=0d0
        !mpartner=0d0
        !
        !do it=1,T
        !    do it2=1,nexp
        !        do ik=1,8
        !            mpartner(ik,it2,:,:,it,:)=1d0/(8*it*na*nu*nfcm)
        !        end do
        !    end do
        !end do
        !
        !do it=1,T
        !    do it2=1,nexp
        !        do ik=1,8
        !            fpartner(ik,it2,:,:,it,:)=1d0/(8*it*na*nu*nfc)
        !        end do
        !    end do
        !end do

        trans_u = 0d0
        prob_u=0d0
        trans_a = 0d0
        prob_a=0d0
        gamma(1,:) = (/ 0.0605927d0, -0.0010648d0, 0.0000093d0 /)
        gamma(2,:) = (/ 0.0784408d0, -0.0025596d0, 0.0000256d0 /)
        !gamma(1,:) = (/ 0.0690d0, -0.00129d0, 0.0d0 /)
        !gamma(2,:) = (/ 0.0430d0, -0.00078d0, 0.0d0 /)
        gamma0=-0.1070d0
        gamma0f=-0.1450d0

        theta(:) = (/ 0.975d0*tax_level_scale, 0.149d0*tax_prog_scale /)
        thetas(:) = (/ 0.895d0*tax_level_scale, 0.140d0*tax_prog_scale /)
        AE = 1d0
        Unemp_benefit=0.201795*AE
        !Unemp_benefit=0d0
        r=alpha*ratio**(alpha-1d0)-delta
        w=(1d0-alpha)*ratio**alpha
        call MakeGrid(nk,k_grid,0d0,100d0,3d0)
        call MakeGrid(nc,c_grid,0.01d0,100d0,3d0)

        exp_grid=0d0
        do it2=2,T
            call MakeGrid(nexp,exp_grid(:,it2),0d0,1d0*(it2-1),1d0)
        end do

        do it2=T+1,T+Tret
            call MakeGrid(nexp,exp_grid(:,it2),0d0,1d0*T,1d0)
        end do

        exp_grid(:,1)=exp_grid(:,2)


        !STOP

        !Filling in US divorce and marriage probabilities

        open(1, file='divprob.txt')

        do it2=1,T
            read (1, *) probd(it2)
        end do

        close(1)

        open(1, file='marprob.txt')

        do it2=1,T
            read (1, *) probm(it2)
        end do

        close(1)



        !probm=0d0
        !probd=0d0

        OmegaActive=1d0

        OmegaRet(1)=1d0-0.014319d0
        OmegaRet(2)=1d0-0.015540d0
        OmegaRet(3)=1d0-0.016920d0
        OmegaRet(4)=1d0-0.018448d0
        OmegaRet(5)=1d0-0.020170d0
        if (testing == 0) then
            OmegaRet(6)=1d0-0.022022d0
            OmegaRet(7)=1d0-0.023973d0
            OmegaRet(8)=1d0-0.026203d0
            OmegaRet(9)=1d0-0.028771d0
            OmegaRet(10)=1d0-0.031629d0
            OmegaRet(11)=1d0-0.034611d0
            OmegaRet(12)=1d0-0.037710d0
            OmegaRet(13)=1d0-0.041264d0
            OmegaRet(14)=1d0-0.045405d0
            OmegaRet(15)=1d0-0.050128d0
            OmegaRet(16)=1d0-0.055339d0
            OmegaRet(17)=1d0-0.061005d0
            OmegaRet(18)=1d0-0.067396d0
            OmegaRet(19)=1d0-0.074476d0
            OmegaRet(20)=1d0-0.082272d0
            OmegaRet(21)=1d0-0.091816d0
            OmegaRet(22)=1d0-0.101898d0
            OmegaRet(23)=1d0-0.112870d0
            OmegaRet(24)=1d0-0.124763d0
            OmegaRet(25)=1d0-0.137597d0
            OmegaRet(26)=1d0-0.151383d0
            OmegaRet(27)=1d0-0.166117d0
            OmegaRet(28)=1d0-0.181778d0
            OmegaRet(29)=1d0-0.198331d0
            OmegaRet(30)=1d0-0.215721d0
            OmegaRet(31)=1d0-0.233874d0
            OmegaRet(32)=1d0-0.252699d0
            OmegaRet(33)=1d0-0.272086d0
            OmegaRet(34)=1d0-0.291912d0
            OmegaRet(35)=1d0-0.312040d0
            OmegaRet(36)=1d0-1d0
        end if

        !OmegaRet=1d0

        OmegaRet2(1)=1d0
        do it2=1,Tret-1
            OmegaRet2(it2+1)=OmegaRet(it2)
        end do



        call tauchen_hans(sigma_am,rho_am,na,a(1,:),trans_a(1,:,:),prob_a(1,:))
        call tauchen_hans2(sigma_um,rho_um,nu,u(1,:),trans_u(1,:,:),prob_u(1,:))

        call tauchen_hans(sigma_af,rho_af,na,a(2,:),trans_a(2,:,:),prob_a(2,:))
        call tauchen_hans2(sigma_uf,rho_uf,nu,u(2,:),trans_u(2,:,:),prob_u(2,:))

        call tauchen_hans(sigma_fcm,rho_fcm,nfc,fc(1,:),trans_fc(1,:,:),prob_fc(1,:))
        call tauchen_hans(sigma_fcs,rho_fcs,nfc,fc(2,:),trans_fc(2,:,:),prob_fc(2,:))

        call tauchen_hans(sigma_fcmm,rho_fcmm,nfc,fcm(1,:),trans_fcm(1,:,:),prob_fcm(1,:))
        call tauchen_hans(sigma_fcsm,rho_fcsm,nfc,fcm(2,:),trans_fcm(2,:,:),prob_fcm(2,:))

        fc(1,:)=fc(1,:)+mu_fcm
        fc(2,:)=fc(2,:)+mu_fcs
        fcm(1,:)=fcm(1,:)+mu_fcmm
        fcm(2,:)=fcm(2,:)+mu_fcsm

        fcage(1,1)=mu_fcm1
        fcage(1,2)=mu_fcm2

        !fc=0d0
        !fcm=0d0
        !
        !fc(1,:)=fc(1,:)+mu_fcm
        !fc(2,:)=fc(2,:)+mu_fcs
        !fcm(1,:)=fcm(1,:)+mu_fcmm
        !fcm(2,:)=fcm(2,:)+mu_fcsm

        !Print *,prob_fc(1,:)
        !Print *,prob_fc(2,:)
        !Print *,prob_fcm(1,:)
        !Print *,prob_fcm(2,:)
        !STOP

        !Print *,fc(1,:)
        !Print *,fc(2,:)
        !Print *,fcm(1,:)
        !Print *,fcm(2,:)
        !STOP
        call random_number(random1m)

        call random_number(random2m)

        call random_number(random3m)

        call random_number(random1f)

        call random_number(random2f)

        call random_number(random3f)

        call random_number(marstatm)

        call random_number(marstatf)

        call random_number(marstatm_init)

        call random_number(marstatf_init)

        call random_number(partshock)

        call random_number(partshock2)



    end subroutine Initialize  

    subroutine init_Tax_ss(test)
        !real(8) :: x1, x2
        real(8) :: f1, f2
        real(8) :: df1, df2
        real(8) :: A(4,4)
        real(8) :: B(4)
        integer :: n, nrhs, lda, ipiv, ldb, info        
        integer, parameter :: ntest = 100
        real(8) :: ytest(ntest), ttest(ntest), dttest(ntest)
        real(8) :: ylo, yhi, dy
        integer :: i
        real(8) :: yy
        logical :: test

        !eps = 0.01d0
        !x1_ss = SocSecCap - ss_eps
        !x2_ss = SocSecCap + ss_eps

        f1 = t_employee1*x1_ss
        f2 = t_employee1*SocSecCap + t_employee2*ss_eps

        df1 = t_employee1
        df2 = t_employee2

        A(1,:) = [1d0, x1_ss, x1_ss**2d0, x1_ss**3d0]
        A(2,:) = [1d0, x2_ss, x2_ss**2d0, x2_ss**3d0]
        A(3,:) = [0d0, 1d0, 2d0*x1_ss, 3d0*x1_ss**2d0]
        A(4,:) = [0d0, 1d0, 2d0*x2_ss**2d0, 3d0*x2_ss**2d0]
        B = [f1, f2, df1, df2]        

        n = 4
        nrhs = 1
        lda = n
        ldb = n
        call DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        if (info == 0) then
            ss_coefs = B
        else
            print *, 'WARNING: failed to find coefficients of the smoothing part of the t_ss function'
        end if

        if (test == .true.) then
            ylo = 0.0d0
            yhi = 3.5d0
            dy = 1.0d0/(ntest-1)
            ytest = [(dy*i, i = 0, ntest-1)]
            ytest = ytest*(yhi - ylo) + ylo

            open(11, file='tss_test.txt')
            do i = 1, ntest
                yy = ytest(i) 
                ttest(i) = tSS_employee(yy)
                dttest(i) = DtSS_employee(yy)
                write(11, '(3f12.6)') yy, ttest(i), dttest(i)
            end do
            close(11)
        end if

    end subroutine init_Tax_ss


end program Laffer
