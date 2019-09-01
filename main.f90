 
 include 'link_fnl_shared.h'
 
!
! old:
! include 'link_f90_dll_smp.h'
      
! Add this source to your project to correct the problem
! that, for a Debug build for certain projects using the
! IMSL Fortran libraries, the linker complains of an
! unresolved external symbol __imp__MPIPRIV
BLOCK DATA MPIPRIV_DEF
COMMON /MPIPRIV/ DUMMY
!DEC$ ATTRIBUTES ALIAS:"__imp__MPIPRIV" :: /MPIPRIV/
END BLOCK DATA MPIPRIV_DEF


 program Laffer
        
        use Utilities
        use Model_Parameters
        use PolicyFunctions
        use Tauchen
        USE CSINT_INT
        USE CSVAL_INT
        USE ANORDF_INT
        use hybrd_wrapper, only: setHybrParams
        implicit none
        integer :: ik,tprint,it2,it3,it4,ium,iam,iuf,iaf,ix,j,iu2,ik2,ifc,counter
        real(8) :: dum,dum2,dum3
        EXTERNAL labor2
        EXTERNAL labor1
        EXTERNAL labor3
        EXTERNAL labors
        !print *, "hello"
        !call TestLinInterp 
        call OMP_SET_NUM_THREADS(40)
        
        call Initialize
        call setHybrParams(2)
        
        !do while(epsilon>0.001)
        
        !Compute optimal policies in retirement

        
        do while(abs(epsilon4)>0.003d0)

            epsilon=1d0
            epsilon2=1d0
            epsilon3=1d0
            epsilon5=1d0

            iter=0

            do while((abs(epsilon)>0.001d0).OR.(abs(epsilon2)>0.001d0).OR.(abs(epsilon3)>0.001d0).OR.(abs(epsilon5)>0.001d0).OR.(abs(epsilon6)>0.01d0))

                iter=iter+1

                do it=1,Tret

                    !$OMP PARALLEL PRIVATE(ik)
                    !$OMP DO SCHEDULE(DYNAMIC)

                    do ik=1,nk 
                        call SolveInRetirement(ik)
                    end do

                    !$OMP END DO    
                    !$OMP END PARALLEL

                    ! Compute spline coefficients:

                    !call d_csint(k_grid, v_ret(:,Tret-it+1), BREAK, ev_spln_coefs_ret(:,:,Tret-it+1))
                    call ppp_csint(k_grid, v_ret(:,Tret-it+1), ev_spln_coefs_ret(:,:,Tret-it+1), nk)

                    !call d_csint(k_grid, vs_ret(:,Tret-it+1), BREAK, evs_spln_coefs_ret(:,:,Tret-it+1))
                    call ppp_csint(k_grid, vs_ret(:,Tret-it+1), evs_spln_coefs_ret(:,:,Tret-it+1), nk)

                end do
 
                dum2= 1.5d0*wage(1,a(1,5),dble(T),u(1,5))/(1d0+t_employer)
                call MakeGrid(nw,wage_grid,0.01d0,dum2,2d0)

                !$OMP PARALLEL PRIVATE(ik)
                !$OMP DO SCHEDULE(DYNAMIC)

                do ik=1,nc
                    call lsupply(ik)
                end do

                !$OMP END DO    
                !$OMP END PARALLEL

                
                !Compute optimal policies at age 64

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
        
                do counter=1,nk*na*nu
                    call Solvelastactive(counter)
                end do
        
                !$OMP END DO    
                !$OMP END PARALLEL
                
                it=0
                
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

                do counter=1,nk*na*nu
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
                
                !!$OMP PARALLEL PRIVATE(counter)
                !!$OMP DO SCHEDULE(DYNAMIC)
                !
                !do counter = 1, nu*na*nfc
                !    call partest4(counter)
                !end do
                !
                !!$OMP END DO    
                !!$OMP END PARALLEL
                
                
                evm=0d0
                

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)

                do counter=1,nk*na*nu
                    call partest(counter)
                end do

                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)

                do counter=1,nk*na*nu
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

                
                !tprint=45
                !
                !open(1,file='csingle.txt')
                !    open(2,file='ksingle.txt')
                !    open(3,file='nsingle.txt')
                !    open(4,file='vsingle.txt')
                !    open(5,file='vsinglem.txt')
                !    do it2 = 1, nk
                !        write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,6,3,3,tprint,20),cs(2,it2,6,3,3,tprint,20),cs(2,it2,6,5,5,tprint,20)
                !        write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,6,3,3,tprint,20),ks(2,it2,6,3,3,tprint,20),ks(2,it2,6,5,5,tprint,20)
                !        write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,6,3,3,tprint,20),ns(2,it2,6,3,3,tprint,20),ns(2,it2,6,5,5,tprint,20)
                !        write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,6,5,5,tprint,20),evs(2,it2,6,3,3,tprint,20),evs(2,it2,6,5,5,tprint,20)
                !        write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,6,5,5,tprint,20),evm(2,it2,6,3,3,tprint,20),evm(2,it2,6,5,5,tprint,20)
                !    end do
                !    close(1)
                !    close(2)
                !    close(3)
                !    close(4)
                !    close(5)
                
                !open(1,file='vtest1.txt')
                !        open(2,file='vtest2.txt')
                !        open(3,file='vtest3.txt')
                !        open(4,file='vtest4.txt')
                !    do it2 = 1, nk
                !        write(1,'(4f12.6)'), k_grid(it2), evm_spln_coefs(1,1,it2,44,3,3,tprint,20),evm_spln_coefs(2,1,it2,44,3,3,tprint,20),evm_spln_coefs(2,1,it2,44,5,5,tprint,20)
                !        write(2,'(4f12.6)'), k_grid(it2), evm_spln_coefs(1,2,it2,44,3,3,tprint,20),evm_spln_coefs(2,2,it2,44,3,3,tprint,20),evm_spln_coefs(2,2,it2,44,5,5,tprint,20)
                !        write(3,'(4f12.6)'), k_grid(it2), evm_spln_coefs(1,3,it2,44,3,3,tprint,20),evm_spln_coefs(2,3,it2,44,3,3,tprint,20),evm_spln_coefs(2,3,it2,44,5,5,tprint,20)
                !        write(4,'(4f12.6)'), k_grid(it2), evm_spln_coefs(1,4,it2,44,3,3,tprint,20),evm_spln_coefs(2,4,it2,44,3,3,tprint,20),evm_spln_coefs(2,4,it2,44,5,5,tprint,20)
                !    end do
                !    close(1)
                !    close(2)
                !    close(3)
                !    close(4)

                
                    !open(1,file='ctest.txt')
                    !    open(2,file='ktest.txt')
                    !    open(3,file='vtest.txt')
                    !    open(4,file='nmtest.txt')
                    !    open(5,file='nftest.txt')
                    !    !open(5,file='evtest.txt')
                    !    !open(4,file='ve_vu.txt')
                    !    !open(4,file='tbc.txt')
                    !    do ik = 1, nk
                    !        write (1,'(4f12.6)') k_grid(ik),c(ik,6,1,1,1,1,tprint,20),c(ik,6,1,1,5,5,tprint,20),c(ik,6,3,3,3,3,tprint,20)
                    !        write (2,'(4f12.6)') k_grid(ik), k(ik,6,1,1,1,1,tprint,20), k(ik,6,1,1,5,5,tprint,20), k(ik,6,3,3,3,3,tprint,20)
                    !        write (3,'(4f12.6)') k_grid(ik),eV(ik,6,1,1,1,1,tprint,20),eV(ik,6,1,1,5,5,tprint,20),eV(ik,6,3,3,3,3,tprint,20)
                    !        write (4,'(4f12.6)') k_grid(ik),nm(ik,6,1,1,1,1,tprint,20),nm(ik,6,1,1,5,5,tprint,20),nm(ik,6,3,3,3,3,tprint,20)
                    !        write (5,'(4f12.6)') k_grid(ik),nf(ik,6,1,1,1,1,tprint,20),nf(ik,6,1,1,5,5,tprint,20),nf(ik,6,3,3,3,3,tprint,20)
                    !    end do
                    !    close(1)
                    !    close(2)
                    !    close(3)   
                    !    close(4)
                    !    close(5)
                    
                    
                    !open(1,file='vtest1.txt')
                    !    open(2,file='vtest2.txt')
                    !    open(3,file='vtest3.txt')
                    !    open(4,file='vtest4.txt')
                    !    do ik = 1, nk
                    !        write (1,'(4f12.6)') k_grid(ik),ev_spln_coefs(1,ik,44,1,1,1,1,tprint,20),ev_spln_coefs(1,ik,44,1,1,5,5,tprint,20),ev_spln_coefs(1,ik,44,3,3,3,3,tprint,20)
                    !        write (2,'(4f12.6)') k_grid(ik),ev_spln_coefs(2,ik,44,1,1,1,1,tprint,20),ev_spln_coefs(2,ik,44,1,1,5,5,tprint,20),ev_spln_coefs(2,ik,44,3,3,3,3,tprint,20)
                    !        write (3,'(4f12.6)') k_grid(ik),ev_spln_coefs(3,ik,44,1,1,1,1,tprint,20),ev_spln_coefs(3,ik,44,1,1,5,5,tprint,20),ev_spln_coefs(3,ik,44,3,3,3,3,tprint,20)
                    !        write (4,'(4f12.6)') k_grid(ik),ev_spln_coefs(4,ik,44,1,1,1,1,tprint,20),ev_spln_coefs(4,ik,44,1,1,5,5,tprint,20),ev_spln_coefs(4,ik,44,3,3,3,3,tprint,20)
                    !    end do
                    !    close(1)
                    !    close(2)
                    !    close(3)
                    !    close(4)

                
                !STOP
                
                
                        
                !Compute optimal policies for age 2-63

                do it=1,T-2

                    Print *,'t is',T-it

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)

                    do counter=1,nk*na*nu
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

                    do counter = 1, nk*na*nu
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
                    
                    !!$OMP PARALLEL PRIVATE(counter)
                    !!$OMP DO SCHEDULE(DYNAMIC)
                    !
                    !do counter = 1, nfc*na*nu
                    !    call partest4(counter)
                    !end do
                    !
                    !!$OMP END DO    
                    !!$OMP END PARALLEL


                    evm=0d0
                    

                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)

                    do counter=1,nk*na*nu
                        call partest(counter)
                    end do

                    !$OMP END DO    
                    !$OMP END PARALLEL


                    !$OMP PARALLEL PRIVATE(counter)
                    !$OMP DO SCHEDULE(DYNAMIC)

                    do counter=1,nk*na*nu
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
                    
            !if(it==25) then     
            !        tprint=20
            !    
            !    open(1,file='csingle.txt')
            !        open(2,file='ksingle.txt')
            !        open(3,file='nsingle.txt')
            !        open(4,file='vsingle.txt')
            !        open(5,file='vsinglem.txt')
            !        do it2 = 1, nk
            !            write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,5,3,3,tprint,20),cs(2,it2,5,3,3,tprint,20),cs(2,it2,5,5,5,tprint,20)
            !            write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,5,3,3,tprint,20),ks(2,it2,5,3,3,tprint,20),ks(2,it2,5,5,5,tprint,20)
            !            write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,5,3,3,tprint,20),ns(2,it2,5,3,3,tprint,20),ns(2,it2,5,5,5,tprint,20)
            !            write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,5,5,5,tprint,20),evs(2,it2,5,3,3,tprint,20),evs(2,it2,5,5,5,tprint,20)
            !            write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,5,5,5,tprint,20),evm(2,it2,5,3,3,tprint,20),evm(2,it2,5,5,5,tprint,20)
            !        end do
            !        close(1)
            !        close(2)
            !        close(3)
            !        close(4)
            !        close(5)
            !    
            !        open(1,file='ctest.txt')
            !            open(2,file='ktest.txt')
            !            open(3,file='vtest.txt')
            !            open(4,file='nmtest.txt')
            !            open(5,file='nftest.txt')
            !            !open(5,file='evtest.txt')
            !            !open(4,file='ve_vu.txt')
            !            !open(4,file='tbc.txt')
            !            do ik = 1, nk
            !                write (1,'(4f12.6)') k_grid(ik),c(ik,5,1,1,1,1,tprint,20),c(ik,5,1,1,5,5,tprint,20),c(ik,5,3,3,3,3,tprint,20)
            !                write (2,'(4f12.6)') k_grid(ik), k(ik,5,1,1,1,1,tprint,20), k(ik,5,1,1,5,5,tprint,20), k(ik,5,3,3,3,3,tprint,20)
            !                write (3,'(4f12.6)') k_grid(ik),eV(ik,5,1,1,1,1,tprint,20),eV(ik,5,1,1,5,5,tprint,20),eV(ik,5,3,3,3,3,tprint,20)
            !                write (4,'(4f12.6)') k_grid(ik),nm(ik,5,1,1,1,1,tprint,20),nm(ik,5,1,1,5,5,tprint,20),nm(ik,5,3,3,3,3,tprint,20)
            !                write (5,'(4f12.6)') k_grid(ik),nf(ik,5,1,1,1,1,tprint,20),nf(ik,5,1,1,5,5,tprint,20),nf(ik,5,3,3,3,3,tprint,20)
            !            end do
            !            close(1)
            !            close(2)
            !            close(3)   
            !            close(4)
            !            close(5)
            !    
            !      STOP
            !      
            !      end if
                    
end do

                
it=T-1               

Print *,'t is',T-it

!$OMP PARALLEL PRIVATE(counter)
!$OMP DO SCHEDULE(DYNAMIC)

    do counter=1,nk*na*nu
        call Solvefirstactive(counter)
    end do

!$OMP END DO    
!$OMP END PARALLEL

                !tprint=1
                !
                !open(1,file='csingle.txt')
                !    open(2,file='ksingle.txt')
                !    open(3,file='nsingle.txt')
                !    open(4,file='vsingle.txt')
                !    open(5,file='vsinglem.txt')
                !    do it2 = 1, nk
                !        write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,1,3,3,tprint,21),cs(2,it2,1,3,3,tprint,21),cs(2,it2,1,5,5,tprint,21)
                !        write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,1,3,3,tprint,21),ks(2,it2,1,3,3,tprint,21),ks(2,it2,1,5,5,tprint,21)
                !        write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,1,3,3,tprint,21),ns(2,it2,1,3,3,tprint,21),ns(2,it2,1,5,5,tprint,21)
                !        write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,1,3,3,tprint,21),evs(2,it2,1,3,3,tprint,21),evs(2,it2,1,5,5,tprint,21)
                !        write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,1,3,3,tprint,21),evm(2,it2,1,3,3,tprint,21),evm(2,it2,1,5,5,tprint,21)
                !    end do
                !    close(1)
                !    close(2)
                !    close(3)
                !    close(4)
                !    close(5)
                !
                !    open(1,file='ctest.txt')
                !        open(2,file='ktest.txt')
                !        open(3,file='vtest.txt')
                !        open(4,file='nmtest.txt')
                !        open(5,file='nftest.txt')
                !        !open(5,file='evtest.txt')
                !        !open(4,file='ve_vu.txt')
                !        !open(4,file='tbc.txt')
                !        do ik = 1, nk
                !            write (1,'(4f12.6)') k_grid(ik),c(ik,1,1,1,1,1,tprint,21),c(ik,1,1,1,5,5,tprint,21),c(ik,1,3,3,3,3,tprint,21)
                !            write (2,'(4f12.6)') k_grid(ik), k(ik,1,1,1,1,1,tprint,21), k(ik,1,1,1,5,5,tprint,21), k(ik,1,3,3,3,3,tprint,21)
                !            write (3,'(4f12.6)') k_grid(ik),eV(ik,1,1,1,1,1,tprint,21),eV(ik,1,1,1,5,5,tprint,21),eV(ik,1,3,3,3,3,tprint,21)
                !            write (4,'(4f12.6)') k_grid(ik),nm(ik,1,1,1,1,1,tprint,21),nm(ik,1,1,1,5,5,tprint,21),nm(ik,1,3,3,3,3,tprint,21)
                !            write (5,'(4f12.6)') k_grid(ik),nf(ik,1,1,1,1,1,tprint,21),nf(ik,1,1,1,5,5,tprint,21),nf(ik,1,3,3,3,3,tprint,21)
                !        end do
                !        close(1)
                !        close(2)
                !        close(3)
                !        close(4)
                !        close(5)
                !
                !   open(1,file='csingle2.txt')
                !    open(2,file='ksingle2.txt')
                !    open(3,file='nsingle2.txt')
                !    open(4,file='vsingle2.txt')
                !    open(5,file='vsinglem2.txt')
                !    do it2 = 1, nk
                !        write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,1,3,3,tprint,11),cs(2,it2,1,3,3,tprint,11),cs(2,it2,1,5,5,tprint,11)
                !        write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,1,3,3,tprint,11),ks(2,it2,1,3,3,tprint,11),ks(2,it2,1,5,5,tprint,11)
                !        write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,1,3,3,tprint,11),ns(2,it2,1,3,3,tprint,11),ns(2,it2,1,5,5,tprint,11)
                !        write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,1,3,3,tprint,11),evs(2,it2,1,3,3,tprint,11),evs(2,it2,1,5,5,tprint,11)
                !        write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,1,3,3,tprint,11),evm(2,it2,1,3,3,tprint,11),evm(2,it2,1,5,5,tprint,11)
                !    end do
                !    close(1)
                !    close(2)
                !    close(3)
                !    close(4)
                !    close(5)
                !
                !    open(1,file='ctest2.txt')
                !        open(2,file='ktest2.txt')
                !        open(3,file='vtest2.txt')
                !        open(4,file='nmtest2.txt')
                !        open(5,file='nftest2.txt')
                !        !open(5,file='evtest.txt')
                !        !open(4,file='ve_vu.txt')
                !        !open(4,file='tbc.txt')
                !        do ik = 1, nk
                !            write (1,'(4f12.6)') k_grid(ik),c(ik,1,1,1,1,1,tprint,11),c(ik,1,1,1,5,5,tprint,11),c(ik,1,3,3,3,3,tprint,11)
                !            write (2,'(4f12.6)') k_grid(ik), k(ik,1,1,1,1,1,tprint,11), k(ik,1,1,1,5,5,tprint,11), k(ik,1,3,3,3,3,tprint,11)
                !            write (3,'(4f12.6)') k_grid(ik),eV(ik,1,1,1,1,1,tprint,11),eV(ik,1,1,1,5,5,tprint,11),eV(ik,1,3,3,3,3,tprint,11)
                !            write (4,'(4f12.6)') k_grid(ik),nm(ik,1,1,1,1,1,tprint,11),nm(ik,1,1,1,5,5,tprint,11),nm(ik,1,3,3,3,3,tprint,11)
                !            write (5,'(4f12.6)') k_grid(ik),nf(ik,1,1,1,1,1,tprint,11),nf(ik,1,1,1,5,5,tprint,11),nf(ik,1,3,3,3,3,tprint,11)
                !        end do
                !        close(1)
                !        close(2)
                !        close(3)
                !        close(4)
                !        close(5)     
                !
                !open(1,file='csingle3.txt')
                !    open(2,file='ksingle3.txt')
                !    open(3,file='nsingle3.txt')
                !    open(4,file='vsingle3.txt')
                !    open(5,file='vsinglem3.txt')
                !    do it2 = 1, nk
                !        write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,1,3,3,tprint,1),cs(2,it2,1,3,3,tprint,1),cs(2,it2,1,5,5,tprint,1)
                !        write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,1,3,3,tprint,1),ks(2,it2,1,3,3,tprint,1),ks(2,it2,1,5,5,tprint,1)
                !        write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,1,3,3,tprint,1),ns(2,it2,1,3,3,tprint,1),ns(2,it2,1,5,5,tprint,1)
                !        write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,1,3,3,tprint,1),evs(2,it2,1,3,3,tprint,1),evs(2,it2,1,5,5,tprint,1)
                !        write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,1,3,3,tprint,1),evm(2,it2,1,3,3,tprint,1),evm(2,it2,1,5,5,tprint,1)
                !    end do
                !    close(1)
                !    close(2)
                !    close(3)
                !    close(4)
                !    close(5)
                !
                !    open(1,file='ctest3.txt')
                !        open(2,file='ktest3.txt')
                !        open(3,file='vtest3.txt')
                !        open(4,file='nmtest3.txt')
                !        open(5,file='nftest3.txt')
                !        !open(5,file='evtest.txt')
                !        !open(4,file='ve_vu.txt')
                !        !open(4,file='tbc.txt')
                !        do ik = 1, nk
                !            write (1,'(4f12.6)') k_grid(ik),c(ik,1,1,1,1,1,tprint,1),c(ik,1,1,1,5,5,tprint,1),c(ik,1,3,3,3,3,tprint,1)
                !            write (2,'(4f12.6)') k_grid(ik), k(ik,1,1,1,1,1,tprint,1), k(ik,1,1,1,5,5,tprint,1), k(ik,1,3,3,3,3,tprint,1)
                !            write (3,'(4f12.6)') k_grid(ik),eV(ik,1,1,1,1,1,tprint,1),eV(ik,1,1,1,5,5,tprint,1),eV(ik,1,3,3,3,3,tprint,1)
                !            write (4,'(4f12.6)') k_grid(ik),nm(ik,1,1,1,1,1,tprint,1),nm(ik,1,1,1,5,5,tprint,1),nm(ik,1,3,3,3,3,tprint,1)
                !            write (5,'(4f12.6)') k_grid(ik),nf(ik,1,1,1,1,1,tprint,1),nf(ik,1,1,1,5,5,tprint,1),nf(ik,1,3,3,3,3,tprint,1)
                !        end do
                !        close(1)
                !        close(2)
                !        close(3)
                !        close(4)
                !        close(5)     
                !        
                !  STOP                      
                
                !open(1, file='policies.txt')
                !
                !read (1, *) c,k,nm,nf,cs,ks,ns,c_ret,k_ret,cs_ret,ks_ret
                !
                !close(1)
                !STOP
                !$OMP PARALLEL PRIVATE(ik)
                !$OMP DO SCHEDULE(DYNAMIC)

                do ik=1,nsim2
                    call Simulation(ik)
                end do

                !$OMP END DO    
                !$OMP END PARALLEL

                call Statistics

            end do

            epsilon4=ratio-ratiodum

            Print *,'epsilon4 is',epsilon4

            ratio=ratio-0.2d0*(ratio-ratiodum)

            w=(1d0-alpha)*ratio**alpha
            r=alpha*ratio**(alpha-1d0)-delta

        end do

        open(1, file='singledist.txt')
        
            write (1, *) fpartner,mpartner
        
        close(1)
        
        open(1, file='abilityprob.txt')
            
            write (1, *) ability_prob
            
        close(1)

        !tprint=1
        !        
        !        open(1,file='csingle.txt')
        !            open(2,file='ksingle.txt')
        !            open(3,file='nsingle.txt')
        !            open(4,file='vsingle.txt')
        !            open(5,file='vsinglem.txt')
        !            do it2 = 1, nk
        !                write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,1,3,3,tprint),cs(2,it2,1,3,3,tprint),cs(2,it2,1,5,5,tprint)
        !                write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,1,3,3,tprint),ks(2,it2,1,3,3,tprint),ks(2,it2,1,5,5,tprint)
        !                write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,1,3,3,tprint),ns(2,it2,1,3,3,tprint),ns(2,it2,1,5,5,tprint)
        !                write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,1,3,3,tprint),evs(2,it2,1,3,3,tprint),evs(2,it2,1,5,5,tprint)
        !                write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,1,3,3,tprint),evm(2,it2,1,3,3,tprint),evm(2,it2,1,5,5,tprint)
        !            end do
        !            close(1)
        !            close(2)
        !            close(3)
        !            close(4)
        !            close(5)
        !        
        !            open(1,file='ctest.txt')
        !                open(2,file='ktest.txt')
        !                open(3,file='vtest.txt')
        !                open(4,file='nmtest.txt')
        !                open(5,file='nftest.txt')
        !                !open(5,file='evtest.txt')
        !                !open(4,file='ve_vu.txt')
        !                !open(4,file='tbc.txt')
        !                do ik = 1, nk
        !                    write (1,'(4f12.6)') k_grid(ik),c(ik,1,1,1,1,1,tprint,3),c(ik,1,1,1,5,5,tprint,3),c(ik,1,3,3,3,3,tprint,3)
        !                    write (2,'(4f12.6)') k_grid(ik), k(ik,1,1,1,1,1,tprint,3), k(ik,1,1,1,5,5,tprint,3), k(ik,1,3,3,3,3,tprint,3)
        !                    write (3,'(4f12.6)') k_grid(ik),eV(ik,1,1,1,1,1,tprint,3),eV(ik,1,1,1,5,5,tprint,3),eV(ik,1,3,3,3,3,tprint,3)
        !                    write (4,'(4f12.6)') k_grid(ik),nm(ik,1,1,1,1,1,tprint,3),nm(ik,1,1,1,5,5,tprint,3),nm(ik,1,3,3,3,3,tprint,3)
        !                    write (5,'(4f12.6)') k_grid(ik),nf(ik,1,1,1,1,1,tprint,3),nf(ik,1,1,1,5,5,tprint,3),nf(ik,1,3,3,3,3,tprint,3)
        !                end do
        !                close(1)
        !                close(2)
        !                close(3)   
        !                close(4)
        !                close(5)
        
        
         !Variables are age, gender, ID number, weight, marital status, asset holdings, household labor income, Household_Labor_Income_Tax_Paid,  Household_consumption_Tax_Paid

        !open(1, file='Simulation_output.txt')
        !
        !    do it2=1,T
        !        do it3=1,nsim2
        !            do it4=1,nsim
        !                write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') it2*1d0, 1d0, nsim*(it3-1)*1d0+it4*1d0, 1d0, Sim1m(it3,it4,it2,10), Sim1m(it3,it4,it2,1), Sim1m(it3,it4,it2,6), Sim1m(it3,it4,it2,7), Sim1m(it3,it4,it2,8)
        !                
        !                write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') it2*1d0, 2d0, nsim*(it3-1)*1d0+it4*1d0, 1d0, Sim1f(it3,it4,it2,10), Sim1f(it3,it4,it2,1), Sim1f(it3,it4,it2,6), Sim1f(it3,it4,it2,7), Sim1f(it3,it4,it2,8)
        !            end do
        !        end do
        !    end do
        
            !do it2=1,Tret
            !    ium=int(weightret(it2)*nsim)
            !    do it3=1,nsim2
            !        do it4=1,ium
            !            write (1,'(F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4,F12.4)') (45+it2)*1d0, ium*(it3-1)*1d0+it4*1d0, 1d0, SimR1(it3,it4,it2,1), SimR1(it3,it4,it2,2), 0d0, 0d0, 0d0, 0d0
            !        end do
            !    end do
            !end do
            
        !close(1)
        
        open(1,file='cpathm.txt')
        
        do it2=1,T
            dum2=0.0
            do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1m(it4,it3,it2,2)
            end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        do it2=1,Tret
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+SimR1m(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        open(1,file='cpathf.txt')
        
        do it2=1,T
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+Sim1f(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        do it2=1,Tret
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+SimR1f(it4,it3,it2,2)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        open(1,file='kpathm.txt')
        
        do it2=1,T
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+Sim1m(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        do it2=1,Tret
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+SimR1m(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        open(1,file='kpathf.txt')
        
        do it2=1,T
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+Sim1f(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        do it2=1,Tret
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+SimR1f(it4,it3,it2,1)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        open(1,file='npathm.txt')
        
        do it2=1,T
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+Sim1m(it4,it3,it2,4)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        open(1,file='npathf.txt')
        
        do it2=1,T
            dum2=0.0
            do it4=1,16
                do it3=1,10000
                    dum2=dum2+Sim1f(it4,it3,it2,4)
                end do
            end do
            dum2=dum2/160000d0
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        open(1,file='lfppathsingle.txt')
        
        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,16
                do it3=1,10000
                if(Sim1f(it4,it3,it2,10)<0.5d0) then
                        dum3=dum3+1d0
                    if(Sim1f(it4,it3,it2,4)>0.001d0) then
                        dum2=dum2+1d0
                    end if
                end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        open(1,file='lfppathmarried.txt')
        
        do it2=1,T
            dum2=0.0d0
            dum3=0.0d0
            do it4=1,16
                do it3=1,10000
                if(Sim1f(it4,it3,it2,10)>0.5d0) then
                        dum3=dum3+1d0
                    if(Sim1f(it4,it3,it2,4)>0.001d0) then
                        dum2=dum2+1d0
                    end if
                end if
                end do
            end do
            dum2=dum2/dum3
            write (1,'(F8.3,F8.3)') it2*1d0, dum2
        end do
        
        close(1)
        
        contains
        
        subroutine Initialize()
            
            !USE ANORDF_INT
            implicit none
            
            allocate(v(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(ev(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(c(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(k(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(nm(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(nf(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(ev_spln_coefs(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(v_spln_coefs(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(vdum(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(cdum(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(gkdum(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(nmdum(nk,nexp,na,nu,na,nu,T,nfc))
            allocate(nfdum(nk,nexp,na,nu,na,nu,T,nfc))
            
            allocate(vs(2,nk,nexp,na,nu,T,nfc))
            allocate(evs(2,nk,nexp,na,nu,T,nfc))
            allocate(evm(2,nk,nexp,na,nu,T,nfc))
            allocate(cs(2,nk,nexp,na,nu,T,nfc))
            allocate(ks(2,nk,nexp,na,nu,T,nfc))
            allocate(ns(2,nk,nexp,na,nu,T,nfc))
            allocate(evs_spln_coefs(2,nk,nexp,na,nu,T,nfc))
            allocate(vs_spln_coefs(2,nk,nexp,na,nu,T,nfc))
            allocate(evm_spln_coefs(2,nk,nexp,na,nu,T,nfc))
            !allocate(vsdum(2,nk,nexp,na,nu,T,nfc))
            !allocate(csdum(2,nk,nexp,na,nu,T,nfc))
            !allocate(gksdum(2,nk,nexp,na,nu,T,nfc))
            !allocate(nsdum(2,nk,nexp,na,nu,T,nfc))
            
            allocate(Sim1m(nsim2,nsim,46,10))
            allocate(Sim1f(nsim2,nsim,46,10))
            allocate(exp1m(nsim2,nsim,46,6))
            allocate(exp1f(nsim2,nsim,46,6))
            allocate(exp2m(nsim2,nsim,46,6))
            allocate(exp2f(nsim2,nsim,46,6))
            allocate(SimR1m(nsim2,nsim,37,4))
            allocate(SimR1f(nsim2,nsim,37,4))
            allocate(Random3m(nsim2,nsim,45))
            allocate(Random3f(nsim2,nsim,45))
            allocate(marstatm(nsim2,nsim,45))
            allocate(marstatf(nsim2,nsim,45))
            allocate(marstatm_init(nsim2,nsim))
            allocate(marstatf_init(nsim2,nsim))
            allocate(partshock(nsim2,nsim,1))
            allocate(Random1m(nsim2,nsim))
            allocate(Random1f(nsim2,nsim))
            allocate(Random2m(nsim2,nsim))
            allocate(Random2f(nsim2,nsim))
            
            allocate(a(2,na))
            allocate(Prob_a(2,na))
            allocate(fc(2,nfc))
            allocate(Prob_fc(2,nfc))
            allocate(u(2,nu))
            allocate(Prob_u(2,nu))
            allocate(trans_u(2,nu,nu))
            allocate(trans_a(2,na,na))
            allocate(trans_fc(2,nfc,nfc))
            allocate(OmegaRet(Tret))
            allocate(OmegaActive(T))
            allocate(Probm(T))
            allocate(Probd(T))
            allocate(WeightRet(Tret))
            allocate(WeightActive(T))
            
            allocate(fpartner(nk,T,na,nu,T,nfc))
            allocate(mpartner(nk,na,nu,T))
            allocate(fpartnerdum(nk,T,na,nu,T,nfc))
            allocate(mpartnerdum(nk,na,nu,T))
            allocate(fpartnerdum2(nk,T,na,nu,T,nfc))
            allocate(mpartnerdum2(nk,na,nu,T))
            allocate(ability_prob(na,na))
            allocate(laborm(nc,nw,nw))
            allocate(laborf(nc,nw,nw))
            allocate(labormwork(nc,nw))
            allocate(laborsinglem(nc,nw))
            allocate(laborsinglef(nc,nw))
            
            allocate(c_grid(nc))
            allocate(wage_grid(nw))
            allocate(k_grid(nk))
            allocate(exp_grid(nexp,T))
            allocate(INTERP2D(nk,nexp))
            allocate(exp_grid_dum(nexp))
            allocate(K_KNOT(nk+KORDER))
            allocate(EXP_KNOT(nexp+EXPORDER,T))
            allocate(ev_spln_coefs_ret(4,nk,Tret))
            ! for testing only
            allocate(p_ev_spln_coefs_ret(4,nk,Tret))

            allocate(evs_spln_coefs_ret(4,nk,Tret))
            ! for testing only
            allocate(p_evs_spln_coefs_ret(4,nk,Tret))

            allocate(c_ret(nk,Tret))
            allocate(v_ret(nk,Tret))
            allocate(k_ret(nk,Tret))
            allocate(vs_ret(nk,Tret))
            allocate(cs_ret(nk,Tret))
            allocate(ks_ret(nk,Tret))
            allocate(break(nk))
            
            open(1, file='singledist.txt')
            
            read (1, *) fpartner,mpartner
            
            close(1)
            
            open(1, file='abilityprob.txt')
            
                read (1, *) ability_prob
            
            close(1)
            
            
            !Print *,ability_prob(5,:)
            !STOP
            !ability_prob=1d0/5d0
            !
            !fpartner=0d0
            !mpartner=0d0
            !
            !do it=1,T
            !    do ik=1,8
            !        mpartner(ik,:,:,it)=1d0/(8*na*nu)
            !    end do
            !end do
            !
            !do it=1,T
            !    do it2=1,it
            !        do ik=1,8
            !            fpartner(ik,it2,:,:,it,:)=1d0/(8*it*na*nu*nfc)
            !        end do
            !    end do
            !end do
            
            trans_u = 0d0
            prob_u=0d0
            trans_a = 0d0
            prob_a=0d0
            gamma(1,:) = (/ 0.1088524d0, -0.001473d0, 0.00000634d0 /)
            gamma(2,:) = (/ 0.0784408d0, -0.0025596d0, 0.0000256d0 /)
            gamma0=0.4938d0
            
            theta(:) = (/ 0.93124354*tax_level_scale, 0.15002363*tax_prog_scale /)
            thetas(:) = (/ 0.81773322*tax_level_scale, 0.11060017*tax_prog_scale /)
            AE = 1.04354601398297d0
            Unemp_benefit=0.201795*AE
            r=alpha*ratio**(alpha-1d0)-delta
            w=(1d0-alpha)*ratio**alpha
            call MakeGrid(nk,k_grid,0d0,60d0,3d0)
            call MakeGrid(nc,c_grid,0.01d0,60d0,3d0)
            
            exp_grid=0d0
            do it2=2,T
                call MakeGrid(nexp,exp_grid(:,it2),0d0,1d0*(it2-1),2d0)
            end do
            
            CALL d_BSNAK(nk, k_grid, KORDER, K_KNOT)
            
            do it2=2,T
                
                CALL d_BSNAK(nexp, exp_grid(:,it2), EXPORDER, EXP_KNOT(:,it2))
                
            end do
            
            exp_grid(:,1)=exp_grid(:,2)
            EXP_KNOT(:,1)=EXP_KNOT(:,2)
            
            !Print *,exp_grid(:,2)
            
            !STOP
            
    !Filling in US divorce and marriage probabilities

    probd(1)=0.11908938
    probd(2)=0.09947980
    probd(3)=0.08335590
    probd(4)=0.07024181
    probd(5)=0.05970394
    probd(6)=0.05134912
    probd(7)=0.04482276
    probd(8)=0.03980704
    probd(9)=0.03601902
    probd(10)=0.03320882
    probd(11)=0.03115778
    probd(12)=0.02967661
    probd(13)=0.02860353
    probd(14)=0.02780249
    probd(15)=0.02716124
    probd(16)=0.02658955
    probd(17)=0.02601734
    probd(18)=0.02539286
    probd(19)=0.02468080
    probd(20)=0.02386051
    probd(21)=0.02292411
    probd(22)=0.02187467
    probd(23)=0.02072434
    probd(24)=0.01949255
    probd(25)=0.01820413
    probd(26)=0.01688749
    probd(27)=0.01557275
    probd(28)=0.01428994
    probd(29)=0.01306712
    probd(30)=0.01192853
    probd(31)=0.01089280
    probd(32)=0.00997105
    probd(33)=0.00916507
    probd(34)=0.00846549
    probd(35)=0.00784992
    probd(36)=0.00728109
    probd(37)=0.00670507
    probd(38)=0.00604934
    probd(39)=0.00522102
    probd(40)=0.00410500
    probd(41)=0.00256208
    probd(42)=0.00042715
    probd(43)=0.00000000
    probd(44)=0.00000000
    probd(45)=0.00000000
    
    probm(1)=0.08301004
    probm(2)=0.09522406
    probm(3)=0.10423554
    probm(4)=0.11041682
    probm(5)=0.11411954
    probm(6)=0.11567485
    probm(7)=0.11539366
    probm(8)=0.11356688
    probm(9)=0.11046570
    probm(10)=0.10634176
    probm(11)=0.10142749
    probm(12)=0.09593630
    probm(13)=0.09006280
    probm(14)=0.08398313
    probm(15)=0.07785512
    probm(16)=0.07181858
    probm(17)=0.06599554
    probm(18)=0.06049049
    probm(19)=0.05539063
    probm(20)=0.05076610
    probm(21)=0.04667024
    probm(22)=0.04313985
    probm(23)=0.04019541
    probm(24)=0.03784131
    probm(25)=0.03606616
    probm(26)=0.03484296
    probm(27)=0.03412940
    probm(28)=0.03386807
    probm(29)=0.03398674
    probm(30)=0.03439857
    probm(31)=0.03500238
    probm(32)=0.03568288
    probm(33)=0.03631093
    probm(34)=0.03674376
    probm(35)=0.03682526
    probm(36)=0.03638617
    probm(37)=0.03524438
    probm(38)=0.03320513
    probm(39)=0.03006129
    probm(40)=0.02559356
    probm(41)=0.01957079
    probm(42)=0.01175015
    probm(43)=0.00187740
    probm(44)=0.00000000
    probm(45)=0.00000000
    
            OmegaActive=1d0
            
            OmegaRet(1)=1d0-0.014319d0
            OmegaRet(2)=1d0-0.015540d0
            OmegaRet(3)=1d0-0.016920d0
            OmegaRet(4)=1d0-0.018448d0
            OmegaRet(5)=1d0-0.020170d0
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
            
    U(1,1)=-(3D0/2D0)*SQRT((sigma_um**2)/(1.0-(rho_um**2)))
    U(1,2)=-(3D0/4D0)*SQRT((sigma_um**2)/(1.0-(rho_um**2)))
    U(1,3)=0D0
    U(1,4)=(3D0/4D0)*SQRT((sigma_um**2)/(1.0-(rho_um**2)))
    U(1,5)=(3D0/2D0)*SQRT((sigma_um**2)/(1.0-(rho_um**2)))
    dum=(5D-1)*(3D0/4D0)*SQRT((sigma_um**2)/(1.0-(rho_um**2)))
    
    trans_u(1,1,1)=D_ANORDF((u(1,1)-rho_um*u(1,1)+dum)/sigma_um)
    trans_u(1,2,1)=D_ANORDF((u(1,1)-rho_um*u(1,2)+dum)/sigma_um)
    trans_u(1,3,1)=D_ANORDF((u(1,1)-rho_um*u(1,3)+dum)/sigma_um)
    trans_u(1,4,1)=D_ANORDF((u(1,1)-rho_um*u(1,4)+dum)/sigma_um)
    trans_u(1,5,1)=D_ANORDF((u(1,1)-rho_um*u(1,5)+dum)/sigma_um)
    trans_u(1,1,2)=D_ANORDF((u(1,2)-rho_um*u(1,1)+dum)/sigma_um)-D_ANORDF((u(1,2)-rho_um*u(1,1)-dum)/sigma_um)
    trans_u(1,2,2)=D_ANORDF((u(1,2)-rho_um*u(1,2)+dum)/sigma_um)-D_ANORDF((u(1,2)-rho_um*u(1,2)-dum)/sigma_um)
    trans_u(1,3,2)=D_ANORDF((u(1,2)-rho_um*u(1,3)+dum)/sigma_um)-D_ANORDF((u(1,2)-rho_um*u(1,3)-dum)/sigma_um)
    trans_u(1,4,2)=D_ANORDF((u(1,2)-rho_um*u(1,4)+dum)/sigma_um)-D_ANORDF((u(1,2)-rho_um*u(1,4)-dum)/sigma_um)
    trans_u(1,5,2)=D_ANORDF((u(1,2)-rho_um*u(1,5)+dum)/sigma_um)-D_ANORDF((u(1,2)-rho_um*u(1,5)-dum)/sigma_um)
    trans_u(1,1,3)=D_ANORDF((u(1,3)-rho_um*u(1,1)+dum)/sigma_um)-D_ANORDF((u(1,3)-rho_um*u(1,1)-dum)/sigma_um)
    trans_u(1,2,3)=D_ANORDF((u(1,3)-rho_um*u(1,2)+dum)/sigma_um)-D_ANORDF((u(1,3)-rho_um*u(1,2)-dum)/sigma_um)
    trans_u(1,3,3)=D_ANORDF((u(1,3)-rho_um*u(1,3)+dum)/sigma_um)-D_ANORDF((u(1,3)-rho_um*u(1,3)-dum)/sigma_um)
    trans_u(1,4,3)=D_ANORDF((u(1,3)-rho_um*u(1,4)+dum)/sigma_um)-D_ANORDF((u(1,3)-rho_um*u(1,4)-dum)/sigma_um)
    trans_u(1,5,3)=D_ANORDF((u(1,3)-rho_um*u(1,5)+dum)/sigma_um)-D_ANORDF((u(1,3)-rho_um*u(1,5)-dum)/sigma_um)
    trans_u(1,1,4)=D_ANORDF((u(1,4)-rho_um*u(1,1)+dum)/sigma_um)-D_ANORDF((u(1,4)-rho_um*u(1,1)-dum)/sigma_um)
    trans_u(1,2,4)=D_ANORDF((u(1,4)-rho_um*u(1,2)+dum)/sigma_um)-D_ANORDF((u(1,4)-rho_um*u(1,2)-dum)/sigma_um)
    trans_u(1,3,4)=D_ANORDF((u(1,4)-rho_um*u(1,3)+dum)/sigma_um)-D_ANORDF((u(1,4)-rho_um*u(1,3)-dum)/sigma_um)
    trans_u(1,4,4)=D_ANORDF((u(1,4)-rho_um*u(1,4)+dum)/sigma_um)-D_ANORDF((u(1,4)-rho_um*u(1,4)-dum)/sigma_um)
    trans_u(1,5,4)=D_ANORDF((u(1,4)-rho_um*u(1,5)+dum)/sigma_um)-D_ANORDF((u(1,4)-rho_um*u(1,5)-dum)/sigma_um)
    trans_u(1,1,5)=1D0-D_ANORDF((u(1,5)-rho_um*u(1,1)-dum)/sigma_um)
    trans_u(1,2,5)=1D0-D_ANORDF((u(1,5)-rho_um*u(1,2)-dum)/sigma_um)
    trans_u(1,3,5)=1D0-D_ANORDF((u(1,5)-rho_um*u(1,3)-dum)/sigma_um)
    trans_u(1,4,5)=1D0-D_ANORDF((u(1,5)-rho_um*u(1,4)-dum)/sigma_um)
    trans_u(1,5,5)=1D0-D_ANORDF((u(1,5)-rho_um*u(1,5)-dum)/sigma_um)
    
    U(2,1)=-(3D0/2D0)*SQRT((sigma_uf**2)/(1.0-(rho_uf**2)))
    U(2,2)=-(3D0/4D0)*SQRT((sigma_uf**2)/(1.0-(rho_uf**2)))
    U(2,3)=0D0
    U(2,4)=(3D0/4D0)*SQRT((sigma_uf**2)/(1.0-(rho_uf**2)))
    U(2,5)=(3D0/2D0)*SQRT((sigma_uf**2)/(1.0-(rho_uf**2)))
    dum=(5D-1)*(3D0/4D0)*SQRT((sigma_uf**2)/(1.0-(rho_uf**2)))
    
    trans_u(2,1,1)=D_ANORDF((u(2,1)-rho_uf*u(2,1)+dum)/sigma_uf)
    trans_u(2,2,1)=D_ANORDF((u(2,1)-rho_uf*u(2,2)+dum)/sigma_uf)
    trans_u(2,3,1)=D_ANORDF((u(2,1)-rho_uf*u(2,3)+dum)/sigma_uf)
    trans_u(2,4,1)=D_ANORDF((u(2,1)-rho_uf*u(2,4)+dum)/sigma_uf)
    trans_u(2,5,1)=D_ANORDF((u(2,1)-rho_uf*u(2,5)+dum)/sigma_uf)
    trans_u(2,1,2)=D_ANORDF((u(2,2)-rho_uf*u(2,1)+dum)/sigma_uf)-D_ANORDF((u(2,2)-rho_uf*u(2,1)-dum)/sigma_uf)
    trans_u(2,2,2)=D_ANORDF((u(2,2)-rho_uf*u(2,2)+dum)/sigma_uf)-D_ANORDF((u(2,2)-rho_uf*u(2,2)-dum)/sigma_uf)
    trans_u(2,3,2)=D_ANORDF((u(2,2)-rho_uf*u(2,3)+dum)/sigma_uf)-D_ANORDF((u(2,2)-rho_uf*u(2,3)-dum)/sigma_uf)
    trans_u(2,4,2)=D_ANORDF((u(2,2)-rho_uf*u(2,4)+dum)/sigma_uf)-D_ANORDF((u(2,2)-rho_uf*u(2,4)-dum)/sigma_uf)
    trans_u(2,5,2)=D_ANORDF((u(2,2)-rho_uf*u(2,5)+dum)/sigma_uf)-D_ANORDF((u(2,2)-rho_uf*u(2,5)-dum)/sigma_uf)
    trans_u(2,1,3)=D_ANORDF((u(2,3)-rho_uf*u(2,1)+dum)/sigma_uf)-D_ANORDF((u(2,3)-rho_uf*u(2,1)-dum)/sigma_uf)
    trans_u(2,2,3)=D_ANORDF((u(2,3)-rho_uf*u(2,2)+dum)/sigma_uf)-D_ANORDF((u(2,3)-rho_uf*u(2,2)-dum)/sigma_uf)
    trans_u(2,3,3)=D_ANORDF((u(2,3)-rho_uf*u(2,3)+dum)/sigma_uf)-D_ANORDF((u(2,3)-rho_uf*u(2,3)-dum)/sigma_uf)
    trans_u(2,4,3)=D_ANORDF((u(2,3)-rho_uf*u(2,4)+dum)/sigma_uf)-D_ANORDF((u(2,3)-rho_uf*u(2,4)-dum)/sigma_uf)
    trans_u(2,5,3)=D_ANORDF((u(2,3)-rho_uf*u(2,5)+dum)/sigma_uf)-D_ANORDF((u(2,3)-rho_uf*u(2,5)-dum)/sigma_uf)
    trans_u(2,1,4)=D_ANORDF((u(2,4)-rho_uf*u(2,1)+dum)/sigma_uf)-D_ANORDF((u(2,4)-rho_uf*u(2,1)-dum)/sigma_uf)
    trans_u(2,2,4)=D_ANORDF((u(2,4)-rho_uf*u(2,2)+dum)/sigma_uf)-D_ANORDF((u(2,4)-rho_uf*u(2,2)-dum)/sigma_uf)
    trans_u(2,3,4)=D_ANORDF((u(2,4)-rho_uf*u(2,3)+dum)/sigma_uf)-D_ANORDF((u(2,4)-rho_uf*u(2,3)-dum)/sigma_uf)
    trans_u(2,4,4)=D_ANORDF((u(2,4)-rho_uf*u(2,4)+dum)/sigma_uf)-D_ANORDF((u(2,4)-rho_uf*u(2,4)-dum)/sigma_uf)
    trans_u(2,5,4)=D_ANORDF((u(2,4)-rho_uf*u(2,5)+dum)/sigma_uf)-D_ANORDF((u(2,4)-rho_uf*u(2,5)-dum)/sigma_uf)
    trans_u(2,1,5)=1D0-D_ANORDF((u(2,5)-rho_uf*u(2,1)-dum)/sigma_uf)
    trans_u(2,2,5)=1D0-D_ANORDF((u(2,5)-rho_uf*u(2,2)-dum)/sigma_uf)
    trans_u(2,3,5)=1D0-D_ANORDF((u(2,5)-rho_uf*u(2,3)-dum)/sigma_uf)
    trans_u(2,4,5)=1D0-D_ANORDF((u(2,5)-rho_uf*u(2,4)-dum)/sigma_uf)
    trans_u(2,5,5)=1D0-D_ANORDF((u(2,5)-rho_uf*u(2,5)-dum)/sigma_uf)
    
    prob_u=1D0/5D0
    
    do it=1,100
        prob_u(1,:)=MATMUL(prob_u(1,:),trans_u(1,:,:))
    end do
    
    do it=1,100
        prob_u(2,:)=MATMUL(prob_u(2,:),trans_u(2,:,:))
    end do

    A(1,1)=-(3D0/2D0)*SQRT((sigma_am**2)/(1.0-(rho_am**2)))
    A(1,2)=-(3D0/4D0)*SQRT((sigma_am**2)/(1.0-(rho_am**2)))
    A(1,3)=0D0
    A(1,4)=(3D0/4D0)*SQRT((sigma_am**2)/(1.0-(rho_am**2)))
    A(1,5)=(3D0/2D0)*SQRT((sigma_am**2)/(1.0-(rho_am**2)))
    dum=(5D-1)*(3D0/4D0)*SQRT((sigma_am**2)/(1.0-(rho_am**2)))
    
    trans_a(1,1,1)=D_ANORDF((a(1,1)-rho_am*a(1,1)+dum)/sigma_am)
    trans_a(1,2,1)=D_ANORDF((a(1,1)-rho_am*a(1,2)+dum)/sigma_am)
    trans_a(1,3,1)=D_ANORDF((a(1,1)-rho_am*a(1,3)+dum)/sigma_am)
    trans_a(1,4,1)=D_ANORDF((a(1,1)-rho_am*a(1,4)+dum)/sigma_am)
    trans_a(1,5,1)=D_ANORDF((a(1,1)-rho_am*a(1,5)+dum)/sigma_am)
    trans_a(1,1,2)=D_ANORDF((a(1,2)-rho_am*a(1,1)+dum)/sigma_am)-D_ANORDF((a(1,2)-rho_am*a(1,1)-dum)/sigma_am)
    trans_a(1,2,2)=D_ANORDF((a(1,2)-rho_am*a(1,2)+dum)/sigma_am)-D_ANORDF((a(1,2)-rho_am*a(1,2)-dum)/sigma_am)
    trans_a(1,3,2)=D_ANORDF((a(1,2)-rho_am*a(1,3)+dum)/sigma_am)-D_ANORDF((a(1,2)-rho_am*a(1,3)-dum)/sigma_am)
    trans_a(1,4,2)=D_ANORDF((a(1,2)-rho_am*a(1,4)+dum)/sigma_am)-D_ANORDF((a(1,2)-rho_am*a(1,4)-dum)/sigma_am)
    trans_a(1,5,2)=D_ANORDF((a(1,2)-rho_am*a(1,5)+dum)/sigma_am)-D_ANORDF((a(1,2)-rho_am*a(1,5)-dum)/sigma_am)
    trans_a(1,1,3)=D_ANORDF((a(1,3)-rho_am*a(1,1)+dum)/sigma_am)-D_ANORDF((a(1,3)-rho_am*a(1,1)-dum)/sigma_am)
    trans_a(1,2,3)=D_ANORDF((a(1,3)-rho_am*a(1,2)+dum)/sigma_am)-D_ANORDF((a(1,3)-rho_am*a(1,2)-dum)/sigma_am)
    trans_a(1,3,3)=D_ANORDF((a(1,3)-rho_am*a(1,3)+dum)/sigma_am)-D_ANORDF((a(1,3)-rho_am*a(1,3)-dum)/sigma_am)
    trans_a(1,4,3)=D_ANORDF((a(1,3)-rho_am*a(1,4)+dum)/sigma_am)-D_ANORDF((a(1,3)-rho_am*a(1,4)-dum)/sigma_am)
    trans_a(1,5,3)=D_ANORDF((a(1,3)-rho_am*a(1,5)+dum)/sigma_am)-D_ANORDF((a(1,3)-rho_am*a(1,5)-dum)/sigma_am)
    trans_a(1,1,4)=D_ANORDF((a(1,4)-rho_am*a(1,1)+dum)/sigma_am)-D_ANORDF((a(1,4)-rho_am*a(1,1)-dum)/sigma_am)
    trans_a(1,2,4)=D_ANORDF((a(1,4)-rho_am*a(1,2)+dum)/sigma_am)-D_ANORDF((a(1,4)-rho_am*a(1,2)-dum)/sigma_am)
    trans_a(1,3,4)=D_ANORDF((a(1,4)-rho_am*a(1,3)+dum)/sigma_am)-D_ANORDF((a(1,4)-rho_am*a(1,3)-dum)/sigma_am)
    trans_a(1,4,4)=D_ANORDF((a(1,4)-rho_am*a(1,4)+dum)/sigma_am)-D_ANORDF((a(1,4)-rho_am*a(1,4)-dum)/sigma_am)
    trans_a(1,5,4)=D_ANORDF((a(1,4)-rho_am*a(1,5)+dum)/sigma_am)-D_ANORDF((a(1,4)-rho_am*a(1,5)-dum)/sigma_am)
    trans_a(1,1,5)=1D0-D_ANORDF((a(1,5)-rho_am*a(1,1)-dum)/sigma_am)
    trans_a(1,2,5)=1D0-D_ANORDF((a(1,5)-rho_am*a(1,2)-dum)/sigma_am)
    trans_a(1,3,5)=1D0-D_ANORDF((a(1,5)-rho_am*a(1,3)-dum)/sigma_am)
    trans_a(1,4,5)=1D0-D_ANORDF((a(1,5)-rho_am*a(1,4)-dum)/sigma_am)
    trans_a(1,5,5)=1D0-D_ANORDF((a(1,5)-rho_am*a(1,5)-dum)/sigma_am)

    a(2,1)=-(3D0/2D0)*SQRT((sigma_af**2)/(1.0-(rho_af**2)))
    a(2,2)=-(3D0/4D0)*SQRT((sigma_af**2)/(1.0-(rho_af**2)))
    a(2,3)=0D0
    a(2,4)=(3D0/4D0)*SQRT((sigma_af**2)/(1.0-(rho_af**2)))
    a(2,5)=(3D0/2D0)*SQRT((sigma_af**2)/(1.0-(rho_af**2)))
    dum=(5D-1)*(3D0/4D0)*SQRT((sigma_af**2)/(1.0-(rho_af**2)))
    
    trans_a(2,1,1)=D_ANORDF((a(2,1)-rho_af*a(2,1)+dum)/sigma_af)
    trans_a(2,2,1)=D_ANORDF((a(2,1)-rho_af*a(2,2)+dum)/sigma_af)
    trans_a(2,3,1)=D_ANORDF((a(2,1)-rho_af*a(2,3)+dum)/sigma_af)
    trans_a(2,4,1)=D_ANORDF((a(2,1)-rho_af*a(2,4)+dum)/sigma_af)
    trans_a(2,5,1)=D_ANORDF((a(2,1)-rho_af*a(2,5)+dum)/sigma_af)
    trans_a(2,1,2)=D_ANORDF((a(2,2)-rho_af*a(2,1)+dum)/sigma_af)-D_ANORDF((a(2,2)-rho_af*a(2,1)-dum)/sigma_af)
    trans_a(2,2,2)=D_ANORDF((a(2,2)-rho_af*a(2,2)+dum)/sigma_af)-D_ANORDF((a(2,2)-rho_af*a(2,2)-dum)/sigma_af)
    trans_a(2,3,2)=D_ANORDF((a(2,2)-rho_af*a(2,3)+dum)/sigma_af)-D_ANORDF((a(2,2)-rho_af*a(2,3)-dum)/sigma_af)
    trans_a(2,4,2)=D_ANORDF((a(2,2)-rho_af*a(2,4)+dum)/sigma_af)-D_ANORDF((a(2,2)-rho_af*a(2,4)-dum)/sigma_af)
    trans_a(2,5,2)=D_ANORDF((a(2,2)-rho_af*a(2,5)+dum)/sigma_af)-D_ANORDF((a(2,2)-rho_af*a(2,5)-dum)/sigma_af)
    trans_a(2,1,3)=D_ANORDF((a(2,3)-rho_af*a(2,1)+dum)/sigma_af)-D_ANORDF((a(2,3)-rho_af*a(2,1)-dum)/sigma_af)
    trans_a(2,2,3)=D_ANORDF((a(2,3)-rho_af*a(2,2)+dum)/sigma_af)-D_ANORDF((a(2,3)-rho_af*a(2,2)-dum)/sigma_af)
    trans_a(2,3,3)=D_ANORDF((a(2,3)-rho_af*a(2,3)+dum)/sigma_af)-D_ANORDF((a(2,3)-rho_af*a(2,3)-dum)/sigma_af)
    trans_a(2,4,3)=D_ANORDF((a(2,3)-rho_af*a(2,4)+dum)/sigma_af)-D_ANORDF((a(2,3)-rho_af*a(2,4)-dum)/sigma_af)
    trans_a(2,5,3)=D_ANORDF((a(2,3)-rho_af*a(2,5)+dum)/sigma_af)-D_ANORDF((a(2,3)-rho_af*a(2,5)-dum)/sigma_af)
    trans_a(2,1,4)=D_ANORDF((a(2,4)-rho_af*a(2,1)+dum)/sigma_af)-D_ANORDF((a(2,4)-rho_af*a(2,1)-dum)/sigma_af)
    trans_a(2,2,4)=D_ANORDF((a(2,4)-rho_af*a(2,2)+dum)/sigma_af)-D_ANORDF((a(2,4)-rho_af*a(2,2)-dum)/sigma_af)
    trans_a(2,3,4)=D_ANORDF((a(2,4)-rho_af*a(2,3)+dum)/sigma_af)-D_ANORDF((a(2,4)-rho_af*a(2,3)-dum)/sigma_af)
    trans_a(2,4,4)=D_ANORDF((a(2,4)-rho_af*a(2,4)+dum)/sigma_af)-D_ANORDF((a(2,4)-rho_af*a(2,4)-dum)/sigma_af)
    trans_a(2,5,4)=D_ANORDF((a(2,4)-rho_af*a(2,5)+dum)/sigma_af)-D_ANORDF((a(2,4)-rho_af*a(2,5)-dum)/sigma_af)
    trans_a(2,1,5)=1D0-D_ANORDF((a(2,5)-rho_af*a(2,1)-dum)/sigma_af)
    trans_a(2,2,5)=1D0-D_ANORDF((a(2,5)-rho_af*a(2,2)-dum)/sigma_af)
    trans_a(2,3,5)=1D0-D_ANORDF((a(2,5)-rho_af*a(2,3)-dum)/sigma_af)
    trans_a(2,4,5)=1D0-D_ANORDF((a(2,5)-rho_af*a(2,4)-dum)/sigma_af)
    trans_a(2,5,5)=1D0-D_ANORDF((a(2,5)-rho_af*a(2,5)-dum)/sigma_af)
    
    prob_a=1D0/5D0
    
    do it=1,100
        prob_a(1,:)=MATMUL(prob_a(1,:),trans_a(1,:,:))
    end do
    
    do it=1,100
        prob_a(2,:)=MATMUL(prob_a(2,:),trans_a(2,:,:))
    end do
    
    !call tauchen_hans(sigma_am,rho_am,na,a(1,:),trans_a(1,:,:),prob_a(1,:))
    !call tauchen_hans(sigma_um,rho_um,nu,u(1,:),trans_u(1,:,:),prob_u(1,:))
    !
    !call tauchen_hans(sigma_af,rho_af,na,a(2,:),trans_a(2,:,:),prob_a(2,:))
    !call tauchen_hans(sigma_uf,rho_uf,nu,u(2,:),trans_u(2,:,:),prob_u(2,:))
    
    call tauchen_hans(sigma_fcm,rho_fcm,nfc,fc(1,:),trans_fc(1,:,:),prob_fc(1,:))
    call tauchen_hans(sigma_fcs,rho_fcs,nfc,fc(2,:),trans_fc(2,:,:),prob_fc(2,:))
    
    !fc(1,1)=-(14D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,2)=-(12D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,3)=-(10D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,4)=-(8D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,5)=-(6D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,6)=-(4D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,7)=-(2D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,8)=0D0
    !fc(1,9)=(2D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,10)=(4D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,11)=(6D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,12)=(8D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,13)=(10D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,14)=(12D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !fc(1,15)=(14D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !dum=(5D-1)*(2D0/7D0)*SQRT((sigma_fcm**2)/(1.0-(rho_fcm**2)))
    !
    !trans_fc(1,1,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,1)+dum)/sigma_fcm)
    !trans_fc(1,2,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,2)+dum)/sigma_fcm)
    !trans_fc(1,3,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,3)+dum)/sigma_fcm)
    !trans_fc(1,4,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,4)+dum)/sigma_fcm)
    !trans_fc(1,5,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,5)+dum)/sigma_fcm)
    !trans_fc(1,6,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,6)+dum)/sigma_fcm)
    !trans_fc(1,7,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,7)+dum)/sigma_fcm)
    !trans_fc(1,8,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,8)+dum)/sigma_fcm)
    !trans_fc(1,9,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,9)+dum)/sigma_fcm)
    !trans_fc(1,10,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,10)+dum)/sigma_fcm)
    !trans_fc(1,11,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,11)+dum)/sigma_fcm)
    !trans_fc(1,12,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,12)+dum)/sigma_fcm)
    !trans_fc(1,13,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,13)+dum)/sigma_fcm)
    !trans_fc(1,14,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,14)+dum)/sigma_fcm)
    !trans_fc(1,15,1)=D_ANORDF((fc(1,1)-rho_fcm*fc(1,15)+dum)/sigma_fcm)
    !
    !trans_fc(1,1,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,2)=D_ANORDF((fc(1,2)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,2)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,3)=D_ANORDF((fc(1,3)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,3)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,4)=D_ANORDF((fc(1,4)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,4)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,5)=D_ANORDF((fc(1,5)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,5)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,6)=D_ANORDF((fc(1,6)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,6)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,7)=D_ANORDF((fc(1,7)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,7)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,8)=D_ANORDF((fc(1,8)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,8)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,9)=D_ANORDF((fc(1,9)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,9)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,10)=D_ANORDF((fc(1,10)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,10)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,11)=D_ANORDF((fc(1,11)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,11)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,12)=D_ANORDF((fc(1,12)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,12)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,13)=D_ANORDF((fc(1,13)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,13)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,1)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,2)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,3)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,4)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,5)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,6)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,7)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,8)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,9)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,10)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,11)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,12)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,13)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,14)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,14)=D_ANORDF((fc(1,14)-rho_fcm*fc(1,15)+dum)/sigma_fcm)-D_ANORDF((fc(1,14)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !trans_fc(1,1,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,1)-dum)/sigma_fcm)
    !trans_fc(1,2,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,2)-dum)/sigma_fcm)
    !trans_fc(1,3,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,3)-dum)/sigma_fcm)
    !trans_fc(1,4,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,4)-dum)/sigma_fcm)
    !trans_fc(1,5,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,5)-dum)/sigma_fcm)
    !trans_fc(1,6,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,6)-dum)/sigma_fcm)
    !trans_fc(1,7,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,7)-dum)/sigma_fcm)
    !trans_fc(1,8,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,8)-dum)/sigma_fcm)
    !trans_fc(1,9,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,9)-dum)/sigma_fcm)
    !trans_fc(1,10,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,10)-dum)/sigma_fcm)
    !trans_fc(1,11,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,11)-dum)/sigma_fcm)
    !trans_fc(1,12,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,12)-dum)/sigma_fcm)
    !trans_fc(1,13,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,13)-dum)/sigma_fcm)
    !trans_fc(1,14,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,14)-dum)/sigma_fcm)
    !trans_fc(1,15,15)=1D0-D_ANORDF((fc(1,15)-rho_fcm*fc(1,15)-dum)/sigma_fcm)
    !
    !
    !fc(2,1)=-(14D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,2)=-(12D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,3)=-(10D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,4)=-(8D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,5)=-(6D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,6)=-(4D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,7)=-(2D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,8)=0D0
    !fc(2,9)=(2D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,10)=(4D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,11)=(6D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,12)=(8D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,13)=(10D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,14)=(12D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !fc(2,15)=(14D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !dum=(5D-1)*(2D0/7D0)*SQRT((sigma_fcs**2)/(1.0-(rho_fcs**2)))
    !
    !trans_fc(2,1,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,1)+dum)/sigma_fcs)
    !trans_fc(2,2,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,2)+dum)/sigma_fcs)
    !trans_fc(2,3,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,3)+dum)/sigma_fcs)
    !trans_fc(2,4,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,4)+dum)/sigma_fcs)
    !trans_fc(2,5,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,5)+dum)/sigma_fcs)
    !trans_fc(2,6,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,6)+dum)/sigma_fcs)
    !trans_fc(2,7,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,7)+dum)/sigma_fcs)
    !trans_fc(2,8,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,8)+dum)/sigma_fcs)
    !trans_fc(2,9,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,9)+dum)/sigma_fcs)
    !trans_fc(2,10,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,10)+dum)/sigma_fcs)
    !trans_fc(2,11,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,11)+dum)/sigma_fcs)
    !trans_fc(2,12,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,12)+dum)/sigma_fcs)
    !trans_fc(2,13,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,13)+dum)/sigma_fcs)
    !trans_fc(2,14,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,14)+dum)/sigma_fcs)
    !trans_fc(2,15,1)=D_ANORDF((fc(2,1)-rho_fcs*fc(2,15)+dum)/sigma_fcs)
    !trans_fc(2,1,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,2)=D_ANORDF((fc(2,2)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,2)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,3)=D_ANORDF((fc(2,3)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,3)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,4)=D_ANORDF((fc(2,4)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,4)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,5)=D_ANORDF((fc(2,5)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,5)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,6)=D_ANORDF((fc(2,6)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,6)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,7)=D_ANORDF((fc(2,7)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,7)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,8)=D_ANORDF((fc(2,8)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,8)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,9)=D_ANORDF((fc(2,9)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,9)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,10)=D_ANORDF((fc(2,10)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,10)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,11)=D_ANORDF((fc(2,11)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,11)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,12)=D_ANORDF((fc(2,12)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,12)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,13)=D_ANORDF((fc(2,13)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,13)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,1)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,2)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,3)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,4)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,5)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,6)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,7)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,8)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,9)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,10)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,11)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,12)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,13)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,14)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,14)=D_ANORDF((fc(2,14)-rho_fcs*fc(2,15)+dum)/sigma_fcs)-D_ANORDF((fc(2,14)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !trans_fc(2,1,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,1)-dum)/sigma_fcs)
    !trans_fc(2,2,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,2)-dum)/sigma_fcs)
    !trans_fc(2,3,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,3)-dum)/sigma_fcs)
    !trans_fc(2,4,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,4)-dum)/sigma_fcs)
    !trans_fc(2,5,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,5)-dum)/sigma_fcs)
    !trans_fc(2,6,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,6)-dum)/sigma_fcs)
    !trans_fc(2,7,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,7)-dum)/sigma_fcs)
    !trans_fc(2,8,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,8)-dum)/sigma_fcs)
    !trans_fc(2,9,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,9)-dum)/sigma_fcs)
    !trans_fc(2,10,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,10)-dum)/sigma_fcs)
    !trans_fc(2,11,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,11)-dum)/sigma_fcs)
    !trans_fc(2,12,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,12)-dum)/sigma_fcs)
    !trans_fc(2,13,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,13)-dum)/sigma_fcs)
    !trans_fc(2,14,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,14)-dum)/sigma_fcs)
    !trans_fc(2,15,15)=1D0-D_ANORDF((fc(2,15)-rho_fcs*fc(2,15)-dum)/sigma_fcs)
    !
    !prob_fc=1D0/15D0
    !
    !do it=1,100
    !    prob_fc(1,:)=MATMUL(prob_fc(1,:),trans_fc(1,:,:))
    !end do
    !
    !do it=1,100
    !    prob_fc(2,:)=MATMUL(prob_fc(2,:),trans_fc(2,:,:))
    !end do
    
    fc(1,:)=fc(1,:)+mu_fcm
    fc(2,:)=fc(2,:)+mu_fcs
    
    !open(1, file='transprob.txt')
    !    
    !    write (1, *) trans_u, Prob_u, trans_a, prob_a, trans_fc, prob_fc
    !    
    !close(1)
    !STOP
    !Print *,prob_u(1,:)
    !
    !Print *,trans_fc(1,1,1)
    !Print *,trans_fc(1,1,2)
    !Print *,trans_fc(1,1,3)
    !Print *,trans_fc(1,1,4)
    !Print *,trans_fc(1,1,5)
    !Print *,trans_fc(1,1,6)
    !Print *,trans_fc(1,1,7)
    !Print *,trans_fc(1,1,8)
    !Print *,trans_fc(1,1,9)
    !Print *,trans_fc(1,1,10)
    !Print *,trans_fc(1,1,11)
    !Print *,trans_fc(1,1,12)
    !Print *,trans_fc(1,1,13)
    !Print *,trans_fc(1,1,14)
    !Print *,trans_fc(1,1,15)
    
    !Print *,trans_u(1,4,:)
    !Print *,trans_u(1,5,:)
    !
    !dum=0D0
    !do it=1,15
    !    dum=dum+trans_fc(1,it,2)
    !end do
    
    !dum=0D0
    !do it=1,15
    !    dum=dum+prob_fc(2,it)
    !end do
    
    !dum=probfam(1,15)+probfam(2,15)+probfam(3,15)
    !
    !Print *,dum
    
    !Print *,fc(1,:)
    
    !Print *,fc(2,:)
    
    !Print *,prob_fc(1,:)
    
    !Print *,prob_fc(2,:)
    
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
    
    
!    open(1, file='random1.txt')
!    
!        read (1, *) random1
!        
!    close(1)
!    
!    open(2, file='random2.txt')
!    
!        read (2, *) random2
!        
!    close(2)
!    
!    open(3, file='random3.txt')
!    
!        read (3, *) random3
!        
!    close(3)
    
            
        end subroutine Initialize  
             
    
    end program Laffer
