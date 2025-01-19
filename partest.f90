subroutine partest(counter)
    !This subroutine computes the expected value function of getting married for a single female
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8) :: dum3 
    real(8) :: vnext
    real(8) :: vnext_test
    integer :: iflag
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))  
    real(8) :: exp_m_prime, exp_f_prime
    
    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1      

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ix=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ix-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ix-1)*na*nu-(ium-1)*na

    j=2
    exp_f_prime = exp_grid(ix,T-it)
    evm(j,ik,ix,iam,ium,T-it,:)=0d0
    do ifc=1,nfc
        do iu2=1,nu
            do ik2=1,nk
                do ixm = 1, nexp
                    exp_m_prime = exp_grid(ixm,T-it)
                    do iaf=1,na
                        do iuf=1,nu
                            do ifcm=1,nfcm
                                dum=k_grid(ik)+k_grid(ik2)
                                vnext = pol_v_aux(iaf,iuf,iam,iu2,ifc,ifcm,WOMEN)%eval([dum,exp_f_prime,exp_m_prime])
                                evm(j,ik,ix,iam,ium,T-it,ifc) = evm(j,ik,ix,iam,ium,T-it,ifc) + &
                                                                trans_u(WOMEN,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*vnext                                
                                
                                
                                !if(dum<k_grid(nk)-0.001d0) then
                                !    call db3val(dum,exp_grid(ix,T-it),exp_grid(ixm,T-it),idx,idy,idz,&
                                !            tx,ty(:,T-it),tz(:,T-it),&
                                !            nk,nexp,nexp,kx,ky,kz,&
                                !            v_bspl(iaf,iuf,iam,iu2,ifc,ifcm)%coefs,vnext,iflag,&
                                !            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)                                    
                                !    evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*vnext
                                !else
                                !    evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*LinInterp(dum,k_grid,v(:,ix,ixm,iaf,iuf,iam,iu2,T-it,ifc,ifcm),nk)
                                !end if
                                
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine partest
    
subroutine partest2(counter)
    !This subroutine computes the expected value function of getting married for a single male
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext
    real(8) :: vnext_test
    integer :: iflag
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))  
    real(8) :: exp_m_prime, exp_f_prime

    
    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1       

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ixm=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ixm-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ixm-1)*na*nu-(ium-1)*na

    
    j=1
    exp_m_prime = exp_grid(ixm,T-it)
    evm(j,ik,ixm,iam,ium,T-it,:)=0d0
    do ifcm=1,nfcm
        do iu2=1,nu
            do ik2=1,nk
                do ix=1,nexp
                    exp_f_prime = exp_grid(ix,T-it)
                    do iaf=1,na
                        do iuf=1,nu
                            do ifc=1,nfc
                                dum=k_grid(ik)+k_grid(ik2)
                                vnext = pol_v_aux(iam,iu2,iaf,iuf,ifc,ifcm,MEN)%eval([dum,exp_f_prime,exp_m_prime])
                                evm(j,ik,ixm,iam,ium,T-it,ifcm) = evm(j,ik,ixm,iam,ium,T-it,ifcm)&
                                         + trans_u(MEN,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext   
                                
                                
                                !if(dum<k_grid(nk)-0.001d0) then
                                !    call db3val(dum,exp_grid(ix,T-it),exp_grid(ixm,T-it),idx,idy,idz,&
                                !    tx,ty(:,T-it),tz(:,T-it),&
                                !    nk,nexp,nexp,kx,ky,kz,&
                                !    v_bspl(iam,iu2,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                !    inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)   
                                !    
                                !    evm(j,ik,ixm,iam,ium,T-it,ifcm)=evm(j,ik,ixm,iam,ium,T-it,ifcm)+trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext
                                !else
                                !    evm(j,ik,ixm,iam,ium,T-it,ifcm)=evm(j,ik,ixm,iam,ium,T-it,ifcm)+trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*LinInterp(dum,k_grid,v(:,ix,ixm,iam,iu2,iaf,iuf,T-it,ifc,ifcm),nk)
                                !end if
                            end do
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
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ium,ifcm
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    integer :: iflag
    integer :: iknot
    real(8), pointer :: exp_grid_ptr(:)    
    real(8), pointer :: v_mar_tmp(:,:,:)
    
    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na
    exp_grid_ptr => exp_grid(:,T-it)
    iknot = 0    

    do ifcm=1,nfcm    
        do iaf=1,na
            do iuf=1,nu
                call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                    v(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm),&
                    kx,ky,kz,iknot,tx,ty(:,T-it),tz(:,T-it),&
                    v_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag) 
                call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                    ev(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm),&
                    kx,ky,kz,iknot,tx,ty(:,T-it),tz(:,T-it),&
                    ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag)  
                
                do j = 1, 2
                    v_mar_tmp => v_aux(:,:,:,iam,ium,iaf,iuf,ifc,ifcm,j)
                    call pol_v_aux(iam, ium, iaf, iuf, ifc, ifcm, j)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)   
                end do                 
                
            end do
        end do
    end do

end subroutine partest3
   
subroutine partest5(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    !USE BS2IN_INT
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ium
    real(8) :: dum3 !c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    !real(8) :: P1,P2,P3,P4,V2,V3,dum2
    integer :: iknot
    integer :: iflag
    real(8), pointer:: exp_grid_ptr(:)    
    
    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na
    
    exp_grid_ptr => exp_grid(:,T-it)
    iknot = 0    
    
    do j=1,2
        !CALL D_BS2IN(k_grid, exp_grid(:,T-it), evs(j,:,:,iam,ium,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), evs_spln_coefs(j,:,:,iam,ium,T-it,ifc), nk, nexp)
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, evs(j,:,:,iam,ium,T-it,ifc), &
            kx, ky, iknot, tx, ty(:,T-it), evs_bspl(j,iam,ium,ifc)%coefs, iflag)          
    end do
            
    do j=1,2
        !CALL D_BS2IN(k_grid, exp_grid(:,T-it), vs(j,:,:,iam,ium,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), vs_spln_coefs(j,:,:,iam,ium,T-it,ifc), nk, nexp)
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, vs(j,:,:,iam,ium,T-it,ifc), &
            kx, ky, iknot, tx, ty(:,T-it), vs_bspl(j,iam,ium,ifc)%coefs, iflag)        
    end do

end subroutine partest5   
    
subroutine partest6(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    !USE BS2IN_INT
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ium
    real(8) :: dum3 
    integer :: iknot
    integer :: iflag
    real(8), pointer:: exp_grid_ptr(:)    
    

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na
    
    iknot = 0
    exp_grid_ptr => exp_grid(:,T-it)    

    do j=1,2
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, evm(j,:,:,iam,ium,T-it,ifc), &
            kx, ky, iknot, tx, ty(:,T-it), evm_bspl(j,iam,ium,ifc)%coefs, iflag)           
    end do

end subroutine partest6
    
subroutine partest7(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8) :: dum3 !,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext, vnext_test
    integer :: iflag
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz

    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1    

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na)*1d0)/(nu*na*1d0))-0.00001d0
    ix=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*nu*na-(ix-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*nexp*na*nu-(ix-1)*na*nu-(ium-1)*na

    if(T-it>1) then

        ev(ik,ix,:,iam,ium,:,:,T-it,:,:)=0d0
        ev_aux(ik,ix,:,iam,ium,:,:,:,:,:) = 0d0
        dum3=K_grid(ik)/2d0

        do ifc=1,nfc
            do ifcm=1,nfcm
                do ixm = 1, nexp
                    do iuf = 1, nu
                        do iaf = 1, na

                            do iu2=1,nu
                                do iu3=1,nu
                                    ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V(ik,ix,ixm,iam,iu2,iaf,iu3,T-it,ifc,ifcm)
                                    ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN) = ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,MEN)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V_aux(ik,ix,ixm,iam,iu2,iaf,iu3,ifc,ifcm,MEN)
                                    ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN) = ev_aux(ik,ix,ixm,iam,ium,iaf,iuf,ifc,ifcm,WOMEN)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V_aux(ik,ix,ixm,iam,iu2,iaf,iu3,ifc,ifcm,WOMEN)
                                end do
                            end do

                            do iu3=1,nu
                                call db2val(dum3,exp_grid(ix,T-it),idx,idy,&
                                    tx,ty(:,T-it),nk,nexp,kx,ky,&
                                    vs_bspl(2,iaf,iu3,ifc)%coefs,vnext,iflag,&
                                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)  
                                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+Probd(T-it-1)*trans_u(2,iuf,iu3)*0.5d0*vnext
                            end do

                            do iu2=1,nu
                                call db2val(dum3,exp_grid(ixm,T-it),idx,idy,&
                                    tx,ty(:,T-it),nk,nexp,kx,ky,&
                                    vs_bspl(1,iam,iu2,ifcm)%coefs,vnext,iflag,&
                                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)                                          
                                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+Probd(T-it-1)*trans_u(1,ium,iu2)*0.5d0*vnext
                            end do

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

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ik
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*na
    
    evs(:,ik,:,iam,ium,T-it,:)=0d0
            
    
    do j=1,2
        do ifc=1,nfc
            do ix = 1, nexp
                    do iu2=1,nu
                        evs(j,ik,ix,iam,ium,T-it,ifc)=evs(j,ik,ix,iam,ium,T-it,ifc)+trans_u(j,ium,iu2)*Vs(j,ik,ix,iam,iu2,T-it,ifc)
                    end do
            end do
        end do
    end do        
            
            
end subroutine partest8
    
        
subroutine partest10(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use bspline_sub_module
    
    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,ifc,ifcm,ium
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8), pointer:: exp_grid_ptr(:)
    integer :: iknot
    integer :: iflag     
    
    !Assigning the grid points
    iam=counter
    exp_grid_ptr => exp_grid(:,T+Tret+1-it)
    iknot = 0
    
    do j=1,2
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, evs_ret(j,:,:,iam,Tret-it+1), &
                    kx, ky, iknot, tx, ty(:,T+Tret+1-it), evs_ret_bspl(j,iam)%coefs, iflag)
        call db2ink(k_grid, nk, exp_grid_ptr, nexp, edcs_ret(j,:,:,iam,Tret-it+1), &
                    kx, ky, iknot, tx, ty(:,T+Tret+1-it), edcs_ret_bspl(j,iam)%coefs, iflag)        
    end do
    
    do iaf=1,na
        call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                    ev_ret(:,:,:,iam,iaf,Tret-it+1),&
                    kx,ky,kz,iknot,tx,ty(:,T+Tret+1-it),tz(:,T+Tret+1-it),&
                    ev_ret_bspl(iam,iaf)%coefs,iflag) 
        call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                    edc_ret(:,:,:,iam,iaf,Tret-it+1),&
                    kx,ky,kz,iknot,tx,ty(:,T+Tret+1-it),tz(:,T+Tret+1-it),&
                    edc_ret_bspl(iam,iaf)%coefs,iflag)     
    end do
            
end subroutine partest10

subroutine update_ev_aux(i_age)
    use PolicyFunctions, only: exp_grid, ev_aux, v_aux, k_grid, pol_ev_aux
    use Model_Parameters, only: MEN, WOMEN, na, nexp, nfc, nk, nu
    implicit none
    integer :: i_age
    integer :: i, ia, iu, ifc
    real(8), pointer :: exp_grid_ptr(:) 
    integer :: iam, iaf, ium, iuf, ifcm, ifcf
    integer :: iexp1, iexp2
    real(8) :: ev_sing_tmp(nk,nexp), ev_mar_tmp(nk,nexp,nexp)
    real(8), pointer :: v_mar_tmp(:,:,:)
    integer :: i_m, i_f, jj
    
    exp_grid_ptr => exp_grid(:,i_age)     
    
    do iam = 1, na
        do iaf = 1, na
            do ium = 1, nu
                do iuf = 1, nu
                    do ifcm = 1, nfc
                        do ifcf = 1, nfc                                
                            do jj = MEN, WOMEN
                                v_mar_tmp => ev_aux(:,:,:,iam,ium,iaf,iuf,ifcm,ifcf,jj)
                                call pol_ev_aux(iam, ium, iaf, iuf, ifcm, ifcf, jj)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                !v_mar_tmp => v_aux(:,:,:,iam,ium,iaf,iuf,ifcm,ifcf,jj)
                                !call pol_v_aux(iam, ium, iaf, iuf, ifcm, ifcf, jj)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                            end do             
                        end do
                    end do
                end do
            end do
        end do
    end do
            
end subroutine update_ev_aux   

subroutine update_lfp_policies(i_age)
    use PolicyFunctions, only: exp_grid, evs, evm, ev, k_grid, v, vs, v_lfp, vs_lfp, pol_v_mar_lfp, pol_v_sing_lfp
    !use GlobParams, only: 
    use Model_Parameters, only: nexp, nk, LFP_M1, LFP_M0, LFP_F1, LFP_F0, LFP_0, LFP_1, MEN, WOMEN
    !use pyplot_module, only : pyplot, wp => pyplot_wp

    integer :: i_age
    integer :: i, ia, iu, ifc
    real(8), pointer :: exp_grid_ptr(:) 
    integer :: iam, iaf, ium, iuf, ifcm, ifcf
    integer :: iexp1, iexp2
    real(8) :: ev_sing_tmp(nk,nexp), ev_mar_tmp(nk,nexp,nexp)
    real(8) :: v_mar_tmp(nk,nexp,nexp), v_sing_tmp(nk,nexp)
    !type(pyplot) :: plt   !! pytplot handler
    character(len=*), parameter :: testdir = "Plots/"
    integer :: istat !! status code
    character(len=2) :: age_ch
    integer :: i_m, i_f, ii

    write(age_ch, '(i2)') i_age

    exp_grid_ptr => exp_grid(:,i_age)
    do i = 1, 2
        do ia = 1, na
            do iu = 1, nu
                do ifc = 1, nfc
                    !ev_sing_tmp = beta*OmegaActive(i_age-1)*( (1d0-Probm(i_age-1))*evs(i,:,:,ia,iu,i_age,ifc) + Probm(i_age-1)*evm(i,:,:,ia,iu,i_age,ifc) )
                    !call pol_ev_single(i, ia, iu, ifc)%set(ev_sing_tmp, k_grid, exp_grid_ptr)
                    !v_sing_tmp = Vs(i,:,:,ia,iu,i_age,ifc)
                    !call pol_v_sing(i, ia, iu, i_age, ifc)%set(v_sing_tmp, k_grid, exp_grid_ptr)
                    do ii = LFP_1, LFP_0
                        v_sing_tmp = Vs_lfp(i,:,:,ia,iu,i_age,ifc,ii)
                        call pol_v_sing_lfp(i, ia, iu, i_age, ifc, ii)%set(v_sing_tmp, k_grid, exp_grid_ptr)                            
                    end do
                end do
            end do
        end do
    end do          


    do iam = 1, na
        do iaf = 1, na
            do ium = 1, nu
                do iuf = 1, nu
                    do ifcm = 1, nfc
                        do ifcf = 1, nfc
                            !ev_mar_tmp = beta*OmegaActive(i_age-1)*ev(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf)
                            !ev_mar_tmp = ev(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf)
                            !call pol_ev(iam, ium, iaf, iuf, ifcm, ifcf)%set(ev_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                            !v_mar_tmp = v(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf)
                            !call pol_v_mar(iam, ium, iaf, iuf, i_age, ifcm, ifcf)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)

                            do i_m = LFP_M1, LFP_M0
                                do i_f = LFP_F1, LFP_F0
                                    do ii = MEN, WOMEN
                                        v_mar_tmp = v_lfp(:,:,:,iam,ium,iaf,iuf,i_age,ifcm,ifcf,i_m,i_f,ii)
                                        call pol_v_mar_lfp(iam, ium, iaf, iuf, i_age, ifcm, ifcf, i_m, i_f,ii)%set(v_mar_tmp, k_grid, exp_grid_ptr, exp_grid_ptr)
                                    end do
                                end do
                            end do

                        end do
                    end do
                end do
            end do
        end do
    end do

    !call plt%initialize(grid=.true.,xlabel='Savings',figsize=[20,10],&
    !                    title='Married',legend=.true.,axis_equal=.true.,&
    !                    tight_layout=.true.)  
    !iam = 2
    !ium = 2
    !iaf = 2
    !iuf = 2
    !ifcm = 1
    !ifcf = 1
    !call plt%add_plot(k_grid,v(:,iexp1,iexp2,iam,ium,iaf,iuf,i_age,ifcm,ifcf),label='iexp=1',linestyle='b-o',markersize=5,linewidth=2,istat=istat)
    !call plt%savefig(testdir//'Vmarried'//age_ch//'.png', pyfile=testdir//'plottest.py',istat=istat)   

end subroutine update_lfp_policies
