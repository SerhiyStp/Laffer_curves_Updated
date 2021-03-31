subroutine partest(counter)
    !This subroutine computes optimal policies at age 64
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
    evm(j,ik,ix,iam,ium,T-it,:)=0d0
    do ifc=1,nfc
        do iu2=1,nu
            do ik2=1,nk
                do ixm = 1, nexp
                    do iaf=1,na
                        do iuf=1,nu
                            do ifcm=1,nfcm
                                dum=k_grid(ik)+k_grid(ik2)
                                if(dum<k_grid(nk)-0.001d0) then
                                    call db3val(dum,exp_grid(ix,T-it),exp_grid(ixm,T-it),idx,idy,idz,&
                                            tx,ty(:,T-it),tz(:,T-it),&
                                            nk,nexp,nexp,kx,ky,kz,&
                                            v_bspl(iaf,iuf,iam,iu2,ifc,ifcm)%coefs,vnext,iflag,&
                                            inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)                                    
                                    evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*vnext
                                else
                                    evm(j,ik,ix,iam,ium,T-it,ifc)=evm(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*LinInterp(dum,k_grid,v(:,ix,ixm,iaf,iuf,iam,iu2,T-it,ifc,ifcm),nk)
                                end if
                            end do
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
    evm(j,ik,ixm,iam,ium,T-it,:)=0d0
        do ifcm=1,nfcm
            do iu2=1,nu
                do ik2=1,nk
                    do ix=1,nexp
                        do iaf=1,na
                            do iuf=1,nu
                                    do ifc=1,nfc
                                        dum=k_grid(ik)+k_grid(ik2)
                                        if(dum<k_grid(nk)-0.001d0) then
                                            call db3val(dum,exp_grid(ix,T-it),exp_grid(ixm,T-it),idx,idy,idz,&
                                                    tx,ty(:,T-it),tz(:,T-it),&
                                                    nk,nexp,nexp,kx,ky,kz,&
                                                    v_bspl(iam,iu2,iaf,iuf,ifc,ifcm)%coefs,vnext,iflag,&
                                                    inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.true.)   
                                            
                                            evm(j,ik,ixm,iam,ium,T-it,:)=evm(j,ik,ixm,iam,ium,T-it,:)+trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext
                                        else
                                            evm(j,ik,ixm,iam,ium,T-it,:)=evm(j,ik,ixm,iam,ium,T-it,:)+trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*LinInterp(dum,k_grid,v(:,ix,ixm,iam,iu2,iaf,iuf,T-it,ifc,ifcm),nk)
                                        end if
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
    !USE BS3IN_INT
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ium,ifcm
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    integer :: iflag
    integer :: iknot
    real(8), pointer :: exp_grid_ptr(:)    
    
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
                !CALL D_BS3IN(k_grid, exp_grid(:,T-it), exp_grid(:,T-it), ev(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm), KORDER, EXPORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), EXP_KNOT(:,T-it), ev_spln_coefs(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm), nk, nexp, nexp)
                !CALL D_BS3IN(k_grid, exp_grid(:,T-it), exp_grid(:,T-it), v(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm), KORDER, EXPORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), EXP_KNOT(:,T-it), v_spln_coefs(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm), nk, nexp, nexp)
                call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                    v(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm),&
                    kx,ky,kz,iknot,tx,ty(:,T-it),tz(:,T-it),&
                    v_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag) 
                call db3ink(k_grid,nk,exp_grid_ptr,nexp,exp_grid_ptr,nexp,&
                    ev(:,:,:,iam,ium,iaf,iuf,T-it,ifc,ifcm),&
                    kx,ky,kz,iknot,tx,ty(:,T-it),tz(:,T-it),&
                    ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs,iflag)                 
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
        !CALL D_BS2IN(k_grid, exp_grid(:,T-it), evm(j,:,:,iam,ium,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), evm_spln_coefs(j,:,:,iam,ium,T-it,ifc), nk, nexp)
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
        dum3=K_grid(ik)/2d0

        do ifc=1,nfc
            do ifcm=1,nfcm
                do ixm = 1, nexp
                    do iuf = 1, nu
                        do iaf = 1, na

                            do iu2=1,nu
                                do iu3=1,nu
                                    ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+(1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V(ik,ix,ixm,iam,iu2,iaf,iu3,T-it,ifc,ifcm)
                                end do
                            end do

                            do iu3=1,nu
                                !vnext = D_BS2VL(dum3,exp_grid(ix,T-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it), nk, nexp, vs_spln_coefs(2,:,:,iaf,iu3,T-it,ifc))
                                call db2val(dum3,exp_grid(ix,T-it),idx,idy,&
                                    tx,ty(:,T-it),nk,nexp,kx,ky,&
                                    vs_bspl(2,iaf,iu3,ifc)%coefs,vnext,iflag,&
                                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)  
                                !if (abs(vnext-vnext_test)>1d-10) then
                                !    print *, 'WARNING'
                                !end if
                                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)=ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm)+Probd(T-it-1)*trans_u(2,iuf,iu3)*0.5d0*vnext
                            end do

                            do iu2=1,nu
                                !vnext = D_BS2VL(dum3,exp_grid(ixm,T-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it), nk, nexp, vs_spln_coefs(1,:,:,iam,iu2,T-it,ifcm))
                                call db2val(dum3,exp_grid(ixm,T-it),idx,idy,&
                                    tx,ty(:,T-it),nk,nexp,kx,ky,&
                                    vs_bspl(1,iam,iu2,ifcm)%coefs,vnext,iflag,&
                                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.true.)       
                                !if (abs(vnext-vnext_test)>1d-10) then
                                !    print *, 'WARNING'
                                !end if                                    
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
