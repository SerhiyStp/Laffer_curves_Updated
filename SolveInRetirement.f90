subroutine SolveInRetirement(counter)
    ! Both men and women are retired (absorbing state)
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    USE QDVAL_INT
    USE BS2VL_INT
    USE BS3VL_INT
    use bspline_sub_module

    implicit none

    integer, INTENT(IN) :: counter
    real(8) :: c2, v2,dum2,dum3, Psi_pension, Psi_pensionf
    real(8) :: d1, d2,P1,P2,P3,P4,v3
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ifcm,ik
    integer :: it2
    real(8) :: test1, test2,r_ret, r_ret_next
    real(8) :: vnext_p2, vnext_p3
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp),pnt1(2),INTERP3D(nk,nexp,nexp),pnt2(3)
    integer :: iflag 
    real(8) :: vnext_test
    integer :: idx, idy, idz, iloy, iloz
    integer :: inbvx, inbvy, inbvz
    real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))
    real(8) :: w1_d2(ky) 
    real(8) :: w0_d2(3*max(kx,ky)) 

    idx=0
    idy=0
    idz=0
    inbvx=1
    inbvy=1
    inbvz=1
    iloy=1
    iloz=1

    !Assigning the grid points
    dum3=((counter*1d0)/(nexp*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nexp*na)*1d0)/(na*1d0))-0.00001d0
    ix=int(dum3)+1
    iam=counter-(ik-1)*nexp*na-(ix-1)*na

    r_ret=((1d0+r)/OmegaRet2(Tret+1-it))-1d0
    if(it>1) then
        r_ret_next=((1d0+r)/OmegaRet2(Tret+2-it))-1d0
    end if

    if(it>1) then
        exp_grid_dum=exp_grid(:,T+Tret+2-it)
    end if

    !Married    

    do iaf = 1, na
        do ixm = 1, nexp

            !Pension depends on expected wage conditional on ability and experience
            !Psi_pension=psi0+psi1*av_earnings(1,1,iam)*min(1d0,exp_grid(ixm,T+Tret+1-it)/35d0)
            !Psi_pensionf=psi0+psi1**av_earnings(2,1,iaf)*min(1d0,exp_grid(ix,T+Tret+1-it)/35d0)
            Psi_pension=psi0
            Psi_pensionf=psi0

            if(it==1) then
                ! Solve the very last period problem:
                k_ret(ik,ix,ixm,iam,iaf,Tret) = 0d0
                c_ret(ik,ix,ixm,iam,iaf,Tret) = ((k_grid(ik) + Gamma_redistr)*(1d0 + r_ret*(1d0-tk)) + Psi_pension+Psi_pensionf+lumpsum)/(1d0+tc)
                Uprime_ret(ik,ix,ixm,iam,iaf,Tret) = dUc(((k_grid(ik) + Gamma_redistr)*(1d0 + r_ret*(1d0-tk)) + Psi_pension+Psi_pensionf+lumpsum)/(1d0+tc))
                v_ret(ik,ix,ixm,iam,iaf,Tret) = Uc(c_ret(ik,ix,ixm,iam,iaf,Tret-it+1))
            else
                ! Solve the Tret-1 to 1st period of retirement:
                !Finding optimal capital by golden search
                P1=0.001d0
                P4=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk)) + Psi_pension+Psi_pensionf+lumpsum)/(1d0+tc)
                do
                    P2 = (P4+P1)/2d0

                    dum2=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk)) + Psi_pension+Psi_pensionf+lumpsum-P2*(1d0+tc))/(1d0+mu)

                    if((dum2<k_grid(nk)).AND.(dum2>k_grid(1))) then
                        !vnext_p2 = D_BS3VL(dum2, exp_grid(ix,T+Tret+1-it), exp_grid(ixm,T+Tret+1-it), KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret+2-it), EXP_KNOT(:,T+Tret+2-it), nk, nexp, nexp, edc_ret_spln_coefs(:,:,:,iam,iaf,Tret+2-it))
                        call db3val(dum2,exp_grid(ix,T+Tret+1-it),exp_grid(ixm,T+Tret+1-it),idx,idy,idz,&
                                tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                nk,nexp,nexp,kx,ky,kz,&
                                edc_ret_bspl(iam,iaf)%coefs,vnext_p2,iflag,&
                                inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)
                        V2=dUc(P2)-beta*OmegaRet(Tret-it+1)*((1d0+r_ret_next*(1d0-tk))/(1d0+mu))*vnext_p2
                    else
                        pnt2=(/dum2, exp_grid(ix,T+Tret-it), exp_grid(ixm,T+Tret-it)/)
                        INTERP3D=edc_ret(:,:,:,iam,iaf,Tret+2-it)
                        vnext_p2 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                        V2=dUc(P2)-beta*OmegaRet(Tret-it+1)*((1d0+r_ret_next*(1d0-tk))/(1d0+mu))*vnext_p2
                    end if

                    if (V2 < 0d0) then
                        P4=P2
                    else
                        P1=P2
                    end if
                    if((P4-P1)<1d-8) exit
                end do

                if(dum2<1d-4) then
                    dum2=k_grid(1)
                    P2=((k_grid(ik)+Gamma_redistr)*(1d0+r_ret*(1d0-tk)) + Psi_pension+Psi_pensionf+lumpsum)/(1d0+tc)
                end if

                if((dum2<k_grid(nk)).AND.(dum2>k_grid(1))) then
                    !vnext_p2 = D_BS3VL(dum2, exp_grid(ix,T+Tret+1-it), exp_grid(ixm,T+Tret+1-it), KORDER, EXPORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret+2-it), EXP_KNOT(:,T+Tret+2-it), nk, nexp, nexp, ev_ret_spln_coefs(:,:,:,iam,iaf,Tret+2-it))
                    call db3val(dum2,exp_grid(ix,T+Tret+1-it),exp_grid(ixm,T+Tret+1-it),idx,idy,idz,&
                                tx,ty(:,T+Tret-it+2),tz(:,T+Tret-it+2),&
                                nk,nexp,nexp,kx,ky,kz,&
                                ev_ret_bspl(iam,iaf)%coefs,vnext_p2,iflag,&
                                inbvx,inbvy,inbvz,iloy,iloz,ww2,ww1,ww0,extrap=.false.)     
                    !if (abs(vnext_p2 - vnext_test)>1d-12) then
                    !    print *, 'WARNING'
                    !end if
                else
                    pnt2=(/dum2, exp_grid(ix,T+Tret-it), exp_grid(ixm,T+Tret-it)/)
                    INTERP3D=ev_ret(:,:,:,iam,iaf,Tret+2-it)
                    vnext_p2 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                end if
                V2=Uc(P2)+beta*OmegaRet(Tret-it+1)*vnext_p2

                v_ret(ik,ix,ixm,iam,iaf,Tret-it+1)=V2
                k_ret(ik,ix,ixm,iam,iaf,Tret-it+1)=dum2
                c_ret(ik,ix,ixm,iam,iaf,Tret-it+1)=P2
                Uprime_ret(ik,ix,ixm,iam,iaf,Tret-it+1)=dUc(P2)
            end if
        end do
    end do

    !Single
    do j=1,2

        !Pension depends on expected earnings conditional on ability, marital status and gender
        !Psi_pension=psi0+psi1*av_earnings(j,2,iam)*min(1d0,exp_grid(ix,T+Tret+1-it)/35d0)
        Psi_pension=psi0

        if(it==1) then
            !Solve the very last period problem:        
            ks_ret(j,ik,ix,iam,Tret) = 0d0
            cs_ret(j,ik,ix,iam,Tret) = ((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pension+lumpsum*0.5d0)/(1d0+tc)
            Uprimes_ret(j,ik,ix,iam,Tret) = dUc(((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pension+lumpsum*0.5d0)/(1d0+tc))
            vs_ret(j,ik,ix,iam,Tret) = Uc(cs_ret(j,ik,ix,iam,Tret)) 
        else
            ! Solve the Tret-1 to 1st period of retirement:
            !Finding optimal capital by golden search
            P1=0.001d0
            P4=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pension+lumpsum*0.5d0)/(1d0+tc)
            do
                P2 = (P4+P1)/2d0

                dum2=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pension+lumpsum*0.5d0-P2*(1d0+tc))/(1d0+mu)

                if((dum2<k_grid(nk)).AND.(dum2>k_grid(1))) then
                    !vnext_p2 = D_BS2VL(dum2, exp_grid(ix,T+Tret+1-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret+2-it), nk, nexp,edcs_ret_spln_coefs(j,:,:,iam,Tret+2-it))
                    call db2val(dum2,exp_grid(ix,T+Tret+1-it),idx,idy,&
                        tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                        edcs_ret_bspl(j,iam)%coefs,vnext_p2,iflag,&
                        inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
                    V2=dUc(P2)-beta*OmegaRet(Tret-it+1)*((1d0+r_ret_next*(1d0-tk))/(1d0+mu))*vnext_p2
                else
                    pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
                    INTERP2D=edcs_ret(j,:,:,iam,Tret+2-it)
                    vnext_p2 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=dUc(P2)-beta*OmegaRet(Tret-it+1)*((1d0+r_ret_next*(1d0-tk))/(1d0+mu))*vnext_p2
                end if

                if (V2 < 0d0) then
                    P4=P2
                else
                    P1=P2
                end if
                if((P4-P1)<1d-8) exit
            end do

            if(dum2<1d-4) then
                dum2=k_grid(1)
                P2=((k_grid(ik)+Gamma_redistr*0.5d0)*(1d0+r_ret*(1d0-tk)) + Psi_pension+lumpsum*0.5d0)/(1d0+tc)
            end if

            if((dum2<k_grid(nk)).AND.(dum2>k_grid(1))) then
                !vnext_p2 = D_BS2VL(dum2, exp_grid(ix,T+Tret+1-it), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+Tret+2-it), nk, nexp,evs_ret_spln_coefs(j,:,:,iam,Tret+2-it))
                call db2val(dum2,exp_grid(ix,T+Tret+1-it),idx,idy,&
                    tx,ty(:,T+Tret+2-it),nk,nexp,kx,ky,&
                    evs_ret_bspl(j,iam)%coefs,vnext_p2,iflag,&
                    inbvx,inbvy,iloy,w1_d2,w0_d2,extrap=.false.)
            else
                pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
                INTERP2D=vs_ret(j,:,:,iam,Tret+2-it)
                vnext_p2 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
            end if
            V2=Uc(P2)+beta*OmegaRet(Tret-it+1)*vnext_p2

            vs_ret(j,ik,ix,iam,Tret-it+1)=V2
            ks_ret(j,ik,ix,iam,Tret-it+1)=dum2
            cs_ret(j,ik,ix,iam,Tret-it+1)=P2
            Uprimes_ret(j,ik,ix,iam,Tret-it+1)=dUc(P2)

            pnt1=(/dum2, exp_grid(ix,T+Tret+1-it)/)
            INTERP2D=cs_ret(j,:,:,iam,Tret+2-it)
            vnext_p2 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
            V2=dUc(P2)-beta*OmegaRet(Tret-it+1)*((1d0+r_ret_next*(1d0-tk))/(1d0+mu))*dUc(vnext_p2)
            Eulers_ret(j,ik,:,iam,Tret-it+1)=V2

        end if

    end do

end subroutine SolveInRetirement
