subroutine Solvefirstactive(counter)
    !This subroutine computes optimal policies durin active worklife.
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    USE QDVAL_INT
    USE QD2VL_INT
    USE BS2VL_INT
    !USE LCONF_INT
    !USE CSVAL_INT
    !USE NEQNF_INT
    !USE ERSET_INT

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ik
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext

!Assigning the grid points
dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
ik=int(dum3)+1
dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
ium=int(dum3)+1
iam=counter-(ik-1)*na*nu-(ium-1)*nu

exp_grid_dum=exp_grid(:,T+1-it)

    
    ce=0d0
    cu=0d0
    ke=0d0
    ku=0d0
    nem=0d0
    nef=0d0
    num=0d0
    nuf=0d0
    ve=0d0
    vu=0d0


    ces=0d0
    cus=0d0
    kes=0d0
    kus=0d0
    nes=0d0
    nus=0d0
    ves=0d0
    vus=0d0

    !if both spouses employed

        do ifc=1,nfc
                do iaf = 1, na
                    do iuf = 1, nu
                        ix = 1
                            !Finding optimal capital by golden search
                            wagem = wage(1,a(1,iam),dble(T-it),u(1,ium))/(1d0+t_employer)
                            wagef = wage(2,a(2,iaf),exp_grid(ix,T-it),u(2,iuf))/(1d0+t_employer)
                            P1=0.01d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+(wagem+wagef)*(1d0-tax_labor(wagem+wagef)-tSS_employee(wagem+wagef)))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt2 = (/P2, wagem, wagef/)
                                dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                y=dum4*wagem+dum5*wagef
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)
                                    !pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                                    !INTERP2D=ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc)
                                    !vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                                    vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V2=V2+beta*OmegaActive(T-it)*vnext
                                else
                                    V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)
                                    !dum2=max(dum2,0d0)
                                    !V2=V2+beta*OmegaActive(T-it)*D_CSVAL(dum2,BREAK,ev_spln_coefs(:,:,ix+1,iam,ium,iaf,iuf,T+1-it,j,ifc))
                                    !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    vnext = D_BS2VL(dum2,exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1-it), nk, nexp, ev_spln_coefs(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V2=V2+beta*OmegaActive(T-it)*vnext

                                end if

                                pnt2 = (/P3, wagem, wagef/)
                                dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                y=dum4*wagem+dum5*wagef
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)
                                    !pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                                    !INTERP2D=ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc)
                                    !vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                                    vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V3=V3+beta*OmegaActive(T-it)*vnext
                                else
                                    V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)
                                    !dum2=max(dum2,0d0)
                                    !V3=V3+beta*OmegaActive(T-it)*D_CSVAL(dum2,BREAK,ev_spln_coefs(:,:,ix+1,iam,ium,iaf,iuf,T+1-it,j,ifc))

                                    !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    vnext = D_BS2VL(dum2,exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1-it), nk, nexp, ev_spln_coefs(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V3=V3+beta*OmegaActive(T-it)*vnext
                                end if

                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-6) exit
                            end do

                            Ve=V2
                            ke=dum2
                            nem=dum4
                            nef=dum5
                            ce=P2
                            if(dum4<0d0) then
                                pause
                            end if
                            if(dum5<0d0) then
                                pause
                            end if
                        

    !If female unemployed

        
       
                            !Finding optimal capital by golden search
                            P1=0.01d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+(wagem)*(1d0-tax_labor(wagem)-tSS_employee(wagem)))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt1 = (/P2, wagem/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagem
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V2=Uc(P2)+Ul(dum4,0d0)
                                    !pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                                    !INTERP2D=ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc)
                                    !vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                                    vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V2=V2+beta*OmegaActive(T-it)*vnext
                                else
                                    V2=Uc(P2)+Ul(dum4,0d0)
                                    !dum2=max(dum2,0d0)
                                    !V2=V2+beta*OmegaActive(T-it)*D_CSVAL(dum2,BREAK,ev_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T+1-it,j,ifc))
                                    !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    vnext = D_BS2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1-it), nk, nexp, ev_spln_coefs(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V2=V2+beta*OmegaActive(T-it)*vnext
                                end if

                                pnt1 = (/P3, wagem/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagem
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                elseif(dum2>k_grid(nk)-0.001d0) then
                                    V3=Uc(P3)+Ul(dum4,0d0)
                                    dum2=max(dum2,0d0)
                                    !pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                                    !INTERP2D=ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc)
                                    !vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                                    vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V3=V3+beta*OmegaActive(T-it)*vnext
                                else
                                    V3=Uc(P3)+Ul(dum4,0d0)
                                    !V3=V3+beta*OmegaActive(T-it)*D_CSVAL(dum2,BREAK,ev_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T+1-it,j,ifc))

                                    !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T+1-it),ev(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    vnext = D_BS2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T+1-it), nk, nexp, ev_spln_coefs(:,:,iam,ium,iaf,iuf,T+1-it,ifc))
                                    V3=V3+beta*OmegaActive(T-it)*vnext
                                end if

                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-6) exit
                            end do

                            Vu=V2
                            ku=dum2
                            num=dum4
                            nuf=0d0
                            cu=P2
                            if(dum4<0d0) then
                                pause
                            end if
                        

    
                            if (ve >= vu) then
                                v(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ve
                                c(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ce
                                k(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ke
                                nm(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=nem
                                nf(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=nef
                            else
                                v(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=vu
                                c(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=cu
                                k(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=ku
                                nm(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=num
                                nf(ik,ix,iam,ium,iaf,iuf,T-it,ifc)=nuf
                            end if
                    end do
            end do
        end do


    !Singles

    j=2

    do ifc=1,nfc
            ix = 1
                !Print *,'ix is',ix
                !Finding optimal capital by golden search
                wagef = wage(2,a(2,iam),exp_grid(ix,T-it),u(2,ium))/(1d0+t_employer)
                P1=0.01d0
                P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+wagef*(1d0-tax_labors(wagef)-tSS_employee(wagef)))/(1d0+tc))
                do
                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                    pnt1 = (/P2, wagef/)
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagef
                    dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V2=-999999999d0
                    elseif(dum2>k_grid(nk)-0.001d0) then
                        V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                        !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                        INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    else
                        V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                        !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evs_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evm_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    end if

                    pnt1 = (/P3, wagef/)
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagef
                    dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V3=-999999999d0
                    elseif(dum2>k_grid(nk)-0.001d0) then
                        V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                        !V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                        INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    else
                        V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                        !V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evs_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix+1,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)+1d0,k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evm_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    end if

                    if (V2 < V3) then
                        P1=P2
                    else
                        P4=P3
                    end if
                    if((P4-P1)<1d-6) exit
                end do

                Ves=V2
                kes=dum2
                nes=dum4
                ces=P2

                if(nes<0d0) then
                    Print *,'ik is',ik
                    Print *,'ix is',ix
                    Print *,'iam is',iam
                    Print *,'ium is',ium
                    Print *,'nes is',nes
                    pause
                end if
            

    !If female unemployed

                !Print *,'ix is',ix
                !Finding optimal capital by golden search
                P1=0.01d0
                P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit)/(1d0+tc))
                do
                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                    dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V2=-999999999d0
                    elseif(dum2>k_grid(nk)-0.001d0) then
                        V2=Uc(P2)
                        !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                        INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    else
                        V2=Uc(P2)
                        !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)*(1d0-deltaexp), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evs_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)*(1d0-deltaexp), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evm_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    end if

                    dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V3=-999999999d0
                    elseif(dum2>k_grid(nk)-0.001d0) then
                        V3=Uc(P3)
                        !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        pnt1=(/dum2, exp_grid(ix,T-it)*(1d0-deltaexp)/)
                        INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                        vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                        V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                    else
                        V3=Uc(P3)
                        !V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))

                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evs(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)*(1d0-deltaexp), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evs_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                        !V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,ix,iam,ium,T+1-it,ifc))                      
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T-it)*(1d0-deltaexp),k_grid,exp_grid(:,T-it+1),evm(j,:,:,iam,ium,T+1-it,ifc))
                        vnext = D_BS2VL(dum2, exp_grid(ix,T-it)*(1d0-deltaexp), KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evm_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                        V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext                      
                    end if                    



                    if (V2 < V3) then
                        P1=P2
                    else
                        P4=P3
                    end if
                    if((P4-P1)<1d-6) exit
                end do

                Vus=V2
                kus=dum2
                nus=0d0
                cus=P2

            
                if (ves >= vus) then
                    vs(j,ik,ix,iam,ium,T-it,ifc)=ves
                    cs(j,ik,ix,iam,ium,T-it,ifc)=ces
                    ks(j,ik,ix,iam,ium,T-it,ifc)=kes
                    ns(j,ik,ix,iam,ium,T-it,ifc)=nes
                else
                    vs(j,ik,ix,iam,ium,T-it,ifc)=vus
                    cs(j,ik,ix,iam,ium,T-it,ifc)=cus
                    ks(j,ik,ix,iam,ium,T-it,ifc)=kus
                    ns(j,ik,ix,iam,ium,T-it,ifc)=nus
                end if
    end do


    !do ium = 1, nu
    !    do ixd = 1, ixm
    !        do iam = 1, na
    !                    ix=min((1+(ixd-1)*4),(T-it))
    !                    vsdum(j,ik,ixd,iam,ium,T-it,ifc)=vs(j,ik,ix,iam,ium,T-it,ifc)
    !                    csdum(j,ik,ixd,iam,ium,T-it,ifc)=cs(j,ik,ix,iam,ium,T-it,ifc)
    !                    gksdum(j,ik,ixd,iam,ium,T-it,ifc)=ks(j,ik,ix,iam,ium,T-it,ifc)
    !                    nsdum(j,ik,ixd,iam,ium,T-it,ifc)=ns(j,ik,ix,iam,ium,T-it,ifc)
    !        end do
    !    end do
    !end do

    !Print *,'nmdum(1,:,1,1,1,1,T) is',nmdum(1,:,1,1,1,1,T)

    !do ium = 1, nu
    !    do ix = 1, T-it
    !        do iam = 1, na
    !                    vs(j,ik,ix,iam,ium,T-it,ifc)=LinInterp(ix*1d0,Expdum,vsdum(j,ik,:,iam,ium,T-it,ifc),ixm)
    !                    cs(j,ik,ix,iam,ium,T-it,ifc)=LinInterp(ix*1d0,Expdum,csdum(j,ik,:,iam,ium,T-it,ifc),ixm)
    !                    ks(j,ik,ix,iam,ium,T-it,ifc)=LinInterp(ix*1d0,Expdum,gksdum(j,ik,:,iam,ium,T-it,ifc),ixm)
    !                    ns(j,ik,ix,iam,ium,T-it,ifc)=LinInterp(ix*1d0,Expdum,nsdum(j,ik,:,iam,ium,T-it,ifc),ixm)
    !                if(ns(j,ik,ix,iam,ium,T-it,ifc)<0d0) then
    !                print *, 'interpolation'
    !                Print *,'ik is',ik
    !                Print *,'ix is',ix
    !                Print *,'iam is',iam
    !                Print *,'ium is',ium
    !                Print *,'iaf is',iaf
    !                Print *,'iuf is',iuf
    !                Print *,'ns(j,ik,ix,iam,ium,T-it,ifc) is',ns(j,ik,ix,iam,ium,T-it,ifc)
    !                Print *,'Expdum is',Expdum
    !                Print *,'nsdum(j,ik,ixd,iam,ium,T-it,ifc) is',nsdum(j,ik,ixd,iam,ium,T-it,ifc)
    !                pause
    !                end if
    !        end do
    !    end do
    !end do

    !Men   
    j=1
    ifc=1
    
            !Print *,'ix is',ix
            !Finding optimal capital by golden search
            wagem = wage(1,a(1,iam),dble(T-it),u(1,ium))/(1d0+t_employer)
            P1=0.01d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+wagem*(1d0-tax_labors(wagem)-tSS_employee(wagem)))/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                pnt1 = (/P2, wagem/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagem
                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)
                    !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))

                    !vnext = D_QDVAL(dum2,k_grid,evs(j,:,1,iam,ium,T+1-it,ifc))
                    pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                    INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                    !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))
                    !vnext = D_QDVAL(dum2,k_grid,evm(j,:,1,iam,ium,T+1-it,ifc))
                    INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                else
                    V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)
                    !V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))

                    !vnext = D_QDVAL(dum2,k_grid,evs(j,:,1,iam,ium,T+1-it,ifc))
                    vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evs_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                    V2=V2+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                    !V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))
                    !vnext = D_QDVAL(dum2,k_grid,evm(j,:,1,iam,ium,T+1-it,ifc))
                    vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evm_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                    V2=V2+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                end if

                pnt1 = (/P3, wagem/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagem
                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                elseif(dum2>k_grid(nk)-0.001d0) then
                    V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)
                    !V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))

                    !vnext = D_QDVAL(dum2,k_grid,evs(j,:,1,iam,ium,T+1-it,ifc))
                    pnt1=(/dum2, exp_grid(ix,T-it)+1d0/)
                    INTERP2D=evs(j,:,:,iam,ium,T+1-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                    !V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))
                    !vnext = D_QDVAL(dum2,k_grid,evm(j,:,1,iam,ium,T+1-it,ifc))
                    INTERP2D=evm(j,:,:,iam,ium,T+1-it,ifc)
                    vnext = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                    V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                else
                    V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)
                    !V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*D_CSVAL(dum2,BREAK,evs_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))

                    !vnext = D_QDVAL(dum2,k_grid,evs(j,:,1,iam,ium,T+1-it,ifc))
                    vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evs_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                    V3=V3+beta*OmegaActive(T-it)*(1d0-Probm(T-it))*vnext

                    !V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*D_CSVAL(dum2,BREAK,evm_spln_coefs(j,:,:,1,iam,ium,T+1-it,ifc))
                    !vnext = D_QDVAL(dum2,k_grid,evm(j,:,1,iam,ium,T+1-it,ifc))
                    vnext = D_BS2VL(dum2, exp_grid(ix,T-it)+1d0, KORDER, EXPORDER, K_KNOT,EXP_KNOT(:,T-it+1), nk, nexp,evm_spln_coefs(j,:,:,iam,ium,T+1-it,ifc))
                    V3=V3+beta*OmegaActive(T-it)*Probm(T-it)*vnext
                end if

                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-6) exit
            end do

            Vs(j,ik,:,iam,ium,T-it,:)=V2
            ks(j,ik,:,iam,ium,T-it,:)=dum2
            ns(j,ik,:,iam,ium,T-it,:)=dum4
            cs(j,ik,:,iam,ium,T-it,:)=P2

            if(ns(j,ik,1,iam,ium,T,1)<0d0) then
                Print *,'ik is',ik
                Print *,'iam is',iam
                Print *,'ium is',ium
                Print *,'ns(j,ik,1,iam,ium,T,1) is',ns(j,ik,1,iam,ium,T-it,1)
                pause
            end if


end subroutine Solvefirstactive
