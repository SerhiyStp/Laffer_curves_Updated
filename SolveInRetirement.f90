    subroutine SolveInRetirement(ik)
        use Model_Parameters
        use PolicyFunctions
        use Utilities
        !USE LCONF_INT
        !USE CSVAL_INT
        implicit none
        integer, INTENT(IN) :: ik
        real(8) :: c2, v2
        real(8) :: d1, d2,P1,P2,P3,P4,v3
        integer :: it2
        
        real(8) :: test1, test2
        real(8) :: vnext_p2, vnext_p3

        !Married    

        if(it==1) then
            ! Solve the very last period problem:           
            k_ret(ik,Tret) = 0d0
            c_ret(ik,Tret) = ((k_grid(ik) + Gamma_redistr)*(1d0 + r*(1d0-tk)) + Psi_pension+lumpsum)/(1d0+tc)
            v_ret(ik,Tret) = Uc(c_ret(ik,Tret)) 
        else
            ! Solve the Tret-1 to 1st period of retirement:
            !Finding optimal capital by golden search
            P1=k_grid(1)
            P4=k_grid(ik)
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
                
                !test1 = D_CSVAL(P2,BREAK,ev_spln_coefs_ret(:,:,Tret-it+2))
                !test2 = ppp_csval(P2, k_grid, nk, p_ev_spln_coefs_ret(:,:,Tret-it+2))
                vnext_p2 = ppp_csval(P2, k_grid, nk, ev_spln_coefs_ret(:,:,Tret-it+2))
                vnext_p3 = ppp_csval(P3, k_grid, nk, ev_spln_coefs_ret(:,:,Tret-it+2))
                
                V2=Uc(((k_grid(ik) + Gamma_redistr)*(1d0+r*(1-tk))+Psi_pension+lumpsum-P2*(1d0+mu))/(1.0+tc))                       
                !V2=V2+beta*OmegaRet(Tret-it+1)*LinInterp(P2,k_grid,v_ret(:,Tret-it+2),nk)
                !V2=V2+beta*OmegaRet(Tret-it+1)*D_CSVAL(P2,BREAK,ev_spln_coefs_ret(:,:,Tret-it+2))
                V2=V2+beta*OmegaRet(Tret-it+1)*vnext_p2
                if(k_grid(ik)*(1d0+r*(1-tk))+Psi_pension+lumpsum-P2*(1d0+mu)<0.01d0) then
                    V2=-999999999.0
                end if

                !test1 = D_CSVAL(P3,BREAK,ev_spln_coefs_ret(:,:,Tret-it+2))
                !test2 = ppp_csval(P3, k_grid, nk, p_ev_spln_coefs_ret(:,:,Tret-it+2))

                V3=Uc(((k_grid(ik) + Gamma_redistr)*(1d0+r*(1-tk))+Psi_pension+lumpsum-P3*(1d0+mu))/(1.0+tc))                       
                !V3=V3+beta*OmegaRet(Tret-it+1)*LinInterp(P3,k_grid,v_ret(:,Tret-it+2),nk)
                !V3=V3+beta*OmegaRet(Tret-it+1)*D_CSVAL(P3,BREAK,ev_spln_coefs_ret(:,:,Tret-it+2))
                V3=V3+beta*OmegaRet(Tret-it+1)*vnext_p3
                if(k_grid(ik)*(1d0+r*(1-tk))+Psi_pension+lumpsum-P3*(1d0+mu)<0.01d0) then
                    V3=-999999999.0
                end if       
                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-8) exit
            end do
            v_ret(ik,Tret-it+1)=V2
            k_ret(ik,Tret-it+1)=P2
            c_ret(ik,Tret-it+1)=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1-tk))+Psi_pension+lumpsum-P2*(1d0+mu))/(1d0+tc)
        end if

        !Single

        if(it==1) then
            ! Solve the very last period problem:           
            ks_ret(ik,Tret) = 0d0
            cs_ret(ik,Tret) = ((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0 + r*(1d0-tk)) + Psi_pension*0.5d0+lumpsum*0.5d0)/(1d0+tc)
            vs_ret(ik,Tret) = Uc(cs_ret(ik,Tret)) 
        else
            ! Solve the Tret-1 to 1st period of retirement:
            !Finding optimal capital by golden search
            P1=k_grid(1)
            P4=k_grid(ik)
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                !test1 = D_CSVAL(P2,BREAK,evs_spln_coefs_ret(:,:,Tret-it+2))
                !test2 = ppp_csval(P2, k_grid, nk, p_evs_spln_coefs_ret(:,:,Tret-it+2))
                vnext_p2 = ppp_csval(P2, k_grid, nk, evs_spln_coefs_ret(:,:,Tret-it+2))
                vnext_p3 = ppp_csval(P3, k_grid, nk, evs_spln_coefs_ret(:,:,Tret-it+2))


                V2=Uc(((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P2*(1d0+mu))/(1.0+tc))                       
                !V2=V2+beta*OmegaRet(Tret-it+1)*LinInterp(P2,k_grid,v_ret(:,Tret-it+2),nk)
                !V2=V2+beta*OmegaRet(Tret-it+1)*D_CSVAL(P2,BREAK,evs_spln_coefs_ret(:,:,Tret-it+2))
                V2=V2+beta*OmegaRet(Tret-it+1)*vnext_p2
                if((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P2*(1d0+mu)<0.01d0) then
                    V2=-999999999.0
                end if
                V3=Uc(((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P3*(1d0+mu))/(1.0+tc))                       
                !V3=V3+beta*OmegaRet(Tret-it+1)*LinInterp(P3,k_grid,v_ret(:,Tret-it+2),nk)
                !V3=V3+beta*OmegaRet(Tret-it+1)*D_CSVAL(P3,BREAK,evs_spln_coefs_ret(:,:,Tret-it+2))
                V3=V3+beta*OmegaRet(Tret-it+1)*vnext_p3
                if((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P3*(1d0+mu)<0.01d0) then
                    V3=-999999999.0
                end if       
                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-8) exit
            end do
            vs_ret(ik,Tret-it+1)=V2
            ks_ret(ik,Tret-it+1)=P2
            cs_ret(ik,Tret-it+1)=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P2*(1d0+mu))/(1d0+tc)
        end if

    end subroutine SolveInRetirement
