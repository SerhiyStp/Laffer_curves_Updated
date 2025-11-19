module ValueFunctions
    
    implicit none
    
contains
    
    function val_mar_T(kp, expm, expf, val_fn_ptr) result(val)
        use PolicyFunctions
        real(8), intent(in) :: kp
        real(8), intent(in) :: expm
        real(8), intent(in) :: expf        
        type(policy_fn_3d), pointer, intent(in) :: val_fn_ptr
        real(8) :: val_tmp
        real(8) :: val_tmp2
        real(8) :: val
        
        if (kp > k_grid(nk)-0.001d0) then
            val_tmp = val_fn_ptr%lin_interp([kp, expf, expm])  
        else
            val_tmp = val_fn_ptr%eval([kp, expf, expm])
        end if
        
        val = val_tmp
    end function val_mar_T    

    function val_mar(kp, expm, expf, iage, iam, ium, iaf, iuf, ifc, ifcm)
        use bspline_sub_module, only: db3val
        use PolicyFunctions
        real(8), intent(in) :: kp
        real(8), intent(in) :: expm
        real(8), intent(in) :: expf
        integer, intent(in) :: iage, iam, ium, iaf, iuf, ifc, ifcm
        integer :: idx, idy, idz, iloy, iloz
        integer :: inbvx, inbvy, inbvz
        integer :: iflag
        real(8) :: ww2(ky,kz),ww1(kz),ww0(3*max(kx,ky,kz))      
        real(8) :: val_mar
        
        call db3val(kp, expm, expf, idx, idy, idz, &
                    tx, ty(:,iage+1), tz(:,iage+1), &
                    nk, nexp, nexp, kx, ky, kz, &
                    ev_bspl(iam,ium,iaf,iuf,ifc,ifcm)%coefs, val_mar, iflag, &
                    inbvx, inbvy, inbvz, iloy, iloz, ww2, ww1, ww0, extrap=.true.)          
        
    end function val_mar
    
end module ValueFunctions
