module GoldenSearch_mod
    
    use PolicyFunctions, only: policy_fn_3d
    
    implicit none


    real(8), parameter :: c_gs = (sqrt(5.0) - 1.0)/2.0
    real(8), parameter :: r_gs = 1.0 - c_gs
    
    type context
        ! This would be adjusted according to the problem to be solved.
        real(8), pointer :: xaux(:)
        real(8), pointer :: yaux(:)
        type(policy_fn_3d), pointer :: val_fn_ptr  
    end type context


contains

    subroutine goldensearch(f, lo, hi, xmax, fmax, data_aux, tol_opt, verbose_opt)
        real(8), intent(in) :: lo, hi
        interface
            function func(x, data_aux) result(y)
                import context
                real(8), intent(in) :: x
                type(context) :: data_aux
                real(8) :: y
            end function func
        end interface
        procedure(func) :: f
        type(context) :: data_aux
        real(8), intent(out) :: xmax
        real(8), optional, intent(in) :: tol_opt
        logical, optional, intent(in) :: verbose_opt
        logical :: verbose
        real(8) :: a, b
        real(8) :: tol
        real(8) :: x1, x2
        real(8) :: fx1, fx2
        real(8) :: fa, fb
        real(8) :: fmax


        if ( .not. present(tol_opt) ) then
            tol = 1d-5
        else
            tol = tol_opt
        end if
        if ( .not. present(verbose_opt) ) then
            verbose = .false.
        else
            verbose = verbose_opt
        end if


        !print *, f(a)
        !print *, f(b)
        a = lo
        b = hi
        ! fa = f(a)
        ! fb = f(b)
        
        x1 = c_gs*a + r_gs*b
        x2 = r_gs*a + c_gs*b
        fx1 = f(x1, data_aux)
        fx2 = f(x2, data_aux)

        if (verbose) then
            print *, '========================================='
            print *, '   x1     x2     f(x1)   f(x2)     b-a '    
            print *, '========================================='
        end if
        do while (abs(b-a) > tol)
            if (fx1 > fx2) then
                b = x2
                x2 = x1
                fx2 = fx1
                x1 = c_gs*a + r_gs*b
                fx1 = f(x1, data_aux)
            else
                a = x1
                x1 = x2
                fx1 = fx2
                x2 = r_gs*a + c_gs*b
                fx2 = f(x2, data_aux)
            end if
            if (verbose) then
                write(*,'(5f8.5)'), x1, x2, fx1, fx2, b-a
            end if
        end do
        
        xmax = (a + b)/2.0d0
        fmax = f(xmax, data_aux)

        ! Check for corner solution:
        !if (fa > fsol .and. fa > fb) then
        !    soln = a
        !else if (fb > fsol .and. fb > fa) then
        !    soln = b
        !end if

        
    end subroutine goldensearch    
    
    
    subroutine goldensearch_h(f, lo, hi, xmax, fmax, data_aux, tol_opt)
        real(8), intent(in) :: lo, hi
        interface
            function func(x, data_aux) result(y)
                import context
                real(8), intent(in) :: x
                type(context) :: data_aux
                real(8) :: y
            end function func
        end interface
        procedure(func) :: f
        type(context) :: data_aux
        real(8), intent(out) :: xmax
        real(8), optional, intent(in) :: tol_opt
        logical :: verbose
        real(8) :: p1, p2, p3, p4 
        real(8) :: tol
        real(8) :: fp2, fp3
        real(8) :: fmax


        if ( .not. present(tol_opt) ) then
            tol = 1d-5
        else
            tol = tol_opt
        end if


        p1 = lo
        p4 = hi

        ! do while (abs(p4-p1) > tol)
        do
            p2 = c_gs*p1 + r_gs*p4
            p3 = r_gs*p1 + c_gs*p4
            fp2 = f(p2, data_aux)
            fp3 = f(p3, data_aux)
            if (fp2 < fp3) then
                p1 = p2
            else
                p4 = p3
            end if
            if (abs(p4-p1) < tol) exit
        end do
        
        xmax = p3 
        fmax = fp2

    end subroutine goldensearch_h       
    
    
end module GoldenSearch_mod
