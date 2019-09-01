module hybrd_wrapper
    
    ! This module is supposed to make using HYBRD non-linear solver earier
    ! Two options - either use setHybrParams and then getHybrSoln,
    ! or use the Hybrd_problem type
    
    implicit none
    
    integer :: ldfjac, lr
    integer :: ml, muh
    integer :: maxfev 
    integer :: mode 
    integer :: nfev
    real(8) :: epsfcn 
    real(8) :: factor
    real(8) :: tol
    
    type Hybrd_problem 
                
        integer :: neqs 
        integer :: ldfjac, lr
        integer :: ml, mu
        integer :: maxfev 
        integer :: mode 
        integer :: nfev
        real(8) :: epsfcn 
        real(8) :: factor
        real(8) :: tol
        
    contains
        procedure :: init
        procedure :: getSoln
        
    end type Hybrd_problem

contains
    
    subroutine setHybrParams(n)
        integer :: n
        
        ldfjac = n
        lr = n*(n+1)/2
        ml = n - 1
        muh = n - 1
        maxfev = 100
        mode = 1
        epsfcn = 1d-6
        factor = 100d0
        tol = 1d-6
    
    end subroutine setHybrParams

    subroutine getHybrSoln(guess, n, sub1, sol, fnorm)
        !use EqCondit
        integer :: n
        integer :: nprint, info, nfev
        real(8) :: guess(n), sol(n)
        real(8), dimension(n) :: x, fvec, diag
        real(8) :: fjac(ldfjac,n), r(lr)
        real(8), dimension(n) :: qtf, wa1, wa2, wa3, wa4
        real(8) :: fnorm
        external :: sub1
        !external fcn
        
        x = guess
        nprint = 0
        call hybrd(sub1, n, x, fvec, tol, maxfev, ml, muh, epsfcn, &
                   diag, mode, factor, nprint, info, nfev, fjac, &
                   ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4)	
        sol = x
        fnorm = sqrt(sum(fvec**2d0))		
    end subroutine getHybrSoln
    
    subroutine init(this, n)
        class(Hybrd_problem) :: this
        integer :: n
        
        this%neqs = n
        this%ldfjac = n
        this%lr = n*(n+1)/2
        this%ml = n - 1
        this%mu = n - 1
        this%maxfev = 100
        this%mode = 1
        this%epsfcn = 1d-8
        this%factor = 100d0
        this%tol = 1d-8
        
    end subroutine init

    subroutine getSoln(this, guess, n, sub1, sol, fnorm)
        class(Hybrd_problem) :: this
        integer :: n
        integer :: nprint, info, nfev
        real(8) :: guess(n), sol(n)
        real(8), dimension(n) :: x(n), fvec(n), diag(n)
        real(8) :: fjac(this%ldfjac,n), r(this%lr)
        real(8), dimension(n) :: qtf(n), wa1(n), wa2(n), wa3(n), wa4(n)
        real(8) :: fnorm
        external :: sub1
        
        x = guess
        nprint = 0
        call hybrd(sub1, n, x, fvec, this%tol, this%maxfev, this%ml, this%mu, this%epsfcn, &
                   diag, this%mode, this%factor, nprint, info, nfev, fjac, &
                   this%ldfjac, r, this%lr, qtf, wa1, wa2, wa3, wa4)	
        sol = x
        fnorm = sqrt(sum(fvec**2d0))		
    end subroutine getSoln
    
end module hybrd_wrapper
