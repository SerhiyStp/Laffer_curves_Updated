SUBROUTINE labor1 (X, F, N)
    !This subroutine contains the maximization problem in labor for a married couple when one spouses work
    USE Model_Parameters
    USE glob0
    USE Utilities
    INTEGER N
    real(8) X, F
    real(8) :: dum2,Y1
    
    if(ind2==1) then
        !Y1 is consumption to be used in the FOC
        Y1=dum
        X=max(X,0.0001d0)
        F=wagem*(1d0-(tax_labor(wagem*X)+DtSS_employee(wagem*X)))-(1d0/AE)*(wagem**2)*X*theta(1)*theta(2)*((wagem*X)*(1d0/AE))**(-1d0-theta(2))-chim*(X**etam)*Y1*(1d0+tc)
    else
        Y1=dum
        X=max(X,0.0001d0)
        F=wagem*(1d0-(tax_labor(wagem*X)+DtSS_employee(wagem*X)))-(1d0/AE)*(wagem**2)*X*theta(1)*theta(2)*((wagem*X)*(1d0/AE))**(-1d0-theta(2))-chif*(X**etaf)*Y1*(1d0+tc)
    end if
        
end subroutine labor1

subroutine labor1_hybrd(n, x, fvec, iflag)
    integer :: n, iflag
    real(8) :: x(n), fvec(n)

    call labor1(x(1), fvec(1), n)
    
end subroutine labor1_hybrd
