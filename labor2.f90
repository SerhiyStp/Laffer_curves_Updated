SUBROUTINE labor2 (X, F, N)
    !This subroutine contains the maximization problem in labor for a married couple when both spouses work
    USE Model_Parameters
    USE glob0
    USE Utilities
    INTEGER N
    real(8) X(N), F(N)
    real(8) :: dum2,Y1
       
    !Y1 is consumption to be used in the FOC
    Y1=dum
    X(1)=max(X(1),0.0001d0)
    X(2)=max(X(2),0.0001d0)
    F(1)=wagem*(1d0-(tax_labor(wagem*X(1)+wagef*X(2))+tSS_employee(wagem*X(1)+wagef*X(2))))-(1d0/AE)*wagem*(wagem*X(1)+wagef*X(2))*theta(1)*theta(2)*((wagem*X(1)+wagef*X(2))*(1d0/AE))**(-1d0-theta(2))-chim*(X(1)**etam)*Y1*(1d0+tc)
    F(2)=wagef*(1d0-(tax_labor(wagem*X(1)+wagef*X(2))+tSS_employee(wagem*X(1)+wagef*X(2))))-(1d0/AE)*wagef*(wagem*X(1)+wagef*X(2))*theta(1)*theta(2)*((wagem*X(1)+wagef*X(2))*(1d0/AE))**(-1d0-theta(2))-chif*(X(2)**etaf)*Y1*(1d0+tc)
    
end subroutine labor2

subroutine labor2_hybrd(n, x, fvec, iflag)
    integer :: n, iflag
    real(8) :: x(n), fvec(n)

    call labor2(x, fvec, n)
    
end subroutine labor2_hybrd
