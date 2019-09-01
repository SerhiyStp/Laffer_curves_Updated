SUBROUTINE labors (X, F, N)
    !This subroutine contains the maximization problem in labor for singles
    USE Model_Parameters
    USE glob0
    USE Utilities
    INTEGER N
    real(8) X(N), F(N)
    real(8) :: dum2,Y1

    !Y1 is consumption to be used in the FOC
    if(ind2==1) then
        Y1=dum
        X(1)=max(X(1),0.0001d0)
        F(1)=wagem*(1d0-(tax_labors(wagem*X(1))+tSS_employee(wagem*X(1))))-(1d0/AE)*(wagem**2)*X(1)*thetas(1)*thetas(2)*((wagem*X(1))*(1d0/AE))**(-1d0-thetas(2))-chims*(X(1)**etam)*Y1*(1d0+tc)
    else
        Y1=dum
        X(1)=max(X(1),0.0001d0)
        F(1)=wagem*(1d0-(tax_labors(wagem*X(1))+tSS_employee(wagem*X(1))))-(1d0/AE)*(wagem**2)*X(1)*thetas(1)*thetas(2)*((wagem*X(1))*(1d0/AE))**(-1d0-thetas(2))-chifs*(X(1)**etaf)*Y1*(1d0+tc)
    end if
    
end subroutine labors


subroutine labors_hybrd(n, x, fvec, iflag)
    integer :: n, iflag
    real(8) :: x(n), fvec(n)

    call labors(x, fvec, n)
    
end subroutine labors_hybrd
