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
    F(1)=wagem*(1d0-(tax_labor(wagem*X(1)+wagef*X(2))+DtSS_employee(wagem*X(1))))-(1d0/AE)*wagem*(wagem*X(1)+wagef*X(2))*theta(1)*theta(2)*((wagem*X(1)+wagef*X(2))*(1d0/AE))**(-1d0-theta(2))-chim*(X(1)**etam)*Y1*(1d0+tc)
    F(2)=wagef*(1d0-(tax_labor(wagem*X(1)+wagef*X(2))+DtSS_employee(wagef*X(2))))-(1d0/AE)*wagef*(wagem*X(1)+wagef*X(2))*theta(1)*theta(2)*((wagem*X(1)+wagef*X(2))*(1d0/AE))**(-1d0-theta(2))-chif*(X(2)**etaf)*Y1*(1d0+tc)
    
end subroutine labor2

SUBROUTINE labor2_gz (X, F, N)
    !This subroutine contains the maximization problem in labor for a married couple when both spouses work
    USE Model_Parameters
    USE glob0
    USE Utilities
    INTEGER N
    real(8) X(N), F(N)
    real(8) :: dum2,Y1
    real(8) :: hm, hf, am, af, am_kt, af_kt, ytot
       
    !Y1 is consumption to be used in the FOC
    Y1=dum
    am = x(1)
    af = x(2)
    hm = 1.0d0 - max(0.0d0,-am)
    !hm = max(hm,0.0001d0)
    hf = 1.0d0 - max(0.0d0,-af)
    !hf = max(hf,0.0001d0)
    am_kt = max(0.0d0,am)
    af_kt = max(0.0d0,af)
    ytot = wagem*hm+wagef*hf
    F(1)=wagem*(1d0-(tax_labor(ytot)+DtSS_employee(wagem*hm)))-(1d0/AE)*wagem*(ytot)*theta(1)*theta(2)*((ytot)*(1d0/AE))**(-1d0-theta(2))-chim*(hm**etam)*Y1*(1d0+tc)-am_kt*Y1*(1d0+tc) + 100.0d0*max(0.0d0,-hm) 
    F(2)=wagef*(1d0-(tax_labor(ytot)+DtSS_employee(wagef*hf)))-(1d0/AE)*wagef*(ytot)*theta(1)*theta(2)*((ytot)*(1d0/AE))**(-1d0-theta(2))-chif*(hf**etaf)*Y1*(1d0+tc)-af_kt*Y1*(1d0+tc) + 100.0d0*max(0.0d0,-hf)
    
end subroutine labor2_gz

subroutine labor2_hybrd(n, x, fvec, iflag)
    integer :: n, iflag
    real(8) :: x(n), fvec(n)

    call labor2(x, fvec, n)
    
end subroutine labor2_hybrd

subroutine labor2_gz_hybrd(n, x, fvec, iflag)
    integer :: n, iflag
    real(8) :: x(n), fvec(n)

    call labor2_gz(x, fvec, n)
    
end subroutine labor2_gz_hybrd
