module defineparameters
    implicit none
    !n is number of grid points for capital, p is the number of grid points for experience, q is the number of grid points for education
    
    integer, parameter :: n=16, p=45, t=75, q=2, l=5
    
    !The F.. are fixed costs, eta is the inverse Frisch elasticity, sigma the coefficient of relative risk aversion
    !beta is the discount factor
    !r is the interest rate
    !chim and chif determines the intertemporal elasticity of substitution between consumption and leisure for men and women
    
    real(8) :: Fsm=5.0,Fsf=7.0,Fmm=5.0,Fmf=7.0,deltam=0.11,deltaf=0.11,etam=1.0/0.3,etaf=1.0/0.6,psim=0.4,psif=0.4,chim=169.0,chif=169.0
    real(8) :: beta=1.04, r=0.011, delta=1.0, eta=3.0
    
    !These variables hold the returns to experience (originally there were 3 education groups, now we only use the first 2 for men and women
    
    real(8) :: gam0=-1.4, gam1=0.0533247, gam2=-0.001456, gam3=0.0000142, marr=0.7372, initmarr=0.1259
    real(8) :: gamsc0=-1.2, gamsc1=0.0720771, gamsc2=-0.0020941, gamsc3=0.0000214
    real(8) :: gamf0=-2.18, gamf1=0.0556236, gamf2=-0.00165, gamf3=0.0000161
    real(8) :: gamscf0=-1.4, gamscf1=0.0713539, gamscf2=-0.0020398, gamscf3=0.0000185
    real(8), parameter :: pi=3.14159265
    
    !AE is the number used in the tax functions
    
    real(8) :: AE=0.37503, AEdum=0.37503, lumpsum=0.015869,lumpsumdum=0.015869, epsilon=0.1, epsilon2=0.1, unempbenefit=0.084
    
    
    !Value and policy functions for married, sinle and retired people
    
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable :: gkM,gnM,gnfM,cM,VM,temp3M,BREAKM
    real(8), dimension (:,:,:,:,:,:,:,:,:), allocatable :: CSCOEFM
    real(8), dimension (1:2,1:45,1:n,1:p,1:q,1:l) :: gkS,cS,VS,temp3S,BREAKS,gnS
    real(8), dimension (1:2,1:45,1:n,1:p,1:q,1:l,1:4) :: CSCOEFS
    real(8), dimension (1:3,1:30,1:n) :: gkr,cr,Vr,temp3r,BREAKr
    real(8), dimension (1:3,1:30,1:n,1:4) :: CSCOEFr
    real(8), dimension (1:2,1:2,1:l) :: U, probu
      
    real(8), dimension (1:2,1:2,1:l,1:l) :: transu
    real(8), dimension (1:q,1:q) :: transa
    
    !K holds the capital grid
    
    real(8), dimension (1:n) :: K
    
    real(8), dimension (1:2,1:q) :: H
    real(8), dimension (1:p) :: Theta
    real(8), dimension (1:q) :: A, proba
    
    !Random1 and Random2 holds the shocks used in the simulation
    
    real(8), dimension (1:10,1:10000,1:45) :: Random1, Random2, Random3, Random4
    real(8), dimension (1:10,1:10000) :: Random5, Random6
    
    !These vectors holds different simulation statistics
    
    real(8), dimension (1:10000,1:46,1:6) :: Malesim1, Femalesim1
    real(8), dimension (1:10000,1:46,1:5) :: XM1,XF1
    integer, dimension (1:10000,1:46,1:4) :: expM1,expF1
    real(8), dimension (1:10000,1:31,1:2) :: XMR1,XFR1
    real(8), dimension (1:100000,1:31,1:2) :: XMR2,XFR2
    real(8), dimension (1:100000,1:46,1:6) :: Malesim2, Femalesim2
    real(8), dimension (1:100000,1:46,1:5) :: XM2,XF2
    integer, dimension (1:100000,1:46,1:4) :: expM2,expF2
    
    
    real(8) :: w=1.0, psi, psi2, zet, epsa, epsu, ma, alpha, siga, rhoa,Knot(n+3),Hmknot(p+3),Hfknot(p+3)
    
    !kmin and kmax are used to assign min and max value for capital. Probm and Probd holds probabilities of marriage and divorce, pprob holds the distribution of single people
    
    real(8) :: kmin, kmax, hmin, hmax, inck, Probm(45), Probd(45), mix(10000,2),mix2(10000,2),mix3(10000,2),mix4(10000,2),tp(1,2),tp2(1,2),pprob(45,46,8,2,l,2),pprob2(45,46,8,2,l,2),pprob3(45,46,8,2,l,2)
    integer, dimension (1:4) :: sim
    integer :: i,d2, iter, count5
    
    !US tax

    !real(8), parameter :: t0=1.727408, t1=-6.44973, t2=8.994808, t3=-4.999817, t4=0.9875019, tc=0.084
    !real(8), parameter :: t0M=2.16239, t1M=-7.301506, t2M=9.221961, t3M=-4.736035, t4M=0.8718943

    !German tax

    !real(8), parameter :: t0=-6.582745, t1=19.08046, t2=-19.22463, t3=8.580912, t4=-1.430125,tc=0.16
    !real(8), parameter :: t0M=-0.5409343, t1M=-0.9886915, t2M=4.474231, t3M=-3.421762, t4M=0.7909097

    !Irish tax

    !real(8), parameter :: t0=-1.75284, t1=2.625375, t2=0.1463597, t3=-1.13193, t4=0.3456357, tc=0.21
    !real(8), parameter :: t0M=1.612143, t1M=-6.871639, t2M=9.391285, t3M=-4.898055, t4M=0.8901651

    !Norwegian tax

    !real(8), parameter :: t0=2.335783, t1=-8.6315, t2=11.83152, t3=-6.471281, t4=1.25354, tc=0.24
    !real(8), parameter :: t0M=-5.335858, t1M=14.96881, t2M=-15.43612, t3M=7.362051, t4M=-1.335945

    !Swiss tax

    !real(8), parameter :: t0=-1.4185, t1=5.181097, t2=-6.488006, t3=3.771889, t4=-0.8035895, tc=0.076
    !real(8), parameter :: t0M=-16.09581, t1M=48.2164, t2M=-53.35435, t3M=26.20165, t4M=-4.78368

    !French tax

    !real(8), parameter :: t0=0.7157418, t1=-2.514716, t2=3.64648, t3=-1.88936, t4=0.3320441, tc=0.196
    !real(8), parameter :: t0M=-0.4677592, t1M=2.062677, t2M=-2.743411, t3M=1.820481, t4M=-0.4305004

    !Spanish tax

    !real(8), parameter :: t0=-2.640157, t1=7.853874, t2=-8.641411, t3=4.527437, t4=-0.9025463, tc=0.16
    !real(8), parameter :: t0M=-2.811092, t1M=8.034616, t2M=-8.401096, t3M=4.023208, t4M=-0.7058137

    !Portugese tax

    real(8), parameter :: t0=2.604929, t1=-9.655736, t2=12.78917, t3=-6.821912, t4=1.293703, tc=0.17
    real(8), parameter :: t0M=3.907341, t1M=-12.23614, t2M=13.88106, t3M=-6.514196, t4M=1.101643

    !Luxembourg tax

    !real(8), parameter :: t0=0.0866169, t1=-2.91607, t2=6.525497, t3=-4.37144, t4=0.9543883, tc=0.15
    !real(8), parameter :: t0M=-0.0840795, t1M=-2.859591, t2M=6.036954, t3M=-3.722134, t4M=0.7483014

    !Italian tax

    !real(8), parameter :: t0=-1.555522, t1=2.965259, t2=-0.9916236, t3=-0.3076185, t4=0.1599916, tc=0.2
    !real(8), parameter :: t0M=-4.143618, t1M=11.07723, t2M=-10.77931, t3M=4.893096, t4M=-0.8552848

    !Austrian tax

    !real(8), parameter :: t0=-5.626168, t1=16.19854, t2=-16.39948, t3=7.397988, t4=-1.250442, tc=0.2
    !real(8), parameter :: t0M=5.591343, t1M=-19.17492, t2M=24.16844, t3M=-12.8056, t4M=2.451535

    !Belgian tax

    !real(8), parameter :: t0=-4.587984, t1=13.62661, t2=-14.19084, t3=6.823648, t4=-1.24974, tc=0.21
    !real(8), parameter :: t0M=-6.1645, t1M=18.35908, t2M=-19.74126, t3M=9.674913, t4M=-1.783321

    !UK tax

    !real(8), parameter :: t0=-0.3775787, t1=0.2900424, t2=1.07663, t3=-0.9579886, t4=0.2236049, tc=0.175
    !real(8), parameter :: t0M=-4.01828, t1M=11.29697, t2M=-11.56235, t3M=5.448592, t4M=-0.9772443

    !Swedish tax

    !real(8), parameter :: t0=5.645098, t1=-18.75109, t2=23.36599, t3=-12.24517, t4=2.322895, tc=0.25
    !real(8), parameter :: t0M=-3.314906, t1M=9.808722, t2M=-10.54196, t3M=5.343565, t4M=-1.032559

    !Finish tax

    !real(8), parameter :: t0=-1.387284, t1=2.706099, t2=-0.9767094, t3=-0.0860593, t4=0.0717587, tc=0.22
    !real(8), parameter :: t0M=-5.062344, t1M=13.81237, t2M=-13.90515, t3M=6.559641, t4M=-1.186956

    !Danish tax

    !real(8), parameter :: t0=0.1422833, t1=-2.357568, t2=5.737164, t3=-3.968169, t4=0.8855884, tc=0.25
    !real(8), parameter :: t0M=-28.1151, t1M=82.6305, t2M=-89.35836, t3M=42.80431, t4M=-7.638321

    !Dutch tax

    !real(8), parameter :: t0=1.126893, t1=-4.322011, t2=6.331867, t3=-3.487033, t4=0.6651015, tc=0.19
    !real(8), parameter :: t0M=-10.87501, t1M=32.46464, t2M=-35.65958, t3M=17.52148, t4M=-3.214915

    !Greek tax

    !real(8), parameter :: t0=-5.55185, t1=14.76655, t2=-14.7313, t3=6.887032, t4=-1.237959, tc=0.18
    !real(8), parameter :: t0M=-15.38484, t1M=48.03587, t2M=-55.50611, t3M=28.30343, t4M=-5.32562

end module defineparameters