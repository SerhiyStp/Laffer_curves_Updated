module Model_Parameters

    implicit none
    
    ! Tax policy options
    integer, parameter :: opt_G = 1
    integer, parameter :: opt_Lumpsum = 2
    
    ! Labels
    integer, parameter :: MEN=1
    integer, parameter :: WOMEN=2  
    integer, parameter :: LFP_1=1
    integer, parameter :: LFP_0=2
    integer, parameter :: LFP_M1=LFP_1
    integer, parameter :: LFP_M0=LFP_0
    integer, parameter :: LFP_F1=LFP_1
    integer, parameter :: LFP_F0=LFP_0    

    real(8), parameter :: etam    = 1d0/0.3d0  ! Inverse Frisch elasticity men
    real(8), parameter :: etaf    = 1d0/0.6d0  ! Inverse Frisch elasticity women
    real(8), parameter :: mu   = 0.02d0       ! The growth rate of the economy
    !real(8), parameter :: beta   = 0.97d0    !Discount factor
    real(8), parameter :: beta   = 1.0019935d0    !Discount factor
    
    integer, parameter :: testing = 1
    integer, parameter :: Tret   = 5*testing + 36*(1-testing) !36 !5         ! Years in retirement  (65-100)
    integer, parameter :: T      = 6*testing + 45*(1-testing) !45 !6         ! Years of active life (20-64)
    integer, parameter :: nfc   = 5         ! Number of fixed costs
    integer, parameter :: nfcm   = 5         ! Number of fixed costs men
    real(8), parameter :: tc    = 0.05d0     ! Consumption tax
    !real(8), parameter :: t_employer=0.0765d0, t_employee=0.0765d0 !Social security tax paid by employer and employee

    ! Social security:
    real(8), parameter :: t_employer=0d0
    real(8), parameter :: t_employee1 = 1d0-(1d0-0.0765d0)/(1d0+0.0765d0)
    real(8), parameter :: t_employee2 = 1d0-(1d0-0.0145d0)/(1d0+0.0145d0)
    real(8), parameter :: SocSecCap = 2.323967451d0
    real(8), parameter :: ss_eps = 0.01d0
    real(8), parameter :: x1_ss = SocSecCap - ss_eps
    real(8), parameter :: x2_ss = SocSecCap + ss_eps
    real(8) :: ss_coefs(4)


    real(8) :: minhours=0.001d0
    real(8), parameter :: tk     = 0.36d0     ! Capital tax
    real(8), parameter :: sigma  = 4d0        ! Risk aversion parameter
    real(8), parameter :: alpha  = 1d0/3d0     ! Capital share
    real(8), parameter :: delta  = 0.0988d0-mu     ! Capital depreciation
    real(8), parameter :: deltaexp  = 0.000d0     !Depreciation of experience
    real(8) :: chim= 31.40d0,  chims=75.50d0       ! Disutility from work men
    real(8) :: chif= 7.42d0,   chifs= 14.77d0       ! Disutility from work women
    real(8), parameter :: sigma_um=0.32228727D0, rho_um=0.3959915D0 ! Parameters governing the process for the  idiosyncratic shock, men
    real(8), parameter :: sigma_am=0.31469361d0, rho_am=0D0     ! Stdev of ability, men
    real(8), parameter :: sigma_uf=0.31004311d0, rho_uf=0.339295 ! Parameters governing the process for the  idiosyncratic shock, women
    real(8), parameter :: sigma_af=0.38475527d0, rho_af=0d0    ! Stdev of ability, women
    real(8), parameter :: sigma_fcm=0.1707d0, rho_fcm=0d0  !Stdev and persistence of fixed costs, married women
    real(8), parameter :: sigma_fcs=0.5590d0, rho_fcs=0d0    !Stdev and persistence of fixed costs, single women
    real(8), parameter :: mu_fcm=0.2400d0, mu_fcs=0.0105d0   ! Mean fixed cost of LFP, women
    real(8) :: mu_fcm1=-0.022782d0, mu_fcm2=0.0004034d0          !First and 2nd order term of married female cost of working
    real(8), parameter :: sigma_fcmm=0.0301d0, rho_fcmm=0d0  !Stdev and persistence of fixed costs, married men
    real(8), parameter :: sigma_fcsm=0.4119d0, rho_fcsm=0d0    !Stdev of fixed costs, single men
    real(8), parameter :: mu_fcmm=0.2008d0, mu_fcsm=0.5565d0   !Mean fixed cost of LFP, men
    real(8) :: gamma(2,3), gamma0, gamma0f  ! Returns to experience parameters
    real(8) :: AE                           ! Average earnings
    real(8) :: match=0.09102d0                 ! Assortative mating parameter
    real(8) :: theta(2), thetas(2)            ! Labor tax parameters
    real(8) :: tax_level_scale=1.0d0, tax_prog_scale=1.0d0  ! Parameters to scale tax- level and progressivity
    real(8), parameter :: ybar = 1.2d0                ! Upper limit for SS tax
    real(8) :: epsilon=1d0, epsilon2=1d0, epsilon3=1d0, epsilon4=1d0, epsilon5=1d0, epsilon6=1d0
    real(8), parameter :: cons_floor = 1d-9   ! Lower limit for consumption of unemployed
    integer, parameter :: na     = 5 !3          ! Number of permanent ability levels
    integer, parameter :: nu     = 5 !3          ! Number of transitory productivity shocks
    integer, parameter :: nexp   = 6 !4          ! Number of gridpoints for experience
    integer, parameter :: nc     = 100        ! Number of gridpoints in "static" hours function
    integer, parameter :: nw     = 100        ! Number of gridpoints in "static" hours function
    integer, parameter :: nk     = 16 !8         ! Number of gridpoints over capital

    integer, parameter :: nsim   = 10000         ! Number of households used for simulation
    integer, parameter :: nsim2   = 16         ! Number of simulations
    integer, parameter :: KORDER   = 3         ! Order of spline in K-dimension
    integer, parameter :: EXPORDER   = 3       ! Order of spline in experience

    real(8), dimension (:,:), allocatable :: a, Prob_a      !Vectors with the value of ability and the probability of each ability
    real(8), dimension (:,:), allocatable :: Fc, Prob_fc, Fcage, Fcm,  Prob_fcm    !Vectors with the value of fixed costs and the probability of each fixed cost
    real(8), dimension (:,:), allocatable :: u, Prob_u      ! Vectors with the value of the idiosyncratic shock and the probability of each shock
    real(8), dimension (:,:,:), allocatable :: trans_u     ! Transition matrix for the shock u
    real(8), dimension (:,:,:), allocatable :: trans_a     ! Dummy transition matrix for a (used to compute the unconditional probability of each a)
    real(8), dimension (:,:,:), allocatable :: trans_fc, trans_fcm ! Dummy transition matrix for fc (used to compute the unconditional probability of each a)

    real(8) :: Gamma_redistr = 2d0*0.01737d0
    real(8) :: Psi0=0.2396d0, w05=1.70423761746912d0, w2=4.69173016308841d0
    real(8) :: Unemp_benefit, lumpsum=2d0*0.0390d0, lumpsumdum=2d0*0.09184923463949926d0

    real(8), dimension (:), allocatable :: OmegaRet, OmegaRet2, OmegaActive, Probm, Probd
    real(8), dimension (:), allocatable :: WeightRet, WeightActive
    real(8) :: ratio=4.39439635920971d0, ratiodum=4.28978442524046d0  !Ratio between capital and labor
    real(8) :: r                                    !Real interest rate
    !real(8) :: debttoGDP=0.6185043538d0             ! Government debt as % of GDP
    real(8) :: debttoGDP=0d0             ! Government debt as % of GDP
    real(8) :: milspendtoGDP=0.03892d0              ! Military spanding as % of GDP
    real(8) :: w                                    ! Wages per efficiency unit

    real(8), dimension (:,:,:,:), allocatable :: Sim1m,Sim1f,exp2m,exp2f ! Array to hold simulation results

    integer, dimension (:,:,:,:), allocatable :: exp1m,exp1f,expR1f,expR1m ! Array to hold simulation results

    real(8), dimension (:,:,:,:), allocatable :: SimR1m, SimR1f ! Array to hold simulation results for retired
    real(8), dimension (:,:,:), allocatable :: Random3m, Random3f, marstatm, marstatf, partshock, partshock2
    real(8), dimension (:,:), allocatable :: Random1m, Random2m, Random1f, Random2f, marstatm_init, marstatf_init
    integer :: it,iter

end module Model_Parameters
