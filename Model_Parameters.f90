module Model_Parameters
        
        implicit none
        
        real(8), parameter :: etam    = 1d0/0.4d0  ! Inverse Frisch elasticity men
        real(8), parameter :: etaf    = 1d0/0.8d0  ! Inverse Frisch elasticity women
        real(8), parameter :: mu   = 0.02d0       ! The growth rate of the economy
        !real(8), parameter :: beta   = 0.97d0    !Discount factor
        real(8), parameter :: beta   = 1.009040d0    !Discount factor
        integer, parameter :: Tret   = 36         ! Years in retirement  (65-100)
        integer, parameter :: T      = 45         ! Years of active life (20-64)
        integer, parameter :: nfc   = 23         ! Number of fixed costs
        real(8), parameter :: tc    = 0.05d0     ! Consumption tax
        real(8), parameter :: t_employer=0.0765d0, t_employee=0.0765d0 !Social security tax paid by employer and employee
        real(8), parameter :: tk     = 0.36d0     ! Capital tax 
        real(8), parameter :: sigma  = 4d0        ! Risk aversion parameter
        real(8), parameter :: alpha  = 1d0/3d0     ! Capital share
        real(8), parameter :: delta  = 0.0988d0-mu     ! Capital depreciation
        real(8), parameter :: deltaexp  = 0.000d0     !Depreciation of experience
        real(8) :: chim= 21.45d0,    chims=60.8d0           ! Disutility from work men
        real(8) :: chif= 3.935d0,   chifs= 8.82d0              ! Disutility from work women
        real(8), parameter :: sigma_um=0.31925155D0, rho_um=0.3962346D0 ! Parameters governing the process for the  idiosyncratic shock, men
        real(8), parameter :: sigma_am=0.33828463d0, rho_am=0D0     ! Stdev of ability, men
        real(8), parameter :: sigma_uf=0.31004311d0, rho_uf=0.339295 ! Parameters governing the process for the  idiosyncratic shock, women
        real(8), parameter :: sigma_af=0.38475527d0, rho_af=0d0    ! Stdev of ability, women
        real(8), parameter :: sigma_fcm=0.1680d0, rho_fcm=0d0 
        real(8), parameter :: sigma_fcs=0.2832d0, rho_fcs=0d0    ! Stdev of fixed costs, women
        real(8), parameter :: mu_fcm=-0.0536d0, mu_fcs=-0.0679d0   ! Mean fixed cost of LFP, women
        real(8) :: gamma(2,3), gamma0           ! Returns to experience parameters
        real(8) :: AE                           ! Average earnings
        real(8) :: match=0.884d0                 ! Assortative mating parameter
        real(8) :: theta(2), thetas(2)            ! Labor tax parameters
        real(8) :: tax_level_scale=1.0, tax_prog_scale=1.0d0  ! Parameters to scale tax- level and progressivity
        real(8), parameter :: ybar = 1.2d0                ! Upper limit for SS tax
        real(8) :: epsilon=1d0, epsilon2=1d0, epsilon3=1d0, epsilon4=1d0, epsilon5=1d0, epsilon6=1d0
        real(8), parameter :: cons_floor = 1d-9   ! Lower limit for consumption of unemployed
        integer, parameter :: na     = 5          ! Number of permanent ability levels
        integer, parameter :: nu     = 5          ! Number of transitory productivity shocks
        integer, parameter :: nexp   = 6          ! Number of gridpoints for experience
        integer, parameter :: nc     = 100        ! Number of gridpoints in "static" hours function
        integer, parameter :: nw     = 100        ! Number of gridpoints in "static" hours function
        integer, parameter :: nk     = 16         ! Number of gridpoints over capital
    
        integer, parameter :: nsim   = 10000         ! Number of households used for simulation
        integer, parameter :: nsim2   = 16         ! Number of simulations
        integer, parameter :: KORDER   = 3         ! Order of spline in K-dimension
        integer, parameter :: EXPORDER   = 3       ! Order of spline in experience
        
        real(8), dimension (:,:), allocatable :: a, Prob_a      ! Vectors with the value of ability and the probability of each ability
        real(8), dimension (:,:), allocatable :: Fc, Prob_fc    ! Vectors with the value of fixed costs and the probability of each fixed cost
        real(8), dimension (:,:), allocatable :: u, Prob_u      ! Vectors with the value of the idiosyncratic shock and the probability of each shock
        real(8), dimension (:,:,:), allocatable :: trans_u     ! Transition matrix for the shock u
        real(8), dimension (:,:,:), allocatable :: trans_a     ! Dummy transition matrix for a (used to compute the unconditional probability of each a)
        real(8), dimension (:,:,:), allocatable :: trans_fc     ! Dummy transition matrix for fc (used to compute the unconditional probability of each a)
        
        real(8) :: Gamma_redistr = 2d0*0.02510210114741359d0
        real(8) :: Psi_pension = 2d0*0.326242842747754d0
        real(8) :: Unemp_benefit, lumpsum=2d0*0.08462158531908406d0, lumpsumdum=2d0*0.09184923463949926d0
        
        real(8), dimension (:), allocatable :: OmegaRet, OmegaActive, Probm, Probd
        real(8), dimension (:), allocatable :: WeightRet, WeightActive
        real(8) :: ratio=4.29195095854055d0, ratiodum=4.28978442524046d0  !Ratio between capital and labor
        real(8) :: r                                    !Real interest rate
        real(8) :: debttoGDP=0.6185043538d0             ! Government debt as % of GDP
        real(8) :: milspendtoGDP=0.03625d0              ! Military spanding as % of GDP
        real(8) :: w                                    ! Wages per efficiency unit
        
        real(8), dimension (:,:,:,:), allocatable :: Sim1m,Sim1f,exp2m,exp2f ! Array to hold simulation results

        integer, dimension (:,:,:,:), allocatable :: exp1m,exp1f ! Array to hold simulation results
        
        real(8), dimension (:,:,:,:), allocatable :: SimR1m, SimR1f ! Array to hold simulation results for retired
        real(8), dimension (:,:,:), allocatable :: Random3m, Random3f, marstatm, marstatf, partshock
        real(8), dimension (:,:), allocatable :: Random1m, Random2m, Random1f, Random2f, marstatm_init, marstatf_init
        integer :: it,iter
           
    end module Model_Parameters