! ==================== complex_wave_precise.f90 ====================
! Precise implementation with exact parameters from literature
!
! Parameters:
!   a = 0.5       (center of region 4)
!   z = -0.7      (center of region 1)  
!   δ = 0.005     (offset for both Gaussian and ellipse)
!   α = 10        (parameter for ellipse)
!   β = log(2) / (36δ²)  (parameter for Gaussian)
!
! Formula:
!   u_0(x) = (1/6)[G(x,β,z-δ) + G(x,β,z+δ) + 4G(x,β,z)],  -0.8 ≤ x ≤ -0.6
!            1,                                            -0.4 ≤ x ≤ -0.2
!            1 - |10(x-0.1)|,                               0 ≤ x ≤ 0.2
!            (1/6)[F(x,α,a-δ) + F(x,α,a+δ) + 4F(x,α,a)],   0.4 ≤ x ≤ 0.6
!            0,                                            otherwise
! =====================================================================

module complex_wave_precise_module
    use kinds, only: dp
    implicit none
    
    private
    public :: complex_wave_precise
    public :: init_precise_params, get_precise_params, calculate_beta
    
    ! 精确参数结构
    type , public ::  PreciseParams
        real(dp) :: a = 0.5_dp          ! center of region 4
        real(dp) :: z = -0.7_dp         ! center of region 1
        real(dp) :: delta = 0.005_dp    ! offset
        real(dp) :: alpha = 10.0_dp     ! parameter for ellipse
        real(dp) :: beta                ! parameter for Gaussian, will be calculated
    end type PreciseParams
    
    type(PreciseParams), save :: params
    
    ! 初始化参数
    data params%beta / 0.0_dp /  ! 将在初始化时计算
    
contains
    
    ! ===================================================================
    ! 初始化函数（计算β）
    ! ===================================================================
    
    subroutine init_precise_params()
        ! 计算 β = log(2) / (36δ²)
        real(dp), parameter :: LOG2 = log(2.0_dp)
        
        params%beta = LOG2 / (36.0_dp * params%delta**2)
        
        print *, "=== Precise Complex Wave Parameters ==="
        print *, "  a (ellipse center)     = ", params%a
        print *, "  z (Gaussian center)    = ", params%z
        print *, "  δ (offset)            = ", params%delta
        print *, "  α (ellipse parameter) = ", params%alpha
        print *, "  β (Gaussian parameter)= ", params%beta, " (log2/(36δ²))"
        print *, ""
    end subroutine init_precise_params
    
    ! 获取参数值
    function get_precise_params() result(p)
        type(PreciseParams) :: p
        p = params
    end function
    
    ! 单独计算β的函数
    function calculate_beta(delta) result(beta)
        real(dp), intent(in) :: delta
        real(dp) :: beta
        real(dp), parameter :: LOG2 = log(2.0_dp)
        beta = LOG2 / (36.0_dp * delta**2)
    end function calculate_beta
    
    ! ===================================================================
    ! 精确的数学函数
    ! ===================================================================
    
    pure function gaussian_precise(x) result(g)
        real(dp), intent(in) :: x
        real(dp) :: g
        ! G(x,β,z) = exp(-β(x-z)²)
        g = exp(-params%beta * (x - params%z)**2)
    end function gaussian_precise
    
    pure function half_ellipse_precise(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f, arg
        ! F(x,α,a) = √max(1 - α²(x-a)², 0)
        arg = 1.0_dp - (params%alpha * (x - params%a))**2
        if (arg > 0.0_dp) then
            f = sqrt(arg)
        else
            f = 0.0_dp
        end if
    end function half_ellipse_precise
    
    ! ===================================================================
    ! 主函数：使用精确参数
    ! ===================================================================
    
    function complex_wave_precise(x) result(u0)
        real(dp), intent(in) :: x
        real(dp) :: u0
        
        real(dp) :: g1, g2, g3, f1, f2, f3
        real(dp) :: arg1, arg2, arg3  ! 添加变量声明
        
        ! 区域1: 组合高斯脉冲
        if (-0.8_dp <= x .and. x <= -0.6_dp) then
            ! G(x,β,z-δ)
            g1 = exp(-params%beta * (x - (params%z - params%delta))**2)
            ! G(x,β,z+δ)
            g2 = exp(-params%beta * (x - (params%z + params%delta))**2)
            ! G(x,β,z)
            g3 = exp(-params%beta * (x - params%z)**2)
            ! (1/6)[G1 + G2 + 4G3]
            u0 = (g1 + g2 + 4.0_dp * g3) / 6.0_dp
        
        ! 区域2: 方波
        else if (-0.4_dp <= x .and. x <= -0.2_dp) then
            u0 = 1.0_dp
        
        ! 区域3: 三角波
        else if (0.0_dp <= x .and. x <= 0.2_dp) then
            u0 = 1.0_dp - abs(10.0_dp * (x - 0.1_dp))
            if (u0 < 0.0_dp) u0 = 0.0_dp  ! 处理浮点误差
        
        ! 区域4: 组合半椭圆
        else if (0.4_dp <= x .and. x <= 0.6_dp) then
            ! F(x,α,a-δ)
            arg1 = 1.0_dp - (params%alpha * (x - (params%a - params%delta)))**2
            if (arg1 > 0.0_dp) then
                f1 = sqrt(arg1)
            else
                f1 = 0.0_dp
            end if
            
            ! F(x,α,a+δ)
            arg2 = 1.0_dp - (params%alpha * (x - (params%a + params%delta)))**2
            if (arg2 > 0.0_dp) then
                f2 = sqrt(arg2)
            else
                f2 = 0.0_dp
            end if
            
            ! F(x,α,a)
            arg3 = 1.0_dp - (params%alpha * (x - params%a))**2
            if (arg3 > 0.0_dp) then
                f3 = sqrt(arg3)
            else
                f3 = 0.0_dp
            end if
            
            ! (1/6)[F1 + F2 + 4F3]
            u0 = (f1 + f2 + 4.0_dp * f3) / 6.0_dp
        
        ! 其他区域
        else
            u0 = 0.0_dp
        end if
        
    end function
    
end module complex_wave_precise_module