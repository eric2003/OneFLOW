! ==================== recon_coefficients.f90 ====================
! Reconstruction coefficients for OneFLOW-CFD
!
! Purpose: Initialize coefficients for ENO/WENO reconstruction schemes
! Author: OneFLOW-CFD Team
! Date: Updated to 7th order
! ===============================================================

module recon_coefficients_module
    use kinds, only: dp
    implicit none
    
    private
    public :: init_coef
    
    ! 定义系数数组大小的常量（最大阶数+1）
    integer, parameter :: MAX_ORDER = 7
    
contains
    
    ! Initialize ENO/WENO coefficients up to 7th order
    subroutine init_coef(spatial_order, coef)
        integer, intent(in) :: spatial_order
        real(dp), intent(out) :: coef(:,:)
        
        integer :: i, j
        
        ! 初始化所有系数为0
        coef = 0.0_dp
        
        ! 检查数组大小是否足够
        if (size(coef, 1) < spatial_order + 1 .or. size(coef, 2) < spatial_order) then
            error stop "Coefficient array too small for requested spatial order"
        end if
        
        select case(spatial_order)
        case(1)
            ! 1st order ENO
            coef(1,1) = 1.0_dp
            coef(2,1) = 1.0_dp
            
        case(2)
            ! 2nd order ENO
            coef(1,1:2) = [ 3.0_dp/2.0_dp, -1.0_dp/2.0_dp ]
            coef(2,1:2) = [ 1.0_dp/2.0_dp,  1.0_dp/2.0_dp ]
            coef(3,1:2) = [ -1.0_dp/2.0_dp, 3.0_dp/2.0_dp ]
            
        case(3)
            ! 3rd order ENO
            coef(1,1:3) = [ 11.0_dp/6.0_dp, -7.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
            coef(2,1:3) = [  1.0_dp/3.0_dp,  5.0_dp/6.0_dp, -1.0_dp/6.0_dp ]
            coef(3,1:3) = [ -1.0_dp/6.0_dp,  5.0_dp/6.0_dp,  1.0_dp/3.0_dp ]
            coef(4,1:3) = [  1.0_dp/3.0_dp, -7.0_dp/6.0_dp, 11.0_dp/6.0_dp ]
            
        case(4)
            ! 4th order ENO
            coef(1,1:4) = [ 25.0_dp/12.0_dp, -23.0_dp/12.0_dp, &
                            13.0_dp/12.0_dp, -1.0_dp/4.0_dp ]
            coef(2,1:4) = [ 1.0_dp/4.0_dp, 13.0_dp/12.0_dp, &
                           -5.0_dp/12.0_dp, 1.0_dp/12.0_dp ]
            coef(3,1:4) = [ -1.0_dp/12.0_dp, 7.0_dp/12.0_dp, &
                            7.0_dp/12.0_dp, -1.0_dp/12.0_dp ]
            coef(4,1:4) = [ 1.0_dp/12.0_dp, -5.0_dp/12.0_dp, &
                            13.0_dp/12.0_dp, 1.0_dp/4.0_dp ]
            coef(5,1:4) = [ -1.0_dp/4.0_dp, 13.0_dp/12.0_dp, &
                           -23.0_dp/12.0_dp, 25.0_dp/12.0_dp ]
            
        case(5)
            ! 5th order ENO
            coef(1,1:5) = [ 137.0_dp/60.0_dp, -163.0_dp/60.0_dp, &
                            137.0_dp/60.0_dp, -21.0_dp/20.0_dp, &
                            1.0_dp/5.0_dp ]
            coef(2,1:5) = [ 1.0_dp/5.0_dp, 77.0_dp/60.0_dp, &
                           -43.0_dp/60.0_dp, 17.0_dp/60.0_dp, &
                           -1.0_dp/20.0_dp ]
            coef(3,1:5) = [ -1.0_dp/20.0_dp, 9.0_dp/20.0_dp, &
                            47.0_dp/60.0_dp, -13.0_dp/60.0_dp, &
                            1.0_dp/30.0_dp ]
            coef(4,1:5) = [ 1.0_dp/30.0_dp, -13.0_dp/60.0_dp, &
                            47.0_dp/60.0_dp, 9.0_dp/20.0_dp, &
                           -1.0_dp/20.0_dp ]
            coef(5,1:5) = [ -1.0_dp/20.0_dp, 17.0_dp/60.0_dp, &
                           -43.0_dp/60.0_dp, 77.0_dp/60.0_dp, &
                            1.0_dp/5.0_dp ]
            coef(6,1:5) = [ 1.0_dp/5.0_dp, -21.0_dp/20.0_dp, &
                            137.0_dp/60.0_dp, -163.0_dp/60.0_dp, &
                            137.0_dp/60.0_dp ]
            
        case(6)
            ! 6th order ENO
            coef(1,1:6) = [ 49.0_dp/20.0_dp, -71.0_dp/20.0_dp, &
                            79.0_dp/20.0_dp, -163.0_dp/60.0_dp, &
                            31.0_dp/30.0_dp, -1.0_dp/6.0_dp ]
            coef(2,1:6) = [ 1.0_dp/6.0_dp, 29.0_dp/20.0_dp, &
                           -21.0_dp/20.0_dp, 37.0_dp/60.0_dp, &
                           -13.0_dp/60.0_dp, 1.0_dp/30.0_dp ]
            coef(3,1:6) = [ -1.0_dp/30.0_dp, 11.0_dp/30.0_dp, &
                            19.0_dp/20.0_dp, -23.0_dp/60.0_dp, &
                            7.0_dp/60.0_dp, -1.0_dp/60.0_dp ]
            coef(4,1:6) = [ 1.0_dp/60.0_dp, -2.0_dp/15.0_dp, &
                            37.0_dp/60.0_dp, 37.0_dp/60.0_dp, &
                           -2.0_dp/15.0_dp, 1.0_dp/60.0_dp ]
            coef(5,1:6) = [ -1.0_dp/60.0_dp, 7.0_dp/60.0_dp, &
                           -23.0_dp/60.0_dp, 19.0_dp/20.0_dp, &
                            11.0_dp/30.0_dp, -1.0_dp/30.0_dp ]
            coef(6,1:6) = [ 1.0_dp/30.0_dp, -13.0_dp/60.0_dp, &
                            37.0_dp/60.0_dp, -21.0_dp/20.0_dp, &
                            29.0_dp/20.0_dp, 1.0_dp/6.0_dp ]
            coef(7,1:6) = [ -1.0_dp/6.0_dp, 31.0_dp/30.0_dp, &
                           -163.0_dp/60.0_dp, 79.0_dp/20.0_dp, &
                           -71.0_dp/20.0_dp, 49.0_dp/20.0_dp ]
            
        case(7)
            ! 7th order ENO
            coef(1,1:7) = [ 363.0_dp/140.0_dp, -617.0_dp/140.0_dp, &
                            853.0_dp/140.0_dp, -2341.0_dp/420.0_dp, &
                            667.0_dp/210.0_dp, -43.0_dp/42.0_dp, &
                            1.0_dp/7.0_dp ]
            coef(2,1:7) = [ 1.0_dp/7.0_dp, 223.0_dp/140.0_dp, &
                           -197.0_dp/140.0_dp, 153.0_dp/140.0_dp, &
                           -241.0_dp/420.0_dp, 37.0_dp/210.0_dp, &
                           -1.0_dp/42.0_dp ]
            coef(3,1:7) = [ -1.0_dp/42.0_dp, 13.0_dp/42.0_dp, &
                            153.0_dp/140.0_dp, -241.0_dp/420.0_dp, &
                            109.0_dp/420.0_dp, -31.0_dp/420.0_dp, &
                            1.0_dp/105.0_dp ]
            coef(4,1:7) = [ 1.0_dp/105.0_dp, -19.0_dp/210.0_dp, &
                            107.0_dp/210.0_dp, 319.0_dp/420.0_dp, &
                           -101.0_dp/420.0_dp, 5.0_dp/84.0_dp, &
                           -1.0_dp/140.0_dp ]
            coef(5,1:7) = [ -1.0_dp/140.0_dp, 5.0_dp/84.0_dp, &
                           -101.0_dp/420.0_dp, 319.0_dp/420.0_dp, &
                            107.0_dp/210.0_dp, -19.0_dp/210.0_dp, &
                            1.0_dp/105.0_dp ]
            coef(6,1:7) = [ 1.0_dp/105.0_dp, -31.0_dp/420.0_dp, &
                            109.0_dp/420.0_dp, -241.0_dp/420.0_dp, &
                            153.0_dp/140.0_dp, 13.0_dp/42.0_dp, &
                           -1.0_dp/42.0_dp ]
            coef(7,1:7) = [ -1.0_dp/42.0_dp, 37.0_dp/210.0_dp, &
                           -241.0_dp/420.0_dp, 153.0_dp/140.0_dp, &
                           -197.0_dp/140.0_dp, 223.0_dp/140.0_dp, &
                            1.0_dp/7.0_dp ]
            coef(8,1:7) = [ 1.0_dp/7.0_dp, -43.0_dp/42.0_dp, &
                            667.0_dp/210.0_dp, -2341.0_dp/420.0_dp, &
                            853.0_dp/140.0_dp, -617.0_dp/140.0_dp, &
                            363.0_dp/140.0_dp ]
            
        case default
            print *, "Warning: Unsupported spatial order = ", spatial_order
            print *, "Supported orders: 1-7"
            error stop "Unsupported spatial order in init_coef"
        end select
        
    end subroutine init_coef
    
    ! 可选：验证系数的辅助函数
    subroutine verify_coef(spatial_order, coef)
        integer, intent(in) :: spatial_order
        real(dp), intent(in) :: coef(:,:)
        
        integer :: i, j
        real(dp) :: sum_row, expected_sum
        
        print *, "Verifying ENO coefficients for order", spatial_order
        
        ! 检查每行的和应该为1（重构的一致性）
        expected_sum = 1.0_dp
        do i = 1, spatial_order + 1
            sum_row = sum(coef(i, 1:spatial_order))
            if (abs(sum_row - expected_sum) > 1.0e-10_dp) then
                print *, "Warning: Row", i, "sum =", sum_row, "expected", expected_sum
            end if
        end do
        
        ! 打印系数以供验证
        print *, "ENO coefficients for order", spatial_order, ":"
        do i = 1, spatial_order + 1
            write(*, '(A,I2,A,*(F10.6,:,", "))') "  coef(", i, ",:) = ", (coef(i,j), j=1,spatial_order)
        end do
        
    end subroutine verify_coef
    
end module recon_coefficients_module