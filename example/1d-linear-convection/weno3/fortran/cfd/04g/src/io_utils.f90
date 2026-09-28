! ==================== io_utils.f90 ====================
module io_utils
	use kinds, only: dp
    implicit none
    private
    
    public :: write_simulation_results
    
contains

    !-----------------------------------------------------------------------
    ! 写入模拟结果到文件
    !-----------------------------------------------------------------------
    subroutine write_simulation_results(filename, xdata, udata, scheme_name)
        character(len=*), intent(in) :: filename      ! 文件名（任意长度）
        real(dp),         intent(in) :: xdata(:)      ! x 坐标数组
        real(dp),         intent(in) :: udata(:)      ! 解数组
        character(len=*), intent(in) :: scheme_name   ! 格式名称
        
        integer :: iunit, i, ncells
        
        ncells = size(xdata)
        
        ! 自动分配文件单元号
        open(newunit=iunit, file=filename, status='replace', action='write')
        
        ! 写入头信息（包含格式名）
        write(iunit, '(3A)') '# x, u (', trim(scheme_name), ')'
        
        ! 写入数据
        do i = 1, ncells
            write(iunit, '(2F12.6)') xdata(i), udata(i)
        end do
        
        close(iunit)
    end subroutine write_simulation_results
    
end module io_utils