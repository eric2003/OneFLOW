! ==================== string_utils.f90 ====================
module string_utils
    implicit none
    private
    
    public :: to_lower, to_upper
    
contains

    ! ===================================================================
    ! 转小写
    ! ===================================================================
    function to_lower(str) result(lower_str)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: lower_str
        integer :: i, code
        
        lower_str = str
        do i = 1, len(lower_str)
            code = iachar(lower_str(i:i))
            ! 大写 A-Z (65-90) 转小写 a-z (97-122)
            if (code >= 65 .and. code <= 90) then
                lower_str(i:i) = achar(code + 32)
            end if
        end do
    end function to_lower

    ! ===================================================================
    ! 转大写
    ! ===================================================================
    function to_upper(str) result(upper_str)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: upper_str
        integer :: i, code
        
        upper_str = str
        do i = 1, len(upper_str)
            code = iachar(upper_str(i:i))
            ! 小写 a-z (97-122) 转大写 A-Z (65-90)
            if (code >= 97 .and. code <= 122) then
                upper_str(i:i) = achar(code - 32)
            end if
        end do
    end function to_upper

end module string_utils