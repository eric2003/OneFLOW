! src/base/precision.f90（简单版本）
module precision_module
    use base_modules, only: wp, ip
    implicit none
    
    ! 重新导出，确保兼容
    public :: wp, ip
    
end module precision_module