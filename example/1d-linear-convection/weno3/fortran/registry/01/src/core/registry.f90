! src/core/registry.f90
! Fortran注册系统模块
module registry_module
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    private  ! 默认所有内容都是私有的
    
    ! ==================== 1. 公开接口 ====================
    ! 只在这里声明一次哪些是公开的
    public :: wp                      ! 精度类型
    public :: component_info          ! 组件信息类型
    public :: component_registry_type ! 注册表类型
    public :: component_registry      ! 全局注册表实例
    public :: register_component      ! 注册函数
    
    ! ==================== 2. 类型定义 ====================
    
    ! 组件信息类型
    type :: component_info
        character(len=50) :: category = ""
        character(len=50) :: name = ""
    contains
        procedure :: print => ci_print
    end type component_info
    
    ! 注册表类型
    type :: component_registry_type
        private
        type(component_info), allocatable :: components(:)
        integer :: count = 0
        integer :: capacity = 10
        logical :: verbose = .true.
    contains
        procedure :: register => cr_register
        procedure :: get => cr_get
        procedure :: list_all => cr_list_all
        procedure :: size => cr_size
        procedure :: clear => cr_clear
    end type component_registry_type
    
    ! ==================== 3. 全局实例 ====================
    ! 这里只声明，不指定属性
    type(component_registry_type) :: component_registry
    
    ! ==================== 4. 接口定义 ====================
    
    interface register_component
        module procedure register_component_simple
    end interface register_component
    
contains
    
    ! ==================== 5. 组件信息方法 ====================
    
    subroutine ci_print(this)
        class(component_info), intent(in) :: this
        print *, "  [", trim(this%category), ".", trim(this%name), "]"
    end subroutine ci_print
    
    ! ==================== 6. 注册表核心方法 ====================
    
    ! 注册组件（简单版）
    subroutine register_component_simple(category, name)
        character(len=*), intent(in) :: category, name
        call component_registry%register(category, name)
    end subroutine register_component_simple
    
    ! 内部注册实现
    subroutine cr_register(this, category, name)
        class(component_registry_type), intent(inout) :: this
        character(len=*), intent(in) :: category, name
        
        type(component_info) :: new_component
        type(component_info), allocatable :: temp(:)
        integer :: i
        
        ! 检查是否已存在
        do i = 1, this%count
            if (trim(this%components(i)%category) == trim(category) .and. &
                trim(this%components(i)%name) == trim(name)) then
                if (this%verbose) then
                    print *, "⚠️  警告: 覆盖注册 ", trim(category), ".", trim(name)
                end if
                ! 更新现有组件
                this%components(i)%category = trim(category)
                this%components(i)%name = trim(name)
                return
            end if
        end do
        
        ! 创建新组件
        new_component%category = trim(category)
        new_component%name = trim(name)
        
        ! 初始化数组（如果需要）
        if (.not. allocated(this%components)) then
            allocate(this%components(this%capacity))
        end if
        
        ! 扩展数组（如果需要）
        if (this%count >= this%capacity) then
            this%capacity = this%capacity * 2
            allocate(temp(this%capacity))
            temp(1:this%count) = this%components(1:this%count)
            call move_alloc(temp, this%components)
        end if
        
        ! 添加组件
        this%count = this%count + 1
        this%components(this%count) = new_component
        
        if (this%verbose) then
            print *, "✅ 已注册: ", trim(category), ".", trim(name)
        end if
    end subroutine cr_register
    
    ! 获取组件
    function cr_get(this, category, name) result(info)
        class(component_registry_type), intent(in) :: this
        character(len=*), intent(in) :: category, name
        type(component_info) :: info
        
        integer :: i
        
        ! 初始化为空
        info%category = ""
        info%name = ""
        
        ! 搜索组件
        do i = 1, this%count
            if (trim(this%components(i)%category) == trim(category) .and. &
                trim(this%components(i)%name) == trim(name)) then
                info = this%components(i)
                return
            end if
        end do
        
        if (this%verbose) then
            print *, "❌ 未找到: ", trim(category), ".", trim(name)
        end if
    end function cr_get
    
    ! 列出所有组件
    subroutine cr_list_all(this)
        class(component_registry_type), intent(in) :: this
        integer :: i
        
        print *, "=== 注册表内容 (" // trim(str(this%count)) // " 个组件) ==="
        if (this%count == 0) then
            print *, "  (空)"
            return
        end if
        
        do i = 1, this%count
            call this%components(i)%print()
        end do
    end subroutine cr_list_all
    
    ! 获取组件数量
    integer function cr_size(this)
        class(component_registry_type), intent(in) :: this
        cr_size = this%count
    end function cr_size
    
    ! 清空注册表
    subroutine cr_clear(this)
        class(component_registry_type), intent(inout) :: this
        if (allocated(this%components)) then
            deallocate(this%components)
        end if
        this%count = 0
        this%capacity = 10
        print *, "🗑️  注册表已清空"
    end subroutine cr_clear
    
    ! ==================== 7. 工具函数 ====================
    
    ! 整数转字符串
    function str(i) result(s)
        integer, intent(in) :: i
        character(len=20) :: s
        write(s, '(I0)') i
        s = trim(adjustl(s))
    end function str
    
end module registry_module