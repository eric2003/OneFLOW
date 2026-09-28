module global
    implicit none
    integer, parameter :: nx = 40
    integer, parameter :: ighost = 10
    integer, parameter :: iorder = 2
    integer, parameter :: isize = iorder * ( iorder + 1 )
    real(8), parameter :: pi = 3.14159265358979323846
    integer il(0:nx),ir(0:nx)
    real(8) :: coef(-1:iorder-1,0:iorder-1)
    real(8) :: dd(0:ighost-1, -ighost:nx+ighost)
    real(8) :: up1_2m(0:nx), up1_2p(0:nx), flux(0:nx)
    real(8) :: res(1:nx)
    real(8) :: dx, dt
end module global

module mesh_module
    use global, only: nx, ighost
    implicit none
    real(8) :: x(-ighost:nx+ighost)
endmodule  mesh_module
    
module field_module
    use global, only: nx, ighost
    implicit none
    real(8) :: pu(-ighost:nx+ighost), un(-ighost:nx+ighost)
endmodule  field_module    
    
subroutine residual(u)
    use global
    implicit none
    real(8) :: u(-ighost:nx+ighost)
    integer i
    
    call reconstruction(u)
    call engquist_osher_flux(up1_2m,up1_2p,flux,nx)
    do i = 1, nx   
        res(i) = - ( flux(i) - flux(i-1) ) / dx
    enddo    
    
end subroutine residual
    
subroutine reconstruction(u)
    use global
    implicit none
    real(8) :: u(-ighost:nx+ighost)
    integer :: i, j, m, k1, k2, l1, l2
      
    !chose the stencil by ENO method
    do j = -ighost, nx + ighost
        dd(0,j)=u(j)
    enddo
    do i=1,iorder-1
        do j=-ighost,nx+ighost-1
            dd(i,j)=dd(i-1,j+1)-dd(i-1,j)
        enddo
    enddo
    
    do j = 0, nx 
        il(j) = j
        ir(j) = j + 1
        do i=1,iorder-1  
            if( abs(dd(i,il(j)-1)) <= abs(dd(i,il(j))) ) then
                il(j)=il(j)-1 
            endif
            if( abs(dd(i,ir(j)-1)) <= abs(dd(i,ir(j))) ) then
                ir(j)=ir(j)-1 
            endif
        enddo
    enddo
    !   reconstruction u(j+1_2)
    do j=0,nx  
        k1=il(j)       
        k2=ir(j)
        l1=j-k1
        l2=j-k2
        up1_2m(j)=0
        up1_2p(j)=0
        do m=0,iorder-1 
            up1_2m(j)=up1_2m(j)+u(k1+m)*coef(l1,m)
            up1_2p(j)=up1_2p(j)+u(k2+m)*coef(l2,m)
        enddo  
    enddo        
end subroutine reconstruction
    
!calculate  numerical flux   
subroutine engquist_osher_flux(up1_2m,up1_2p,flux,nx)
    implicit none
    real(8) :: up1_2m(0:nx),up1_2p(0:nx),flux(0:nx)
    integer :: i, nx
      
    do i = 0, nx
        if ( up1_2m(i) >= 0 ) then
            if ( up1_2p(i) >= 0 )  then
                flux(i) = 0.5 * up1_2m(i) * up1_2m(i)
            else 
                flux(i) = 0.5 * ( up1_2m(i) * up1_2m(i) + up1_2p(i) * up1_2p(i) )
            endif    
        else
            if ( up1_2p(i) >= 0 )  then
                flux(i) = 0
            else
                flux(i) = 0.5 * up1_2p(i) * up1_2p(i)
            endif
        endif    
    enddo            
end subroutine engquist_osher_flux
    
subroutine boundary( u, nx, ighost )
    implicit none
    integer :: nx, ighost
    real(8) :: u(-ighost:nx+ighost)
    integer :: i
         
    do i = 0, - ighost, - 1
        u( i ) = u( i + nx )
    enddo
         
    do i = nx + 1, nx + ighost
        u( i ) = u( i - nx )
    enddo
end subroutine boundary
    
subroutine update_oldfield(un, pu, nx, ighost)
    implicit none
    integer :: nx, ighost
    real(8) :: un(-ighost:nx+ighost), pu(-ighost:nx+ighost)
    integer :: i
    
    do i = -ighost, nx + ighost
        un( i ) = pu( i )
    enddo
end subroutine update_oldfield
    
subroutine init_coef
    use global
    implicit none
    real(8) :: values(isize) = [1.5d0, -0.5d0, 0.5d0, 0.5d0, -0.5d0, 1.5d0]
    integer :: i, j, icount
    
    icount = 1
    do i = -1, iorder-1
        do j = 0, iorder-1
            coef(i, j) = values(icount)
            icount = icount + 1
        end do
    end do    
    
end subroutine init_coef
    
subroutine init_mesh
    use global
    use mesh_module
    implicit none
    integer :: i
    
    do i = -ighost, nx + ighost
        x(i) = ( i - 1 ) * dx + dx/2 - 1.0
    enddo
   
end subroutine init_mesh

subroutine init_field()
    use global
    use mesh_module
    use field_module
    implicit none
    integer :: i
    
    do i = 1, nx
        pu(i) = 0.25 + 0.5 * sin( pi * x(i) )
    enddo
    
    call boundary( pu, nx, ighost )
    call update_oldfield(un, pu, nx, ighost)    

end subroutine init_field
    
subroutine runge_kutta_3()
    use global
    use field_module
    implicit none
    integer :: i
    real(8) :: c1, c2, c3
    
    call residual(pu)
    do i = 1, nx
        pu(i) = pu(i) + dt * res(i)
    enddo
    call boundary( pu, nx, ighost )
    
    call residual(pu)
    
    do i = 1, nx
        pu(i) = 0.75 * un(i) + 0.25 * pu(i) + 0.25 * dt * res(i)
    enddo
    
    call boundary( pu, nx, ighost )
    
    call residual(pu)
    
    c1 = 1.0 / 3.0
    c2 = 2.0 / 3.0
    c3 = 2.0 / 3.0
    
    do i = 1, nx                               
        pu(i) = c1 * un(i) + c2 * pu(i) + c3 * dt * res(i)
    enddo
    call boundary( pu, nx, ighost )
    
    call update_oldfield(un, pu, nx, ighost)

end subroutine runge_kutta_3
    
subroutine visualize
    use global
    use mesh_module
    use field_module
    implicit none
    integer :: i
    open(1,file='solution_total.plt',status='unknown')
    do i = -ighost, nx + ighost
        write(1,101) x(i),pu(i)
    enddo
    close(1)
    
    open(2,file='solution.plt',status='unknown')
    do i = 1, nx
        write(2,101) x(i),pu(i)
    enddo
    close(2)
    101 format(1x,e20.10,e20.10)
end subroutine visualize
    
program main  
    use global
    use mesh_module
    use field_module
    implicit none
    integer :: i
    real(8) :: t, simu_time
      
    dx = 2.0 / nx
    dt = dx * 0.5
    write(*,*) 'Input T:'
    read(*,*) simu_time
    
    call init_coef()
    call init_mesh()
    call init_field()
 
    t = 0
    do while( t < simu_time ) 
        call runge_kutta_3()
       
        t = t + dt
        if ( t + dt > simu_time ) then
            dt = simu_time - t
        endif
    enddo
       
    write(*,*) t
    
    call visualize()

end program main