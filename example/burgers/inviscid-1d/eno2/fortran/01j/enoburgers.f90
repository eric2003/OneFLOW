module global
    implicit none
    integer, parameter :: nx = 40
    integer, parameter :: ighost = 10
    integer, parameter :: iorder = 2
    integer, parameter :: ishift = ighost + 1
    integer, parameter :: ist = 1 + ishift
    integer, parameter :: ied = nx + ishift
    integer, parameter :: ntcell = nx + ishift + ighost
    integer, parameter :: isize = iorder * ( iorder + 1 )
    real(8), parameter :: pi = 3.14159265358979323846
    integer :: il(0:nx), ir(0:nx)
    real(8) :: coef(0:iorder,0:iorder-1)
    real(8) :: dd(0:ighost-1, 1:ntcell)
    real(8) :: up1_2m(0:nx), up1_2p(0:nx), flux(0:nx)
    real(8) :: res(1:nx)
    real(8) :: dx, dt
end module global

module mesh_module
    use global, only: ntcell
    implicit none
    real(8) :: x(1:ntcell)
endmodule  mesh_module
    
module field_module
    use global, only: ntcell
    implicit none
    real(8) :: u(1:ntcell), un(1:ntcell)
endmodule  field_module    
    
subroutine residual(q)
    use global
    implicit none
    real(8) :: q(1:ntcell)
    integer :: i
    
    call reconstruction(q)
    call engquist_osher_flux(up1_2m,up1_2p,flux)
    do i = 1, nx
        res(i) = - ( flux(i) - flux(i-1) ) / dx
    enddo    
    
end subroutine residual
    
subroutine reconstruction(q)
    use global
    implicit none
    real(8) :: q(1:ntcell)
    integer :: i, j, m, k1, k2, l1, l2
    integer :: jj, kk, k
      
    !chose the stencil by ENO method
    do j = 1, ntcell
        dd(0,j) = q(j)
    enddo
    
    do k = 1, iorder - 1
        do j = 1, ntcell - 1
            dd(k,j) = dd(k-1,j+1)-dd(k-1,j)
        enddo
    enddo
    
    do i = 0, nx
        j = i + ist - 1
        il(i) = j
        ir(i) = j + 1
        do k=1,iorder-1  
            if (  abs(dd(k,il(i)-1)) <= abs(dd(k,il(i))) ) then
                il(i) = il(i) - 1 
            endif
            if ( abs(dd(k,ir(i)-1)) <= abs(dd(k,ir(i))) ) then
                ir(i) = ir(i) - 1 
            endif
        enddo
    enddo
    !   reconstruction u(j+1_2)
    do i = 0, nx
        j = i + ist - 1
        k1 = il(i)
        k2 = ir(i)
        l1 = j - k1 + 1
        l2 = j - k2 + 1
        up1_2m(i) = 0
        up1_2p(i) = 0
        do m=0,iorder-1 
            up1_2m(i) = up1_2m(i) + q(k1+m) * coef(l1,m)
            up1_2p(i) = up1_2p(i) + q(k2+m) * coef(l2,m)
        enddo  
    enddo        
end subroutine reconstruction
    
!calculate  numerical flux   
subroutine engquist_osher_flux(up1_2m,up1_2p,flux)
    use global, only: ist, ied, nx
    implicit none
    real(8) :: up1_2m(0:nx), up1_2p(0:nx), flux(0:nx)
    integer :: i
      
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
    
subroutine boundary( u )
    use global
    implicit none
    real(8) :: u(1:ntcell)
    integer :: i
         
    do i = - ighost, 0
        u( ist - 1 + i ) = u( ied + i )
    enddo
    
    do i = 1, ighost
        u( ied + i ) = u( ist - 1 + i )
    enddo
    
end subroutine boundary
    
subroutine update_oldfield(qn, q)
    use global, only: ntcell
    implicit none
    real(8),dimension(1:ntcell) :: qn, q
    qn = q
end subroutine update_oldfield
    
subroutine init_coef
    use global
    implicit none
    real(8) :: values(isize) = [1.5d0, -0.5d0, 0.5d0, 0.5d0, -0.5d0, 1.5d0]
    integer :: i, j, icount
    
    icount = 1
    do i = 0, iorder
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
    integer :: i, ii
    
    do i = -ighost, nx + ighost
        ii = i + ighost + 1
        x(ii) = ( i - 1 ) * dx + dx/2 - 1.0
        !write(*,*)'i=',i,'ii=',ii
    enddo
    
    i = 0
   
end subroutine init_mesh

subroutine init_field()
    use global
    use mesh_module
    use field_module
    implicit none
    integer :: i
    
    do i = ist, ied
        u(i) = 0.25 + 0.5 * sin( pi * x(i) )
    enddo
    
    call boundary( u )
    call update_oldfield(un, u)    

end subroutine init_field
    
subroutine runge_kutta_3()
    use global
    use field_module
    implicit none
    integer :: i, ii, j
    real(8) :: c1, c2, c3
    
    call residual(u)
    do i = 1, nx
        j = ist + i - 1
        ii = i - ist + 1
        u(j) = u(j) + dt * res(i)
    enddo
    call boundary( u )
    
    call residual(u)
    
    do i = 1, nx
        j = ist + i - 1
        u(j) = 0.75 * un(j) + 0.25 * u(j) + 0.25 * dt * res(i)
    enddo
    
    call boundary( u )
    
    call residual( u )
    
    c1 = 1.0 / 3.0
    c2 = 2.0 / 3.0
    c3 = 2.0 / 3.0
    
    do i = 1, nx
        j = ist + i - 1
        u(j) = c1 * un(j) + c2 * u(j) + c3 * dt * res(i)
    enddo
    call boundary( u )
    
    call update_oldfield(un, u)

end subroutine runge_kutta_3
    
subroutine visualize
    use global
    use mesh_module
    use field_module
    implicit none
    integer :: i
    open(1,file='solution_total.plt',status='unknown')
    do i = 1, ntcell
        write(1,101) x(i),u(i)
    enddo
    close(1)
    
    open(2,file='solution.plt',status='unknown')
    do i = ist, ied
        write(2,101) x(i),u(i)
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