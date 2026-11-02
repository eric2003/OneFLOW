subroutine reconstruction(u,nx,up1_2m,up1_2p,dd,il,ir,coef,iorder,ighost)
    implicit none
    real(8) :: u(-ighost:nx+ighost)
    real(8) :: up1_2m(0:nx),up1_2p(0:nx)            
    real(8) :: dd(0:ighost-1,-ighost:nx+ighost)
    real(8) :: coef(-1:iorder-1,0:iorder-1)
    integer :: il(0:nx), ir(0:nx)
    integer :: i, j, k1, k2, l1, l2, m, nx
    integer :: iorder, ighost
      
    !chose the stencil by ENO method
    do j=-ighost,nx+ighost
        dd(0,j)=u(j)
    enddo
    do i=1,iorder-1
        do j=-ighost,nx+ighost-1
            dd(i,j)=dd(i-1,j+1)-dd(i-1,j)
        enddo
    enddo 
    
    do j=0,nx 
        il(j)=j
        ir(j)=j+1
        do i=1,iorder-1  
            if( abs(dd(i,il(j)-1)) <= abs(dd(i,il(j))) ) then
                il(j)=il(j)-1 
            endif
            if( abs(dd(i,ir(j)-1)) <= abs(dd(i,ir(j))) ) then
                ir(j)=ir(j)-1 
            endif
        enddo
    enddo
    
    !reconstruction u(j+1_2)
    do j = 0, nx
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
subroutine getflux(up1_2m,up1_2p,flux,nx)
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
end subroutine getflux
    
subroutine rhs(up1_2m,up1_2p,nx,dx,res)
    implicit none
    real(8) :: up1_2m(0:nx), up1_2p(0:nx), flux(0:nx)
    real(8) :: res(1:nx)
    real(8) :: dx
    integer :: nx, i
    
    call getflux(up1_2m,up1_2p,flux,nx)
    
    do i = 1, nx
        res(i) = - ( flux(i) - flux(i-1) ) / dx
    enddo    
end subroutine rhs
    
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
end subroutine
      
program main  
    implicit none
    integer, parameter :: nx = 40
    integer, parameter :: ighost = 10
    integer, parameter :: iorder = 2
    integer, parameter :: isize = iorder * ( iorder + 1 )
    real(8), parameter :: pi = 3.14159265358979323846
    real(8) :: pu(-ighost:nx+ighost), su(-ighost:nx+ighost)
    real(8) :: u1(-ighost:nx+ighost), u2(-ighost:nx+ighost)
    real(8) :: u0(1:nx), x(-ighost:nx+ighost)
    real(8) :: up1_2m(0:nx), up1_2p(0:nx), flux(0:nx)
    real(8) :: res(1:nx)
    real(8) :: dd(0:ighost-1, -ighost:nx+ighost)
    real(8) :: coef(-1:iorder-1,0:iorder-1) 
    integer :: il(0:nx), ir(0:nx)
    integer :: i, j, icount
    real(8) :: supt, t, temp, t1, t2, error, it
    real(8) :: dx, dt
    real(8) :: values(isize) = [1.5d0, -0.5d0, 0.5d0, 0.5d0, -0.5d0, 1.5d0]
            
    dx = 2.0 / nx
    dt = dx * 0.5
    
    write(*,*) 'dx = ', dx, 'dt = ', dt
    write(*,*) 'Input T:'
    read(*,*) supt 
      
!   2nd-order coefficients       
    icount = 1
    do i = -1, iorder-1
        do j = 0, iorder-1
            coef(i, j) = values(icount)
            icount = icount + 1
        end do
    end do    
!   Initialize grid and initial conditions
    do i= - ighost, nx + ighost
        x(i) = (i-1) * dx + dx/2 - 1.0
    enddo
      
!   Initial mean value 1      
    do i = 1, nx 
        pu(i) = 0.25 + 0.5 * sin( pi * x(i) )
    enddo
    
    call boundary( pu, nx, ighost )
    
    do i = 0, nx
       u0(i) = pu(i)
    enddo                           

!   Time stepping
    t = 0                      
    do while ( t < supt ) 
        call reconstruction(pu,nx,up1_2m,up1_2p,dd,il,ir,coef,iorder,ighost)
        call rhs(up1_2m,up1_2p,nx,dx,res)
        do i = 1, nx
            !ress = - ( flux(i) - flux(i-1) ) / dx
            u1(i) = pu(i) + dt * res(i)
        enddo
        
        call boundary( u1, nx, ighost )
        
        call reconstruction(u1,nx,up1_2m,up1_2p,dd,il,ir,coef,iorder,ighost)
        call rhs(up1_2m,up1_2p,nx,dx,res)

        do i = 1, nx
            u2(i) = 3.0/4.0 * pu(i) + 1.0/4.0 * u1(i) + 1.0/4.0 * dt * res(i)
        enddo
        
        call boundary( u2, nx, ighost )
        
        call reconstruction(u2,nx,up1_2m,up1_2p,dd,il,ir,coef,iorder,ighost)
        call rhs(up1_2m,up1_2p,nx,dx,res)

        do i = 1, nx                                
            t1 = 1.0/3
            t2 = 2.0/3
            su(i) = t1 * pu(i) + t2 * u2(i) + t2 * dt * res(i)
        enddo
        
        call boundary( su, nx, ighost )
           
        do i=-ighost,nx+ighost
            pu(i)=su(i)
        enddo
        
        t = t + dt
        if ( t + dt > supt ) then
            dt = supt - t
        endif 
    enddo       

     
    do i = 1, nx
        !if( x(i) > 0.2 .and. x(i) < 0.8 ) then
        error = error + abs( u0(i) - pu(i) )
        it = it + 1
        !endif  
    enddo
    error = error / it
     
    write(*, *) 'Final time:', t
    write(*, *) 'Error:', error
    
    open(1,file='temp.plt',status='unknown')
    do i=-ighost,nx+ighost
        write(1,101) x(i),pu(i)
    enddo
    close(1)
    open(2,file='solution.plt',status='unknown')
    do i=1,nx
        write(2,101) x(i),pu(i)
    enddo
    close(2)
    101 format(1x,e20.10,e20.10)
end program main