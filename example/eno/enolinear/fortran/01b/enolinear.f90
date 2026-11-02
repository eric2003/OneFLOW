subroutine reconstruction(u,nx,up1_2,dd,ir,coef,iorder,ighost)
    implicit none
    real(8) :: u(-ighost:nx+ighost)
    real(8) :: up1_2(0:nx)
    real(8) :: dd(0:ighost-1,-ighost:nx+ighost)
    real(8) :: coef(0:iorder-1,0:iorder-1)
    integer :: ir(0:nx)
    integer :: i, j, k, l, m
    integer :: nx, iorder, ighost
      
!     Choose the stencil by ENO method
    do j=-ighost,nx+ighost
        dd(0,j) = u(j)
    enddo
    do i=1,iorder-1
        do j=-ighost,nx+ighost-1
            dd(i,j) = dd(i-1,j+1) - dd(i-1,j)
        enddo
    enddo 
    do j=0,nx
        ir(j) = j
        do i=1,iorder-1
            if ( abs(dd(i,ir(j)-1)).le.abs(dd(i,ir(j))) ) then
                ir(j)=ir(j)-1 
            endif
        enddo
    enddo
    
!   Reconstruction u(j+1/2)
    do j=0,nx         
        k = ir(j)
        l = j - k  
        up1_2(j) = 0
        do m=0,iorder-1
            up1_2(j) = up1_2(j) + u(k+m) * coef(l,m)
        enddo  
    enddo        
end subroutine reconstruction       
      
program main  
    implicit none
    integer, parameter :: nx = 40
    integer, parameter :: ighost = 10
    integer, parameter :: iorder = 2
    real(8), parameter :: pi = 3.14159265358979323846
    real(8) :: pu(-ighost:nx+ighost), su(-ighost:nx+ighost)
    real(8) :: u1(-ighost:nx+ighost), u2(-ighost:nx+ighost)
    real(8) :: u0(1:nx), up1_2(0:nx), x(-ighost:nx+ighost)
    real(8) :: dd(0:ighost-1, -ighost:nx+ighost)
    real(8) :: coef(0:iorder-1, 0:iorder-1)
    integer :: ir(0:nx)
    integer :: i, j, count
    real(8) :: supt, t, temp, t1, t2, error, it
    real(8) :: dx, dt
            
      open(unit=1,file='\temp2.plt',status='unknown')

      dx=1.0/nx
      dt=dx*0.5
      write(*,*) 'Input T:'
      read(*,*) supt 
      
!     2nd-order coefficients       
      data ((coef(i,j),j=0,iorder-1),i=0,iorder-1) &
     &     /0.5,0.5,-0.5,1.5/
!     Initialize grid and initial conditions
      do i=-ighost,nx+ighost
          x(i)=(i-1)*dx+dx/2
      enddo
      
    count = 0  ! 初始化计数器

    do i=-ighost,nx+ighost
        write(*,'(1x,f10.6)',advance='no') x(i)
        count = count + 1
        if (count == 5) then
            write(*,*)  ! 换行
            count = 0  ! 重置计数器
        endif
    enddo
      
      
!     Initial mean value 1      
    do i=1,nx 
       pu(i)=-(cos(2.0*pi*(x(i)+dx/2))-cos(2.0*pi*(x(i)-dx/2))) &
                /(dx*2.0*pi)
    enddo    
      do i=0,-ighost,-1
          pu(i)=pu(i+nx)
      enddo 
      do i=nx+1,nx+ighost
          pu(i)=pu(i-nx)
      enddo  
        
      do i=0,nx
        u0(i)=pu(i)
      enddo                           

!     Time stepping
      t=0                      
      do while(t.lt.supt) 
        call reconstruction(pu,nx,up1_2,dd,ir,coef,iorder,ighost) 
        do i=1,nx        
          temp=up1_2(i)-up1_2(i-1)
          u1(i)=pu(i)-dt*temp/dx
        enddo
        do i=0,-ighost,-1
          u1(i)=u1(i+nx)
        enddo 
        do i=nx+1,nx+ighost
          u1(i)=u1(i-nx)
        enddo  
        
        call reconstruction(u1,nx,up1_2,dd,ir,coef,iorder,ighost)
        do i=1,nx   
          temp=up1_2(i)-up1_2(i-1)
          u2(i)=3.0/4.0*pu(i)+1.0/4.0*(u1(i)-dt*temp/dx)
        enddo  
        do i=0,-ighost,-1
          u2(i)=u2(i+nx)
        enddo 
        do i=nx+1,nx+ighost
          u2(i)=u2(i-nx)
        enddo  
        
        call reconstruction(u2,nx,up1_2,dd,ir,coef,iorder,ighost)
        do i=1,nx                                
          temp=up1_2(i)-up1_2(i-1)   
          t1=1.0/3 
          t2=2.0/3
          su(i)=t1*pu(i)+t2*(u2(i)-dt*temp/dx)
        enddo
        do i=0,-ighost,-1
          su(i)=su(i+nx)
        enddo 
        do i=nx+1,nx+ighost
          su(i)=su(i-nx)
        enddo 
           
        do i=-ighost,nx+ighost
          pu(i)=su(i)
        enddo
        t=t+dt
        if(t+dt.gt.supt) then
          dt=supt-t
        endif 
      enddo       

!     Output results to Tecplot file
      open(unit=2, file='solution.plt', status='unknown')
      write(2, *) 'TITLE = "Numerical and Exact Solutions"'
      write(2, *) 'VARIABLES = "x", "Numerical", "Exact"'
      write(2, *) 'ZONE T="Solution", I=', nx + 1, ', F=POINT'
      do i = 0, nx
          write(2, *) x(i), pu(i), sin(2 * pi * (x(i) - t))
      end do
      close(2)
      
      do i=1,nx
!        if(x(i).gt.0.2.and.x(i).lt.0.8) then
          error=error+abs(u0(i)-pu(i))      
          it=it+1
!        endif  
      enddo
      error=error/it
     
      write(*, *) 'Final time:', t
      write(*, *) 'Error:', error
	
end program main