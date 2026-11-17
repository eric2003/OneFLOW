program main  
	implicit none
	integer :: i, j
	integer, parameter :: iorder = 2
    double precision :: coef(0:iorder-1, 0:iorder-1)
    data ((coef(i,j),j=0,iorder-1),i=0,iorder-1) &
         /0.5,0.5,-0.5,1.5/
    do i = 0, iorder-1
        do j = 0, iorder-1
            print *, 'coef(', i, ',', j, ') = ', coef(i,j)
        end do
    end do	
end program main 