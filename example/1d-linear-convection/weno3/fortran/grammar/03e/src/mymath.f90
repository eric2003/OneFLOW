module mymath_module  
	implicit none 
	integer, parameter :: N = 1024 
   
contains      
   subroutine show_N()          
      print*, "N = ", N          
   end subroutine show_N 
   
end module mymath_module 