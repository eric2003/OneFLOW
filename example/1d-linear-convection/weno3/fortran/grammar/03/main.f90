module mymod1  
implicit none 

   integer, parameter :: N = 1024 
   
contains      
   subroutine show_N()          
      print*, "N = ", N          
   end subroutine show_N 
   
end module mymod1 

module mymod2
implicit none 

   integer, parameter :: M = 256
   
contains      
   subroutine show_M()          
      print*, "M = ", M          
   end subroutine show_M 
   
end module mymod2

program main_prog     
use mymod1
use mymod2
implicit none     
   
   call show_N() 
   call show_M()
   
end program main_prog