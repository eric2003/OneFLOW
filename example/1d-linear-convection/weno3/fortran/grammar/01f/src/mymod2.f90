module mymod2
implicit none 

   integer, parameter :: M = 256
   
contains      
   subroutine show_M()          
      print*, "M = ", M          
   end subroutine show_M 
   
end module mymod2