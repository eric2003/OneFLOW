program main_prog
implicit none

    call sub1
    call sub2
   
end program main_prog


subroutine sub1() 
implicit none

    print*,"haha1"
    
end subroutine sub1

subroutine sub2() 
implicit none

    print*,"haha2"
    
end subroutine sub2