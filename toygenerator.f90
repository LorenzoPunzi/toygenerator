program toygenerator

    implicit none
    character(len = 100) :: clarg
    integer :: nargs

    nargs = command_argument_count()
    if (nargs == 0 .or. nargs >= 2) then
        print *, 'Wrong input! Usage:'
        print *, '$toygenerator path/to/inputcard.dat' 
        print *, '$toygenerator search path/to/searchdirectory' 
        stop
    end if


end program toygenerator