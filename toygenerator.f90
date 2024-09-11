module numeric

    implicit none
    integer, parameter :: r14 = selected_real_kind(14,99)
    integer, parameter :: i10 = selected_int_kind(10)
    real(r14) :: pi = 4.0 * atan(1.0)  ! pi = 4 * arctan(1) gives the value of pi
    real(r14) :: alpha = 1/(137.035999174)  ! alpha QED
    real(r14) :: mmu = 0.10565837  ! muon mass [GeV]
    real(r14) :: me = 0.000510998950  ! electron mass [GeV]

    


end module numeric    



module inputs

    implicit none
    integer, parameter :: r14 = selected_real_kind(14,99)
    integer, parameter :: i10 = selected_int_kind(10)
    character(len = 100) :: cardname, readline, opt, optval, evsave = '', histsave = ''  !!! extend it to be flexible
    integer :: nargs, iu, ios, seed
    integer(i10) :: ngen = -999, nmax = -999
    logical :: exists, wghtopt = .false., isropt = .false.
    real(r14) :: cme = -999, muthcutmin = 0, muthcutmax = 180, qqcutmin = 0, qqcutmax = -999

    contains
        subroutine loadinput() 
            nargs = command_argument_count()
            if (nargs == 0 .or. nargs >= 2) then
                print *, "Wrong input! Usage:"
                print *, "$toygenerator path/to/inputcard.dat"
                print *, "$toygenerator search path/to/searchdirectory"
                stop
            end if

            call get_command_argument(1,cardname)

            inquire(file=cardname, exist=exists)
            if(exists) then
                print *, "Reading '" , trim(cardname), "' as input card"
                print*, ""
            else
                print *, "File '", trim(cardname), "' does not exist, aborting!"
                print*, ""
                stop
            end if

            open(newunit=iu, file=cardname, iostat = ios, iomsg=readline, action='read')
            if (ios /= 0) then
                print *, readline
                stop
            end if
            do
                read(iu, '(a)', iostat=ios) readline
                if (ios == 0) then
                    if (readline ==' ') exit
                    read(readline, *) opt, optval

                    !!! Might not be optimal, for each line we check if it's cme, seed, etc
                    if (opt == 'cme') then
                        read(optval, *, iostat=ios) cme
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'cme' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (cme<=0) then
                            print*, "Input error! Center of mass energy 'cme' MUST be positive real number [GeV]! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'weight') then
                        if (trim(optval) == 'yes' .or. trim(optval) == 'y') then
                            wghtopt = .true.
                        else if (trim(optval) == 'no' .or. trim(optval) == 'n') then
                            wghtopt = .false.
                        else
                            print *, "Input error! Invalid value for 'weight' in input card. &
                            Only acceptable options are 'yes'/'y' and 'no'/'n'. Aborting!"
                            print*, ""
                            stop
                        endif

                    else if (opt == 'seed') then
                        read(optval, *, iostat=ios) seed
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'seed' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (seed<0) then
                            print*, "Input error! Seed of the generation 'seed' MUST be positive! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'ngen') then
                        read(optval, *, iostat=ios) ngen
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'ngen' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (ngen<=0) then
                            print*, "Input error! Number of generated events 'ngen' MUST be positive! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'nmax') then
                        read(optval, *, iostat=ios) nmax
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'nmax' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (nmax<=0) then
                            print*, "Input error! Number of events to find maximum 'nmax' MUST be positive! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'isr') then
                        if (trim(optval) == 'yes' .or. trim(optval) == 'y') then
                            isropt = .true.
                        else if (trim(optval) == 'no' .or. trim(optval) == 'n') then
                            isropt = .false.
                        else
                            print *, "Input error! Invalid value for 'isr' in input card.&
                             Only acceptable options are 'yes'/'y' and 'no'/'n'. Aborting!"
                            print*, ""
                            stop
                        endif

                    else if (opt == 'evsave') then
                        evsave = trim(optval)

                    else if (opt == 'histsave') then
                        histsave = trim(optval)

                    else if (opt == 'muthcutmin') then
                        read(optval, *, iostat=ios) muthcutmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'muthcutmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (muthcutmin < 0 .or. muthcutmin > 180) then
                            print*, "Input error! Polar angle cut minumum for muons 'muthcutmin'&
                             (",muthcutmin, ") MUST be a real number between 0 and 180! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'muthcutmax') then
                        read(optval, *, iostat=ios) muthcutmax
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'muthcutmax' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (muthcutmax < 0 .or. muthcutmax > 180) then
                            print*, "Input error! Polar angle cut maximum for muons 'muthcutmax'&
                             (",muthcutmax, ") MUST be a real number between 0 and 180! Aborting..."
                            print *, ""
                            stop
                        endif
                    else if (opt == 'qqcutmin') then
                        read(optval, *, iostat=ios) qqcutmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'qqcutmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (qqcutmin < 0) then
                            print*, "Input error! Minimum of invariant mass squared of outgoing"&
                            " muons 'qqcutmin' (",qqcutmin, ") MUST be positive [GeV^2]! Aborting..."
                            print *, ""
                            stop
                        endif
                    else if (opt == 'qqcutmax') then
                        read(optval, *, iostat=ios) qqcutmax
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'qqcutmax' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (qqcutmax < 0) then
                            print*, "Input error! Maximum of invariant mass squared of outgoing"&
                            " muons 'qqcutmax' (",qqcutmax, ") MUST be positive [GeV^2]! Aborting..."
                            print *, ""
                            stop
                        endif

                    endif
                else
                    exit
                end if
            end do

            if (qqcutmax == -999) qqcutmax = cme**2 !!! Numerical problem?

            if (muthcutmin >= muthcutmax) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'muthcutmin' value (", muthcutmin, ")&
                must be lower than 'muthcutmax' value (", muthcutmin, ")"
                print *, ""
                stop    
            endif

            if (qqcutmin >= qqcutmax .and. qqcutmax >= 0) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'qqcutmin' value (", qqcutmin, ") must &
                be lower than 'qqcutmax' value (", qqcutmax, ")"
                print *, ""
                stop    
            endif

            if (qqcutmin >= cme**2) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'qqcutmin' value &
                (", qqcutmin, ") must be lower than 'cme' value squared (", cme**2, ")"

                print *, ""
                stop    
            endif
            
        end subroutine loadinput

end module inputs



program toygenerator

    use numeric
    use inputs

    implicit none



    call loadinput()

    

    print '(f0.2)', cme
    print *, wghtopt
    print '(i0.2)', seed
    print *, evsave
    print '(f0.2)', muthcutmin
    print '(f0.2)', muthcutmax
    print '(f0.2)', qqcutmin
    print '(f0.2)', qqcutmax





end program toygenerator