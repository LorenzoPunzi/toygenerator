module numeric

    implicit none
    integer, parameter :: r14 = selected_real_kind(14,99)
    integer, parameter :: i10 = selected_int_kind(10)
    real(r14), parameter :: pi = 4.0 * atan(1.0)  ! pi = 4 * arctan(1) gives the value of pi
    real(r14), parameter :: alpha = 1/(137.035999174)  ! alpha QED
    real(r14), parameter :: mmu = 0.10565837  ! muon mass [GeV]
    real(r14), parameter :: me = 0.000510998950  ! electron mass [GeV]
    real(r14), parameter :: gev2nbarn = 389379.292  ! Conversion from GeV^2 to nbarn
    real(r14), parameter :: degtorad = pi/180.  ! Conversion from degrees to radians
    real(r14), parameter :: radtodeg = 180./pi  ! Conversion from radians to degrees



end module numeric    



module inputs

    use numeric
    implicit none
    character(len = 100) :: cardname, readline, opt, optval, evsave = '', histsave = ''  !!! extend it to be flexible
    integer :: nargs, iu, ios, seed(8)
    integer(i10) :: ngen = -999, nmax = -999
    logical :: exists, wghtopt = .false., isr = .false.
    real(r14) :: cme = -999, thmucutmin = 0, thmucutmax = 180, qqcutmin = 0, qqcutmax = -999, sinv

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
                        read(optval, *, iostat=ios) seed(1)
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'seed' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (seed(1)<0) then
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
                            isr = .true.
                        else if (trim(optval) == 'no' .or. trim(optval) == 'n') then
                            isr = .false.
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

                    else if (opt == 'thmucutmin') then
                        read(optval, *, iostat=ios) thmucutmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'thmucutmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (thmucutmin < 0 .or. thmucutmin > 180) then
                            print*, "Input error! Polar angle cut minumum for muons 'thmucutmin'&
                             (",thmucutmin, ") MUST be a real number between 0 and 180! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'thmucutmax') then
                        read(optval, *, iostat=ios) thmucutmax
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'thmucutmax' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (thmucutmax < 0 .or. thmucutmax > 180) then
                            print*, "Input error! Polar angle cut maximum for muons 'thmucutmax'&
                             (",thmucutmax, ") MUST be a real number between 0 and 180! Aborting..."
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


            ! Processing of loaded input

            sinv = cme**2

            if (qqcutmax == -999) qqcutmax = sinv !!! Numerical problem?

            if (thmucutmin >= thmucutmax) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'thmucutmin' value (", thmucutmin, ")&
                must be lower than 'thmucutmax' value (", thmucutmin, ")"
                print *, ""
                stop    
            endif

            if (qqcutmin >= qqcutmax .and. qqcutmax >= 0) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'qqcutmin' value (", qqcutmin, ") must &
                be lower than 'qqcutmax' value (", qqcutmax, ")"
                print *, ""
                stop    
            endif

            if (qqcutmin >= sinv) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'qqcutmin' value &
                (", qqcutmin, ") must be lower than 'cme' value squared (", sinv, ")"

                print *, ""
                stop    
            endif

            
            
        end subroutine loadinput

end module inputs

real(r14) function inverseCDF(x) !!! Putting it inside the module eventgen doesn't work

            use numeric
            real(r14), intent(in) :: x
            
            inverseCDF = 1/(2 - 4*x + sqrt(5 - 16*x + 16* x**2))**(1./3.) - (2 - 4*x +&
             sqrt(5 - 16*x + 16*x**2))**(1./3.)
    
end function inverseCDF

module eventgen

    use numeric
    use inputs
    
    implicit none
    real(r14) :: intemax = 1., tmpmin = 0., tmpmax = 0., sinthmu, phimu, cosphimu, sinphimu
    real(r14) :: pmod1, pmod2
    real(r14) :: rndm(2)
    real(r14), allocatable :: pmu1(:,:), pmu2(:,:), pgamma(:,:), costhmu(:), qq(:)
    integer(i10)  :: naccpt = 0, iev
    integer :: run
    logical :: accepted



    contains

        subroutine genborn()

            call bornthetmu()
            call bornphimu()
            
            pmu1(iev,4) = cme/2

            pmod1 = sqrt( pmu1(iev,4)**2 - mmu**2 )

            pmu1(iev,1) = pmod1 * sinthmu * cosphimu
            pmu1(iev,2) = pmod1 * sinthmu * sinphimu
            pmu1(iev,3) = pmod1 * costhmu(iev)

            pmu2(iev,1) = -pmu1(iev,1)
            pmu2(iev,2) = -pmu1(iev,2)
            pmu2(iev,3) = -pmu1(iev,3)
            
            
        end subroutine genborn
        
        subroutine bornthetmu
            real(r14) :: inverseCDF
            costhmu(iev) = inverseCDF(rndm(1))
            sinthmu = sqrt(1.0_r14-costhmu(iev)**2)
            
        end subroutine bornthetmu

        subroutine bornphimu
            phimu = 2*pi*rndm(2)
            cosphimu = cos(phimu)
            sinphimu = sin(phimu)
            
        end subroutine bornphimu

        subroutine testcuts

            if(isr .eqv. .false.) then

                if ( costhmu(iev) <= cos(thmucutmin * degtorad) .and. costhmu(iev) >= cos(thmucutmax * degtorad)) then
                    accepted = .true.
                else
                    accepted = .false.
                end if

            endif
            
        end subroutine testcuts
    

end module eventgen

module histogram

    use numeric
    use inputs
    implicit none
    real(r14), allocatable :: h_pmod1(:,:), h_thmu1(:,:)
    real(r14) :: ene1_max, pmod1_max 
    integer :: allerr
    contains
        subroutine inithistosborn()

            ene1_max = cme/2
            pmod1_max = sqrt(ene1_max**2-mmu**2)
            allocate(h_pmod1(nbins,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu1(nbins,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu- histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            
        end subroutine inithistosborn
    

end module histogram



program toygenerator

    use numeric
    use inputs
    use eventgen
    implicit none
    integer :: err

    call loadinput()

    call random_seed(put=seed)

    if (isr .eqv. .false.) then ! BORN CASE

        allocate(pmu1(ngen,4), stat=err)
        if (err /= 0) then
            print *, "pmu1 array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        allocate(pmu2(ngen,4), stat=err)
        if (err /= 0) then
            print *, "pmu2 array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        allocate(costhmu(ngen), stat=err)
        if (err /= 0) then
            print *, "costhmu array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        

        do iev = 1, ngen !!! THINK ABOUT PARALLELISATION

            call random_number(rndm)
            call genborn()
            call testcuts()
            
            if (accepted) then
            
                print*, iev, radtodeg*acos(costhmu(iev))
                naccpt = naccpt + 1 


            endif

            
        end do

    endif

!        run = 1
 !        do iev = 1, nmax !!! Think about parallelisation
! 
!             call random_number(rndm)
!             call genborn()
            
 !        end do

!         run = 2
!         intemax = tmpmax

 ! !        do iev = 1, ngen

 !            call random_number(rndm)
 !            call genborn()
            
 !        end do
 !    endif
    

    print '(f0.2)', cme
    print *, wghtopt
    print '(i0.2)', seed(1)
    print *, evsave
    print '(f0.2)', thmucutmin
    print '(f0.2)', thmucutmax
    print '(f0.2)', qqcutmin
    print '(f0.2)', qqcutmax





end program toygenerator