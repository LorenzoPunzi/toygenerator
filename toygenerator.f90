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
    real(r14), parameter :: qqmin = 4.*mmu*mmu ! Minimum invariant mass squared of the mu-mu system



end module numeric    



module inputs

    use numeric
    implicit none
    character(len = 100) :: cardname, readline, opt, optval, evsave = '', histsave = ''  !!! extend it to be flexible
    integer :: nargs, iu, ios, seed(8), nbins = -999 !!! Add messages for when the necessary arguments are not initialised (-999)
    integer(i10) :: ngen = -999, nmax = -999
    logical :: exists, wghtopt = .false., isr = .false.
    real(r14) :: cme = -999, thmucutmin = 0, thmucutmax = 180, qqcutmin = 0, qqcutmax = -999, sinv, gmin = 0

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

                    else if (opt == 'nbins') then
                        read(optval, *, iostat=ios) nbins
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'nbins' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (seed(1)<=0) then
                            print*, "Input error! Number of histogram bins 'nbins' MUST be at least 1! Aborting..."
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
                    else if (opt == 'gmin') then
                        read(optval, *, iostat=ios) gmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'gmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (gmin < 0) then
                            print*, "Input error! Minimum energy of isr photon"&
                            " 'gmin' (",gmin, ") MUST be positive [GeV^2]! Aborting..."
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
            qqmax = sinv-2.*sqrt(sinv)*gmin ! Maximum invariant mass squared of the mu-mu system

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


module histogram

    use numeric
    use inputs
    implicit none
    real(r14), allocatable :: h_pmod1(:,:), h_thmu1(:,:)
    real(r14) :: ene1_max, pmod1_max 
    integer :: allerr
    contains
        subroutine inithistsborn()

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
            
        end subroutine inithistsborn

        subroutine updatehist(hist,histmin,histmax,val)
            real(r14), intent(inout) :: hist(:,:) 
            real(r14), intent(in) :: histmin,histmax,val
            integer :: bin

            bin = val/(histmax-histmin)*nbins
            hist(bin,1) = hist(bin,1) + 1
            
        end subroutine updatehist
    

end module histogram

module eventgen

    use numeric
    use inputs
    use histogram
    implicit none
    real(r14) :: intemax = 1., tmpmin = 0., tmpmax = 0., sinthmu, phimu, cosphimu, sinphimu
    real(r14) :: rndm(3), jacqq, jacgam
    real(r14), allocatable :: pmu1(:,:), pmu2(:,:), pmod1(:), pmod2(:), pgam(:,:), costhmu1(:), &
    costhmu2(:), costhgam(:), qq(:), accpt(:), wght(:) ! Reduce to fewer higher dimensional arrays 
    integer(i10)  :: naccpt = 0, iev
    integer :: run
    logical :: accepted



    contains

        subroutine genborn()

            call bornmuangles()
            
            pmu1(iev,4) = cme/2

            pmod1(iev) = sqrt( pmu1(iev,4)**2 - mmu**2 )

            pmu1(iev,1) = pmod1(iev) * sinthmu * cosphimu
            pmu1(iev,2) = pmod1(iev) * sinthmu * sinphimu
            pmu1(iev,3) = pmod1(iev) * costhmu1(iev)

            if (histsave /= '') then
            
                call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)))
            
            endif
            
            
        end subroutine genborn

        subroutine bornmuangles
            real(r14) :: inverseCDF
            costhmu1(iev) = inverseCDF(rndm(1))
            sinthmu = sqrt(1.0_r14-costhmu1(iev)**2)
            phimu = 2*pi*rndm(2)
            cosphimu = cos(phimu)
            sinphimu = sin(phimu)
            
        end subroutine bornmuangles

        subroutine genisr()

            call isrqq()
            call isrgamma()
            
            pmu1(iev,4) = cme/2

            pmod1(iev) = sqrt( pmu1(iev,4)**2 - mmu**2 )

            pmu1(iev,1) = pmod1(iev) * sinthmu * cosphimu
            pmu1(iev,2) = pmod1(iev) * sinthmu * sinphimu
            pmu1(iev,3) = pmod1(iev) * costhmu1(iev)

            if (histsave /= '') then
            
                call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)))
            
            endif
            
            
        end subroutine genisr

        subroutine isrqq

            real(r14) :: fak1, amin, amax, a, bmin, b, p, ppp, y
            fak1 = -1./sinv
            amin = fak1*log(sinv-qqmin)
            amax = fak1*log(sinv-qqmax)
            a = amax-amin
            bmin = log(qqmin/sinv)/sinv
            b    = log(qqmax/qqmin)/sinv
            p = rndm(1)   
            ppp  = a/(a+b)
            if (p < ppp) then
                y  = amin+a*x
                qq = sinv-dExp(y/fak1)                                       
            else
                y  = bmin+b*x
                qq = sinv*exp(sinv*y)
            endif
            jacqq = (a+b)/(1./(sinv*(sinv-qq)) + 1./sinv/qq)
            
        end subroutine isrqq

        subroutine isrgamma
            real(r14) :: x, phigam, b, cmin, cmax, y
            x = rndm(2)
            phigam = 2.*pi*rndm(3)
            b = sqrt(1.-4.*me*me/sinv)
            cmin = log((1.+b*cosmin)/(1.-b*cosmin))/(2.*b)
            cmax = log((1.+b*cosmax)/(1.-b*cosmax))/(2.*b)
            y = cmin+x*(cmax-cmin)
            costhgam(iev) = tanh(b*y)/b
            jacgam = 2.*pi*(1.-b*b*costheta*costheta)*(cmax-cmin)
            return
            
        end subroutine isrgamma

        subroutine testcuts

            if(isr .eqv. .false.) then

                if ( costhmu1(iev) <= cos(thmucutmin * degtorad) .and. costhmu1(iev) >= cos(thmucutmax * degtorad)) then
                    accepted = .true.
                else
                    accepted = .false.
                end if

            endif
            
        end subroutine testcuts
    

end module eventgen

module output

    use numeric
    use inputs
    use histogram
    use eventgen
    implicit none
    integer :: idx

    contains
        subroutine writevents
            if (evsave /= '') then
                open(newunit=iu, file=evsave, iostat=ios)
                if (ios == 0) then
                    print*, "Saving events to output file ", evsave
                    if ( isr .eqv. .false. ) then
                        write(iu, *) 'Output file of generation...' !!! More details on the generation
                        write(iu, '(*(A6, 3x))') 'px-',   'py-',  'pz-',    'E-',    'pmod-',    'th-',    'px+',   &
                        'py+', 'pz+', 'E+',  'pmod+',  'th+', 'qq' !!! More details on the generation
                        do idx = 1, ngen
                            if (accpt(idx) == 1) then
                                write (iu,'(*(f6.2, 3x))') pmu1(idx,1), pmu1(idx,2), pmu1(idx,3), pmu1(idx,4), pmod1(idx), &
                                radtodeg*acos(costhmu1(idx)), -pmu1(idx,1), -pmu1(idx,2), -pmu1(idx,3), &
                                pmu1(idx,4), pmod1(idx), 180-radtodeg*acos(costhmu1(idx)), sinv
                            endif
                        end do
                        close(unit=iu)
                    end if
                    if ( isr .eqv. .true. ) then
                        write(iu, *) 'Output file of generation...' !!! More details on the generation
                        write(iu, '(*(A6, 3x))') 'px-',   'py-',  'pz-',    'E-',    'pmod-',    'th-',    'px+',   &
                        'py+', 'pz+', 'E+',  'pmod+',  'th+', 'qq', 'pgamx', 'pgamy', 'pgamz', 'Egam' !!! More details on the generation
                        do idx = 1, ngen
                            if (accpt(idx) == 1) then
                                write (iu,'(*(f6.2, 3x))') pmu1(idx,1), pmu1(idx,2), pmu1(idx,3), pmu1(idx,4), pmod1(idx), &
                                radtodeg*acos(costhmu1(idx)), -pmu1(idx,1), -pmu1(idx,2), -pmu1(idx,3), &
                                pmu1(idx,4), pmod1(idx), 180-radtodeg*acos(costhmu1(idx)), sinv
                            endif
                        end do
                        close(unit=iu)
                    end if
                else
                    print *, 'ERROR while opening file ', evsave
                end if
            else    
                return
            endif
            
        end subroutine writevents

    

end module output


program toygenerator

    use numeric
    use inputs
    use histogram
    use eventgen
    use output
    implicit none
    integer :: err

    call loadinput()

    call random_seed(put=seed)

    allocate(pmu1(ngen,4), stat=err) !!! Maybe put the ones in common outside of if
    if (err /= 0) then
        print *, "pmu1 array allocation request denied, aborting!"
        print*, ""
        stop
    endif
    allocate(pmod1(ngen), stat=err)
    if (err /= 0) then
        print *, "pmod1 array allocation request denied, aborting!"
        print*, ""
        stop
    endif
    allocate(costhmu1(ngen), stat=err)
    if (err /= 0) then
        print *, "costhmu1 array allocation request denied, aborting!"
        print*, ""
        stop
    endif

    if (isr .eqv. .false.) then ! BORN CASE

        if (histsave /= '') call inithistsborn()

        allocate(accpt(ngen), stat=err)
        if (err /= 0) then
            print *, "accpt array allocation request denied, aborting!"
            print*, ""
            stop
        endif        

        do iev = 1, ngen !!! THINK ABOUT PARALLELISATION

            call random_number(rndm)
            call genborn()
            call testcuts()
            
            if (accepted) then

                accpt(iev) = 1
                naccpt = naccpt + 1 
            else
                accpt(iev) = 0
            endif

        end do

        call writevents()

    else

        allocate(pmu2(ngen,4), stat=err) 
        if (err /= 0) then
            print *, "pmu1 array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        allocate(pmod2(ngen), stat=err)
        if (err /= 0) then
            print *, "pmod1 array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        allocate(costhmu2(ngen), stat=err)
        if (err /= 0) then
            print *, "costhmu1 array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        allocate(pgam(ngen), stat=err)
        if (err /= 0) then
            print *, "pgam array allocation request denied, aborting!"
            print*, ""
            stop
        endif
        allocate(costhgam(ngen), stat=err)
        if (err /= 0) then
            print *, "costhgam array allocation request denied, aborting!"
            print*, ""
            stop
        endif

        if (wghtopt .eqv. .false.) then
            allocate(accpt(ngen), stat=err)
            if (err /= 0) then
                print *, "accpt array allocation request denied, aborting!"
                print*, ""
                stop
            endif
        else
            allocate(wght(ngen), stat=err)
            if (err /= 0) then
                print *, "wght array allocation request denied, aborting!"
                print*, ""
                stop
            endif
        endif
        
        
        !if (histsave /= '') call inithistsisr()
        run = 1

        do iev = 1, nmax !!! Think about parallelisation
 
             call random_number(rndm)
             call genisr()
            
        end do

         run = 2
         intemax = tmpmax

        do iev = 1, ngen

            call random_number(rndm)
            call genisr()
            
        end do

    endif
    

    !print '(f0.2)', cme
    !print *, wghtopt
    !print '(i0.2)', seed(1)
    !print *, evsave
    !print '(f0.2)', thmucutmin
    !!print '(f0.2)', thmucutmax
    !print '(f0.2)', qqcutmin
    !print '(f0.2)', qqcutmax





end program toygenerator