module numeric

    implicit none
    integer, parameter :: r14 = selected_real_kind(14,99) 
    integer, parameter :: i10 = selected_int_kind(10)
    real(r14), parameter :: pi = 4.0_r14 * atan(1.0_r14)  ! pi = 4 * arctan(1) gives the value of pi
    real(r14), parameter :: alpha = 1.0_r14/(137.0359991740_r14)  ! alpha QED
    real(r14), parameter :: mmu = 0.105658370_r14  ! muon mass [GeV]
    real(r14), parameter :: me = 0.000510998950_r14 ! electron mass [GeV]
    real(r14), parameter :: gev2nbarn = 389379.2920_r14  ! Conversion from GeV^2 to nbarn
    real(r14), parameter :: degtorad = pi/180.0_r14  ! Conversion from degrees to radians
    real(r14), parameter :: radtodeg = 180.0_r14/pi  ! Conversion from radians to degrees
    real(r14), parameter :: qqmin = 4.0_r14*mmu*mmu ! Minimum invariant mass squared of the mu-mu system

    contains
        function cubicroot(x) 
            real(r14), intent(in) :: x
            real(r14) :: cubicroot
            
            cubicroot = sign(1.0_r14, x)* (abs(x))**(1.0_r14/3.)
    
        end function cubicroot

        function compareals(x,y)
            real(r14), intent(in) :: x, y
            logical :: compareals

            if ( abs(x-y) < epsilon(x) ) then
                compareals = .true.
            else 
                compareals = .false.
            end if

        end function compareals




end module numeric    



module inputs

    use numeric
    implicit none
    character(len = 100) :: cardname, readline, opt, optval, evsave = '', histsave = ''
    integer :: nargs, iu, ios, seed(8), nbins
    integer(i10) :: ngen, nmax
    logical :: exists, wghtopt = .false., isr = .false., cmeflg = .false., ngenflg = .false., &
    nmaxflg = .false., isrflg = .false., gminflg = .false., seedflg = .false., nbinflg = .false., qqcutmaxflg = .false.
    real(r14) :: cme, thmucutmin = 0.0_r14, thmucutmax = 180.0_r14, qqcutmin = 0.0_r14, qqcutmax, &
    sinv, gmin, qqmax, thgamcutmin = 0.0_r14, thgamcutmax = 180.0_r14, costhgammin, &
    costhgammax, costhmumin, costhmumax, intemax = 1.0_r14, bornsigma

    contains
        subroutine loadinput() 
            nargs = command_argument_count()
            if (nargs == 0 .or. nargs >= 2) then
                print*, ""
                print *, "Wrong input!"
                print*, ""
                print*, "MC generation usage:"
                print*, ""
                print*, "$make"
                print *, "$./toygenerator path/to/inputcard"
                print*, ""
                print*, "For help mode use"
                print*, ""
                print*, "$make"
                print *, "$./toygenerator --help"
                print*, "or"
                print*, "$make"
                print *, "$./toygenerator -h"
                print*, ""
                stop
            end if

            print*, ""
            print*, "***********************************************************************************"
            print*, "***********************************************************************************"
            print*, ""
            print*, ""
            print*, "       WELCOME TO THE e+ e- ==> mu+ mu- (gamma) TOY MONTE CARLO GENERATOR!"
            print*, ""
            print*, ""
            print*, "***********************************************************************************"
            print*, "***********************************************************************************"
            print*, ""
            print*, ""

            call get_command_argument(1,cardname)

            if ( cardname == '--help' .or. cardname == '-h') then
                print*, "============================================================================================"
                print*, ""
                print*, "                              TOYGENERATOR HELP MODE"
                print*, ""
                print*, "The Monte Carlo generator simulates e+ e- ==> mu+ mu- (gamma) interactions"
                print*, "where an ISR gamma can be present in the final state. It generates the relevant"
                print*, "kinematical quantities event by event, and calculates the total cross section"
                print*, "with the selected cuts. The cross section, the seed used and other important"
                print*, "information are stored in a file named 'info.out' in the same directory the "
                print*, "program is run in. The generation can be performed in a 'weighted' manner, "
                print*, "where all events are accepted but stored witha  weight < 1 and equal to their"
                print*, "acceptance probability."
                print*, ""
                print*, "The Monte Carlo generator requires a configuration input card to be"
                print*, "given as input (e.g. '$./toygenerator path/to/inputcard'). The input card"
                print*, "should be a column file with no header. Each variable is separated by an"
                print*, "arbitrary number of spaces from its variable. Nothing is read after the"
                print*, "first empty line. Before running the generator, compile with '$make'."
                print*, ""
                print*, "Required input:"
                print*, ""
                print*, "* seed : generation seed used (positive integer)"
                print*, "* ngen : number of events the program should attempt to generate (positive integer)"
                print*, "* nmax : number of events to find the integrand maximum (ISR only, positive integer)"
                print*, "* isr : generation mode, whether BORN ('n','no') or ISR ('y','yes')"
                print*, "* cme : center of mass energy of the initial state (positive real)[GeV]"
                print*, "* gmin : minimum energy of ISR photon (positive real)[GeV]"
                print*, "* nbins : number of histogram bins (positive integer), mandatory if 'histsave' is given"
                print*, ""
                print*, "Optional input:"
                print*, ""
                print*, "* evsave : path to file where events should be written"
                print*, "* histsave : path to file where histograms should be written"
                print*, "* weight : generation mode, whether weighted ('y','yes') or not ('n','no'). When enabled, "
                print*, "           ALL events are stored, with a weight equal to their acceptance probability."
                print*, "           This mode is only available in ISR mode. By default the generation is not weighted."
                print*, "-------------------------------------------CUTS--------------------------------------------------"
                print*, "* qqcutmin : minimum invariant mass of the muon system (positive real)[GeV^2]"
                print*, "* qqcutmax : maximum invariant mass of the muon system (positive real)[GeV^2]"
                print*, "* thmucutmin : minimum polar angle of muons in the center of mass system (positive real)[degrees]"
                print*, "* thmucutmax : maximum polar angle of muons in the center of mass system (positive real)[degrees]"
                print*, "* thgamcutmin : minimum polar angle of ISR photon in the center of mass system"
                print*, "                (positive real)[degrees]. This cut is only available in ISR mode."
                print*, "* thgamcutmax : maximum polar angle of ISR photon in the center of mass system"
                print*, "                (positive real)[degrees]. This cut is only available in ISR mode."
                print*, ""
                print*, "============================================================================================"
                print*, ""
                stop
            end if

            inquire(file=cardname, exist=exists)
            if(exists) then
                print *, "Reading '" , trim(cardname), "' as input card..."
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

                    if (opt == 'cme') then
                        cmeflg = .true.
                        read(optval, *, iostat=ios) cme
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'cme' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (cme<= 0.0_r14) then
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
                            &Only acceptable options are 'yes'/'y' and 'no'/'n'. Aborting!"
                            print*, ""
                            stop
                        endif

                    else if (opt == 'seed') then
                        seedflg = .true.
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
                        nbinflg = .true.
                        read(optval, *, iostat=ios) nbins
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'nbins' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (nbins<=0) then
                            print*, "Input error! Number of histogram bins 'nbins' MUST be at least 1! Aborting..."
                            print *, ""
                            stop
                        endif

                    else if (opt == 'ngen') then
                        ngenflg = .true.
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
                        nmaxflg = .true.
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
                        isrflg = .true.
                        if (trim(optval) == 'yes' .or. trim(optval) == 'y') then
                            isr = .true.
                        else if (trim(optval) == 'no' .or. trim(optval) == 'n') then
                            isr = .false.
                        else
                            print *, "Input error! Invalid value for 'isr' in input card.&
                             &Only acceptable options are 'yes'/'y' and 'no'/'n'. Aborting!"
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
                        if (thmucutmin < 0.0_r14 .or. thmucutmin > 180.0_r14) then
                            print*, "Input error! Polar angle cut minumum for muons 'thmucutmin'&
                             &(",thmucutmin, ") MUST be a real number between 0 and 180! Aborting..."
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
                        if (thmucutmax < 0.0_r14 .or. thmucutmax > 180.0_r14) then
                            print*, "Input error! Polar angle cut maximum for muons 'thmucutmax'&
                             &(",thmucutmax, ") MUST be a real number between 0 and 180! Aborting..."
                            print *, ""
                            stop
                        endif
                    else if (opt == 'thgamcutmin') then
                        read(optval, *, iostat=ios) thgamcutmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'thgamcutmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (thgamcutmin < 0.0_r14 .or. thgamcutmin > 180.0_r14) then
                            print*, "Input error! Polar angle cut minumum for ISR photon 'thgamcutmin'&
                            & (",thgamcutmin, ") MUST be a real number between 0 and 180! Aborting..."
                            print *, ""
                            stop
                        endif
                    else if (opt == 'thgamcutmax') then
                        read(optval, *, iostat=ios) thgamcutmax
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'thgamcutmax' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (thgamcutmax < 0.0_r14 .or. thgamcutmax > 180.0_r14) then
                            print*, "Input error! Polar angle cut maximum for ISR photon 'thgamcutmax'&
                            & (",thgamcutmax, ") MUST be a real number between 0 and 180! Aborting..."
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
                        if (qqcutmin < 0.0_r14) then
                            print*, "Input error! Minimum of invariant mass squared of outgoing"&
                            " muons 'qqcutmin' (",qqcutmin, ") MUST be positive [GeV^2]! Aborting..."
                            print *, ""
                            stop
                        endif
                    else if (opt == 'qqcutmax') then
                        qqcutmaxflg = .true.
                        read(optval, *, iostat=ios) qqcutmax
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'qqcutmax' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (qqcutmax < 0.0_r14) then
                            print*, "Input error! Maximum of invariant mass squared of outgoing"&
                            " muons 'qqcutmax' (",qqcutmax, ") MUST be positive [GeV^2]! Aborting..."
                            print *, ""
                            stop
                        endif
                    else if (opt == 'gmin') then
                        gminflg = .true.
                        read(optval, *, iostat=ios) gmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'gmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (gmin < 0.0_r14) then
                            print*, "Input error! Minimum energy of isr photon"&
                            " 'gmin' (",gmin, ") MUST be positive [GeV]! Aborting..."
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
            if (isr .eqv. .true.) qqmax = sinv-2.0_r14*sqrt(sinv)*gmin ! Maximum invariant mass squared of the mu-mu system
            costhgammin = cos(radtodeg*thgamcutmax)
            costhgammax = cos(radtodeg*thgamcutmin)
            costhmumin = cos(radtodeg*thmucutmax)
            costhmumax = cos(radtodeg*thmucutmin)

            if (qqcutmaxflg .eqv. .false.) qqcutmax = sinv 

            if (thmucutmin >= thmucutmax) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'thmucutmin' value (", thmucutmin, ")&
                &must be lower than 'thmucutmax' value (", thmucutmin, ")"
                print *, ""
                stop    
            endif

            if (thgamcutmin >= thgamcutmax) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'thgamcutmin' value (", thgamcutmin, ")&
                &must be lower than 'thgamcutmax' value (", thgamcutmax, ")"
                print *, ""
                stop    
            endif

            if (qqcutmin >= qqcutmax .and. qqcutmax >= 0) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'qqcutmin' value (", qqcutmin, ") must &
                &be lower than 'qqcutmax' value (", qqcutmax, ")"
                print *, ""
                stop    
            endif

            if (qqcutmin >= sinv) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'qqcutmin' value &
                &(", qqcutmin, ") must be lower than 'cme' value squared (", sinv, ")"
                print *, ""
                stop    
            endif

            if ( seedflg .eqv. .false. ) then
                print*, "No value given for generation seed 'seed' in input card! Aborting..."
                print*, ""
                stop
            end if
            if ( isrflg .eqv. .false. ) then
                print*, "No option given ('n'/'no','y'/'yes') for ISR generation 'isr' in input card! Aborting..."
                print*, ""
                stop
            end if

            if ( cmeflg .eqv. .false. ) then
                print*, "No value given for center of mass energy [GeV] 'cme' in input card! Aborting..."
                print*, ""
                stop
            end if

            if ( ngenflg .eqv. .false. ) then
                print*, "No value given for number of events to generate 'ngen' in input card! Aborting..."
                print*, ""
                stop
            end if
            if ( (isr .eqv. .true.) .and. (nmaxflg .eqv. .false.) ) then
                print*, "No value given for number of events to generate in order to find integrand&
                & maximum 'nmax' in input card! Aborting..."
                print*, ""
                stop
            end if
            if ( (isr .eqv. .true.) .and. (gminflg .eqv. .false.) ) then
                print*, "No value given for minimum ISR photon energy [GeV] 'gmin' in input card! Aborting..."
                print*, ""
                stop
            end if
            if ( (histsave /= '') .and. (nbinflg .eqv. .false.) ) then
                print*, "No bin number for the histograms 'nbin' given in input card! Aborting..."
                print*, ""
                stop
            end if

        close(unit=iu)

        print*, "Loading of input card completed!"
        print*, ""
        print*, "============================================"
        print*, ""
        print '(A, i0)', "    Seed : ", seed(1)
        print*, ""
        if ( isr .eqv. .false. ) then
            print '(A)', "    Generation mode : BORN"
            print*, ""
        else
            print '(A)', "    Generation mode : ISR"
            print*, ""
            if (wghtopt .eqv. .true. ) then
                print '(A)', "    Weighted generation : ON"
                print*, ""
            else 
                print '(A)', "    Weighted generation : OFF"
                print*, ""
            end if
        end if
        print '(A, i0)', "    Number of events to generate : ", ngen
        print*, ""
        if ( isr .eqv. .true. ) then
            print '(A, i0)', "    Number of events&
            & to find maximum of integrand : ", nmax
            print*, ""
        endif 
        print '(A, f0.2, A)', "    Center of mass energy = ", cme, " GeV"
        print*, "" 
        if ( compareals(thmucutmin,0.0_r14) .and. compareals(thmucutmax,180.0_r14) ) then
            print '(A)', "    No angular cuts on muons"
        else 
            print '(A, f0.2, A, f0.2)', "    Angular cuts on muons : theta min = ", thmucutmin, " theta max =  ", thmucutmax
        end if 
        print*, ""
        if ( isr .eqv. .true. ) then
            if ( compareals(thgamcutmin,0.0_r14) .and. compareals(thgamcutmax,180.0_r14) ) then
                print '(A)',  "    No angular cuts on ISR photon"
            else 
                print '(A, f0.2, A, f0.2)', "    Angular cuts on ISR photon :&
                & theta min = ", thgamcutmin, " theta max =  ", thgamcutmax
            end if 
            print*, ""
        endif
        if ( compareals(qqcutmin,0.0_r14) .and. compareals(qqcutmax,sinv) ) then
            print '(A)', "    No cuts on muon system invariant mass"
        else 
            print '(A, f0.2, A, f0.2)', "    Angular cuts on muons' invariant mass: qq min = ", qqcutmin, " qq max =  ", qqcutmax
        end if
        print*, ""
        if ( isr .eqv. .true. ) then
            print '(A, f0.2, A)', "    Minimum ISR photon energy = ", gmin, " GeV"
            print*, ""
        endif 
        print*, "============================================"
        print*, ""
        
            
            
        end subroutine loadinput

end module inputs

module utilities
    use numeric
    implicit none

    contains
        function inverseCDF(x)
            real(r14), intent(in) :: x
            real(r14) :: inverseCDF
            
            inverseCDF = 1/cubicroot(2 - 4*x + sqrt(5 - 16*x + 16* x**2)) - cubicroot(2 - 4*x +&
            sqrt(5 - 16*x + 16*x**2))
    
        end function inverseCDF

        subroutine deboost(boosted, boostvector, deboosted, ii)

            real(r14), intent(inout) :: boosted(0:)
            real(r14), intent(inout) :: boostvector(0:)
            real(r14), intent(inout) :: deboosted(:,0:)
            integer(i10), intent(in) :: ii
            real(r14) :: lormtrx(0:3,0:3), boostE, boostp, beta, gamma, boostcosth, boostsinth, &
            boostcosphi, boostsinphi
            integer :: k, j

            boostE = boostvector(0)
            boostp = sqrt(boostvector(1)**2+boostvector(2)**2+ boostvector(3)**2)
            beta  = boostp/boostE
            gamma = 1.0_r14/sqrt(1.0_r14-beta**2)
            ! Invert the spatial components since we are DEBOOSTING (use - sign)
            boostcosth = -boostvector(3)/boostp
            boostsinth = sqrt(1.0_r14-boostcosth**2)
            boostcosphi   = -boostvector(1)/(boostp*boostsinth)
            boostsinphi   = -boostvector(2)/(boostp*boostsinth)


   
            ! Lorentz boost matrix elements for a boost in the direction of (theta, phi)
            !  lormtrx = 
            !  | gamma    -gamma*beta*v_x  -gamma*beta*v_y  -gamma*beta*v_z  |
            !  | -gamma*beta*v_x ...                                         |
            !  | -gamma*beta*v_y ...             (spatial part)              |
            !  | -gamma*beta*v_z ...                                         |

            ! Time-space components
            lormtrx(0,0) = gamma
            lormtrx(0,1) = -gamma * beta * boostsinth * boostcosphi
            lormtrx(0,2) = -gamma * beta * boostsinth * boostsinphi
            lormtrx(0,3) = -gamma * beta * boostcosth

            ! Space-time components
            lormtrx(1,0) = -gamma * beta * boostsinth * boostcosphi
            lormtrx(2,0) = -gamma * beta * boostsinth * boostsinphi
            lormtrx(3,0) = -gamma * beta * boostcosth

            ! Spatial-spatial components
            lormtrx(1,1) = 1.0_r14+ (gamma-1.)*(boostsinth**2 * boostcosphi**2)
            lormtrx(1,2) = (gamma-1.)*(boostsinth**2 * boostcosphi * boostsinphi)
            lormtrx(1,3) = (gamma-1.)*(boostsinth * boostcosth * boostcosphi)

            lormtrx(2,1) = (gamma-1.)*(boostsinth**2 * boostcosphi * boostsinphi)
            lormtrx(2,2) = 1.0_r14+ (gamma-1.)*(boostsinth**2 * boostsinphi**2)
            lormtrx(2,3) = (gamma-1.)*(boostsinth * boostcosth * boostsinphi)

            lormtrx(3,1) = (gamma-1.)*(boostsinth * boostcosth * boostcosphi)
            lormtrx(3,2) = (gamma-1.)*(boostsinth * boostcosth * boostsinphi)
            lormtrx(3,3) = 1.0_r14+ (gamma-1.)*(boostcosth**2)

            do k = 0, 3
                deboosted(ii,k) = 0.
                do j = 0, 3
                    deboosted(ii,k) = deboosted(ii,k)+lormtrx(k,j)*boosted(j)
                enddo
            enddo
            
            
        end subroutine deboost

        real(r14) function metric(mu,nu)
            integer :: mu,nu
            if (mu /= nu) then
                metric = 0.0_r14
            else if (mu == 0) then
                metric = 1.0_r14
            else
                metric = -1.0_r14
            endif


        end function metric




end module utilities



module histogram

    use numeric
    use inputs
    implicit none
    real(r14), allocatable :: h_pmod1(:,:), h_emu1(:,:), h_thmu1(:,:), h_pmod2(:,:), &
    h_thmu2(:,:), h_emu2(:,:),  h_pgam(:,:), h_thgam(:,:), h_qq(:,:)
    real(r14) :: enemu_max, pmodmu_max, pgam_max
    integer :: allerr
    contains
        subroutine inithistsborn()

            enemu_max = cme/2
            pmodmu_max = sqrt(enemu_max**2-mmu**2)
            allocate(h_pmod1(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_emu1(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu1(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu- histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif

            h_pmod1 = 0.0_r14
            h_emu1 = 0.0_r14
            h_thmu1 = 0.0_r14
            
        end subroutine inithistsborn

        subroutine inithistsisr()

            enemu_max = cme/2
            pmodmu_max = sqrt(enemu_max**2-mmu**2)
            pgam_max = cme
            allocate(h_pmod1(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_emu1(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu1(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu- histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_pmod2(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu+ p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_emu2(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu+ energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu2(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu+ histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_qq(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Muons invariant mass histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_pgam(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "ISR gamma energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thgam(0:nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "ISR gamma theta histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif

            h_pmod1 = 0.0_r14
            h_emu1 = 0.0_r14
            h_thmu1 = 0.0_r14
            h_pmod2 = 0.0_r14
            h_emu2 = 0.0_r14
            h_thmu2 = 0.0_r14
            h_qq = 0.0_r14
            h_pgam = 0.0_r14
            h_thgam = 0.0_r14
            
        end subroutine inithistsisr

        subroutine updatehist(hist,histmin,histmax,val,weight)
            real(r14), intent(inout) :: hist(0:,:) 
            real(r14), intent(in) :: histmin,histmax,val,weight
            integer :: bin

            if (val >= histmax) then
                hist(nbins+1,1)= hist(nbins+1,1) + weight
                hist(nbins+1,2)= hist(nbins+1,2) + weight**2
            else if (val< histmin) then
                hist(0,1) = hist(0,1) + weight
                hist(0,2) = hist(0,2) + weight**2
            else
                bin = int((val-histmin)/(histmax-histmin)*nbins) +1
                hist(bin,1)= hist(bin,1) + weight
                hist(bin,2)= hist(bin,2) + weight**2
            endif 
            
            
        end subroutine updatehist

        subroutine endhistisr(hist,histmin,histmax)
            real(r14), intent(inout) :: hist(0:,:) 
            real(r14), intent(in) :: histmin, histmax
            integer :: bb

            do bb = 0, nbins+1
                hist(bb,2) = sqrt(hist(bb,2))
                hist(bb,1) = intemax/dble(ngen)*hist(bb,1)
                hist(bb,2) = intemax/dble(ngen)*hist(bb,2)
                hist(bb,1)=hist(bb,1)*dble(nbins)/(histmax-histmin)
                hist(bb,2)=hist(bb,2)*dble(nbins)/(histmax-histmin)
            enddo
        end subroutine endhistisr

        subroutine endhistborn(hist,histmin,histmax)
            real(r14), intent(inout) :: hist(0:,:) 
            real(r14), intent(in) :: histmin, histmax
            integer :: bb

            do bb = 0, nbins+1
                hist(bb,2) = sqrt(hist(bb,2))
                hist(bb,1)=hist(bb,1)*bornsigma/dble(ngen)
                hist(bb,2)=hist(bb,2)*bornsigma/dble(ngen)
                hist(bb,2)=hist(bb,2)*dble(nbins)/(histmax-histmin)
                hist(bb,1)=hist(bb,1)*dble(nbins)/(histmax-histmin)
                hist(bb,2)=hist(bb,2)*dble(nbins)/(histmax-histmin)
            enddo

        end subroutine endhistborn

        subroutine endbornhistos()

            call endhistborn(h_pmod1, 0.0_r14, pmodmu_max)
            call endhistborn(h_emu1, 0.0_r14, enemu_max)
            call endhistborn(h_thmu1, 0.0_r14, 180.0_r14)
            
        end subroutine endbornhistos

        subroutine endisrhistos()

            call endhistisr(h_pmod1, 0.0_r14, pmodmu_max)
            call endhistisr(h_emu1, 0.0_r14, enemu_max)
            call endhistisr(h_thmu1, 0.0_r14, 180.0_r14)
            call endhistisr(h_pmod2, 0.0_r14, pmodmu_max)
            call endhistisr(h_emu2, 0.0_r14, enemu_max)
            call endhistisr(h_thmu2, 0.0_r14, 180.0_r14)
            call endhistisr(h_pgam, 0.0_r14, pgam_max)
            call endhistisr(h_thgam, 0.0_r14, 180.0_r14)
            call endhistisr(h_qq, 0.0_r14, qqmax)
            
        end subroutine endisrhistos
    

end module histogram

module eventgen

    use numeric
    use inputs
    use histogram
    use utilities
    implicit none
    real(r14) :: inte, tmpintemin = 0.0_r14, tmpintemax = 0.0_r14, sinthmu, phimu, cosphimu, &
    sinphimu, sinthgam, Muontensor(0:3,0:3), Electrontensor(0:3,0:3), invampl, naccpt = 0.0_r14
    real(r14) :: rndm(7), jacqq, jacgam, jacmuang, pgamvirt(0:3), z, ppos(0:3), pel(0:3), &
    eff, sigma, dsigma 
    real(r14), allocatable :: pmu1(:,:), pmu2(:,:), pmod1(:), pmod2(:), pgam(:,:), costhmu1(:), &
    costhmu2(:), costhgam(:), qq(:), wght(:)  ! Reduce to fewer higher dimensional arrays 
    integer(i10) :: iev, arraylen
    integer(i10), allocatable :: accpt(:)
    integer :: run, err
    logical :: accepted, warnflag = .false.



    contains
        subroutine allocarrays
            arraylen = max(ngen,nmax)

            allocate(pmu1(arraylen,0:3), stat=err) !!! Maybe only allocate the big event arrays if evsave is called
            if (err /= 0) then
                print *, "pmu1 array allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(pmod1(arraylen), stat=err)
            if (err /= 0) then
                print *, "pmod1 array allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(costhmu1(arraylen), stat=err)
            if (err /= 0) then
                print *, "costhmu1 array allocation request denied, aborting!"
                print*, ""
                stop
            endif

            if (isr .eqv. .false.) then ! BORN CASE

                allocate(accpt(arraylen), stat=err)
                if (err /= 0) then
                    print *, "accpt array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif        

            else ! ISR CASE

                allocate(pmu2(arraylen,0:3), stat=err) 
                if (err /= 0) then
                    print *, "pmu1 array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif

                allocate(pmod2(arraylen), stat=err)
                if (err /= 0) then
                    print *, "pmod1 array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif
                allocate(costhmu2(arraylen), stat=err)
                if (err /= 0) then
                    print *, "costhmu1 array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif
                allocate(pgam(arraylen,0:3), stat=err)
                if (err /= 0) then
                    print *, "pgam array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif
                allocate(costhgam(arraylen), stat=err)
                if (err /= 0) then
                    print *, "costhgam array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif
                allocate(qq(arraylen), stat=err)
                if (err /= 0) then
                    print *, "qq array allocation request denied, aborting!"
                    print*, ""
                    stop
                endif

                if (wghtopt .eqv. .false.) then
                
                    allocate(accpt(arraylen), stat=err)
                    if (err /= 0) then
                        print *, "accpt array allocation request denied, aborting!"
                        print*, ""
                        stop
                    endif
                else
                    allocate(wght(arraylen), stat=err)
                    if (err /= 0) then
                        print *, "wght array allocation request denied, aborting!"
                        print*, ""
                        stop
                    endif
                endif
                
            endif


            
        end subroutine allocarrays

        subroutine genborn()

            call bornmuangles()
            
            pmu1(iev,0) = cme/2

            pmod1(iev) = sqrt( pmu1(iev,0)**2 - mmu**2 )

            pmu1(iev,1) = pmod1(iev) * sinthmu * cosphimu
            pmu1(iev,2) = pmod1(iev) * sinthmu * sinphimu
            pmu1(iev,3) = pmod1(iev) * costhmu1(iev)

            call testcuts()
            
            if (accepted) then

                accpt(iev) = 1
                naccpt = naccpt + 1.0_r14 

                if (histsave /= '') then
            
                    call updatehist(h_pmod1, 0.0_r14, pmodmu_max, pmod1(iev),1.0_r14)
                    call updatehist(h_emu1, 0.0_r14, enemu_max, pmu1(iev,0),1.0_r14)
                    call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)),1.0_r14)
            
                endif

            else
                accpt(iev) = 0
            endif

            
        end subroutine genborn

        subroutine bornmuangles()
            costhmu1(iev) = inverseCDF(rndm(1))
            sinthmu = sqrt(1.0_r14-costhmu1(iev)**2)
            phimu = 2*pi*rndm(2)
            cosphimu = cos(phimu)
            sinphimu = sin(phimu)
            
        end subroutine bornmuangles

        subroutine genisr()

            call elepos()
            call isrmuqq()
            call isrgammas()
            call isrmumom()
            call buildampl()

            call testcuts()

            if (accepted) then

                inte = gev2nbarn*jacqq*jacgam*jacmuang/(4.0_r14*pi*sinv)*invampl

                if (run == 1) then
                    if(inte > tmpintemax) tmpintemax = inte
                    if(inte < tmpintemin) tmpintemin = inte
                endif

                if (run == 2) then
                    
                    if (inte > intemax .and. (warnflag .eqv. .false.)) then
                        print*, "WARNING : Found integrand greater than estimated maximum!&
                        & Try increasing 'nmax'..."
                        print*, ""
                        warnflag = .true.
                    endif
                    if (wghtopt .eqv. .false.) then

                        z = rndm(7)
                        if ( z <= inte/intemax ) then
                            accpt(iev) = 1
                            naccpt = naccpt + 1.

                            if (histsave /= '') then
                        
                            call updatehist(h_pmod1, 0.0_r14, pmodmu_max, pmod1(iev),1.0_r14)
                            call updatehist(h_emu1, 0.0_r14, enemu_max, pmu1(iev,0),1.0_r14)
                            call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)),1.0_r14)
                            call updatehist(h_pmod2, 0.0_r14, pmodmu_max, pmod2(iev),1.0_r14)
                            call updatehist(h_emu2, 0.0_r14, enemu_max, pmu2(iev,0),1.0_r14)
                            call updatehist(h_thmu2, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu2(iev)),1.0_r14)
                            call updatehist(h_pgam, 0.0_r14, pgam_max, pgam(iev,0),1.0_r14)
                            call updatehist(h_thgam, 0.0_r14, 180.0_r14, radtodeg*acos(costhgam(iev)),1.0_r14)
                            call updatehist(h_qq, 0.0_r14, qqmax, qq(iev),1.0_r14)
                        
                            endif
                        
                        end if
                    endif
                    if (wghtopt .eqv. .true.) then ! WEIGHTED GENERATION
                        wght(iev) = inte/intemax
                        naccpt = naccpt + wght(iev) 

                        if (histsave /= '') then
                    
                        call updatehist(h_pmod1, 0.0_r14, pmodmu_max, pmod1(iev),wght(iev))
                        call updatehist(h_emu1, 0.0_r14, enemu_max, pmu1(iev,0),wght(iev))
                        call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)),wght(iev))
                        call updatehist(h_pmod2, 0.0_r14, pmodmu_max, pmod2(iev),wght(iev))
                        call updatehist(h_emu2, 0.0_r14, enemu_max, pmu2(iev,0),wght(iev))
                        call updatehist(h_thmu2, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu2(iev)),wght(iev))
                        call updatehist(h_pgam, 0.0_r14, pgam_max, pgam(iev,0),wght(iev))
                        call updatehist(h_thgam, 0.0_r14, 180.0_r14, radtodeg*acos(costhgam(iev)),wght(iev))
                        call updatehist(h_qq, 0.0_r14, qqmax, qq(iev),wght(iev))
                    
                        endif
                    endif
                endif       
            else
                if (wghtopt .eqv. .false.) accpt(iev) = 0
                if (wghtopt .eqv. .true.) wght(iev) = -999
            endif            
            
        end subroutine genisr

        subroutine elepos()

            pel(0) = cme/2.
            pel(1) = 0.
            pel(2) = 0.
            pel(3) = sqrt(pel(0)**2-me**2)

            ppos(0) = cme/2.
            ppos(1) = 0.
            ppos(2) = 0.
            ppos(3) = -sqrt(pel(0)**2-me**2)
            
        end subroutine elepos

        subroutine isrmuqq()

            real(r14) :: x, fak1, amin, amax, a, bmin, b, p, ppp, y
            x = rndm(1)
            fak1 = -1.0_r14/sinv
            amin = fak1*log(sinv-qqmin)
            amax = fak1*log(sinv-qqmax)
            a = amax-amin
            bmin = log(qqmin/sinv)/sinv
            b    = log(qqmax/qqmin)/sinv
            p = rndm(2)   
            ppp  = a/(a+b)
            if (p < ppp) then
                y  = amin+a*x
                qq(iev) = sinv-exp(y/fak1)                                       
            else
                y  = bmin+b*x
                qq(iev) = sinv*exp(sinv*y)

            endif
            jacqq = (a+b)/(1.0_r14/(sinv*(sinv-qq(iev))) + 1.0_r14/sinv/qq(iev))

            !print*, " a = ", a, " b = ", b, " qq = ", qq(iev), " y = ", y, " amax = ", amax, " ppp = ", ppp

            
        end subroutine isrmuqq

        subroutine isrgammas()

            real(r14) :: x, phigam, b, cmin, cmax, y, sinthgam
            integer :: ii
            x = rndm(3)
            phigam = 2.0_r14*pi*rndm(4)
            b = sqrt(1.0_r14 - 4.0_r14*me**2/sinv)
            cmin = log((1.0_r14+b*costhgammin)/(1.0_r14-b*costhgammin))/(2.0_r14*b)
            cmax = log((1.0_r14+b*costhgammax)/(1.0_r14-b*costhgammax))/(2.0_r14*b)
            y = cmin+x*(cmax-cmin)
            costhgam(iev) = tanh(b*y)/b
            sinthgam = sqrt(1.0_r14-costhgam(iev)**2)
            jacgam = 2.0_r14*pi*(1.0_r14- b**2 * costhgam(iev)**2)*(cmax-cmin)

            
            pgam(iev,0) = cme/2*(1-qq(iev)/sinv)
            pgam(iev,1) = pgam(iev,0) * sinthgam * cos(phigam)
            pgam(iev,2) = pgam(iev,0) * sinthgam * sin(phigam)
            pgam(iev,3) = pgam(iev,0) * costhgam(iev)

            pgamvirt(0) = cme/2*(1+qq(iev)/sinv)
            do ii = 1,3
                pgamvirt(ii) = -pgam(iev,ii)
            enddo

            
        end subroutine isrgammas

        subroutine isrmumom()
            
            real(r14) :: x, costhmucm, sinthmucm, phimucm, pmodcm, tmpvect(0:3)
            integer :: ii
            x = rndm(5)
            phimucm = 2.0_r14*pi*rndm(6)
            costhmucm = costhmumin+(costhmumax-costhmumin)*x
            jacmuang = 2.0_r14*pi*(costhmumax-costhmumin)
            sinthmucm = sqrt(1.0_r14-costhmucm*costhmucm)
            pmu1(iev,0) = sqrt(qq(iev))/2
            pmodcm = sqrt( pmu1(iev,0)**2 - mmu**2 )
            pmu1(iev,1) = pmodcm * sinthmucm * cos(phimucm)
            pmu1(iev,2) = pmodcm * sinthmucm * sin(phimucm)
            pmu1(iev,3) = pmodcm * costhmucm
            pmu2(iev,0) = pmu1(iev,0) 
            do ii = 1,3
                pmu2(iev,ii) = -pmu1(iev,ii)
            enddo
            do ii = 0,3
                tmpvect(ii) = pmu1(iev,ii)
            enddo
            call deboost(tmpvect,pgamvirt,pmu1,iev)

            do ii = 0,3
                tmpvect(ii) = pmu2(iev,ii)
            enddo
            call deboost(tmpvect,pgamvirt,pmu2,iev)

            pmod1(iev) = sqrt( pmu1(iev,0)**2 - mmu**2 )
            pmod2(iev) = sqrt( pmu2(iev,0)**2 - mmu**2 )
            costhmu1(iev) = pmu1(iev,3)/pmod1(iev)
            costhmu2(iev) = pmu2(iev,3)/pmod2(iev)


        end subroutine isrmumom

        subroutine maketensors
            real(r14) :: dps, a00, a11, a12, a22, m2, q2, uq2, b, x, y1, y2, globfact
            integer :: mu,nu
                
            dps =  sqrt(1.0_r14-4.0_r14*mmu**2/qq(iev))/(32.0_r14*pi**2)  ! Phase space factors
            do mu = 0,3
                do nu = 0,3
                    Muontensor(mu,nu) = 16.0_r14*pi*alpha*(pmu2(iev,mu)*pmu1(iev,nu)+&
                    pmu1(iev,mu)*pmu2(iev,nu)-qq(iev)/2.0_r14*metric(mu,nu))*dps
                enddo
            enddo

            m2 = me*me/sinv
            q2 = qq(iev)/sinv
            uq2 = 1.0_r14-q2
            b = sqrt(1.0_r14-4.0_r14*m2)
            x = b*costhgam(iev)
            y1 = uq2*(1.0_r14-x)/2.
            y2 = uq2*(1.0_r14+x)/2.
            globfact = (4.0_r14*pi*alpha/sinv)**2/(q2**2)


            a00 = globfact*( 2.0_r14*m2*Q2*uq2*uq2/(y1*y2)- (2.0_r14*q2+y1*y1+y2*y2) )/(y1*y2)
            a11 = globfact*(8.0_r14*m2/y2-4.0_r14*q2/y1)/y2
            a22 = globfact*(8.0_r14*m2/y1-4.0_r14*q2/y2)/y1
            a12 = -globfact*8.0_r14*m2/(y1*y2)

            dps = (1.0_r14-qq(iev)/sinv)/(32.0_r14*pi**2)  ! Phase space factor

            do mu = 0,3
                do nu = 0,3
                    Electrontensor(mu,nu) = (a00*metric(mu,nu)+ a11*pel(mu)*pel(nu)/sinv+&
                    a22*ppos(mu)*ppos(nu)/sinv+a12*(ppos(mu)*pel(nu)+&
                    pel(mu)*ppos(nu))/sinv)*dps
                enddo
            enddo
            
        end subroutine maketensors

        subroutine buildampl

            integer :: mu, nu
            real(r14) :: metric1, metric2

            call maketensors()

            invampl = 0.0_r14
            do mu = 0,3
                metric1 = 1.0_r14
                if (mu.eq.0) metric1 = -1.0_r14
                do nu = 0,3
                    metric2 = 1.0_r14
                if (nu.eq.0) metric2 = -1.0_r14               
                invampl = invampl + metric1*metric2*Electrontensor(mu,nu)*Muontensor(mu,nu) 
                enddo
            enddo
            
        end subroutine buildampl

        subroutine testcuts()

            logical :: log1, log2, log3, log4

            if(isr .eqv. .false.) then

                if ( costhmu1(iev) <= cos(thmucutmin * degtorad) .and. costhmu1(iev) >= cos(thmucutmax * degtorad)) then
                    accepted = .true.
                else
                    accepted = .false.
                end if

            endif

            if(isr .eqv. .true.) then
                
                log1 = costhmu1(iev) <= cos(thmucutmin * degtorad) .and. costhmu1(iev) >= cos(thmucutmax * degtorad)
                log2 = qqcutmin <= qq(iev) .and. qq(iev) <= qqcutmax
                log3 = pgam(iev,0) >= gmin
                log4 = costhgam(iev) <= cos(thgamcutmin * degtorad) .and. costhgam(iev) >= cos(thgamcutmax * degtorad)

                if ( log1 .and. log2 .and. log3 .and. log4 ) then
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
    integer(i10) :: idx
    integer :: ou, hu

    contains
        subroutine writevents
            inquire(file=evsave, exist=exists)
            if(exists) then
                print *, "Writing events to file ", trim(evsave), " (overwriting)..."
                print*, ""
            else
                print *, "Writing events to file", trim(evsave)
                print*, ""
            end if
            open(newunit=ou, file=evsave, iostat=ios)
            if (ios == 0) then
                if ( isr .eqv. .false. ) then
                    write(ou,*) "==========================================&
                    &======================================"
                    write(ou,*) "                              TOYGENERATOR EVENTS&
                    &                                       "
                    write(ou,*) "==========================================&
                    &======================================"
                    write(ou,*) ""
                    write(ou,'(A, i0, A, f0.2, A)') "   Born generation using seed = ", seed(1), " at center of mass&
                    & energy = ", cme, " GeV"
                    write(ou,*) ""
                    write(ou, '(*(A6, 3x))') 'px-',   'py-',  'pz-',    'E-',    'pmod-',    'th-',    'px+',   &
                    'py+', 'pz+', 'E+',  'pmod+',  'th+', 'qq'
                    do idx = 1, ngen
                        if (accpt(idx) == 1) then
                            write (ou,'(*(f6.2, 3x))') pmu1(idx,1), pmu1(idx,2), pmu1(idx,3), pmu1(idx,0), pmod1(idx), &
                            radtodeg*acos(costhmu1(idx)), -pmu1(idx,1), -pmu1(idx,2), -pmu1(idx,3), &
                            pmu1(idx,0), pmod1(idx), 180-radtodeg*acos(costhmu1(idx)), sinv
                        endif
                    end do
                end if
                if ( isr .eqv. .true. ) then
                    if ( wghtopt .eqv. .false. ) then
                        write(ou,*) "==========================================&
                        &======================================"
                        write(ou,*) "                              TOYGENERATOR EVENTS&
                        &                                       "
                        write(ou,*) "==========================================&
                        &======================================"
                        write(ou,*) ""
                        write(ou,'(A, i0, A, f0.2, A)') "   ISR generation using seed = ", seed(1), " at center of mass&
                        & energy = ", cme, " GeV"
                        write(ou,*) ""
                        write(ou, '(*(A6, 3x))') 'px-',   'py-',  'pz-',    'E-',    'pmod-',    'th-',    'px+',   &
                        'py+', 'pz+', 'E+',  'pmod+',  'th+', 'qq', 'pgamx', 'pgamy', 'pgamz', 'Egam', 'thgam'
                        do idx = 1, ngen
                            if (accpt(idx) == 1) then
                                write (ou,'(*(f6.2, 3x))') pmu1(idx,1), pmu1(idx,2), pmu1(idx,3), pmu1(idx,0), pmod1(idx), &
                                radtodeg*acos(costhmu1(idx)), pmu2(idx,1), pmu2(idx,2), pmu2(idx,3), &
                                pmu2(idx,0), pmod2(idx), radtodeg*acos(costhmu2(idx)), qq(idx), pgam(idx,1)&
                                , pgam(idx,2), pgam(idx,3), pgam(idx,0), radtodeg*acos(costhgam(idx))
                            endif
                        end do
                    end if
                    if ( wghtopt .eqv. .true. ) then
                        write(ou,*) "==========================================&
                        &======================================"
                        write(ou,*) "                              TOYGENERATOR EVENTS&
                        &                                       "
                        write(ou,*) "==========================================&
                        &======================================"
                        write(ou,*) ""
                        write(ou,'(A, i0, A, f0.2, A)') "   Weighted ISR generation using seed = ", seed(1), " at center of mass&
                        & energy = ", cme, " GeV"
                        write(ou,*) ""
                        write(ou, '(*(A6, 3x))') 'px-',   'py-',  'pz-',    'E-',    'pmod-',    'th-',    'px+',   &
                        'py+', 'pz+', 'E+',  'pmod+',  'th+', 'qq', 'pgamx', 'pgamy', 'pgamz', 'Egam', 'thgam', 'wght'
                        do idx = 1, ngen
                            if (wght(idx) >= 0) then
                                write (ou,'(*(f6.2, 3x))') pmu1(idx,1), pmu1(idx,2), pmu1(idx,3), pmu1(idx,0), pmod1(idx), &
                                radtodeg*acos(costhmu1(idx)), pmu2(idx,1), pmu2(idx,2), pmu2(idx,3), &
                                pmu2(idx,0), pmod2(idx), radtodeg*acos(costhmu2(idx)), qq(idx), pgam(idx,1)&
                                , pgam(idx,2), pgam(idx,3), pgam(idx,0), radtodeg*acos(costhgam(idx)), wght(idx)
                            endif
                        end do
                    end if
                end if
                close(unit=ou)
            else
                print *, 'ERROR while opening file ', trim(evsave)
            end if
            
        end subroutine writevents

        subroutine appendhist(outunit, hist, histname, histmin, histmax)
            real(r14), intent(inout) :: hist(0:,:) 
            real(r14), intent(in) :: histmin,histmax
            real(r14) :: xlow, xhigh
            character(len=*) :: histname
            integer :: bb, outunit

            write(outunit, *) histname
            write(outunit, '(*(A12, 4x))') "binlow", "binhigh", "content", "error"
            write (outunit,'(*(f12.3, 4x))') histmin, histmin, hist(0,1), hist(0,2)
            do bb = 1, nbins
                xlow = histmin + (bb-1)* (histmax-histmin)/nbins
                xhigh = xlow + (histmax-histmin)/nbins
                write (outunit,'(*(f12.3, 4x))') xlow, xhigh, hist(bb,1), hist(bb,2)
            enddo          
            write (outunit,'(*(f12.3, 4x))') histmax, histmax, hist(nbins+1,1), hist(nbins+1,2)
            
        end subroutine appendhist

        subroutine writehists

            inquire(file=histsave, exist=exists)
            if(exists) then
                print *, "Writing histograms to file ", trim(histsave), " (overwriting)..."
                print*, ""
            else
                print *, "Writing histograms to file", trim(histsave)
                print*, ""
            end if

            open(newunit=hu, file=histsave, iostat=ios)    
            if (ios == 0) then
                if ( isr .eqv. .false. ) then
                    call endbornhistos()
                    write(hu,*) "==========================================&
                    &======================================"
                    write(hu,*) "                              TOYGENERATOR HISTOGRAMS&
                    &                                       "
                    write(hu,*) "==========================================&
                    &======================================"
                    write(hu,*) ""
                    write(hu,'(A, i0, A, f0.2, A)') "   Born generation using seed = ", seed(1), " at center of mass&
                    & energy = ", cme, " GeV"
                    write(hu,*) ""
                    call appendhist(hu, h_pmod1, "h_pmod1", 0.0_r14, pmodmu_max)
                    call appendhist(hu, h_emu1, "h_emu1", 0.0_r14, enemu_max)
                    call appendhist(hu, h_thmu1, "h_thmu1", 0.0_r14, 180.0_r14)
                end if
                if ( isr .eqv. .true. ) then
                    call endisrhistos()
                    write(hu,*) "==========================================&
                    &======================================"
                    write(hu,*) "                              TOYGENERATOR HISTOGRAMS&
                    &                                       "
                    write(hu,*) "==========================================&
                    &======================================"
                    write(hu,*) ""
                    write(hu,'(A, i0, A, f0.2, A)') "   ISR generation using seed = ", seed(1), " at center of mass&
                    & energy = ", cme, " GeV"
                    write(hu,*) ""
                    call appendhist(hu, h_pmod1, "h_pmod1", 0.0_r14, pmodmu_max)
                    call appendhist(hu, h_emu1, "h_emu1", 0.0_r14, enemu_max)
                    call appendhist(hu, h_thmu1, "h_thmu1", 0.0_r14, 180.0_r14)
                    call appendhist(hu, h_pmod2, "h_pmod2", 0.0_r14, pmodmu_max)
                    call appendhist(hu, h_emu2, "h_emu2", 0.0_r14, enemu_max)
                    call appendhist(hu, h_thmu2, "h_thmu2", 0.0_r14, 180.0_r14)
                    call appendhist(hu, h_pgam, "h_pgam", 0.0_r14, pgam_max)
                    call appendhist(hu, h_thgam, "h_thgam", 0.0_r14, 180.0_r14)
                    call appendhist(hu, h_qq, "h_qq", 0.0_r14, qqmax)
                    
                end if
                close(unit=hu)
            else
                print *, 'ERROR while opening file ', evsave
            end if
            
        end subroutine writehists

        subroutine finalout

            inquire(file="info.out", exist=exists)
            if(exists) then
                print *, "Writing this generation output to file 'info.out' (overwriting)..."
                print*, ""
            else
                print *, "Writing this generation output to file 'info.out'..."
                print*, ""
            end if
            open(newunit=ou, file="info.out", iostat=ios)
            if (ios == 0) then
                if ( isr .eqv. .false. ) then ! BORN CASE
                    write(ou,*) "==========================================&
                    &======================================"
                    write(ou,*) "                              TOYGENERATOR OUTPUT&
                    &                                       "
                    write(ou,*) "==========================================&
                    &======================================"
                    write(ou,*) ""
                    write(ou,'(A, i0)') "* Born generation (LO) using seed = ", seed(1)
                    write(ou,*) ""
                    write(ou,'(A, f0.2)') "* Center of mass energy = ", cme, " GeV"
                    write(ou,*) ""
                    if ( compareals(thmucutmin,0.0_r14) .and. compareals(thmucutmax,180.0_r14) ) then
                        write(ou,'(A)') "* No angular cuts on muons"
                    else 
                        write(ou,'(A, f0.2, A, f0.2)') "* Angular cuts on muons :&
                        & theta min = ", thmucutmin, " theta max =  ", thmucutmax
                    end if 
                    write(ou,*) ""
                    if ( compareals(qqcutmin,0.0_r14) .and. compareals(qqcutmax,sinv) ) then
                        write(ou,'(A)') "* No cuts on muon system invariant mass"
                    else 
                        write(ou,'(A, f0.2, A, f0.2)') "* Angular cuts on muons' &
                        &invariant mass: qq min = ", qqcutmin, " qq max =  ", qqcutmax
                    end if
                    write(ou,*) "" 
                    write(ou,'(A, i0, A, i0, A, f5.3)') "* ", int(naccpt), " events accepted of ", ngen, " events generated&
                    & ====> Generation efficiency = ", eff 
                    write(ou,*) ""
                    write(ou,'(A, f0.2, A)') "* Total Born cross section without cuts = &
                    &", bornsigma, " nb"
                    write(ou,*) ""
                    write(ou,'(A, f0.2, A, f0.2, A)') "* Total Born cross section with&
                    & selected cuts= ", sigma, " +- ", dsigma, " nb"
                    write(ou,*) ""
                    write(ou,*) "==========================================&
                    &======================================"
                end if
                if ( isr .eqv. .true. ) then
                    write(ou,*) "==========================================&
                    &======================================"
                    write(ou,*) "                              TOYGENERATOR OUTPUT&
                    &                                       "
                    write(ou,*) "==========================================&
                    &======================================"
                    write(ou,*) ""
                    write(ou,'(A, i0)') "* ISR generation (NLO) using seed = ", seed(1)
                    write(ou,*) ""
                    if ( wghtopt .eqv. .false. ) then
                        write(ou,'(A)') "* Weighted generation : OFF " 
                        write(ou,*) ""
                    end if
                    if ( wghtopt .eqv. .true. ) then
                        write(ou,'(A)') "* Weighted generation : ON " 
                        write(ou,*) ""
                    end if
                    write(ou,'(A, f0.2, A)') "* Center of mass energy = ", cme, " GeV"
                    write(ou,*) ""
                    if ( compareals(thmucutmin,0.0_r14) .and. compareals(thmucutmax,180.0_r14) ) then
                        write(ou,'(A)') "* No angular cuts on muons"
                    else 
                        write(ou,'(A, f0.2, A, f0.2)') "* Angular cuts on muons :&
                        & theta min = ", thmucutmin, " theta max =  ", thmucutmax
                    end if 
                    write(ou,*) ""
                    if ( compareals(thgamcutmin,0.0_r14) .and. compareals(thgamcutmax,180.0_r14) ) then
                        write(ou,'(A)') "* No angular cuts on ISR photon"
                    else 
                        write(ou,'(A, f0.2, A, f0.2)') "* Angular cuts on ISR photon :&
                        & theta min = ", thgamcutmin, " theta max =  ", thgamcutmax
                    end if 
                    write(ou,*) ""
                    if ( compareals(qqcutmin,0.0_r14) .and. compareals(qqcutmax,sinv) ) then
                        write(ou,'(A)') "* No cuts on muon system invariant mass"
                    else 
                        write(ou,'(A, f0.2, A, f0.2)') "* Angular cuts on muons' &
                        &invariant mass: qq min = ", qqcutmin, " qq max =  ", qqcutmax
                    end if
                    write(ou,*) ""
                    write(ou,'(A, f0.2, A)') "* Minimum ISR photon energy = ", gmin, " GeV"
                    write(ou,*) "" 
                    write(ou,'(A, i0, A, i0, A, f5.3)') "* ", int(naccpt), " events accepted of ", ngen, " events generated&
                    & ====> Generation efficiency = ", eff 
                    write(ou,*) ""
                    write(ou,'(A, f0.2, A, f0.2, A)') "* Total ISR cross section with selected cuts= ", sigma, " +- ", dsigma, " nb"
                    write(ou,*) ""
                    write(ou,*) "==========================================&
                    &======================================"
                    
                end if
                close(unit=ou)
            else
                print *, 'ERROR while opening file ', evsave
            end if

            print*, "Finished writing output!"
            
        end subroutine finalout

end module output


program toygenerator
    use numeric
    use inputs
    use histogram
    use eventgen
    use output
    implicit none

    call loadinput()

    call random_seed(put=seed)

    call allocarrays()

    if (isr .eqv. .false.) then ! BORN CASE

        if (histsave /= '') call inithistsborn()    

        do iev = 1, ngen 

            call random_number(rndm)
            call genborn()

        end do
        
        if (evsave /= '') call writevents()
        if (histsave /= '') call writehists()

        bornsigma = gev2nbarn*4.0_r14*pi/3.0_r14*alpha**2/sinv
        eff = dble(naccpt)/ngen
        sigma = bornsigma * eff
        dsigma = bornsigma/ngen * sqrt(ngen * eff * (1-eff))
        


    else ! ISR CASE
        
        if (histsave /= '') call inithistsisr()
        run = 1

        do iev = 1, nmax 
             call random_number(rndm)
             call genisr()
            
        end do

        run = 2
        intemax = tmpintemax

        do iev = 1, ngen

            call random_number(rndm)
            call genisr()
            
        end do
        if (evsave /= '') call writevents()
        if (histsave /= '') call writehists()

        !sigma = integral = fraction of accepted (naccpt/ngen) * max (intemax) * range of random variables (1)  

        eff = naccpt/ngen
        sigma = intemax * eff
        if ( wghtopt .eqv. .true. ) then
            dsigma = 0
            do iev = 1, ngen
                dsigma = dsigma + intemax**2/(ngen*(ngen-1))*(wght(iev)-naccpt/ngen)**2
            end do
            dsigma = sqrt(dsigma)
        else 
            dsigma = intemax/ngen * sqrt(ngen * eff * (1-eff))
        end if

    endif

    call finalout()

    print*, ""
    print*, "***********************************************************************************"
    print*, "***********************************************************************************"
    print*, ""

end program toygenerator