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
    real(r14) :: cme = -999, thmucutmin = 0, thmucutmax = 180, qqcutmin = 0, qqcutmax = -999, &
    sinv, gmin = 0, qqmax, thgamcutmin = 0, thgamcutmax = 180, costhgammin, &
    costhgammax, costhmumin, costhmumax

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
                    else if (opt == 'thgamcutmin') then
                        read(optval, *, iostat=ios) thgamcutmin
                        if (ios /= 0) then
                            print *, "Input error! Value given for 'thgamcutmin' is not valid! Aborting..."
                            print *, ""
                            stop
                        endif
                        if (thgamcutmin < 0 .or. thgamcutmin > 180) then
                            print*, "Input error! Polar angle cut minumum for ISR photon 'thgamcutmin'&
                             (",thgamcutmin, ") MUST be a real number between 0 and 180! Aborting..."
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
                        if (thgamcutmax < 0 .or. thgamcutmax > 180) then
                            print*, "Input error! Polar angle cut maximum for ISR photon 'thgamcutmax'&
                             (",thgamcutmax, ") MUST be a real number between 0 and 180! Aborting..."
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
            costhgammin = cos(radtodeg*thgamcutmax)
            costhgammax = cos(radtodeg*thgamcutmin)
            costhmumin = cos(radtodeg*thmucutmax)
            costhmumax = cos(radtodeg*thmucutmin)

            if (qqcutmax == -999) qqcutmax = sinv !!! Numerical problem?

            if (thmucutmin >= thmucutmax) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'thmucutmin' value (", thmucutmin, ")&
                must be lower than 'thmucutmax' value (", thmucutmin, ")"
                print *, ""
                stop    
            endif

            if (thgamcutmin >= thgamcutmax) then
                print '(A, f0.2, A, f0.2, A)', "Input error! 'thgamcutmin' value (", thgamcutmin, ")&
                must be lower than 'thgamcutmax' value (", thgamcutmax, ")"
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

module utilities
    use numeric
    implicit none

    contains
        function inverseCDF(x) !!! Putting it inside the module eventgen doesn't work
            real(r14), intent(in) :: x
            real(r14) :: inverseCDF
            
            inverseCDF = 1/(2 - 4*x + sqrt(5 - 16*x + 16* x**2))**(1./3.) - (2 - 4*x +&
            sqrt(5 - 16*x + 16*x**2))**(1./3.)
    
        end function inverseCDF

        subroutine deboost(boosted, boostvector, deboosted, ii)

            real(r14), intent(inout) :: boosted(:)
            real(r14), intent(inout) :: boostvector(:)
            real(r14), intent(inout) :: deboosted(:,:)
            integer(i10), intent(in) :: ii
            real(r14) :: lormtrx(4,4), boostE, boostp, beta, gamma, boostcosth, boostsinth, &
            boostcosphi, boostsinphi
            integer :: k, j

            boostE = boostvector(4)
            boostp = sqrt(boostvector(1)**2+boostvector(2)**2+ boostvector(3)**2)
            beta  = boostp/boostE
            gamma = 1./sqrt(1.-beta**2)
            ! Invert the spatial components since we are DEBOOSTING
            do k = 1, 3
               boostvector(k) = -boostvector(k)
            end do
            boostcosth = boostvector(3)/boostp
            boostsinth = sqrt(1.-boostcosth**2)
            boostcosphi   = boostvector(1)/(boostp*boostsinth)
            boostsinphi   = boostvector(2)/(boostp*boostsinth)


   
            ! Fill in the Lorentz boost matrix elements in terms of theta and phi

            lormtrx(1,1) = 1. + (gamma-1.) * boostsinth**2 * boostcosphi**2
            lormtrx(1,2) = (gamma-1.) * boostsinth**2 * boostsinphi * boostcosphi
            lormtrx(1,3) = (gamma-1.) * boostsinth * boostcosth * boostcosphi
            lormtrx(1,4) = -gamma * beta * boostsinth * boostcosphi

            lormtrx(2,1) = (gamma-1.) * boostsinth**2 * boostsinphi * boostcosphi
            lormtrx(2,2) = 1. + (gamma-1.) * boostsinth**2 * boostsinphi**2
            lormtrx(2,3) = (gamma-1.) * boostsinth * boostcosth * boostsinphi
            lormtrx(2,4) = -gamma * beta * boostsinth * boostsinphi

            lormtrx(3,1) = (gamma-1.) * boostsinth * boostcosth * boostcosphi
            lormtrx(3,2) = (gamma-1.) * boostsinth * boostcosth * boostsinphi
            lormtrx(3,3) = 1. + (gamma-1.) * boostcosth**2
            lormtrx(3,4) = -gamma * beta * boostcosth

            lormtrx(4,1) = -gamma * beta * boostsinth * boostcosphi
            lormtrx(4,2) = -gamma * beta * boostsinth * boostsinphi
            lormtrx(4,3) = -gamma * beta * boostcosth
            lormtrx(4,4) = gamma


            do k = 1, 4
                deboosted(ii,k) = 0.
                do j = 1, 4
                    deboosted(ii,k) = deboosted(ii,k)+lormtrx(k,j)*boosted(j)
                enddo
            enddo
            
            
        end subroutine deboost




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
        subroutine inithistsborn() !!! Add other histos

            enemu_max = cme/2
            pmodmu_max = sqrt(enemu_max**2-mmu**2)
            allocate(h_pmod1(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_emu1(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu1(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu- histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            
        end subroutine inithistsborn

        subroutine inithistsisr() !!! Add other histos

            enemu_max = cme/2
            pmodmu_max = sqrt(enemu_max**2-mmu**2)
            pgam_max = cme
            allocate(h_pmod1(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_emu1(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu- energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu1(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu- histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_pmod2(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu+ p modulus histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_emu2(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Mu+ energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thmu2(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Theta mu+ histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_qq(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "Muons invariant mass histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_pgam(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "ISR gamma energy histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif
            allocate(h_thgam(nbins+1,2), stat=allerr)
            if (allerr /= 0) then
                print *, "ISR gamma theta histogram allocation request denied, aborting!"
                print*, ""
                stop
            endif

            
        end subroutine inithistsisr

        subroutine updatehist(hist,histmin,histmax,val)
            real(r14), intent(inout) :: hist(:,:) 
            real(r14), intent(in) :: histmin,histmax,val
            integer :: bin

            if (val >= histmax) then
                hist(nbins+1,1)= hist(nbins+1,1) + 1
            else
                bin = int((val-histmin)/(histmax-histmin)*nbins) +1
                hist(bin,1)= hist(bin,1) + 1
            endif 
            
            
        end subroutine updatehist
    

end module histogram

module eventgen

    use numeric
    use inputs
    use histogram
    use utilities
    implicit none
    real(r14) :: inte, intemax = 1., tmpintemin = 0., tmpintemax = 0., sinthmu, phimu, cosphimu, &
    sinphimu, sinthgam
    real(r14) :: rndm(7), jacqq, jacgam, jacmuang, pgamvirt(4), z
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

            call testcuts()
            
            if (accepted) then

                accpt(iev) = 1
                naccpt = naccpt + 1 

                if (histsave /= '') then
            
                    call updatehist(h_pmod1, 0.0_r14, pmodmu_max, pmod1(iev))
                    call updatehist(h_emu1, 0.0_r14, enemu_max, pmu1(iev,4))
                    call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)))
            
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

            call isrmuqq()
            call isrgammas()
            call isrmumom()


            call testcuts()


            if (accepted) then

                inte = gev2nbarn*jacqq*jacgam*jacmuang/(4.*pi*sinv)!*Matrix(Leptonic,Hadronic)

                if (run == 1) then
                    if(inte > tmpintemax) tmpintemax = inte
                    if(inte < tmpintemin) tmpintemin = inte
                endif

                if (run == 2) then
                    
                    if (inte > intemax) print*, "Warning: Found integrand greater than estimated maximum!"
                    z = rndm(7)

                    if ( z <= inte/intemax ) then
                        accpt(iev) = 1
                        naccpt = naccpt + 1 

                        if (histsave /= '') then
                    
                            call updatehist(h_pmod1, 0.0_r14, pmodmu_max, pmod1(iev))
                            call updatehist(h_emu1, 0.0_r14, enemu_max, pmu1(iev,4))
                            call updatehist(h_thmu1, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu1(iev)))
                            call updatehist(h_pmod2, 0.0_r14, pmodmu_max, pmod2(iev))
                            call updatehist(h_emu2, 0.0_r14, enemu_max, pmu2(iev,4))
                            call updatehist(h_thmu2, 0.0_r14, 180.0_r14, radtodeg*acos(costhmu2(iev)))
                            call updatehist(h_pgam, 0.0_r14, pgam_max, pgam(iev,4))
                            call updatehist(h_thgam, 0.0_r14, 180.0_r14, radtodeg*acos(costhgam(iev)))
                            call updatehist(h_qq, 0.0_r14, qqmax, radtodeg*acos(costhmu1(iev)))
                    
                        endif
                    
                    end if
                endif

            else
                accpt(iev) = 0
            endif            
            
        end subroutine genisr

        subroutine isrmuqq()

            real(r14) :: x, fak1, amin, amax, a, bmin, b, p, ppp, y
            x = rndm(1)
            fak1 = -1./sinv
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
            jacqq = (a+b)/(1./(sinv*(sinv-qq(iev))) + 1./sinv/qq(iev))

            
        end subroutine isrmuqq

        subroutine isrgammas()

            real(r14) :: x, phigam, b, cmin, cmax, y, sinthgam
            integer :: ii
            x = rndm(3)
            phigam = 2.*pi*rndm(4)
            b = sqrt(1. - 4.*me*me/sinv)
            cmin = log((1.+b*costhgammin)/(1.-b*costhgammin))/(2.*b)
            cmax = log((1.+b*costhgammax)/(1.-b*costhgammax))/(2.*b)
            y = cmin+x*(cmax-cmin)
            costhgam(iev) = tanh(b*y)/b
            sinthgam = sqrt(1.-costhgam(iev)**2)
            jacgam = 2.*pi*(1.- b**2 * costhgam(iev)**2)*(cmax-cmin)

            
            pgam(iev,4) = cme/2*(1-qq(iev)/sinv)
            pgam(iev,1) = pgam(iev,4) * sinthgam * cos(phigam)
            pgam(iev,2) = pgam(iev,4) * sinthgam * sin(phigam)
            pgam(iev,3) = pgam(iev,4) * costhgam(iev)

            pgamvirt(4) = cme/2*(1+qq(iev)/sinv)
            do ii = 1,3
                pgamvirt(ii) = -pgam(iev,ii)
            enddo

            
        end subroutine isrgammas

        subroutine isrmumom()
            
            real(r14) :: x, costhmucm, sinthmucm, phimucm, pmodcm, tmpvect(4)
            integer :: ii
            x = rndm(5)
            phimucm = 2.*pi*rndm(6)
            costhmucm = costhmumin+(costhmumax-costhmumin)*x
            jacmuang = 2.*pi*(costhmumax-costhmumin)
            sinthmucm = sqrt(1.-costhmucm*costhmucm)
            pmu1(iev,4) = sqrt(qq(iev))/2
            pmodcm = sqrt( pmu1(iev,4)**2 - mmu**2 )
            pmu1(iev,1) = pmodcm * sinthmucm * cos(phimucm)
            pmu1(iev,2) = pmodcm * sinthmucm * sin(phimucm)
            pmu1(iev,3) = pmodcm * costhmucm
            pmu2(iev,4) = pmu1(iev,4) 
            do ii = 1,3
                pmu2(iev,ii) = -pmu1(iev,ii)
            enddo
            do ii = 1,4
                tmpvect(ii) = pmu1(iev,ii)
            enddo
            call deboost(tmpvect,pgamvirt,pmu1,iev)

            do ii = 1,4
                tmpvect(ii) = pmu2(iev,ii)
            enddo
            call deboost(tmpvect,pgamvirt,pmu2,iev)


            
        end subroutine isrmumom

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
                log3 = pgam(iev,4) >= gmin
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

    allocate(pmu1(ngen,4), stat=err) !!! Maybe only allocate the big event arrays if evsave is called
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
        allocate(pgam(ngen,4), stat=err)
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
        allocate(qq(ngen), stat=err)
        if (err /= 0) then
            print *, "qq array allocation request denied, aborting!"
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
        
        if (histsave /= '') call inithistsisr()
        run = 1

        do iev = 1, nmax !!! Think about parallelisation

             call random_number(rndm)
             call genisr()
            
        end do

        run = 2
        intemax = tmpintemax

        do iev = 1, ngen

            call random_number(rndm)
            call genisr()
            
        end do

        call writevents()

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