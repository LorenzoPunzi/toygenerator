c ======================================================================== c
c                                                                          c
c   Monte Carlo generator for simulating hadron-production with radiated   c
c   photons in e+ e- collisions at high luminosity meson factories         c
c                                                                          c 
c                         PHOKHARA version 1.0                             c
c                                                                          c
c   simulates the production of two charged pions or muons together with   c
c   one or two hard photons emitted from the initial state;  it includes   c
c   virtual and soft photon corrections to the emission  of  one  single   c 
c   real photon.                                                           c
c                                                                          c
c ------------------------------------------------------------------------ c
c     Based on EVA: S.Binner, J.H.Kuehn and K.Melnikov, PLB459(1999)279    c
c     Includes NLO corrections to ISR from                                 c 
c     G.Rodrigo, A.Gehrmann, M.Guilleaume and J.H.Kuehn,                   c 
c         EPJC22(2001)81, hep-ph/0106132                                   c
c     G.Rodrigo, APPB32(2001)3833, hep-ph/0111151                          c
c     G.Rodrigo, H.Czyz, J.H.Kuehn and M.Szopa, hep-ph/0112184             c
c ======================================================================== c
c     version 1.0: (c) November 2001,   http://cern.ch/grodrigo/phokhara/  c
c ======================================================================== c




c *************************************************************************************
c ****************************** MAIN PROGRAM ********************************
c *************************************************************************************




c   The include statement is used in legacy fortran (today we would use modules) in order to load some 
c   variables and potentially subroutines. Here are the contents of "phokhara.inc"
c
c      implicit none 
c --- couplings, masses, momenta ---
c double precision pi,gev2nbarn,gev2pbarn,alpha,me,mmu,mtau,mtop,
c & mrho,mpi,grhoee,gammarho,mrhol,grhol,momega,gomega,al,be,
c & S,ebeam,momenta(10,0:3)
c --- cuts ---     
c double precision Emin,gmin,phot1cut,phot2cut,pi1cut,pi2cut,
c & piphcut,accecut,q2min,w
c --- histograms ---
c character*20 title(10)
c double precision xlow(10),xup(10)
c integer bins(10)
c --- Maxima ---
c double precision Mmax(2),gross(2),klein(2),tr(2),count(2)
c ---
c common /ctes/    pi,gev2nbarn,gev2pbarn,alpha,me,mmu,mtau,mtop,
c & mrho,mpi,grhoee,gammarho,mrhol,grhol,momega,gomega,al,be,
c & S,ebeam,momenta
c common/cuts/ Emin,gmin,phot1cut,phot2cut,pi1cut,pi2cut,
c & piphcut,accecut,q2min,w
c common /histo/   xlow,xup,bins,title
c common /maxima/  Mmax,gross,klein,tr,count
c
c Where the COMMON command makes the variables accessible to any subprogram (see SP2)
c

c ==================================== MAIN PROGRAM ====================================

      include 'phokhara.inc'
      double precision qqmin,qqmax,
     &  cos1min,cos1max,cos2min,cos2max,cos3min,cos3max,
     &  dsigm1,dsigm2,sigma1,sigma2,sigma,dsigm,Ar(14)
      integer nges,nm,i,j,k,nlo,vacua,pion,i1,n1,n2
      character outfile*20,fname*20

c   Reads from input "seed.dat" and initialises RM48 CERN RNG : https://cds.cern.ch/record/2050843/files/v116.pdf
c --- reads the seed for random routine RM48 ------
      open(9,file='seed.dat',type='old')
      read (9,*) i1,n1,n2
      call rm48in(i1,n1,n2)
c      close (9,DISP='delete')
      
c --- input parameters ----------------------------
c This call saves parameters from the input card "input.dat" to either variables ocntained in
c "phokhara.inc" or the following variables explicitly stated
c 
c nges = number of generated events
c nm = number of events to determine the maximum
c outfile = name of output file
c fname = name of output hbook file
c nlo = Born(0) or NLO (1)
c vacua = vacuum polarization: yes(1), no(0)
c pions = muons(0), pions(>0) 

      call input(nges,nm,outfile,fname,nlo,vacua,pion)

c --- open output file for generated momenta ------
      open (8,file=outfile,type='new')

c --- print run data ------------------------------
      write (*,*) '---------------------------------------------------'
      if (pion.eq.0) then 
         write (*,*) '       EVA MC: e^+ e^- -> mu^+ mu^- gamma'
      else
         write (*,*) '       EVA MC: e^+ e^- -> pi^+ pi^- gamma'        
      endif
      write (*,*) '---------------------------------------------------'
      write (*,100) 'cms total energy     = ',dSqrt(S),' GeV'
      write (*,100) 'tagged photon energy = ',gmin,' GeV'
      write (*,*) 'cuts on the photon angle = ',
     &  int(phot1cut),',',int(phot2cut)
      if (pion.eq.0) then 
         write (*,*) 'cuts on the muons  angle = ',
     &     int(pi1cut),',',int(pi2cut)
      else
         write (*,*) 'cuts on the pions  angle = ',
     &     int(pi1cut),',',int(pi2cut)
      endif
      write (*,110) 'minimal hadrons-photon inv mass = ',q2min,' GeV^2' 
      if (nlo.eq.0) write (*,*) 'Born'
      if (nlo.eq.1) write (*,*) 'NLO: soft photon with w = ',w
      write (*,*) 'vacua = ',vacua

c --- book (initialise) histograms -----------------------------
      call inithisto

c --- set cuts ------------------------------------
      cos1min = dCos(phot2cut*pi/180.d0)    ! photon1 angle cuts in the 
      cos1max = dCos(phot1cut*pi/180.d0)    ! LAB rest frame            
      cos2min = -1.d0                       ! photon2 angle limits      
      cos2max = 1.d0                        !                           
      cos3min = -1.d0                       ! pion angle limits in the  
      cos3max = 1.d0                        ! rho rest frame            
      if (pion.eq.0) then                   ! virtual photon energy cut 
         qqmin = 4.d0*mmu*mmu
      else
         qqmin = 4.d0*mpi*mpi
      endif
c     S = sqrt(s)^2 of e+e-   gmin = minimum tagged photon energy
      qqmax = S-2.d0*dSqrt(S)*gmin          ! if only one photon, from 4-momentum conservation
           
c =================================================
c --- finding the maximum -------------------------
      k = nm                        
      do i = 1,2
         Mmax(i) = 1.d0
         gross(i) = 0.d0 ! Temporary values of maximum (gross) and minimum (klein)
         klein(i) = 0.d0      
      enddo 
      if (nlo.eq.0) Mmax(2)=0.d0   ! only 1 photon events generated

c   ************ BEGINNING OF THE TWO GENERATION RUNS **************+
c   Do two runs, first (nm) to find the maximum (Mmax), the second (nges) to actually generate events 
      
      do i = 1,2        ! initializing the MC loop
      tr(1) = 0.d0    ! Number of accepted 1/2 photon events    
      tr(2) = 0.d0
      count(1) = 0.d0   ! Number of generated 1/2 photon events 
      count(2) = 0.d0

c =================================================
c --- beginning the MC loop event generation ------

c ----- First loop to find Mmax
      do j = 1,k

c   Toss a uniform (0,1) to see if one or two photons using the ratio of the maxima (first run it's 50/50)
         call rm48(Ar,1)        

         if (Ar(1).le.(Mmax(1)/(Mmax(1)+Mmax(2)))) then 

            count(1) = count(1)+1.d0
            call gen_1ph(i,S,qqmin,qqmax,cos1min,cos1max,
     &	         cos3min,cos3max,pion,nlo,vacua)
 
         else
	 
            count(2) = count(2)+1.d0
            call gen_2ph(i,S,qqmin,cos1min,cos1max,
     &           cos2min,cos2max,cos3min,cos3max,pion,nlo,vacua)

         endif
       
      enddo
c --- end of the MC loop --------------------------
c =================================================

c --- for the second run ---
      k = nges
      if (i.eq.1) then
c       For actual generation run Mmax is set to the value found in first run (with a bit of room)
         Mmax(1) = gross(1)+.05d0*Sqrt(gross(1)*gross(1))
         Mmax(2) = gross(2)+(.03d0+.02d0*S)*Sqrt(gross(2)*gross(2)) 
      endif

      enddo       
c --- end of the second run -----------------------
c =================================================

c   ************ END OF THE TWO GENERATION RUNS **************+

c --- save histograms -----------------------------
      call endhisto(fname)

c --- value of the cross section ------------------
      if (nlo.eq.0) then 
         sigma = Mmax(1)/count(1)*tr(1)
         dsigm = 
     &    Mmax(1)*dSqrt((tr(1)/count(1)-(tr(1)/count(1))**2)/count(1))
      else

c   sigma = integral = fraction of accepted (tr/count) * max (Mmax) * range of random variables (1)  

         sigma1 = Mmax(1)/count(1)*tr(1)
         sigma2 = Mmax(2)/count(2)*tr(2)
         sigma = sigma1+sigma2
         dsigm1= 
     &    Mmax(1)*dSqrt((tr(1)/count(1)-(tr(1)/count(1))**2)/count(1))
         dsigm2= 
     &    Mmax(2)*dSqrt((tr(2)/count(2)-(tr(2)/count(2))**2)/count(2))
         dsigm = dSqrt(dsigm1**2+dsigm2**2)
      endif

c --- output --------------------------------------
      write (*,*) '---------------------------------------------------'      
      write (*,*) int(tr(1)+tr(2)),   ' total events accepted of '
      write (*,*) nges,' total events generated'
      write (*,*) int(tr(1)),   ' one photon events accepted of '
      write (*,*) int(count(1)),' events generated'
      write (*,*) int(tr(2)),' two photon events accepted of '
      write (*,*) int(count(2)),' events generated'
      write (*,*)
      if (nlo.ne.0) write (*,*) 'sigma1(nbarn) = ',sigma1,' +- ',dsigm1
      if (nlo.ne.0) write (*,*) 'sigma2(nbarn) = ',sigma2,' +- ',dsigm2
      write (*,*) 'sigma (nbarn) = ',sigma, ' +- ',dsigm
      write (*,*)
      write (*,*) 'maximum1 = ',gross(1),'  minimum1 = ',klein(1)
      write (*,*) 'Mmax1    = ',Mmax(1)            
      write (*,*) 'maximum2 = ',gross(2),'  minimum2 = ',klein(2)
      write (*,*) 'Mmax2    = ',Mmax(2)            
      write (*,*) '---------------------------------------------------'
      write (*,*) 'output file: ',fname
 100  format (a24,f8.4,a5)
 110  format (a35,f8.4,a6)

c --- saves the new seed --------------------------
      open(9,file='seed.dat',type='new')
      call rm48ut(i1,n1,n2) ! to pick up the generation where you left off if you run the program again
      write (9,*) i1,n1,n2
      
      
      end

c End of the main program (without explicit name because legacy)




c ======================================================================
c --- input parameters -------------------------------------------------
c ======================================================================
      subroutine input(nges,nm,outfile,fname,nlo,vacua,pion)
      include 'phokhara.inc'
      integer nges,nm,i,nlo,vacua,pion
      character outfile*20,fname*20
      double precision E

      real*8 p1(4),p2(4),dme,el_m2              ! inteface to helicity amp
      common /cp1p2/p1,p2,dme,el_m2

c --- input file ---------------------------------
      open(7,file='input.dat',type='old')
c --- input generation parameters ----------------
      read(7,*)           !                                   
      read(7,*) nges      ! number of generated events        
      read(7,*) nm        ! events to determine the maximum   
      read(7,*) outfile   ! output file                       
      read(7,*) nlo       ! Born(0),NLO(1)
      read(7,*) vacua     ! vacuum polarization: yes(1), no(0)
      read(7,*) w         ! soft photon cutoff                
      read(7,*) pion      ! muons(0), pions(>0)    
c --- input couplings, masses and meson widths ---
      read(7,*)           !                                   
      read(7,*) alpha     ! 1/alpha (QED)                     
      read(7,*) me        ! Electron mass                     
      read(7,*) mmu       ! Muon mass                         
      read(7,*) mtau      ! Tau mass                          
      read(7,*) mtop      ! Top mass                          
      read(7,*) mpi       ! Pion mass                         
      read(7,*) mrho      ! Rho mass                          
      read(7,*) grhoee    ! Rho->e+e- width                   
      read(7,*) gammarho  ! Total rho width                   
      read(7,*) mrhol     ! Rho' mass                         
      read(7,*) grhol     ! Rho' width                        
      read(7,*) momega    ! Omega mass                        
      read(7,*) gomega    ! Omega width                       
      read(7,*) al        ! Pion form factor parameters a     
      read(7,*) be        ! and b                             
c --- input collider parameters -------------------
      read(7,*)           !                                   
      read(7,*) E         ! CMS-energy                        
c --- input experimental cuts ---------------------
      read(7,*)           !                                   
      read(7,*) q2min     ! minimal pi,pi,gamma-inv mass squared     
      read(7,*) gmin      ! minimal photon energy             
      read(7,*) phot1cut  ! minimal photon angle              
      read(7,*) phot2cut  ! maximal photon angle              
      read(7,*) pi1cut    ! minimal pions/muons angle               
      read(7,*) pi2cut    ! maximal pions/muons angle               
c --- read histogram paremeters -------------------
      read(7,*)           !                                   
      read(7,*) fname     ! hbook output file                 
      do i = 1,10         ! read title, limits and bins       
         read(7,*) title(i)
         read(7,*) xlow(i),xup(i),bins(i)
      enddo
c --- fix constants -------------------------------
      alpha = 1.d0/alpha
      S = E*E                                ! CMS-energy squared
      ebeam = dSqrt(S)/2.d0                  ! beam energy
      pi = 4.d0*dAtan(1.d0)   
      gev2nbarn = .389379292d6               ! from GeV^2 to nbarn
      gev2pbarn = .389379292d9

      dme   = me                             ! inteface to helicity amp
      el_m2 = me**2
      p1(1) = ebeam
      p1(2) = 0.d0
      p1(3) = 0.d0
      p1(4) = sqrt(ebeam**2-el_m2)
      p2(1) = ebeam
      p2(2) = 0.d0
      p2(3) = 0.d0
      p2(4) = -p1(4)
      return
      end







c ******************************************************************************************
c ********************************* BEGIN SUBROUTINES **************************************
c ******************************************************************************************




c ======================================================================
c --- generates one photon ---------------------------------------------
c ======================================================================
      subroutine gen_1ph(i,Sp,qqmin,qqmax,cos1min,cos1max,
     &                   cos3min,cos3max,pion,nlo,vacua)
      include 'phokhara.inc'
      double precision Sp,qqmin,qqmax,qq,jac0,
     &  cos1min,cos1max,cos1,phi1,jac1,
     &  cos3min,cos3max,cos3,phi3,jac4,
     &  Ar(14),z,inte,vacuumpolarization,
     &  Matrix,Leptonic(0:3,0:3),Hadronic(0:3,0:3)
      integer i,pion,nlo,vacua
      logical accepted      

      call rm48(Ar,7)
      
c --- get the variables ---------------------------
c I think jacobians are used because you always generate U(0,1)
      call qquadrat(Sp,qqmin,qqmax,Ar,pion,qq,jac0) 
      call photonangles1(Sp,cos1min,cos1max,Ar,cos1,phi1,jac1)
      call pionangles(cos3min,cos3max,Ar,cos3,phi3,jac4)
      z = Mmax(1)*Ar(6)  
             
c --- These subroutines fill the momenta(10,0:3) array with the 4-momenta in the CMS frame ---
      call leptonmomenta1(Sp,qq,cos1,phi1) ! Momenta of the phisical and virtual photon and e+,e-
      call LeptonicTensor1ISR(Sp,qq,cos1,nlo,Leptonic) ! Build the leptonic tensor (physics inside)
      call hadronmomenta(Sp,qq,pion,cos3,phi3) ! Momenta of the two mu/π
      call HadronicTensorISR(qq,pion,Hadronic) ! Build the hadronic tensor (doesn't have pion FF if muons)

c --- Check if the event passes the required cuts ---
      call testcuts(1,accepted,pion)
c --- value of the integrand ---
      if (accepted) then
c ---
         inte = gev2nbarn*Matrix(Leptonic,Hadronic)*
     &       jac0*jac1*jac4/(4.d0*pi*Sp) ! cfr eq 108 "quest for precision..."
         if (vacua.eq.1) inte = inte/(1.d0-vacuumpolarization(qq))   !vp
c ---
         if (inte.gt.gross(1)) gross(1) = inte ! Update maximum
         if (inte.lt.klein(1)) klein(1) = inte ! Update minimum
c --- in the second run ---
         if (i.eq.2) then
            if (Mmax(1).lt.inte) write(*,*) 'Warning! Max(1) too small!' ! shouldn't happen
            if (inte.lt.0.d0)    write(*,*) 'Warning! negative weight'
c --- event accepted? i.e. event uniformly between 0 and Mmax (i.e. z) is < f(x)?---
            if (z.le.inte) then
               tr(1) = tr(1) + 1.d0
               call writeevent(pion)
c --- add to the histogrammes ---
               call addiere(1.d0,qq,1)
            endif 
         endif
      else
         inte = 0.d0
      endif
      return
      end


c ======================================================================
c --- generates two photons --------------------------------------------
c ======================================================================
     subroutine gen_2ph(i,Sp,qqmin,cos1min,cos1max,
     &   cos2min,cos2max,cos3min,cos3max,pion,nlo,vacua)
      include 'phokhara.inc'
      double precision Sp,qqmin,qqmax,qq,jac0,
     &  cos1min,cos1max,cos1,phi1,jac1,sin1,cos12,
     &  cos2min,cos2max,cos2,phi2,jac2,sin2,
     &  w1,w2,w1min,jac3,
     &  cos3min,cos3max,cos3,phi3,jac4,
     &  Ar(14),z,inte,vacuumpolarization,
     &  Matrix,Leptonic(0:3,0:3),Hadronic(0:3,0:3),helicityamp
      integer i,pion,nlo,vacua
      logical accepted      

      call rm48(Ar,12)
      
c --- get the variables -----------------------------------------------
c --- one of the photons is generated inside the angular cuts and -----
c --- the other is generated everywhere -------------------------------
      if (Ar(12).lt.0.5d0) then ! randomly choose which photon is going to have which angle
         call photonangles1(Sp,cos1min,cos1max,Ar,cos1,phi1,jac1)
         call photonangles2(Sp,cos2min,cos2max,Ar,cos2,phi2,jac2)  
         if (cos2.lt.cos1min.or.cos2.gt.cos1max) jac1=2.d0*jac1 ! ????
       else
         call photonangles1(Sp,cos2min,cos2max,Ar,cos1,phi1,jac1)
         call photonangles2(Sp,cos1min,cos1max,Ar,cos2,phi2,jac2)  
         if (cos1.lt.cos1min.or.cos1.gt.cos1max) jac2=2.d0*jac2
      endif  
         sin1 = dSqrt(1.d0-cos1*cos1)
         sin2 = dSqrt(1.d0-cos2*cos2)
         cos12 = sin1*sin2*dCos(phi1-phi2)+cos1*cos2
         w1min = gmin/dSqrt(Sp)
         qqmax = Sp*(1.d0-2.d0*(w1min+w)+2.d0*w1min*w*(1.d0-cos12))
      call qquadrat(Sp,qqmin,qqmax,Ar,pion,qq,jac0)
      call photonenergy2(Sp,qq,cos1min,cos1max,cos1,cos2,cos12
     &                        ,Ar,w1,w2,jac3)
      call pionangles(cos3min,cos3max,Ar,cos3,phi3,jac4)
      z = Mmax(2)*Ar(6)         
      
c --- 4-momenta in the CMS frame ------
      call leptonmomenta2(Sp,qq,w1,w2,cos1,phi1,cos2,phi2)
      call LeptonicTensor2ISR(Sp,qq,w1,w2,cos1,phi1,cos2,phi2,Leptonic)
      call hadronmomenta(Sp,qq,pion,cos3,phi3) ! This assumes the angles of the muons are flat like π ????
      call HadronicTensorISR(qq,pion,Hadronic)

c --- tests cuts ---
      call testcuts(2,accepted,pion)
c --- value of the integrand ---
      if (accepted) then
c --- leptonic tensor, or ---
c         inte = gev2nbarn*Matrix(Leptonic,Hadronic)*
c     &       jac0*jac1*jac2*jac3*jac4/(4.d0*pi*Sp)
c --- helicity amplitudes ---
         inte = gev2nbarn*helicityamp(Sp,qq,pion)*
     &       jac0*jac1*jac2*jac3*jac4/(4.d0*pi*Sp)
c ---     
         if (vacua.eq.1) inte = inte/(1.d0-vacuumpolarization(qq))   !vp
c ---
         if (inte.gt.gross(2)) gross(2) = inte
         if (inte.lt.klein(2)) klein(2) = inte
c --- in the second rund ---
         if (i.eq.2) then
            if (Mmax(2).lt.inte) write(*,*) 'Warning! Max(2) too small!'
            if (inte.lt.0.d0)    write(*,*) 'Warning! negative weight'
c --- event accepted? ---
            if (z.le.inte) then
               tr(2) = tr(2) + 1.d0
               call writeevent(pion)
c --- add to the histogrammes ---
               call addiere(1.d0,qq,2)
            endif
         endif
      else
         inte = 0.d0
      endif
      return
      end


c =================================================
c --- print the momenta of the generated event ----
c =================================================
      subroutine writeevent(pion)
      include 'phokhara.inc'
      integer pion

      write(8,*),'-------------------------------------------'
      write(8,*),'Photon1:',momenta(3,0),momenta(3,1),
     &   momenta(3,2),momenta(3,3)
      write(8,*),'Photon2:',momenta(4,0),momenta(4,1),
     &   momenta(4,2),momenta(4,3)
      if (pion.eq.0) then 
         write(8,*),'Mu+:    ',momenta(6,0),momenta(6,1),
     &      momenta(6,2),momenta(6,3)
         write(8,*),'Mu-:    ',momenta(7,0),momenta(7,1),
     &      momenta(7,2),momenta(7,3)
      else
         write(8,*),'Pi+:    ',momenta(6,0),momenta(6,1),
     &      momenta(6,2),momenta(6,3)
         write(8,*),'Pi-:    ',momenta(7,0),momenta(7,1),
     &      momenta(7,2),momenta(7,3)
      endif 
      return
      end


c ========================================================================
c ------------------------------------------------------------------------
c --- Generates: photon virtuality Q2 using a two-component substitution -
c ------------------------------------------------------------------------
      subroutine qquadrat(Sp,qqmin,qqmax,Ar,pion,qq,jacobian)
      include 'phokhara.inc'
      double precision Sp,qqmin,qqmax,Ar(14),qq,jacobian,
     &  x,a,b,amin,amax,bmin,bmax,fak1,fak2,p,y,ppp
      integer pion

      x = Ar(1)      
c --- flat ---
c      a = qqmax-qqmin
c      qq = qqmin+a*x
c      jacobian = a

c --- muons -----------------------------------------------------
c --- the Q2 distribution is peaked at threshold and Q2->Sp (soft photon). 
      if (pion.eq.0) then
      fak1 = -1.d0/Sp
      amin = fak1*dLog(Sp-qqmin)
      amax = fak1*dLog(Sp-qqmax)
      a = amax-amin
      bmin = dLog(qqmin/Sp)/Sp
      b    = dLog(qqmax/qqmin)/Sp
c --- which substitution? ---
      p = Ar(7)      
      ppp  = a/(a+b)
      if(p.lt.ppp)then
         y  = amin+a*x
         qq = Sp-dExp(y/fak1)                                       
      else
         y  = bmin+b*x
         qq = Sp*exp(Sp*y)
      endif
      jacobian = (a+b)/(1.d0/(Sp*(Sp-qq)) + 1.d0/Sp/qq)  

c --- pions -----------------------------------------------------
c --- the Q2 distribution is peaked at Q2=rho and w mass^2       
c --- resonances (pion form factor) and at Q2->Sp (soft photon). 
      else
      fak1 = -1.d0/Sp
      amin = fak1*dLog(Sp-qqmin)
      amax = fak1*dLog(Sp-qqmax)
      a = amax-amin
      fak2 = 1.d0/gammarho/mrho
      bmin = fak2*dAtan((qqmin-mrho**2)*fak2)
      bmax = fak2*dAtan((qqmax-mrho**2)*fak2)
      b = bmax-bmin
c --- two channels ---
      p = Ar(7)      
      ppp  = a/(a+b)
      if(p.lt.ppp)then
         y  = amin+a*x
         qq = Sp-dExp(y/fak1)                                       
      else
         y  = bmin+b*x
         qq = mrho*(mrho+gammarho*dTan(y/fak2))                    
      endif
      jacobian = (a+b)/( (1.d0/(Sp*(Sp-qq)) +
     &     1.d0/(((qq-mrho**2)**2+(gammarho*mrho)**2))) )
      endif
      return
      end


c ========================================================================
c ------------------------------------------------------------------------
c --- Generates: real photon costheta and phi angles in the e^+e^- CMS ---
c ------------------------------------------------------------------------
      subroutine photonangles1(Sp,cosmin,cosmax,Ar,
     &           costheta,phi,jacobian)
      include 'phokhara.inc'
      double precision Sp,cosmin,cosmax,Ar(14),
     &  costheta,phi,jacobian,x,b,cmin,cmax,y

      x   = Ar(2)
      phi = 2.d0*pi*Ar(3)

c --- flat ---
c      costheta = cosmin+(cosmax-cosmin)*x
c      jacobian = 2.d0*pi*(cosmax-cosmin)
c --- peaked at costheta = +-1 ---
      b = dSqrt(1.d0-4.d0*me*me/Sp)
      cmin = dLog((1.d0+b*cosmin)/(1.d0-b*cosmin))/(2.d0*b)
      cmax = dLog((1.d0+b*cosmax)/(1.d0-b*cosmax))/(2.d0*b)
      y = cmin+x*(cmax-cmin)
      costheta = dTanh(b*y)/b
      jacobian = 2.d0*pi*(1.d0-b*b*costheta*costheta)*(cmax-cmin)
      return
      end
      
c ========================================================================
c ------------------------------------------------------------------------
c --- Generates: real photon costheta and phi angles in the e^+e^- CMS ---
c Difference with photonangles1 should just be the element of Ar used for RNG
c ------------------------------------------------------------------------
      subroutine photonangles2(Sp,cosmin,cosmax,Ar,
     &           costheta,phi,jacobian)
      include 'phokhara.inc'
      double precision Sp,cosmin,cosmax,Ar(14),
     &  costheta,phi,jacobian,x,b,cmin,cmax,y

      x   = Ar(10)
      phi = 2.d0*pi*Ar(11)

c --- flat ---
c      costheta = cosmin+(cosmax-cosmin)*x
c      jacobian = 2.d0*pi*(cosmax-cosmin)
c --- peaked at costheta = +-1 ---
      b = dSqrt(1.d0-4.d0*me*me/Sp)
      cmin = dLog((1.d0+b*cosmin)/(1.d0-b*cosmin))/(2.d0*b)
      cmax = dLog((1.d0+b*cosmax)/(1.d0-b*cosmax))/(2.d0*b)
      y = cmin+x*(cmax-cmin)
      costheta = dTanh(b*y)/b
      jacobian = 2.d0*pi*(1.d0-b*b*costheta*costheta)*(cmax-cmin)
      return
      end      
      

c ========================================================================
c ------------------------------------------------------------------------
c --- Generates: real photon energy in the e^+e^- CMS normalized to the --
c --- CMS energy ---------------------------------------------------------
c ------------------------------------------------------------------------
      subroutine photonenergy2(Sp,qq,c1min,c1max,cos1,cos2,cos12,Ar
     &                        ,w1,w2,jacobian)
      include 'phokhara.inc'
      double precision Sp,qq,c1min,c1max,cos1,cos2,cos12,Ar(14),
     &  w1,w2,jacobian,x,wmin,q2,qqb,w1max,w2max,u,umin,umax

      x = Ar(9)
      
c --- flat ---
c      w2 = w+(w2max-w)*x
c      jacobian = (w2max-w)
c --- peaked at w1, w2 = w ---

      wmin = gmin/dSqrt(Sp)
      q2 = qq/Sp
      qqb = Sp*(1.d0-4.d0*wmin+2.d0*wmin*wmin*(1.d0-cos12))

c --- photon 2 inside the angular cuts, photon 1 outside ---
c --- then w2 > wmin and w1 > w ----------------------------
      if(((cos1.gt.c1max).or.(cos1.lt.c1min)).and.
     &   ((cos2.le.c1max).and.(cos2.ge.c1min))) then

         w2max = (1.d0-q2-2.d0*w)/(2.d0*(1.d0-w*(1.d0-cos12)))         
         umin = dLog(wmin/(1.d0-q2-2.d0*wmin))
         umax = dLog(w2max/(1.d0-q2-2.d0*w2max))
         u  = umin+x*(umax-umin)
         w2 = (1.d0-q2)/(2.d0+dExp(-u))
         w1 = (1.d0-q2-2.d0*w2)/(2.d0*(1.d0-w2*(1.d0-cos12)))
	 jacobian = 1.d0/2.d0 * w2*w2*w1*w1 * (umax-umin)/(1.d0-q2)
            
c --- photon 1 inside the angular cuts, photon 2 outside ---
c --- then w1 > wmin and w2 > w ----------------------------
      elseif(((cos2.gt.c1max).or.(cos2.lt.c1min)).and.
     &   ((cos1.le.c1max).and.(cos1.ge.c1min))) then

         w1max = (1.d0-q2-2.d0*w)/(2.d0*(1.d0-w*(1.d0-cos12)))         
         umin = dLog(wmin/(1.d0-q2-2.d0*wmin))
         umax = dLog(w1max/(1.d0-q2-2.d0*w1max))
         u  = umin+x*(umax-umin)
         w1 = (1.d0-q2)/(2.d0+dExp(-u))
         w2 = (1.d0-q2-2.d0*w1)/(2.d0*(1.d0-w1*(1.d0-cos12)))            
         jacobian = 1.d0/2.d0 * w2*w2*w1*w1 * (umax-umin)/(1.d0-q2)

c --- both photons pass the angular cuts ---
      elseif(((cos2.le.c1max).and.(cos2.ge.c1min)).and.
     &   ((cos1.le.c1max).and.(cos1.ge.c1min))) then

         if(qq.le.qqb)then
	   
            w1max = (1.d0-q2-2.d0*w)/(2.d0*(1.d0-w*(1.d0-cos12)))         
            umin = dLog(w/(1.d0-q2-2.d0*w))
            umax = dLog(w1max/(1.d0-q2-2.d0*w1max))
            u  = umin+x*(umax-umin)
            w1 = (1.d0-q2)/(2.d0+dExp(-u))
            w2 = (1.d0-q2-2.d0*w1)/(2.d0*(1.d0-w1*(1.d0-cos12)))            
            jacobian = 1.d0/2.d0 * w2*w2*w1*w1 * (umax-umin)/(1.d0-q2)

         else

            if(Ar(8).lt.0.5d0)then ! Random choice I think at this point ????

            w2max = (1.d0-q2-2.d0*w)/(2.d0*(1.d0-w*(1.d0-cos12)))         
            umin = dLog(wmin/(1.d0-q2-2.d0*wmin))
            umax = dLog(w2max/(1.d0-q2-2.d0*w2max))
            u  = umin+x*(umax-umin)
            w2 = (1.d0-q2)/(2.d0+dExp(-u))
            w1 = (1.d0-q2-2.d0*w2)/(2.d0*(1.d0-w2*(1.d0-cos12)))            
            jacobian = w2*w2*w1*w1* (umax-umin)/(1.d0-q2)
	      
            else

            w1max = (1.d0-q2-2.d0*w)/(2.d0*(1.d0-w*(1.d0-cos12)))         
            umin = dLog(wmin/(1.d0-q2-2.d0*wmin))
            umax = dLog(w1max/(1.d0-q2-2.d0*w1max))
            u  = umin+x*(umax-umin)
            w1 = (1.d0-q2)/(2.d0+dExp(-u))
            w2 = (1.d0-q2-2.d0*w1)/(2.d0*(1.d0-w1*(1.d0-cos12)))            
            jacobian = w2*w2*w1*w1* (umax-umin)/(1.d0-q2)

            endif

         endif

      endif
      return
      end

c ------------------------------------------------------------------
c    lepton four-momenta: one real photon                           
c ------------------------------------------------------------------
      subroutine leptonmomenta1(Sp,qq,costheta,phi)
      include 'phokhara.inc'
      double precision Sp,qq,q2,E,costheta,sintheta,phi
      integer i

      q2 = qq/Sp
      sintheta = dSqrt(1.d0-costheta*costheta)
c --- positron ---
      momenta(1,0) = dSqrt(Sp)/2.d0
      momenta(1,1) = 0.d0
      momenta(1,2) = 0.d0
      momenta(1,3) = momenta(1,0)*dSqrt(1.d0-4.d0*me*me/Sp)
c --- electron ---
      momenta(2,0) = momenta(1,0)
      do i = 1,3
         momenta(2,i) = -momenta(1,i)
      enddo
c --- real photon1 ---
      E = (1.d0-q2)*dSqrt(Sp)/2.d0
      momenta(3,0) = E
      momenta(3,1) = E*sintheta*dCos(phi)
      momenta(3,2) = E*sintheta*dSin(phi)
      momenta(3,3) = E*costheta
c --- real photon2 ---
      momenta(4,0) = 0.d0
      momenta(4,1) = 0.d0
      momenta(4,2) = 0.d0
      momenta(4,3) = 0.d0
c --- virtual photon ---
      momenta(5,0) = (1.d0+q2)*dSqrt(Sp)/2.d0
      do i = 1,3
        momenta(5,i) = -momenta(3,i)
      enddo
      return
      end

c ------------------------------------------------------------------
c    lepton four-momenta: two real photons                          
c ------------------------------------------------------------------
      subroutine leptonmomenta2(Sp,qq,w1,w2,cos1,phi1,cos2,phi2)
      include 'phokhara.inc'
      double precision Sp,qq,w1,w2,cos1,sin1,phi1,cos2,sin2,phi2,
     &  E1,E2
      integer i

      sin1 = dSqrt(1.d0-cos1*cos1)
      sin2 = dSqrt(1.d0-cos2*cos2)
c --- positron ---
      momenta(1,0) = dSqrt(Sp)/2.d0
      momenta(1,1) = 0.d0
      momenta(1,2) = 0.d0
      momenta(1,3) = momenta(1,0)*dSqrt(1.d0-4.d0*me*me/Sp)
c --- electron ---
      momenta(2,0) = momenta(1,0)
      do i = 1,3
         momenta(2,i) = -momenta(1,i)
      enddo
c --- real photon1 ---
      E1 = w1*dSqrt(Sp)
      momenta(3,0) = E1
      momenta(3,1) = E1*sin1*dCos(phi1)
      momenta(3,2) = E1*sin1*dSin(phi1)
      momenta(3,3) = E1*cos1
c --- real photon2 ---
      E2 = w2*dSqrt(Sp)      
      momenta(4,0) = E2
      momenta(4,1) = E2*sin2*dCos(phi2)
      momenta(4,2) = E2*sin2*dSin(phi2)
      momenta(4,3) = E2*cos2
c --- virtual photon ---
      momenta(5,0) = dSqrt(Sp)-E1-E2
      do i = 1,3
        momenta(5,i) = -momenta(3,i)-momenta(4,i)
      enddo
      return
      end


c ========================================================================
c ------------------------------------------------------------------------
c --- Generates: pion costheta and phi angles in the Q2 CMS system, ------
c afterwards boosted into the e^+ e^- CMS --------------------------------
c ------------------------------------------------------------------------
      subroutine pionangles(cosmin,cosmax,Ar,costheta,phi,jacobian)
      include 'phokhara.inc'
      double precision cosmin,cosmax,Ar(14),costheta,phi,jacobian,x 

      x   = Ar(4)
      phi = 2.d0*pi*Ar(5)

c --- flat ---
      costheta = cosmin+(cosmax-cosmin)*x
      jacobian = 2.d0*pi*(cosmax-cosmin)
      return
      end

c ------------------------------------------------------------------
c     hadron four momenta                                            
c ------------------------------------------------------------------
      subroutine hadronmomenta(Sp,qq,pion,costheta,phi)
      include 'phokhara.inc'
      double precision Sp,qq,p,costheta,sintheta,phi,m2,
     &        cmsvector(0:3),boostvector(0:3),labvector(0:3)
      integer pion,i
      
      sintheta = dSqrt(1.d0-costheta*costheta)
      if (pion.eq.0) then 
         m2 = mmu*mmu
      else
         m2 = mpi*mpi
      endif
c --- pions/muon in the qq-rest frame ---
c --- pi^+/mu^+ ---
      momenta(6,0) = dSqrt(qq)/2.d0
      p            = momenta(6,0)*dSqrt(1.d0-4.d0*m2/qq)
      momenta(6,1) = p*sintheta*dCos(phi)
      momenta(6,2) = p*sintheta*dSin(phi)
      momenta(6,3) = p*costheta
c --- pi^-/mu^- ---
      momenta(7,0) = momenta(6,0) 
      do i = 1,3
        momenta(7,i) = -momenta(6,i)
      enddo
c --- boost the hadron momenta into the e^+ e^- CMS ---
c     positive charge particle
      do i =0,3
         boostvector(i) = momenta(5,i)
         cmsvector(i)   = momenta(6,i)
      enddo   
      call boost(cmsvector,boostvector,labvector)

c     negative charged particle
      do i =0,3
         momenta(6,i) = labvector(i) ! Save positive charge momenta from previous boost
         cmsvector(i) = momenta(7,i) ! Prepare for new boost of negative charge
      enddo   
      call boost(cmsvector,boostvector,labvector)
      do i =0,3
        momenta(7,i) = labvector(i) 
      enddo         
      return
      end

c -----------------------------------------------------------
c     boost cmsvector by boostvector into labvector          
c -----------------------------------------------------------
      subroutine boost(cmsvector,boostvector,labvector)
      implicit none
      double precision cmsvector(0:3),boostvector(0:3),
     &  labvector(0:3),m(0:3,0:3),
     &  E,p,beta,gamma,costheta,sintheta,cosphi,sinphi
      integer i,j
      
      E = boostvector(0)
      p = dSqrt(boostvector(1)**2+boostvector(2)**2+
     &    boostvector(3)**2)
      beta  = p/E
      gamma = 1.d0/dSqrt(1.d0-beta*beta)
      costheta = boostvector(3)/p
      sintheta = dSqrt(1.d0-costheta*costheta)
      cosphi   = boostvector(1)/(p*sintheta)
      sinphi   = boostvector(2)/(p*sintheta)

        m(0,0) = gamma
        m(0,1) = 0.d0
        m(0,2) = 0.d0
        m(0,3) = beta*gamma
        m(1,0) = beta*gamma*sintheta*cosphi
        m(1,1) = costheta*cosphi
        m(1,2) = -sinphi
        m(1,3) = gamma*sintheta*cosphi
        m(2,0) = beta*gamma*sintheta*sinphi
        m(2,1) = costheta*sinphi
        m(2,2) = cosphi
        m(2,3) = gamma*sintheta*sinphi
        m(3,0) = beta*gamma*costheta
        m(3,1) = -sintheta
        m(3,2) = 0.d0
        m(3,3) = gamma*costheta

      do i=0,3
         labvector(i) = 0.d0
         do j=0,3
            labvector(i) = labvector(i)+m(i,j)*cmsvector(j)
         enddo
      enddo
      return
      end


c ========================================================================
c ------------------------------------------------------------------------
c --- Test experimental cuts ---------------------------------------------
c ------------------------------------------------------------------------
      subroutine testcuts(n,accepted,pion)
      include 'phokhara.inc'
      double precision piplab,pimlab,m2,phot1,invmom(0:3),invm2
      integer n,pion,i,j
      logical accepted,accept(3:4)

      if (pion.eq.0) then 
         m2 = mmu*mmu
      else
         m2 = mpi*mpi
      endif

c --- one of the photons has energy > gmin ---
      accepted=(momenta(3,0).ge.gmin.or.momenta(4,0).ge.gmin)

c --- angular cuts on the pions/muons ---
      piplab=dacos(momenta(6,3)/dSqrt(momenta(6,0)**2-m2))*180.d0/pi
      pimlab=dacos(momenta(7,3)/dSqrt(momenta(7,0)**2-m2))*180.d0/pi
      accepted = (accepted.and.(piplab.ge.pi1cut.and.piplab.le.pi2cut))
      accepted = (accepted.and.(pimlab.ge.pi1cut.and.pimlab.le.pi2cut))
      
c --- invariant mass of the observed photon and the pions ---
      if (n.eq.2) then 
      do i = 3,4
         accept(i) = .false.
         phot1 = dacos(momenta(i,3)/momenta(i,0))*180.d0/pi

c       Selects only the i corresponding to the observed photon
         if (momenta(i,0).ge.gmin.and.
     &      (phot1.ge.phot1cut.and.phot1.le.phot2cut)) then 
            do j = 0,3
	        invmom(j) = momenta(i,j)+momenta(6,j)+momenta(7,j)
	    enddo 
	    invm2 = invmom(0)**2-invmom(1)**2-invmom(2)**2-invmom(3)**2
            accept(i) = (invm2.ge.q2min)
         endif
      enddo
      accepted = (accepted.and.(accept(3).or.accept(4)))
      endif
      end
      
      
c ========================================================================
c ------------------------------------------------------------------------
c --- vacuum polarization contribution -----------------------------------
      double precision function vacuumpolarization(qq)
      include 'phokhara.inc'
      double precision qq,deltalep

      vacuumpolarization = deltalep(qq,me) + deltalep(qq,mmu)
     &  + deltalep(qq,mtau)
      end

c ------------------------------------------------------------------------

      double precision function deltalep(qq,ml)
      include 'phokhara.inc'
      double precision qq,ml,bet,ml2,rat

      ml2 = 4.d0*ml**2
      rat = ml2/qq

      if(qq.gt.ml2)then
         bet = sqrt(1.d0-rat)
         deltalep = alpha/(3.d0*pi)*(1.d0/3.d0-(1.d0+ml2/2.d0)
     &            *(2.d0+bet*log((1.d0-bet)/(1.d0+bet))))
      elseif(qq.ge.0.d0)then
         bet = sqrt(rat-1.d0)
         deltalep = alpha/(3.d0*pi)*(1.d0/3.d0-(1.d0+ml2/2.d0)
     &            *(2.d0-bet*atan(1.d0/bet)))
      else
         write(6,*)'deltalep:qq<0 not supported'
         stop
      endif
      end
      
      
c ====================================================================
c ====================================================================
c --------------------------------------------------------------------
c --- Matrix element squared: contract Leptonic and Hadronic tensors -
c --------------------------------------------------------------------
      double precision function Matrix(Leptonic,Hadronic)
c Still don't understand the contraction performed
      include 'phokhara.inc'       
      double precision metric1,metric2,
     &  Leptonic(0:3,0:3),Hadronic(0:3,0:3)
      integer mu,nu
      
      Matrix = 0.d0
      do mu = 0,3
         metric1 = 1.d0
         if (mu.eq.0) metric1 = -1.d0
         do nu = 0,3
             metric2 = 1.d0
             if (nu.eq.0) metric2 = -1.d0               
             Matrix = Matrix + metric1*metric2*
     &          Leptonic(mu,nu)*Hadronic(mu,nu) 
         enddo
      enddo
      end


c ========================================================================
c --- Metric tensor ------------------------------------------------------
      double precision function metric(mu,nu)
      integer mu,nu
      if (mu.ne.nu) then 
         metric = 0.d0
      else 
         if (mu.eq.0) then 
             metric = 1.d0
	 else
	     metric = -1.d0
	 endif
      endif
      end 


c ========================================================================
c ------------------------------------------------------------------------
c     Leptonic tensor for the process                                     
c                                                                         
c         e^+(p1) + e^-(p2) ---->  photon^* + photon                      
c                                                                         
c     ISR at the NLO: virtual corrections + soft photons                  
c                                                                         
c     Sp : CMS energy squared
c     qq : virtuality of the photon coupled to the hadronic system        
c     cosphoton: cosinus of the hard photon-positron angle                
c          in the e^+ e^- CMS system                                       
c     w  : soft photon energy normalized to the CMS energy                 
c                                                                         
c ------------------------------------------------------------------------
c                                                                         
c         a00 : coefficient of g(mu,nu)                                   
c         a11 :      "      p1(mu)*p1(nu)/s                               
c         a22 :      "      p2(mu)*p2(nu)/s                               
c         a12 :      "      (p1(mu)*p2(nu)+p2(mu)*p1(nu))/s               
c                                                                         
c ------------------------------------------------------------------------
c (c) German Rodrigo 2000                                                 
c ------------------------------------------------------------------------
      subroutine LeptonicTensor1ISR(Sp,qq,cosphoton,nlo,Leptonic)
      include 'phokhara.inc'       
      double precision Sp,qq,cosphoton,Leptonic(0:3,0:3),
     &  a00,a11,a22,a12,dps,metric
      integer mu,nu,nlo

      dps = (1.d0-qq/Sp)/(32.d0*pi*pi)        ! Phase space factors

c --- ISR ---
      call Bornvirtualsoft(Sp,qq,cosphoton,nlo,a00,a11,a22,a12)
      do mu = 0,3
         do nu = 0,3
           Leptonic(mu,nu) = (a00*metric(mu,nu)+
     &        a11*momenta(1,mu)*momenta(1,nu)/Sp+
     &        a22*momenta(2,mu)*momenta(2,nu)/Sp+
     &        a12*(momenta(1,mu)*momenta(2,nu)+
     &             momenta(2,mu)*momenta(1,nu))/Sp)*dps
         enddo
      enddo
      return
      end
 
c ------------------------------------------------------------------------
 
      subroutine Bornvirtualsoft(Sp,qq,cosphoton,nlo,a00,a11,a22,a12)
      include 'phokhara.inc'
      double precision Sp,qq,cosphoton,a00,a11,a22,a12,all,
     &  m2,q2,uq2,b,x,y1,y2,globalfactor,
     &  a00NLO,a11NLO,a22NLO,a12NLO,api,
     &  t1,t2,t3,t4,t5,t6,t7,t8,t9,s1,t10,t11,t12,t13,
     &  soft,coll,extramass
      integer nlo
      complex*16 cdilog, dcmplx
       
      m2 = me*me/Sp
      q2 = qq/Sp
      uq2 = 1.d0-q2
      b = dSqrt(1.d0-4.d0*m2)
      x = b*cosphoton
      y1 = uq2*(1.d0-x)/2.d0
      y2 = uq2*(1.d0+x)/2.d0
      globalfactor = (4.d0*pi*alpha/Sp)**2/(q2*q2)

c --- LO ---

      a00 = ( 2.d0*m2*Q2*uq2*uq2/(y1*y2)-
     &  (2.d0*q2+y1*y1+y2*y2) )/(y1*y2)
      a11 = (8.d0*m2/y2-4.d0*q2/y1)/y2
      a22 = (8.d0*m2/y1-4.d0*q2/y2)/y1
      a12 = -8.d0*m2/(y1*y2)

c --- NLO ---

      if (nlo.ne.0) then    

         api = alpha/pi      
         t1 = dLog(m2)
         t2 = dLog(m2/q2)      
         t3 = pi*pi/3.d0
         t4 = dLog(y1/q2)
         t5 = dLog(y2/q2)
         t6 = dLog(q2)
         t7 = dreal(cdilog(dcmplx(1.d0-1.d0/q2)))
         t8 = dreal(cdilog(dcmplx(-y1/q2)))-t7+dLog(q2+y1)*dLog(y1/q2)
         t9 = dreal(cdilog(dcmplx(-y2/q2)))-t7+dLog(q2+y2)*dLog(y2/q2) 
         s1 = y1*y1+y2*y2      
         t10 = dLog(y1/m2)
	 t11 = dLog(y2/m2)
         t12 = dreal(cdilog(dcmplx(1.d0-y1/m2)))
	 t13 = dreal(cdilog(dcmplx(1.d0-y2/m2)))

      a00NLO = ( -q2*uq2/2.d0-y1*y2+
     &  y1/2.d0*(4.d0-y1-3.d0*(1.d0+q2)/(1.d0-y2))*t4+
     &  y2/2.d0*(4.d0-y2-3.d0*(1.d0+q2)/(1.d0-y1))*t5-
     &  (q2+2.d0*y1*y2/uq2)*t6-
     &  (1.d0+(1.d0-y2)**2+y1*q2/y2)*t8-
     &  (1.d0+(1.d0-y1)**2+y2*q2/y1)*t9 )/(y1*y2)

      a11NLO = ( (1.d0+q2)**2*(1.d0/(1.d0-y1)-1.d0/uq2)-
     &  4.d0*(1.d0-y2)*y1/uq2-q2*(1.d0+2.d0/y2)*t4-
     &  q2*(2.d0-3.d0*y1)*(1.d0-y2)**2/(y1*(1.d0-y1)**2)*t5-
     &  2.d0*q2/uq2*((1.d0-y2)*(1.d0/y2+q2/y1+2.d0*y1/uq2)+
     &  2.d0*q2/uq2)*t6-2.d0*q2*(1.d0+1.d0/(y2*y2))*t8-
     &  2.d0*q2*(3.d0+2.d0*q2/y1+q2*q2/(y1*y1))*t9 )/(y1*y2)

      a22NLO = ( (1.d0+q2)**2*(1.d0/(1.d0-y2)-1.d0/uq2)-
     &  4.d0*(1.d0-y1)*y2/uq2-q2*(1.d0+2.d0/y1)*t5-
     &  q2*(2.d0-3.d0*y2)*(1.d0-y1)**2/(y2*(1.d0-y2)**2)*t4-
     &  2.d0*q2/uq2*((1.d0-y1)*(1.d0/y1+q2/y2+2.d0*y2/uq2)+
     &  2.d0*q2/uq2)*t6-2.d0*q2*(1.d0+1.d0/(y1*y1))*t9-
     &  2.d0*q2*(3.d0+2.d0*q2/y2+q2*q2/(y2*y2))*t8 )/(y1*y2)

      a12NLO = ( q2/(1.d0-y1)+q2/(1.d0-y2)-(4.d0*q2+(y1-y2)**2)/uq2-
     &  2.d0*q2/(1.d0-y2)*(1.d0-y1+q2/y2-q2/(2.d0*(1.d0-y2)))*t4-
     &  2.d0*q2/(1.d0-y1)*(1.d0-y2+q2/y1-q2/(2.d0*(1.d0-y1)))*t5-
     &  2.d0*q2*(q2/(y1*y2)+(1.d0+q2-2.d0*y1*y2)/(uq2*uq2))*t6-
     &  2.d0*q2*(1.d0+q2/y2+q2/(y2*y2))*t8-
     &  2.d0*q2*(1.d0+q2/y1+q2/(y1*y1))*t9 )/(y1*y2) 

c --- NLO mass corrections       
      
      all = m2 * ( 2.d0*q2*(t6*(t1+1.5d0*t6-2.d0*dLog(y1))+
     &  2.d0*t7+t3/4.d0-t12/2.d0) - uq2*(1.d0-t10+
     &  (m2/y1)*(t12-t3/2.d0)) )/(y1*y1)+ extramass(m2,y1,q2)
     &    + m2 * ( 2.d0*q2*(t6*(t1+1.5d0*t6-2.d0*dLog(y2))+
     &  2.d0*t7+t3/4.d0-t13/2.d0) - uq2*(1.d0-t11+
     &  (m2/y2)*(t13-t3/2.d0)) )/(y2*y2)+ extramass(m2,y2,q2)

      a00NLO = a00NLO - all 
      a11NLO = a11NLO - 4.d0*q2*all/uq2**2 
      a22NLO = a22NLO - 4.d0*q2*all/uq2**2
      a12NLO = a12NLO - 4.d0*q2*all/uq2**2

c --- soft and collinear logs 

      soft = -dLog(4.d0*w*w)*(1.d0+t1)
      coll = (-1.5d0*t2-2.d0+t3)

c --- final result

      a00 = a00 + api*(a00*(soft+coll)+a00NLO)
      a11 = a11 + api*(a11*(soft+coll)+a11NLO)
      a22 = a22 + api*(a22*(soft+coll)+a22NLO)
      a12 = a12 + api*(a12*(soft+coll)+a12NLO)                   
      endif 

      a00 = globalfactor * a00
      a11 = globalfactor * a11
      a22 = globalfactor * a22
      a12 = globalfactor * a12
      return 
      end

c ---

      double precision function extramass(m2,y,q2)
      double precision m2,y,q2
      extramass = (q2+(1.d0-3.d0*q2)*dLog(y/m2))/(2.d0*y)
      if ((m2-y)**2.le.m2*m2/4.d0) then 
         do n = 0,4
            extramass = extramass - .5d0*(q2/(dble(n)+2.d0)+
     &        (1.d0-3.d0*q2)/(dble(n)+1.d0))*(m2-y)**n/(m2**(n+1))
         enddo
      else
         extramass = extramass +
     &      (q2+(1.d0-3.d0*q2)*dLog(y/m2))/(2.d0*(m2-y)) +
     &      m2*q2/(2.d0*(m2-y)**2)*dLog(y/m2)
      endif
      end


c ========================================================================
c ------------------------------------------------------------------------
c     Leptonic tensor for the process                                     
c                                                                         
c         e^+(p1) + e^-(p2) ----> photon^* + photon(1) + photon(2)        
c                                                                         
c     Sp : CMS energy squared
c     qq : virtuality of the photon coupled to the hadronic system        
c     w1,w2: photon energies normalized to the CMS energy                    
c     cos1: cosinus of the photon1-positron angle in the CMS system         
c     cos2: cosinus of the photon2-positron angle in the CMS system         
c     phi1, phi2: azimutal photon angles                                  
c                                                                         
c ------------------------------------------------------------------------
c                                                                         
c         b00 : coefficient of g(mu,nu)                                   
c         b11 :      "      p1(mu)*p1(nu)/s                               
c         b22 :      "      p2(mu)*p2(nu)/s                               
c         b12 :      "      (p1(mu)*p2(nu)+p2(mu)*p1(nu))/s               
c         b33 :      "      k1(mu)*k1(nu)                                 
c         b13 :      "      (k1(mu)*p1(nu)+p1(mu)*k1(nu))/s               
c         b23 :      "      (k1(mu)*p2(nu)+p2(mu)*k1(nu))/s               
c                                                                         
c ------------------------------------------------------------------------
c (c) German Rodrigo 2000                                                 
c ------------------------------------------------------------------------
      subroutine LeptonicTensor2ISR(Sp,qq,w1,w2,
     &  cos1,phi1,cos2,phi2,Leptonic)
      include 'phokhara.inc'       
      double precision Sp,qq,w1,w2,cos1,phi1,cos2,phi2,
     &  Leptonic(0:3,0:3),metric,b00,b11,b22,b12,b33,b13,b23,dps
      integer mu,nu

      dps = Sp/(4.d0*(2.d0*pi)**5)            ! Phase space factors

c --- ISR ---
      call realhard(Sp,qq,w1,w2,cos1,phi1,cos2,phi2,
     &              b00,b11,b22,b12,b33,b13,b23)
      do mu = 0,3
         do nu = 0,3
            Leptonic(mu,nu) = (b00*metric(mu,nu)+
     &         b11*momenta(1,mu)*momenta(1,nu)/Sp+
     &         b22*momenta(2,mu)*momenta(2,nu)/Sp+
     &         b12*(momenta(1,mu)*momenta(2,nu)+
     &              momenta(2,mu)*momenta(1,nu))/Sp+
     &         b33*momenta(3,mu)*momenta(3,nu)/Sp+
     &         b13*(momenta(1,mu)*momenta(3,nu)+
     &              momenta(3,mu)*momenta(1,nu))/Sp+
     &         b23*(momenta(2,mu)*momenta(3,nu)+
     &              momenta(3,mu)*momenta(2,nu))/Sp)*dps
         enddo
      enddo
      return
      end

c ------------------------------------------------------------------------

      subroutine realhard(Sp,qq,w1,w2,cos1,phi1,cos2,phi2,
     &  b00,b11,b22,b12,b33,b13,b23)
      implicit double precision (s-t)      
      double precision pi,gev2nbarn,gev2pbarn,alpha,me,mmu,mtau,mtop,
     & mrho,mpi,grhoee,gammarho,mrhol,grhol,momega,gomega,al,be,
     & S,ebeam,momenta(10,0:3)
      common /ctes/    pi,gev2nbarn,gev2pbarn,alpha,me,mmu,mtau,mtop,
     & mrho,mpi,grhoee,gammarho,mrhol,grhol,momega,gomega,al,be,
     & S,ebeam,momenta
      double precision Sp,qq,w1,w2,cos1,cos2,phi1,phi2,
     &  b00,b11,b22,b12,b33,b13,b23,
     &  m2,q2,b,x1,x2,y11,y12,y21,y22,globalfactor

      m2 = me*me/Sp
      q2 = qq/Sp
      b = dSqrt(1.d0-4.d0*m2)
      x1 = b*cos1
      x2 = b*cos2      
      y11 = w1*(1.d0-x1)
      y12 = w1*(1.d0+x1)
      y21 = w2*(1.d0-x2)
      y22 = w2*(1.d0+x2)
      globalfactor = (4.d0*pi*alpha/Sp)**3/(y11*y12*y21*y22*q2*q2)

      t1 = y12**2
      t2 = y11*t1
      t4 = y12*y21
      t5 = y22**2
      t6 = t4*t5
      t7 = y11**2
      t8 = t7*y12
      t14 = y11*y21
      t16 = t7*y21
      t18 = y11*y12
      t19 = y21*t5
      t24 = y21*y22
      t27 = t14*y22
      t29 = y21**2
      t30 = t29*y22
      t32 = t7*y11
      t41 = Q2**2
      t44 = (-1.D0+y21+y22)**2
      t55 = t1*y12
      t70 = -t2*y22-t6-t8*y21+t8*t5-t14*t5+t16*t5+3.D0*t18*t19-t8*y22+3.
     #D0*t2*t24-4.D0*t18*t24+t27+3.D0*t8*t24+3.D0*t18*t30+2.D0*m2*(t32*y
     #12+t7*(2.D0*t1+(2.D0*Q2-y21)*y22+y12*(-2.D0+2.D0*Q2+y22))+y21*(t41
     #*y22+y22*(-t1+t4+t44)+2.D0*Q2*(t1+y22*(-3.D0+y21+y22)+y12*(-2.D0+y
     #21+2.D0*y22)))+y11*(t55+t1*(-2.D0+2.D0*Q2+y21)+y22*(t24+2.D0*Q2*(-
     #2.D0+2.D0*y21+y22))+y12*(1.D0+t41-t29-t5+Q2*(-6.D0+4.D0*y21+4.D0*y
     #22))))
      t73 = y12*t29
      t74 = t73*y22
      t75 = m2**2
      t76 = t75*Q2
      t77 = y11+y21
      t78 = y12+y22
      t84 = t1*y21
      t87 = t18*y21
      t88 = t18*y22
      t89 = t4*y22
      t93 = -1.D0+y22
      t96 = -1.D0+y21
      t97 = t1*t96
      t110 = t1*t29*y22-t74+8.D0*t76*t77*t78-t2*y21-t18*t5-t18*t29-t84*y
     #22+t2*t29+t87+t88+t89-t16*y22-y11*t29*y22+Q2*(2.D0*t7*t93*y22+y21*
     #(2.D0*t97+y12*(2.D0-2.D0*y21-3.D0*y22)+2.D0*y22)-y11*(y22*(-2.D0+3
     #.D0*y21+2.D0*y22)+y12*(-2.D0+3.D0*y21+3.D0*y22)))
      t112 = -1.D0+Q2+y11+y21
      t113 = 1.d0/t112
      t115 = -1.D0+Q2+y12+y22
      t116 = 1.d0/t115
      t118 = m2*y11
      t119 = t118*y12
      t120 = m2*Q2
      t124 = 1.d0/y21
      t127 = m2*y12
      t129 = y22*t115
      t132 = 1.d0/y11
      t135 = y12*y22
      t136 = -1.D0+y12+y22
      t139 = y11*(-1.D0+y12+2.D0*y22)
      t148 = -1.D0+y11+y12+y21+y22
      t157 = t115**2
      t160 = y11*y22
      t163 = y11*t5
      t164 = t7*y22
      t165 = -2.D0*(t70+t110)*t113*t116+2.D0*t4+(-4.D0*t119*y22*(2.D0*t1
     #20+y12*t115)*t124-4.D0*t127*y21*y22*(2.D0*t120+t129)*t132+2.D0*t13
     #5*(-8.D0*t76-2.D0*t41-t136*(-2.D0+t7-2.D0*y12*t96+y21+t29+2.D0*y22
     #-t24-t139)+Q2*(2.D0*y12*(-2.D0+y21)+(-4.D0+y21)*t93+t139)+2.D0*m2*
     #(t41+t136*t148+Q2*(-2.D0+3.D0*y11+2.D0*y12+3.D0*y21+2.D0*y22))))/t
     #157-2.D0*Q2*(2.D0+t4+t160)-4.D0*t135-4.D0*t14-16.D0*t76+2.D0*t160-
     #2.D0*t27-2.D0*t163-2.D0*t84-2.D0*t164-2.D0*t73
      t169 = 1.d0/y22
      t173 = y21*t112
      t176 = 1.d0/y12
      t179 = -2.D0*t24-3.D0*t4-2.D0*t7+y12-2.D0*t29-4.D0*t14-8.D0*t76-3.
     #D0*t160+3.D0*t27-t163-t84-2.D0+2.D0*t164+y22
      t194 = 2.D0*t73+4.D0*y21+3.D0*t87+t1+4.D0*y11+t5-2.D0*t41+t8-t2+t3
     #0-t19-2.D0*t18+Q2*(y12*(-1.D0+2.D0*y21)+t96*(-4.D0+y22)+y11*(-4.D0
     #+y12+2.D0*y22))+2.D0*m2*(t41+(-1.D0+y11+y21)*t148+Q2*(-2.D0+2.D0*y
     #11+3.D0*y12+2.D0*y21+3.D0*y22))
      t198 = t112**2
      t201 = t120*y12
      t206 = t120*y11
      t217 = y12*t93
      t227 = y12*(-1.D0+Q2+y21+y22-4.D0*t24)
      t234 = t93*y22
      t259 = Q2*t1
      t261 = 2.D0*t8*t93+16.D0*t76*t78-2.D0*Q2*y21*y22-2.D0*t6+4.D0*Q2*y
     #12*y22+2.D0*t24-4.D0*t135-2.D0*t84+2.D0*t74+4.D0*t259*y21-2.D0*t89
     #+2.D0*t1
      t263 = t1*y22
      t265 = y12*t5
      t275 = t5*y22
      t288 = -t1+t55-t4+t73+2.D0*t135-2.D0*t263-t24+2.D0*t89+t30-t5-2.D0
     #*t265+t275+t7*t78+Q2*(t1+y12*(-2.D0+y11+y21-2.D0*y22)+y22*(-2.D0+y
     #11+y21+y22))+y11*(t96*y22+y12*(-1.D0+y21+2.D0*y22))
      t298 = 2.D0*t5-2.D0*t55-2.D0*t30-2.D0*t259+2.D0*t55*y21+2.D0*t263-
     #2.D0*Q2*t5+2.D0*t265+2.D0*y11*(t1*(y21-y22)-t227+y22*(y21*t93+y22*
     #(-1.D0+2.D0*Q2+y22)))-2.D0*t275+4.D0*m2*t288+8.D0*t201*(-1.D0+2.D0
     #*m2+y12)*y21*t132+8.D0*t206*y22*(-1.D0+2.D0*m2+y22)*t124
      t303 = 1.D0+t41+t7-2.D0*y12+t1-2.D0*y21+4.D0*t4+t29-2.D0*y22+2.D0*
     #t24+t5+2.D0*Q2*(-3.D0+y11+y12+y21+y22)+2.D0*t139
      s1 = -2.D0*t87-2.D0*t88-2.D0*t89+(-4.D0*t119*y21*(2.D0*t120+y11*t1
     #12)*t169-4.D0*t118*y21*(2.D0*t120+t173)*y22*t176+2.D0*t14*(t179+t1
     #94))/t198+(8.D0*t201*y21*(-1.D0+2.D0*m2+y21)*t169+8.D0*t206*(-1.D0
     #+2.D0*m2+y11)*y22*t176+16.D0*t76*t77+2.D0*t32*t93+2.D0*y21*((-1.D0
     #+y12)*t29+y21*(1.D0+Q2*(-1.D0+2.D0*y12)+t217)-t129)+2.D0*t7*(1.D0+
     #y21-t4-y22+t135+Q2*(-1.D0+2.D0*y22))+2.D0*y11*(t97-t227+y21*(-2.D0
     #+2.D0*Q2+y21-y22-t24+t5))+4.D0*m2*(t32-t7*(1.D0+2.D0*y21)+y21*(t1-
     #y21+t29+t217+t234)+Q2*(t7+y11*(-2.D0+y12-2.D0*y21+y22)+y21*(-2.D0+
     #y12+y21+y22))+y11*(t1-2.D0*t29+t234+2.D0*y21*(1.D0+y22)+y12*(-1.D0
     #+2.D0*y21+y22))))*t113+(t261+t298)*t116
      t323 = s1-4.D0*m2*t303-8.D0*t76*y12*y21*t132*t169-8.D0*t76*y11*y22
     #*t176*t124-4.D0*t127*t173*t169-4.D0*t118*t112*y22*t176-4.D0*t127*y
     #21*t115*t132-4.D0*t118*t129*t124
      b00 = t165+t323

      t1 = m2*y12
      t2 = 1.d0/y11
      t3 = y21*t2
      t5 = m2**2
      t7 = 1.d0/y22
      t10 = y12*y22
      t14 = -3.D0+Q2+2.D0*y12+y21+3.D0*y22
      t17 = y12*y21
      t18 = t17*t7
      t20 = m2*y11
      t25 = y21*y22
      t26 = 1.d0/y12
      t30 = y12**2
      t31 = -1.D0+Q2+y11+y21
      t41 = t31**2
      t44 = -1.D0+y12+y22
      t50 = Q2*y22
      t51 = -1.D0+y21+y22
      t57 = 1.d0/(-1.D0+Q2+y12+y22)
      t75 = y22**2
      t85 = 1.d0/t31
      s1 = -16.D0*t1*t3-32.D0*t5*y12*t3*t7-8.D0*Q2-8.D0*y12-8.D0*y21
      s2 = s1+8.D0*t10-8.D0*m2*t14
      s3 = s2+(-32.D0*t5*y11*t18-16.D0*t20*(1.D0+2.D0*m2-Q2-y11-y21)*t25
     #*t26-8.D0*y11*y21*(8.D0*t5-t30+2.D0*y12*t31+(-1.D0+Q2+y11+y21-y22)
     #*y22-4.D0*m2*(y12+y22)))/t41
      b11 = s3+(16.D0*t1*y21*t44*t2-8.D0*y12*(Q2+m2*(1.D0+Q2+2.D0*y12-3.
     #D0*y21-y22)-t50-t44*t51))*t57+8.D0*(y11*(t30*(1.D0+2.D0*m2-y21-y22
     #)+t50*(-1.D0+2.D0*m2+y22)+y12*(m2*(1.D0+Q2-3.D0*y21-y22)-(-1.D0+y2
     #2)*t51))+y21*(-t30*(-1.D0+2.D0*m2-2.D0*Q2+y21+2.D0*y22)+y12*(-1.D0
     #+y21+m2*(2.D0-4.D0*y22)+4.D0*y22-t25-2.D0*t75)+y22*(-1.D0+y22+m2*(
     #-1.D0+3.D0*Q2+3.D0*y21+y22))))*t85*t57+(32.D0*m2*(2.D0*m2-Q2)*t18-
     #16.D0*t20*(Q2+y21)*y22*t26+64.D0*t5*y21+8.D0*y21*(-1.D0+Q2-Q2*y12-
     #t30+y21-y22)+8.D0*y11*(y12+y21+2.D0*t17+y22-t10+t25-t75)+8.D0*m2*(
     #-2.D0*y21*(1.D0+Q2-y12+y21-y22)+y11*t14))*t85

      t1 = m2*y11
      t2 = 1.d0/y12
      t3 = y22*t2
      t5 = m2**2
      t6 = t5*y11
      t7 = 1.d0/y21
      t14 = -3.D0+Q2+2.D0*y11+3.D0*y21+y22
      t16 = y12*y22
      t19 = m2*y12
      t25 = 1.d0/y11
      t28 = y11**2
      t31 = -1.D0+Q2+y12+y22
      t39 = t31**2
      t42 = -1.D0+y11+y21
      t49 = -1.D0+y21+y22
      t55 = 1.d0/(-1.D0+Q2+y11+y21)
      t66 = y12*y21
      t67 = y21**2
      t70 = Q2*y22
      t72 = y21*y22
      t74 = y22**2
      t82 = 32.D0*m2*(2.D0*m2-Q2)*y11*y22*t7-16.D0*t19*y21*(Q2+y22)*t25+
     #8.D0*t66-8.D0*y12*t67-8.D0*y22+64.D0*t5*y22+8.D0*t70-8.D0*t28*y22+
     #8.D0*t16-8.D0*t72+8.D0*t66*y22+8.D0*t74+8.D0*y11*(y12-t66-t70+2.D0
     #*t16)+8.D0*m2*(-2.D0*y22*(1.D0+Q2-y11-y21+y22)+y12*t14)
      t83 = 1.d0/t31
      b22 = -16.D0*t1*t3-32.D0*t6*t3*t7-8.D0*Q2-8.D0*y11+8.D0*y11*y21-8.
     #D0*y22-8.D0*m2*t14+(-32.D0*t6*t16*t7-16.D0*t19*y21*(1.D0+2.D0*m2-Q
     #2-y12-y22)*y22*t25-8.D0*t16*(8.D0*t5-t28-4.D0*m2*(y11+y21)+2.D0*y1
     #1*t31+y21*(-1.D0+Q2+y12-y21+y22)))/t39+(16.D0*t1*t42*y22*t2-8.D0*y
     #11*(Q2-Q2*y21+m2*(1.D0+Q2+2.D0*y11-y21-3.D0*y22)-t42*t49))*t55+t82
     #*t83+8.D0*(t28*(y12*(1.D0+2.D0*m2-y21-y22)-y22*(-1.D0+2.D0*m2-2.D0
     #*Q2+2.D0*y21+y22))+y11*(y22*(-1.D0+m2*(2.D0-4.D0*y21)+4.D0*y21-2.D
     #0*t67+y22-t72)+y12*(m2*(1.D0+Q2-y21-3.D0*y22)-(-1.D0+y21)*t49))+y2
     #1*(Q2*(y12*(-1.D0+2.D0*m2+y21)+3.D0*m2*y22)+y22*(-1.D0+y21+m2*(-1.
     #D0+y21+3.D0*y22))))*t55*t83

      t1 = m2*y12
      t2 = 1.d0/y11
      t5 = m2*y11
      t6 = 1.d0/y12
      t9 = m2**2
      t15 = y11*y21
      t16 = y12*y22
      t21 = y12**2
      t30 = y21*y22
      t31 = y22**2
      t32 = -1.D0+8.D0*t9+Q2-y11+t15-2.D0*m2*(y11+3.D0*y21-3.D0*y22)-2.D
     #0*y22+2.D0*Q2*y22+2.D0*y11*y22+t30+2.D0*t31
      t43 = 1.d0/(-1.D0+Q2+y12+y22)
      t47 = y12+y22
      t49 = y11**2
      t52 = -1.D0+y21+y22
      t62 = -1.D0+y21+y22-t30+Q2*(1.D0+y21+y22)
      t67 = y21**2
      t94 = 1.d0/(-1.D0+Q2+y11+y21)
      t113 = -1.D0+8.D0*t9+Q2-y12-2.D0*y21+2.D0*Q2*y21+2.D0*y12*y21+2.D0
     #*t67+t16+t30-2.D0*m2*(y12-3.D0*y21+3.D0*y22)
      s1 = -8.D0*t1*y21*t2-8.D0*t5*y22*t6-8.D0+32.D0*t9+8.D0*Q2+4.D0*y11
     #-8.D0*m2*(1.D0+Q2-y11-y12)
      s2 = s1+4.D0*y12+4.D0*y21+8.D0*t15
      s3 = s2+4.D0*y22+8.D0*t16
      s4 = s3+(-8.D0*t1*y21*(1.D0+4.D0*m2-Q2-y12-y22)*t2-4.D0*t21*(1.D0+
     #2.D0*m2-y21+2.D0*y22)-4.D0*y12*t32-4.D0*y22*(1.D0+Q2*(-1.D0+y11)+y
     #21-y22+2.D0*m2*(-1.D0+Q2-y11-y21+y22)))*t43
      b12 = s4+4.D0*(8.D0*t9*(y11+y21)*t47+t49*((-1.D0+y22)*y22+y12*t52)
     #+y21*(t21*(-1.D0+y21)-y22*(y21+y22)-y12*t62)+y11*(t21*t52+y12*(t67
     #-y21*(2.D0+Q2-6.D0*y22)+y22*(-2.D0-Q2+y22))-y22*t62)-2.D0*m2*(t49*
     #t47+y21*(t21+(1.D0+Q2)*y22+y12*(-1.D0+3.D0*Q2+2.D0*y22))+y11*(t21+
     #(-1.D0+3.D0*Q2+2.D0*y21)*y22+y12*(-1.D0+3.D0*Q2+2.D0*y21+2.D0*y22)
     #)))*t94*t43+(-8.D0*t5*(1.D0+4.D0*m2-Q2-y11-y21)*y22*t6-4.D0*t49*(1
     #.D0+2.D0*m2+2.D0*y21-y22)-4.D0*y21*(1.D0+Q2*(-1.D0+y12)-y21+2.D0*m
     #2*(-1.D0+Q2-y12+y21-y22)+y22)-4.D0*y11*t113)*t94

      t1 = m2*y12
      t2 = 1.d0/y11
      t3 = y21*t2
      t5 = 1.d0/y22
      t8 = m2**2
      t15 = m2*y11
      t16 = 1.d0/y12
      t17 = y22*t16
      t19 = 1.d0/y21
      t27 = y12*y21
      t28 = y11*y22
      t32 = y11+y21
      t34 = y12+y22
      t37 = t27+t28
      t38 = -1.D0+y11+y21
      t39 = -1.D0+y12+y22
      t45 = 1.d0/(-1.D0+Q2+y11+y21)
      t48 = 1.d0/(-1.D0+Q2+y12+y22)
      b33 = -16.D0*t1*t3-16.D0*t1*y21*t5-32.D0*t8*y12*t3*t5-16.D0*t15*t1
     #7-16.D0*t15*y22*t19-32.D0*t8*y11*t17*t19-16.D0+64.D0*t8-8.D0*t27-8
     #.D0*t28-16.D0*m2*(-2.D0+Q2+y11+y12+y21+y22)-8.D0*(2.D0*m2*Q2*t32*t
     #34+t37*(-Q2+t38*t39))*t45*t48+(16.D0*t1*y21*t38*t5+16.D0*t15*t38*y
     #22*t16-8.D0*m2*t32*(2.D0+y11-3.D0*y12+y21-3.D0*y22)+8.D0*t38*t37)*
     #t45+(16.D0*t1*y21*t39*t2+16.D0*t15*y22*t39*t19-8.D0*m2*t34*(2.D0-3
     #.D0*y11+y12-3.D0*y21+y22)+8.D0*t39*t37)*t48

      t1 = m2*y12
      t2 = 1.d0/y11
      t3 = y21*t2
      t5 = 1.d0/y22
      t8 = m2**2
      t15 = m2*y11
      t16 = 1.d0/y12
      t19 = y12*y21
      t20 = -1.D0+y22
      t21 = y11*t20
      t30 = y12**2
      t38 = y21*y22
      t45 = 1.d0/(-1.D0+Q2+y12+y22)
      t48 = y22**2
      t50 = y21**2
      t64 = y11**2
      t82 = 1.d0/(-1.D0+Q2+y11+y21)
      t96 = y12*y22
      s1 = 16.D0*t1*t3+8.D0*t1*y21*t5+32.D0*t8*y12*t3*t5+8.D0*t15*y22*t1
     #6+8.D0-32.D0*t8+4.D0*y12
      s2 = s1+4.D0*y21+4.D0*t19+4.D0*t21
      s3 = s2-4.D0*y22+4.D0*m2*(-4.D0+2.D0*Q2+y11+y12+3.D0*y21+3.D0*y22)
      s4 = s3+(-16.D0*t1*y21*(-1.D0+y12+y22)*t2+4.D0*t30*(1.D0+3.D0*m2-2
     #.D0*y21)-4.D0*y22*(-1.D0+Q2+y22+m2*(-1.D0+Q2+3.D0*y21+y22))+4.D0*y
     #12*(-1.D0+Q2+2.D0*y21-2.D0*t38+m2*(3.D0+Q2-3.D0*y11-6.D0*y21+2.D0*
     #y22)))*t45
      b13 = s4+4.D0*(Q2*y11*t48+t30*(2.D0*t50-y11*y22+y21*(-2.D0-Q2+2.D0
     #*y11+y22))+y12*(2.D0*t50*t20-t21*y22+y21*(-2.D0*Q2+t20*(-2.D0+2.D0
     #*y11+y22)))+m2*(3.D0*t64*y12+y21*(2.D0*t30+y22*(2.D0-3.D0*y21+y22)
     #+y12*(-2.D0+4.D0*Q2+3.D0*y22))-y11*(t30+y22*(-2.D0+3.D0*y21+2.D0*y
     #22)+y12*(2.D0-4.D0*Q2-3.D0*y21+3.D0*y22))))*t82*t45+(-8.D0*t1*y21*
     #(-1.D0+4.D0*m2-Q2+y11+y21)*t5+8.D0*t15*(1.D0+4.D0*m2-Q2-y11-y21)*y
     #22*t16-4.D0*t64*(-1.D0+m2+y22)-4.D0*y21*(-1.D0+8.D0*t8+Q2-2.D0*y12
     #-Q2*y12+y21+t19-2.D0*y22+t96+m2*(-1.D0-3.D0*Q2+2.D0*y12-3.D0*y21+3
     #.D0*y22))-4.D0*y11*(1.D0-8.D0*t8+2.D0*y12+t19+Q2*t20-t96+t38+m2*(-
     #3.D0+3.D0*Q2+3.D0*y12-2.D0*y21+4.D0*y22)))*t82

      t1 = m2*y12
      t2 = 1.d0/y11
      t5 = m2*y11
      t6 = 1.d0/y12
      t7 = y22*t6
      t9 = 1.d0/y21
      t15 = m2**2
      t20 = y11*y22
      t34 = y12**2
      t38 = y11*y21
      t43 = -1.D0+y21
      t47 = y21*y22
      t52 = 1.d0/(-1.D0+Q2+y12+y22)
      t55 = y21**2
      t57 = y11**2
      t58 = y21-2.D0*y22
      t90 = 1.d0/(-1.D0+Q2+y11+y21)
      s1 = 8.D0*t1*y21*t2+16.D0*t5*t7+8.D0*t5*y22*t9+32.D0*t15*y11*t7*t9
     #+8.D0-32.D0*t15+4.D0*y11-4.D0*y12
      s2 = s1-4.D0*y21+4.D0*y12*y21+4.D0*y22
      s3 = s2+4.D0*t20+4.D0*m2*(-4.D0+2.D0*Q2+y11+y12+3.D0*y21+3.D0*y22)
      s4 = s3+(8.D0*t1*y21*(1.D0+4.D0*m2-Q2-y12-y22)*t2-8.D0*t5*y22*(-1.
     #D0+4.D0*m2-Q2+y12+y22)*t9-4.D0*t34*(-1.D0+m2+y21)-4.D0*y22*(-1.D0+
     #8.D0*t15+Q2-2.D0*y11-Q2*y11-2.D0*y21+t38+m2*(-1.D0-3.D0*Q2+2.D0*y1
     #1+3.D0*y21-3.D0*y22)+y22+t20)-4.D0*y12*(1.D0-8.D0*t15+2.D0*y11+Q2*
     #t43-t38+m2*(-3.D0+3.D0*Q2+3.D0*y11+4.D0*y21-2.D0*y22)+t20+t47))*t5
     #2
      b23 = s4+4.D0*(Q2*y12*t55+t57*(-y12*t58+y22*(-2.D0-Q2+y21+2.D0*y22
     #))+y11*(-y12*t43*t58+y22*(-2.D0*Q2+t43*(-2.D0+y21+2.D0*y22)))+m2*(
     #-t57*(y12-2.D0*y22)+y21*(y12*(2.D0-2.D0*y21-3.D0*y22)+(2.D0+y21-3.
     #D0*y22)*y22)+y11*(3.D0*t34+(-2.D0+4.D0*Q2+3.D0*y21)*y22+y12*(-2.D0
     #+4.D0*Q2-3.D0*y21+3.D0*y22))))*t90*t52+(-16.D0*t5*(-1.D0+y11+y21)*
     #y22*t6+4.D0*t57*(1.D0+3.D0*m2-2.D0*y22)+4.D0*y11*(-1.D0+Q2+m2*(3.D
     #0+Q2-3.D0*y12+2.D0*y21-6.D0*y22)+2.D0*y22-2.D0*t47)-4.D0*y21*(-1.D
     #0+Q2+y21+m2*(-1.D0+Q2+y21+3.D0*y22)))*t90
      
      b00 = globalfactor * b00
      b11 = globalfactor * b11
      b22 = globalfactor * b22
      b12 = globalfactor * b12
      b33 = globalfactor * b33
      b13 = globalfactor * b13
      b23 = globalfactor * b23      

      return 
      end



c ========================================================================
c ------------------------------------------------------------------------
c     Hadronic Tensor for the process                                     
c                                                                         
c         photon^* (Q2) ----->  mu^+ mu^-, pi^+ pi^-                      
c ------------------------------------------------------------------------

      subroutine HadronicTensorISR(qq,pion,Hadronic)
      include 'phokhara.inc'       
      double precision qq,Hadronic(0:3,0:3),pionFF,
     &   metric,PionFormFactor2,dps
      integer mu,nu,pion

c --- muons ---      
      if (pion.eq.0)  then 
      pionFF = 16.d0*pi*alpha
      dps =  dSqrt(1.d0-4.d0*mmu*mmu/qq)/(32.d0*pi*pi)  ! Phase space factors
      do mu = 0,3
         do nu = 0,3
           Hadronic(mu,nu) = pionFF*(momenta(6,mu)*momenta(7,nu)+
     &        momenta(7,mu)*momenta(6,nu)-qq/2.d0*metric(mu,nu))*dps
         enddo
      enddo
     
      else
      
c --- pions ---        
      pionFF = 4.d0*pi*alpha*PionFormFactor2(qq,qq)
      dps =  dSqrt(1.d0-4.d0*mpi*mpi/qq)/(32.d0*pi*pi)  ! Phase space factors
      do mu = 0,3
         do nu = 0,3
            Hadronic(mu,nu) = pionFF*(momenta(6,mu)-momenta(7,mu))*
     &         (momenta(6,nu)-momenta(7,nu))*dps 
         enddo
      enddo
      endif
      return
      end
      
c --------------------------------------------------------------------

      double precision function PionFormFactor2(a,b)
      include 'phokhara.inc'       
      double precision a,b
      complex*16 f1,f2,BW
       
      f1 = (BW(mrho,gammarho,a,1)*(1.D0+al*BW(momega,gomega,a,1))/
     &     (1.d0+al)+be*BW(mrhol,grhol,a,1))/(1.d0+be)
      f2 = (BW(mrho,gammarho,b,2)*(1.D0+al*BW(momega,gomega,b,2))/
     &     (1.d0+al)+be*BW(mrhol,grhol,b,2))/(1.d0+be)
      PionFormFactor2 = dreal(f1*f2)
      return
      end

c --------------------------------------------------------------------

      complex*16 function BW(m,breite,x,k)
      include 'phokhara.inc'       
      integer k
      double precision m,breite,x,g
      complex *16 i

      if(breite.eq.gomega)then 
         g=breite
      else
         g=breite*m*m/x*(x-4.d0*mpi*mpi)**(1.5d0)/
     &     (m*m-4.d0*mpi*mpi)**(1.5d0)
      endif
      i=(0.d0,1.d0)
      if(k.eq.1)then
         BW=m*m/(m*m-x-i*sqrt(x)*g)
      else
         BW=m*m/(m*m-x+i*sqrt(x)*g)
      endif
      end


c ========================================================================
c ---- complex dilogarithm -----------------------------------------------                       
c ------------------------------------------------------------------------
       complex function cdilog*16(z)                                    
       complex*16 z,zl,coef,dilog1,u,caux                               
       real*8 pi,sign                                                   
       integer n 
       pi=3.141592653589793238462643d0                                  
       zl=z                                                             
       dilog1=dcmplx(pi**2/6.d0)                                        
       if(dreal(zl).eq.1.and.dimag(zl).eq.0.) then                      
          cdilog=dilog1                                                    
          return                                                           
       else if (cdabs(zl).lt.1.d-2) then   
          n=-40./dlog(cdabs(zl))                                           
          caux=(0.d0,0.d0)                                                 
          do i=1,n                                                         
             caux=caux+zl**i/dble(i**2)                                       
          enddo                                                           
          cdilog=caux                                                      
          return                                                           
       else if(cdabs(zl).lt.1.) then                                    
          sign=1.d0                                                       
          coef=dcmplx(dble(0.))                                           
       else                                                            
          coef=-cdlog(-zl)**2/2.d0-dilog1                                 
          sign=-1.d0                                                      
          zl=1.d0/zl                                                      
       endif                                                          
       if(dreal(zl).gt.0.5) then                   
          coef=coef+sign*(dilog1-cdlog(zl)*cdlog(1.d0-zl))                
          sign=-sign                                                      
          zl=1.d0-zl                                                      
       else   
       endif  
       u=-cdlog(1.d0-zl)                                               
       cdilog=u-u**2/4.d0+u**3/36.d0-u**5/3600.d0+u**7/211680.d0       
     &  -u**9/10886400.d0+u**11*5.d0/2634508800.d0                     
       cdilog=cdilog-u**13*691.d0/2730.d0/6227020800.d0                
       cdilog=cdilog+u**15*7.d0/6.d0/1.307674368d12                    
       cdilog=cdilog-u**17*3617.d0/510.d0/3.5568742810d14              
       cdilog=cdilog+u**19*43867.d0/798.d0/1.2164510041d17              
       cdilog=cdilog-u**21*174611.d0/330.d0/5.1090942172d19            
       cdilog=sign*cdilog+coef                                         
       return                                                          
       end                                                             


c ========================================================================
c ------------------------------------------------------------------------
c --- all about the histogrammes -----------------------------------------

      subroutine addiere(wgt,qq,i)
      include 'phokhara.inc'
      double precision wgt,qq
      integer i

      if (i.eq.1) call addhisto(1,qq,wgt)      ! one photon events
      if (i.eq.2) call addhisto(2,qq,wgt)      ! two photon events
      return
      end

c ------------------------------------------------------
c     create histograms                                 
c ------------------------------------------------------
      subroutine inithisto
      include 'phokhara.inc'
      integer i
c --- PAW common variables ---
      double precision h
      integer nwpawc
      parameter (nwpawc = 30000)
      common/pawc/h(nwpawc)
      call hlimit(nwpawc)
c --- book the histograms ----
      do i=1,10
        call hbook1(i,title(i),bins(i),sngl(xlow(i)),sngl(xup(i)),0.)
        call hbarx(i)
      enddo
      return 
      end

c ------------------------------------------------------
c     add value to histo i at x                         
c ------------------------------------------------------
      subroutine addhisto(i,x,value)   
      include 'phokhara.inc'
      double precision x,value,histo(10,100),error(10,100)
      integer i,j
      common/histograms/histo,error
c --- add --------------------
      j = 1+(x-xlow(i))/(xup(i)-xlow(i))*bins(i)
      if (j.ge.1.and.j.le.bins(i)) then 
         histo(i,j)=histo(i,j)+value
      endif
      end

c ------------------------------------------------------
c     save histograms to file fname                     
c ------------------------------------------------------
      subroutine endhisto(fname)
      include 'phokhara.inc'
c --- PAW common variables ---
      double precision h,histo(10,100),error(10,100),x
      integer nwpawc,i,j
      parameter (nwpawc = 30000)
      common/pawc/h(nwpawc)
      common/histograms/histo,error    
c --- input/output file ------
      character*20 fname
c --- fill histograms --------
      do i=1,2                            !=========================
        do j=1,bins(i)
           x = xlow(i)+(j-.5d0)*(xup(i)-xlow(i))/dble(bins(i))
           if (count(i).ne.0.d0) then 
              error(i,j) = Mmax(i)*dSqrt((histo(i,j)/count(i)-
     &	                   (histo(i,j)/count(i))**2)/count(i))
              histo(i,j) = Mmax(i)/count(i)*histo(i,j)
              histo(i,j)=histo(i,j)*dble(bins(i))/(xup(i)-xlow(i))
              error(i,j)=error(i,j)*dble(bins(i))/(xup(i)-xlow(i))
           endif
           call hf1e(i,sngl(x),sngl(histo(i,j)),sngl(error(i,j)))
c	   write (*,*) x,histo(i,j),error(i,j)
        enddo
      enddo
c --- save histograms --------
      call hrput( 0, fname, 'N' )
      return
      end

 

c ========================================================================
c ========================================================================
c ------------------------------------------------------------------------
c ************************************************************************
c from here helicity amplitudes: H.C. 04.07.2001
c muons added : H.C. 04.10.2001
c ************************************************************************
c ------------------------------------------------------------------------

      double precision function helicityamp(Sp,qq,pion)
      include 'phokhara.inc'     
      double precision Sp,qq,rk1(4),rk2(4),q(4),dps,amp_h,amp
      complex*16 gam(4),gammu(4,2,2)
      integer i1,pion,ic1,ic2

      call gam1(gam,gammu,qq,pion)
      do i1=1,4
         rk1(i1) = momenta(3,i1-1)
         rk2(i1) = momenta(4,i1-1)
         q(i1)   = momenta(5,i1-1)
      enddo
c
      dps = Sp/(4.d0*(2.d0*pi)**5)  ! Phase space factors
c
c --- muons ---
      if(pion.eq.0)then
      dps = dps*dSqrt(1.d0-4.d0*mmu*mmu/qq)/(32.d0*pi*pi)
      amp_h = 0.d0
      do ic1=1,2
      do ic2=1,2
        do i1=1,4
         gam(i1)=gammu(i1,ic1,ic2)
        enddo
        amp_h = amp_h + amp(rk1,rk2,gam,q)
      enddo
      enddo
      amp_h = ((4.d0*pi*alpha)**4/qq**2)*amp_h/4.d0
c
c --- pions ---
      else
      dps = dps*dSqrt(1.d0-4.d0*mpi*mpi/qq)/(32.d0*pi*pi) 
      amp_h = ((4.d0*pi*alpha)**4/qq**2)*amp(rk1,rk2,gam,q)/4.d0
      endif
c
      helicityamp = amp_h*dps
      return
      end
      
c ---------------------------------------------------------------------

      subroutine gam1(gam,gammu,qq,pion)
      include 'phokhara.inc'       
      integer mu,pion,ic1,ic2
      real*8 qq,th1,th2,sphi1,cphi1,sphi2,cphi2,cth1d2,sth1d2,
     &       cth2d2,sth2d2,sq1,sq2,em1,em2,pm1,pm2,sth1,sth2
      complex*16 BW,gam(4),f1,gammu(4,2,2)
     &      ,v1(2,2),v2(2,2),up1(2,2),up2(2,2),ex1,ex2

c --- muons ---   
c  four different combinations of \mu^+, \mu^- helicities
c ++,+-,-+,-- ; v1(i,1)==v1(i,+),v1(i,2)==v1(i,-) etc.
c see notes p.4-6   
c
      if (pion.eq.0)then 
        em1 = momenta(7,0)
        em2 = momenta(6,0)
        pm1 = sqrt(em1**2-mmu**2)
        pm2 = sqrt(em2**2-mmu**2)
        sq1 = sqrt(em1+pm1)
        sq2 = sqrt(em2+pm2)
        th1 = acos(momenta(7,3)/pm1)
        th2 = acos(momenta(6,3)/pm2)
        cth1d2 = cos(th1/2.d0)
        cth2d2 = cos(th2/2.d0)
c
         v1(2,1)= -sq2*cth2d2
         v1(1,2)= mmu/sq2*cth2d2
         v2(2,1)= v1(1,2)
         v2(1,2)= v1(2,1)
c
         up1(1,1) = mmu/sq1*cth1d2
         up1(2,2) = sq1*cth1d2
         up2(1,1) = up1(2,2)
         up2(2,2) = up1(1,1)
c
        if((th2.eq.0.d0).or.(th2.eq.pi))then
         v1(1,1)= 0.d0
         v1(2,2)= 0.d0
         v2(1,1)= 0.d0
         v2(2,2)= 0.d0
        else
c
         sth2  = sin(th2)
         sth2d2= sin(th2/2.d0)
         cphi2 = momenta(6,1)/pm2/sth2
         sphi2 = momenta(6,2)/pm2/sth2
         ex2 = cphi2+dcmplx(0.d0,1.d0)*sphi2
c
         v1(1,1)= sq2*sth2d2/ex2
         v1(2,2)= mmu/sq2*ex2*sth2d2
         v2(1,1)= -mmu/sq2/ex2*sth2d2
         v2(2,2)= -sq2*ex2*sth2d2
        endif
c
        if((th1.eq.0.d0).or.(th1.eq.pi))then
         up1(1,2)=0.d0
         up1(2,1)=0.d0
         up2(1,2)=0.d0
         up2(2,1)=0.d0
        else
         sth1  = sin(th1)
         sth1d2= sin(th1/2.d0)
         cphi1 = momenta(7,1)/pm1/sth1
         sphi1 = momenta(7,2)/pm1/sth1
         ex1 = cphi1+dcmplx(0.d0,1.d0)*sphi1
c
         up1(1,2)= -sq1*ex1*sth1d2
         up1(2,1)= mmu/sq1/ex1*sth1d2
         up2(1,2)= -mmu/sq1*ex1*sth1d2
         up2(2,1)= sq1/ex1*sth1d2
        endif
c
       do ic1=1,2
       do ic2=1,2
        gammu(1,ic1,ic2) = up1(1,ic1)*v1(1,ic2)+up1(2,ic1)*v1(2,ic2)
     &                    +up2(1,ic1)*v2(1,ic2)+up2(2,ic1)*v2(2,ic2)
        gammu(2,ic1,ic2) = -up1(1,ic1)*v1(2,ic2)-up1(2,ic1)*v1(1,ic2)
     &                     +up2(1,ic1)*v2(2,ic2)+up2(2,ic1)*v2(1,ic2)
        gammu(3,ic1,ic2) = dcmplx(0.d0,1.d0)* 
     &                    (up1(1,ic1)*v1(2,ic2)-up1(2,ic1)*v1(1,ic2)
     &                    -up2(1,ic1)*v2(2,ic2)+up2(2,ic1)*v2(1,ic2))
        gammu(4,ic1,ic2) = -up1(1,ic1)*v1(1,ic2)+up1(2,ic1)*v1(2,ic2)
     &                     +up2(1,ic1)*v2(1,ic2)-up2(2,ic1)*v2(2,ic2)
       enddo
       enddo
c
      else
c --- pions ---
      f1 = (BW(mrho,gammarho,qq,1)*(1.D0+al*BW(momega,gomega,qq,1))/
     1     (1.d0+al)+be*BW(mrhol,grhol,qq,1))/(1.d0+be)
      do mu = 0,3
         gam(mu+1) = (momenta(6,mu)-momenta(7,mu))*f1
      enddo             
      endif
      return
      end


c --------------------------------------------------------------------- 
c     definicje macierzy, definicje iloczynow skalarnych, mi=minus,
c     pl=plus, eck1=epsilon conjug od k1 itd.
c
      subroutine skalar1(rk1,rk2,gam,q,eck1,eck2)
c      
      implicit real*8 (a-h,o-z)
      complex*16 gam,eck1,eck2,p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      complex*16 qpl(2,2),qmi(2,2),gampl(2,2),gammi(2,2),k1pl(2,2),
     1          k1mi(2,2),k2pl(2,2),k2mi(2,2),eck1pl(2,2),eck1mi(2,2),
     2          eck2pl(2,2),eck2mi(2,2),I(2,2) 
c
      dimension p1(4),p2(4),rk1(4),rk2(4),gam(4),q(4),eck1(4),eck2(4)
c
      common/iloczs1/p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      common/matri/qpl,qmi,gampl,gammi,k1pl,k1mi,k2pl,k2mi,eck1pl,
     1               eck1mi,eck2pl,eck2mi,I
      common/iloczs2/rk1p1,rk1p2,rk2p1,rk2p2,rk1rk2,anaw1,anaw2
      common /cp1p2/p1,p2,dme,el_m2
c
c
      rat1 = el_m2/(p1(1)+p1(4))
c
      cos1 = rk1(4) / rk1(1)
      cos2 = rk2(4) / rk2(1)
c
      rk1p1 = rk1(1) * ( rat1 + p1(4) * (1.d0 - cos1) )
c
      rk1p2 = rk1(1) * ( rat1 + p1(4) * (1.d0 + cos1) )
c
      rk2p1 = rk2(1) * ( rat1 + p1(4) * (1.d0 - cos2) ) 
c
      rk2p2 = rk2(1) * ( rat1 + p1(4) * (1.d0 + cos2) )
c
      rk1rk2 = rk1(1)*rk2(1) - rk1(2)*rk2(2) - rk1(3)*rk2(3) -
     1          rk1(4)*rk2(4)
c
      anaw1 = rk1rk2 - rk1p1 - rk2p1
      anaw2 = rk1rk2 - rk1p2 - rk2p2
c
      I(1,1)=dcmplx(1.d0,0.d0)
      I(1,2)=dcmplx(0.d0,0.d0)
      I(2,1)=dcmplx(0.d0,0.d0)
      I(2,2)=dcmplx(1.d0,0.d0)
c
      qpl(1,1)= dcmplx(q(1)-q(4),0.d0)
      qpl(1,2)=-q(2)+dcmplx(0.d0,1.d0)*q(3)
      qpl(2,1)=-q(2)-dcmplx(0.d0,1.d0)*q(3)
      qpl(2,2)=q(1)+q(4)
c
      qmi(1,1)=q(1)+q(4)
      qmi(1,2)=q(2)-dcmplx(0.d0,1.d0)*q(3)
      qmi(2,1)=q(2)+dcmplx(0.d0,1.d0)*q(3)
      qmi(2,2)=q(1)-q(4)
c
      gampl(1,1)=gam(1)-gam(4)
      gampl(1,2)=-gam(2)+dcmplx(0.d0,1.d0)*gam(3)
      gampl(2,1)=-gam(2)-dcmplx(0.d0,1.d0)*gam(3)
      gampl(2,2)=gam(1)+gam(4)
c
      gammi(1,1)=gam(1)+gam(4)
      gammi(1,2)=gam(2)-dcmplx(0.d0,1.d0)*gam(3)
      gammi(2,1)=gam(2)+dcmplx(0.d0,1.d0)*gam(3)
      gammi(2,2)=gam(1)-gam(4)
c
      k1pl(1,1)=rk1(1)-rk1(4)
      k1pl(1,2)=-rk1(2)+dcmplx(0.d0,1.d0)*rk1(3)
      k1pl(2,1)=-rk1(2)-dcmplx(0.d0,1.d0)*rk1(3)
      k1pl(2,2)=rk1(1)+rk1(4)
c
      k1mi(1,1)=rk1(1)+rk1(4)
      k1mi(1,2)=rk1(2)-dcmplx(0.d0,1.d0)*rk1(3)
      k1mi(2,1)=rk1(2)+dcmplx(0.d0,1.d0)*rk1(3)
      k1mi(2,2)=rk1(1)-rk1(4)
c
      k2pl(1,1)=rk2(1)-rk2(4)
      k2pl(1,2)=-rk2(2)+dcmplx(0.d0,1.d0)*rk2(3)
      k2pl(2,1)=-rk2(2)-dcmplx(0.d0,1.d0)*rk2(3)
      k2pl(2,2)=rk2(1)+rk2(4)
c
      k2mi(1,1)=rk2(1)+rk2(4)
      k2mi(1,2)=rk2(2)-dcmplx(0.d0,1.d0)*rk2(3)
      k2mi(2,1)=rk2(2)+dcmplx(0.d0,1.d0)*rk2(3)
      k2mi(2,2)=rk2(1)-rk2(4)
c
c     
c     iloczyny skalarne sa mnozone przez 2, co nie jest
c     uwzglednione w nazwach tych iloczynow!
c
      p1gam=2.d0*(p1(1)*gam(1)-p1(2)*gam(2)-p1(3)*gam(3)-
     1            p1(4)*gam(4))
      p2gam=2.d0*(p2(1)*gam(1)-p2(2)*gam(2)-p2(3)*gam(3)-
     1            p2(4)*gam(4))      
c
      return
      end
c---------------------------------------------------------------------
c----------------------------------------------------------------------
c     definicje macierzy, definicje iloczynow skalarnych, mi=minus,
c     pl=plus, eck1=epsilon conjug od k1 itd.
c
      subroutine skalar2(rk1,rk2,gam,q,eck1,eck2)
c      
      implicit real*8 (a-h,o-z)
      complex*16 gam,eck1,eck2,p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      complex*16 qpl(2,2),qmi(2,2),gampl(2,2),gammi(2,2),k1pl(2,2),
     1          k1mi(2,2),k2pl(2,2),k2mi(2,2),eck1pl(2,2),eck1mi(2,2),
     2          eck2pl(2,2),eck2mi(2,2),I(2,2) 
c
      dimension p1(4),p2(4),rk1(4),rk2(4),gam(4),q(4),eck1(4),eck2(4)
c
      common/iloczs1/p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      common/matri/qpl,qmi,gampl,gammi,k1pl,k1mi,k2pl,k2mi,eck1pl,
     1               eck1mi,eck2pl,eck2mi,I
      common /cp1p2/p1,p2,dme,el_m2
c
c
      eck1pl(1,1)=eck1(1)-eck1(4)
      eck1pl(1,2)=-eck1(2)+dcmplx(0.d0,1.d0)*eck1(3)
      eck1pl(2,1)=-eck1(2)-dcmplx(0.d0,1.d0)*eck1(3)
      eck1pl(2,2)=eck1(1)+eck1(4)
c
      eck1mi(1,1)=eck1(1)+eck1(4)
      eck1mi(1,2)=eck1(2)-dcmplx(0.d0,1.d0)*eck1(3)
      eck1mi(2,1)=eck1(2)+dcmplx(0.d0,1.d0)*eck1(3)
      eck1mi(2,2)=eck1(1)-eck1(4)     
c
      eck2pl(1,1)=eck2(1)-eck2(4)
      eck2pl(1,2)=-eck2(2)+dcmplx(0.d0,1.d0)*eck2(3)
      eck2pl(2,1)=-eck2(2)-dcmplx(0.d0,1.d0)*eck2(3)
      eck2pl(2,2)=eck2(1)+eck2(4)
c
      eck2mi(1,1)=eck2(1)+eck2(4)
      eck2mi(1,2)=eck2(2)-dcmplx(0.d0,1.d0)*eck2(3)
      eck2mi(2,1)=eck2(2)+dcmplx(0.d0,1.d0)*eck2(3)
      eck2mi(2,2)=eck2(1)-eck2(4)           
c     
c     iloczyny skalarne sa mnozone przez 2, co nie jest
c     uwzglednione w nazwach tych iloczynow!
c
      p1eck1=2.d0*(p1(1)*eck1(1)-p1(2)*eck1(2)-p1(3)*eck1(3)-
     1            p1(4)*eck1(4))
      p1eck2=2.d0*(p1(1)*eck2(1)-p1(2)*eck2(2)-p1(3)*eck2(3)-
     1            p1(4)*eck2(4))
      p2eck1=2.d0*(p2(1)*eck1(1)-p2(2)*eck1(2)-p2(3)*eck1(3)-
     1            p2(4)*eck1(4))
      p2eck2=2.d0*(p2(1)*eck2(1)-p2(2)*eck2(2)-p2(3)*eck2(3)-
     1            p2(4)*eck2(4))
c
      return
      end
c---------------------------------------------------------------------
c     mnozenie macierzy 2x2
      subroutine matr(mat1,mat2,mat3)
c
      complex*16  mat1(2,2),mat2(2,2),mat3(2,2)
c
      do i=1,2
         do j=1,2
            mat3(i,j)=(0.d0,0.d0)
         enddo
      enddo
c   
      do i=1,2
         do j=1,2
            do k=1,2
               mat3(i,j)=mat3(i,j)+mat1(i,k)*mat2(k,j)
            enddo
         enddo
      enddo
c
      end
c---------------------------------------------------------------------
c     mnozenie macierzy przez stala (to do iloczynow skalarnych)
      subroutine conmat(alfa,mat,amat)
c
      complex*16 alfa
      complex*16 mat(2,2),amat(2,2)
c
      do i=1,2
         do j=1,2
            amat(i,j)=alfa*mat(i,j)
         enddo
      enddo
c
      end
c--------------------------------------------------------------------
c     odejmowanie macierzy
      subroutine minmat(mat1,mat2,mat3)
c
      complex*16 mat1(2,2),mat2(2,2),mat3(2,2)
c
      do i=1,2
         do j=1,2
            mat3(i,j)=mat1(i,j)-mat2(i,j)
         enddo
      enddo
c
      end
c--------------------------------------------------------------------
c     dodawanie 6 macierzy: with proper denominators
c
      subroutine plumat(mat1,mat2,mat3,mat4,mat5,mat6,mat7)
c      
      implicit real*8 (a-h,o-z)
c
      complex*16 mat1(2,2),mat2(2,2),mat3(2,2),mat4(2,2),mat5(2,2),
     1          mat6(2,2),mat7(2,2)
c
      common/iloczs2/rk1p1,rk1p2,rk2p1,rk2p2,rk1rk2,anaw1,anaw2
c
      do i=1,2
         do j=1,2
            mat7(i,j)= -0.25d0*mat1(i,j)/rk2p2/anaw2
     1                 -0.25d0*mat2(i,j)/rk1p2/anaw2
     2                 +0.25d0*mat3(i,j)/rk1p2/rk2p1
     3                 +0.25d0*mat4(i,j)/rk2p2/rk1p1
     4                 -0.25d0*mat5(i,j)/rk2p1/anaw1
     5                 -0.25d0*mat6(i,j)/rk1p1/anaw1
c
         enddo
      enddo
c
      end
c--------------------------------------------------------------------
c
      subroutine blocks
c
      implicit real*8 (a,c-h,o-z)
      implicit complex*16 (b,m,n)
      complex*16 p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      complex*16 qpl(2,2),qmi(2,2),gampl(2,2),gammi(2,2),k1pl(2,2),
     1          k1mi(2,2),k2pl(2,2),k2mi(2,2),eck1pl(2,2),eck1mi(2,2),
     2          eck2pl(2,2),eck2mi(2,2),I(2,2)
      complex*16 m1(2,2),m2(2,2),m3(2,2),m4(2,2),m5(2,2),m6(2,2),
     1          n1(2,2),n2(2,2),n3(2,2),n4(2,2),n5(2,2),n6(2,2),n7(2,2),
     2          n8(2,2),n9(2,2),n10(2,2),n11(2,2),n12(2,2)
      complex*16 block1(2,2),block2(2,2),block3(2,2),block4(2,2),
     1          block5(2,2),block6(2,2),block7(2,2),block8(2,2),
     2          block9(2,2),block10(2,2),block11(2,2),block12(2,2)
      complex*16 m1amp1(2,2),mamp1a(2,2),m2amp1(2,2),mamp1b(2,2),
     1           m1amp2(2,2),mamp2a(2,2),m2amp2(2,2),mamp2b(2,2),
     2           m1amp3(2,2),mamp3a(2,2),m2amp3(2,2),mamp3b(2,2), 
     3           m1amp4(2,2),mamp4a(2,2),m2amp4(2,2),mamp4b(2,2),
     4           m1amp5(2,2),mamp5a(2,2),m2amp5(2,2),mamp5b(2,2),
     5           m1amp6(2,2),mamp6a(2,2),m2amp6(2,2),mamp6b(2,2)
      complex*16 ma(2,2),mb(2,2)
c
      common/iloczs1/p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      common/matri/qpl,qmi,gampl,gammi,k1pl,k1mi,k2pl,k2mi,eck1pl,
     1               eck1mi,eck2pl,eck2mi,I
      common/matri1/ma,mb
c
      call conmat(p1eck1,I,m1)
      call conmat(p1eck2,I,m2)
      call conmat(p2eck1,I,m3)
      call conmat(p2eck2,I,m4)
      call conmat(p1gam,I,m5)
      call conmat(p2gam,I,m6)
c
      call matr(qpl,gammi,n1)
      call matr(qmi,gampl,n2)
      call matr(gammi,qpl,n3)
      call matr(gampl,qmi,n4)
      call matr(eck1mi,k1pl,n5)
      call matr(eck1pl,k1mi,n6)
      call matr(eck2mi,k2pl,n7)
      call matr(eck2pl,k2mi,n8)
      call matr(k2pl,eck2mi,n9)
      call matr(k2mi,eck2pl,n10)
      call matr(k1pl,eck1mi,n11)
      call matr(k1mi,eck1pl,n12)
c
      call minmat(n7,m4,block1)
      call minmat(n8,m4,block2)
      call minmat(m5,n1,block3)
      call minmat(m5,n2,block4)
      call minmat(n5,m3,block5)
      call minmat(n6,m3,block6)
      call minmat(m2,n9,block7)
      call minmat(m2,n10,block8)
      call minmat(m1,n11,block9)
      call minmat(m1,n12,block10)
      call minmat(n3,m6,block11)
      call minmat(n4,m6,block12)
c
c     m1amp1=macierz powstala w wyniku mnozenia pierwszych dwoch 
c     macierzy, z pierwszego czlony w ampl1
c     mamp1a=macierz powstala w wyniku wymnozenia wszystkich macierzy
c     w pierwszym czlonie amp1 itd.
c
c     mnozenie macierzy dla amp1
c
      call matr(block1,eck1mi,m1amp1)
      call matr(m1amp1,block3,mamp1a)
      call matr(block2,eck1pl,m2amp1)
      call matr(m2amp1,block4,mamp1b)
c
c     mnozenie macierzy dla amp2
c
      call matr(block5,eck2mi,m1amp2)
      call matr(m1amp2,block3,mamp2a)
      call matr(block6,eck2pl,m2amp2)
      call matr(m2amp2,block4,mamp2b)
c
c     mnozenie macierzy dla amp3
c
      call matr(block5,gammi,m1amp3)
      call matr(m1amp3,block7,mamp3a)
      call matr(block6,gampl,m2amp3)
      call matr(m2amp3,block8,mamp3b)
c
c     mnozenie macierzy dla amp4
c
      call matr(block1,gammi,m1amp4)
      call matr(m1amp4,block9,mamp4a)
      call matr(block2,gampl,m2amp4)
      call matr(m2amp4,block10,mamp4b)
c
c     mnozenie macierzy dla amp5
c
      call matr(block11,eck1mi,m1amp5)
      call matr(m1amp5,block7,mamp5a)
      call matr(block12,eck1pl,m2amp5)
      call matr(m2amp5,block8,mamp5b)
c
c     mnozenie macierzy dla amp6
c
      call matr(block11,eck2mi,m1amp6)
      call matr(m1amp6,block9,mamp6a)
      call matr(block12,eck2pl,m2amp6)
      call matr(m2amp6,block10,mamp6b)
c
c     dodawanie macierzy typu a i typu b
c
      call plumat(mamp1a,mamp2a,mamp3a,mamp4a,mamp5a,mamp6a,ma)
      call plumat(mamp1b,mamp2b,mamp3b,mamp4b,mamp5b,mamp6b,mb)
c
      return
      end
c-----------------------------------------------------------
c     definicje spinorow 
c
      subroutine spinor(rk1,rk2)
c
      implicit real*8 (a-h,o-z)
      real*8 rk1(4),rk2(4)
      complex*16 epsk1(2,4),epsk2(2,4)
c
      common/polari/epsk1,epsk2
c
c
c     ideks 1 odpowiada epsilon(p,1)
c     ideks 2 odpowiada epsilon(p,2)
c
      cth1 = rk1(4)/rk1(1)
      sth1 = sqrt(1.d0-cth1**2)
      cth2 = rk2(4)/rk2(1)
      sth2 = sqrt(1.d0-cth2**2)
c
      if(sth1.ne.0.d0)then
        cphi1 = rk1(2)/sth1/rk1(1)
        sphi1 = rk1(3)/sth1/rk1(1)
      else
        if(cth1.gt.0.d0)then
           cphi1 = 1.d0
           sphi1 = 0.d0
        else
           cphi1 = -1.d0
           sphi1 = 0.d0
        endif
      endif
c
      if(sth2.ne.0.d0)then
        cphi2 = rk2(2)/sth2/rk2(1)
        sphi2 = rk2(3)/sth2/rk2(1)
      else
        if(cth2.gt.0.d0)then
           cphi2 = 1.d0
           sphi2 = 0.d0
        else
           cphi2 = -1.d0
           sphi2 = 0.d0
        endif
      endif
c
c      do i=1,4
c        epsk1(1,i)=rk1(i)
c        epsk1(2,i)=rk1(i)
c      enddo
c
c helicity basis
c
      epsk1(1,1)=dcmplx(0.d0,0.d0)
      epsk1(1,2)=dcmplx(cth1*cphi1,sphi1)/dsqrt(2.D0)
      epsk1(1,3)=dcmplx(cth1*sphi1,-cphi1)/dsqrt(2.D0)
      epsk1(1,4)=dcmplx(-sth1,0.d0)/dsqrt(2.D0)
      epsk1(2,1)=dcmplx(0.d0,0.d0)
      epsk1(2,2)=dcmplx(-cth1*cphi1,sphi1)/dsqrt(2.D0)
      epsk1(2,3)=dcmplx(-cth1*sphi1,-cphi1)/dsqrt(2.D0)
      epsk1(2,4)=dcmplx(sth1,0.d0)/dsqrt(2.D0)
c
      epsk2(1,1)=dcmplx(0.d0,0.d0)
      epsk2(1,2)=dcmplx(cth2*cphi2,sphi2)/dsqrt(2.D0)
      epsk2(1,3)=dcmplx(cth2*sphi2,-cphi2)/dsqrt(2.D0)
      epsk2(1,4)=dcmplx(-sth2,0.d0)/dsqrt(2.D0)
      epsk2(2,1)=dcmplx(0.d0,0.d0)
      epsk2(2,2)=dcmplx(-cth2*cphi2,sphi2)/dsqrt(2.D0)
      epsk2(2,3)=dcmplx(-cth2*sphi2,-cphi2)/dsqrt(2.D0)
      epsk2(2,4)=dcmplx(sth2,0.d0)/dsqrt(2.D0)
c
c cartesian basis
c
c      epsk1(1,1)=dcmplx(0.d0,0.d0)
c      epsk1(1,2)=dcmplx(cth1*cphi1,0.d0)
c      epsk1(1,3)=dcmplx(cth1*sphi1,0.d0)
c      epsk1(1,4)=dcmplx(-sth1,0.d0)
c      epsk1(2,1)=dcmplx(0.d0,0.d0)
c      epsk1(2,2)=dcmplx(-sphi1,0.d0)
c      epsk1(2,3)=dcmplx(cphi1,0.d0)
c      epsk1(2,4)=dcmplx(0.d0,0.d0)
c
c      epsk2(1,1)=dcmplx(0.d0,0.d0)
c      epsk2(1,2)=dcmplx(cth2*cphi2,0.d0)
c      epsk2(1,3)=dcmplx(cth2*sphi2,0.d0)
c      epsk2(1,4)=dcmplx(-sth2,0.d0)
c      epsk2(2,1)=dcmplx(0.d0,0.d0)
c      epsk2(2,2)=dcmplx(-sphi2,0.d0)
c      epsk2(2,3)=dcmplx(cphi2,0.d0)
c      epsk2(2,4)=dcmplx(0.d0,0.d0)
c
      return
      end
c----------------------------------------------------------
c
c
      real*8 function amp(rk1,rk2,gam,q)
c
      implicit real*8 (a-h,o-z)
c
      complex*16 p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      complex*16 epsk1(2,4),epsk2(2,4)
      complex*16 gam(4),eck1(4),eck2(4)
      complex*16 ma(2,2),mb(2,2)
c
      dimension rk1(4),rk2(4),p1(4),p2(4),q(4)
c
      common/polari/epsk1,epsk2
      common/iloczs1/p1eck1,p1eck2,p2eck1,p2eck2,p1gam,p2gam
      common/matri1/ma,mb
      common /cp1p2/p1,p2,dme,el_m2
c     
      amp = 0.d0
      call spinor(rk1,rk2)
      call skalar1(rk1,rk2,gam,q,eck1,eck2)
      ebppb = p1(1)+p1(4)
      do i=1,2
         do j=1,2
                  eck1(1)=epsk1(i,1)
                  eck1(2)=epsk1(i,2)
                  eck1(3)=epsk1(i,3)
                  eck1(4)=epsk1(i,4)
c                  
                  eck2(1)=epsk2(j,1)
                  eck2(2)=epsk2(j,2)
                  eck2(3)=epsk2(j,3)
                  eck2(4)=epsk2(j,4)
c
                  call skalar2(rk1,rk2,gam,q,eck1,eck2)
                  call blocks
      amp = amp+(dme*cdabs(mb(1,1)-ma(1,1)))**2
     1         +(dme*cdabs(mb(2,2)-ma(2,2)))**2
     2         +(cdabs(-ebppb*ma(1,2)+el_m2/ebppb*mb(1,2)))**2
     3         +(cdabs(ebppb*mb(2,1)-el_m2/ebppb*ma(2,1)))**2
c
c
         enddo
      enddo
c
      return
      end
c*************************************************************************
c ========================================================================


c ========================================================================
c === PHOKARA 1.0, (c) November 2001 =====================================








