program main
implicit none

integer,  parameter :: dp = kind(1.0D0)
real(dp), parameter :: pi = 3.14159265358979323846D0
real(dp), parameter :: hbar  = 1.0545726D-34     ! J*s
real(dp), parameter :: ev2j    = 1.60217733D-19  ! J
real(dp), parameter :: kB   = 1.380658D-23       ! J/K
real(dp), parameter :: thz2ev = 0.0041356        ! eV
real(dp), parameter :: ev2cm = 8065.541154       ! cm^-1
real(dp), parameter :: rad2ev   = 6.58551D-4     ! eV
real(dp) cc, vol, nb, t1, factor, scales, vec(3,3), dt, factor1, factor2,gauss
real(dp) dfreq,sigma,minq,delta,maxq,tdos,tdostemp,freq1,tdos1,tdos2,tdos3,tdostemp1,tdostemp2,tdostemp3
integer nq, nmode, iq, imode, aa, bb, i, n, nt, j, k, nfreq, ifreq
real(dp), allocatable::tau(:,:,:),tau0(:,:,:),tau_epc(:,:,:),temp(:)
real(dp), allocatable::vel(:,:,:),freq(:,:),kappal0(:,:,:),kappal(:,:,:)
real(dp), allocatable::kappal0_dos(:,:,:,:,:),kappal_dos(:,:,:,:,:),cv2_dos(:,:,:,:,:)
CHARACTER*50, ALLOCATABLE :: filename1(:),filename2(:)
character*1 nomeaning

open(unit=8,file="input")
read(8,*) nq, nmode, t1, nt, dt
read(8,*) sigma, delta, minq, maxq
read(8,*) scales
read(8,*) vec(1,1:3)
read(8,*) vec(2,1:3)
read(8,*) vec(3,1:3)
allocate(filename1(nt),filename2(nt))
do i = 1, nt
   read(8,*) filename1(i)         ! w file at each T
enddo
do i = 1, nt
   read(8,*) filename2(i)         ! phononlinewidth file at each T and Ef
enddo
close(8)

allocate(tau(nt,nq,nmode),tau0(nt,nq,nmode),tau_epc(nt,nq,nmode))
allocate(vel(nq,nmode,3),freq(nq,nmode),kappal0(nt,3,3),kappal(nt,3,3),temp(nt))
allocate(kappal0_dos(nt,nq,nmode,3,3),kappal_dos(nt,nq,nmode,3,3),cv2_dos(nt,nq,nmode,3,3))
! read frequency
open(unit=9,file="BTE.omega_full")
do iq = 1, nq
   read(9,*) freq(iq,:)  ! in rad/ps
enddo
freq=freq*rad2ev*ev2j    ! in J
close(9) 
! Read group velocity
open(unit=10,file='BTE.v_full')
do iq= 1, nq
   read(10,*) vel(iq,:,:) ! in km/s
enddo
vel=vel*1000.0           ! in m/s
close(10)
! Calculate volume
vec=vec*scales*10.0      ! in A
vol=vec(1,1)*(vec(2,2)*vec(3,3)-vec(3,2)*vec(2,3))+vec(2,1)*(vec(3,2)*vec(1,3)-vec(1,2)*vec(3,3))+vec(3,1)*(vec(1,2)*vec(2,3)-vec(2,2)*vec(1,3))
vol=vol*1.0D-30          ! in m^3
! Read RT of P-P
do i = 1, nt
   temp(i)=t1+dt*(i-1)   ! in K
   open(unit=10+i, file=filename1(i))
   do iq = 1, nq
      read(10+i,*) tau0(i,iq,1:nmode)   ! in ps-1
   enddo
   close(10+i)
enddo
tau0 = tau0*1.0D12       ! in s-1
! Read RT of EPC
do i = 1, nt
   open(unit=20+i, file=filename2(i))
   read(20+i,*)
   read(20+i,*)
   do iq = 1, nq
      do imode = 1, nmode
         read(20+i,*) aa, bb, cc, tau_epc(i,iq,imode)     ! in meV
      enddo
   enddo
   close(20+i)
enddo
tau_epc=(1.0D-3)*ev2j*tau_epc/hbar  ! in s-1
!tau_epc=tau_epc*10.0*(1.0D12)/6.58211
open(unit=30,file='out')
write(30,*) "T(K) w(rad/ps) vx(km/s) vy(km/s) vz(km/s) tau0^-1(ps^-1) tau_epc^-1(ps^-1)"
do i = 1, nt
   do iq = 1, nq
      do imode = 1, nmode
         write(30,"(F6.2,2X,6ES20.10)") temp(i), freq(iq,imode)/rad2ev/ev2j,vel(iq,imode,:)*1.0D-3,tau0(i,iq,imode)*1.0D-12,tau_epc(i,iq,imode)*1.0D-12
      enddo
    enddo
    write(30,*)
enddo
close(30)
! Calculate KappaL
do i = 1, nt
   factor=kB*(temp(i)**2)*vol*nq  ! in J*K*m^3
   do j = 1, 3
      do k = 1, 3
         kappal0(i,j,k)=0.0
         kappal(i,j,k)=0.0
         do iq = 1, nq
            do imode = 1, nmode
               if (freq(iq,imode).le.0.0) then
                  nb = 0.0
               else
                  nb=1.0/(exp(freq(iq,imode)/kB/temp(i))-1.0)
               endif
               tau(i,iq,imode)=tau0(i,iq,imode)+tau_epc(i,iq,imode)
               if (abs(tau(i,iq,imode)).gt. 0.0) then
                  factor2=1.0/tau(i,iq,imode)
               else
                  factor2=0.0
               endif
               if (abs(tau0(i,iq,imode)).gt. 0.0) then
                  factor1=1.0/tau0(i,iq,imode)
               else
                  factor1=0.0
               endif
               kappal(i,j,k)=kappal(i,j,k)+nb*(nb+1.0)*(freq(iq,imode)**2)*vel(iq,imode,j)*vel(iq,imode,k)*factor2
               kappal0(i,j,k)=kappal0(i,j,k)+nb*(nb+1.0)*(freq(iq,imode)**2)*vel(iq,imode,j)*vel(iq,imode,k)*factor1
               kappal_dos(i,iq,imode,j,k)=nb*(nb+1.0)*(freq(iq,imode)**2)*vel(iq,imode,j)*vel(iq,imode,k)*factor2/kB/(temp(i)**2)/vol
               kappal0_dos(i,iq,imode,j,k)=nb*(nb+1.0)*(freq(iq,imode)**2)*vel(iq,imode,j)*vel(iq,imode,k)*factor1/kB/(temp(i)**2)/vol
               cv2_dos(i,iq,imode,j,k)=nb*(nb+1.0)*(freq(iq,imode)**2)*vel(iq,imode,j)*vel(iq,imode,k)/kB/(temp(i)**2)/vol               
            enddo   ! nmode
         enddo      ! nq
         kappal(i,j,k)=kappal(i,j,k)/factor
         kappal0(i,j,k)=kappal0(i,j,k)/factor
      enddo         ! direction
   enddo            ! direction
enddo               ! temp

open(unit=50, file='out_kappal')
do i = 1, nt
   write(50,*) "T(k)                          KappaL(W/mK))"
   write(50,"(F6.2,2x,9ES20.10)") temp(i),kappal0(i,:,:)
   write(50,"(F6.2,2x,9ES20.10)") temp(i),kappal(i,:,:)
enddo
close(50)

! get dos
nfreq=(maxq-minq)/delta+1
open(unit=38,file="totaldos.txt")
write(38,*) "freq(rad/ps)          TDOS      CvV^2        KL0     KL"
do i = 1,nt
   do ifreq=1,nfreq
      freq1=minq+(ifreq-1)*delta
      tdos=0.0
      tdos1=0.0
      tdos2=0.0
      tdos3=0.0
      do iq = 1, nq
         tdostemp=0.0
         tdostemp1=0.0
         tdostemp2=0.0
         tdostemp3=0.0
         do imode=1, nmode
            dfreq=freq(iq,imode)/ev2j/rad2ev-freq1
            gauss=exp(-dfreq*dfreq/2/sigma/sigma)/sqrt(2*pi*sigma*sigma)
            tdostemp=tdostemp+gauss
            tdostemp1=tdostemp1+(cv2_dos(i,iq,imode,1,1)+cv2_dos(i,iq,imode,2,2)+cv2_dos(i,iq,imode,3,3))*gauss
            tdostemp2=tdostemp2+(kappal0_dos(i,iq,imode,1,1)+kappal0_dos(i,iq,imode,2,2)+kappal0_dos(i,iq,imode,3,3))*gauss
            tdostemp3=tdostemp3+(kappal_dos(i,iq,imode,1,1)+kappal_dos(i,iq,imode,2,2)+kappal_dos(i,iq,imode,3,3))*gauss
         enddo
         tdos=tdos+tdostemp/nq
         tdos1=tdos1+tdostemp1/nq
         tdos2=tdos2+tdostemp2/nq
         tdos3=tdos3+tdostemp3/nq
      enddo
   write(38,"(5E20.10)") freq1, tdos, tdos1, tdos2, tdos3
   enddo
   write(38,*)
enddo
close(38)


end         




      
