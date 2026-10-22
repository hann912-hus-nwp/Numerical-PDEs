program test3

implicit none
integer, parameter :: n=25
real :: temp0(n), temp1(n)
real :: K
real :: dt ,dz, time, endtime
real :: E, SW, LWs, LWa
real :: rho_soil, Cp_soil
real :: r
real :: Ta
integer :: istep, nsteps, istep_print
integer :: unit_out

! Khai báo hằng số
K = 3e-7
dt = 30.0
dz = 0.04
endtime = 24.0*3600.0

! Tính chất đất (giá trị lựa chọn đại diện)
rho_soil = 2000.0
Cp_soil = 1000.0

! Tính số bước và bước in
nsteps = int(endtime/ dt)
istep_print = int(3600.0/ dt)

! Tính số Fourier (kiểm tra độ ổn định)
r = K*dt / (dz*dz)

! Gọi subroutine khởi tạo

call init(temp0, temp1, n)

! Ghi tiêu đề vào file TXT
        unit_out = 10
      open(unit=unit_out, file='result_version3.txt', &
     status='replace',action='write')
      write(unit_out,'(A10,2X,A10,2X,A10,2X,A10,2X,A10,2X,A10,2X,A10,2X,A10)') &
'Time(h)','Ts(do C)','SW(W/m2)','LWa(W/m2)', &
'LWs(W/m2)','E(W/m2)','T8cm(do C)','T20cm(do C)'

! In tiêu đề bảng Kết quả
    write(*,'(A12,3X,A12,3X,A12,3X,A12,3X,A12,3X,A12,3X,A12,3X,A12)') &
'Time(h)','Ts(do C)','SW(W/m2)','LWa(W/m2)', &
'LWs(W/m2)','E(W/m2)','T8cm(do C)','T20cm(do C)' 

! Vòng lặp chính - 24 giờ
do istep = 0, nsteps
       time = istep*dt   ! Tính thời gian hiện tại
  call caculate_E_version3(temp0,time, SW, LWs, Lwa, E)
  if(mod(istep, istep_print) == 0) then
    write(*,'(F12.1,3X,F12.2,3X,F12.2,3X,F12.2,3X,F12.2,3X,F12.2,3X,F12.2,3X,F12.2)') &
    time/3600.0, temp0(1)-273.15, SW, LWa, LWs, E, &
    temp0(3)-273.15, temp0(6)-273.15
  end if

    ! Ghi kết quả vào file TXT
if (mod(istep, istep_print) == 0) then
 write(unit_out,'(F10.1,2X,F10.2,2X,F10.2,2X,F10.2,2X,F10.2,2X,F10.2,2X,F10.2,2X,F10.2)') &
time/3600.0, temp0(1)-273.15, SW, LWa, LWs, E, &
temp0(3)-273.15, temp0(6)-273.15
endif

 if(istep >= nsteps) exit

 call update_temprature(temp0, temp1, K, dt, dz, n, E, rho_soil, Cp_soil)

end do

close(unit_out)
contains

subroutine init(temp0,temp1, n)
! Khởi tạo nhiệt độ ban đầu

implicit none
integer, intent(in) :: n
real, intent(out) :: temp0(n), temp1(n)
integer :: i 

  do i = 1,n
    temp0(i) = 25.0 +273.15
    temp1(i) = 25.0 +273.15
  end do
end subroutine init

subroutine caculate_E_version3(temp0, time, SW, LWs, LWa, E)
 ! Tính E (E = SW - LWs + LWa)
implicit none
real, intent(in) :: temp0(:)
real, intent(in) :: time
real, intent(out) :: SW
real, intent(out) :: LWa
real, intent(out) :: LWs
real, intent(out) :: E
real :: Ts
real :: Ta
real :: Ta_K
real ::pi, omega, time_from_6am, amplitude
real, parameter :: sigma = 5.67e-8
real, parameter :: epsilon_s = 0.95
real, parameter :: epsilon_a = 1.0
real, parameter :: SW_max = 1000.0

! Lấy nhiệt độ bề mặt hiện tại 
Ts = temp0(1)

pi = 4.0*atan(1.0)
omega = (2.0*pi) / (24.0*3600.0)
time_from_6am = time - 6.0*3600.0 

! Tính SW
  if(time_from_6am >= 0.0.and.time_from_6am <= 12.0*3600.0) then
      SW = SW_max*sin(omega*time_from_6am)
         if(SW < 0.0) SW = 0.0
        else
         SW = 0.0
  end if

! Tính LWa

amplitude = 5.0   ! Biên độ

    Ta = 23.0 + amplitude * &
         sin(omega * time_from_6am)

Ta_K = Ta + 273.15

    LWa = epsilon_a * sigma * Ta_K**4

! Tính LWs
   LWs = epsilon_s * sigma * (Ts**4.0)

! Tính E
   E = SW - LWs + LWa

end subroutine caculate_E_version3

subroutine update_temprature(temp0, temp1, K, dt, dz, n, E, rho_soil, Cp_soil)
! Cập nhật nhiệt độ
implicit none
integer, intent(in) :: n
real, intent(inout) :: temp0(n), temp1(n)
real, intent(in) :: K
real, intent(in) :: dt
real, intent(in) :: dz
real, intent(in) :: E
real, intent(in) :: rho_soil
real, intent(in) :: Cp_soil
integer :: i
real :: r
real :: rho_cp_dz

! Tính hệ số Fourier
r = K * dt/(dz*dz)

! Tính rho_cp_dz
rho_cp_dz = rho_soil*Cp_soil*dz

! Cập nhật node 1 (bề mặt)
temp1(1) = temp0(1) + dt *(&
           K*(temp0(2) - temp0(1)) / (dz*dz) + &
           E / rho_cp_dz )

! Cập nhật node 2 -> n-1 (trung gian)
do i = 2,n-1 
  temp1(i) = temp0(i) + r*(temp0(i-1) - 2.0*temp0(i) + temp0(i+1))
end do

! Cập nhật node n (biên sâu)
temp1(n) = 25.0 +273.15

! Hoán đổi mảng
do i = 1,n
  temp0(i) = temp1(i)
end do

end subroutine update_temprature

end program test3