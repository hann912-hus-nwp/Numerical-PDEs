program test4

implicit none
integer, parameter :: n = 25
real :: temp0(n), temp1(n)
real :: K, dt, dz, time, endtime
real :: E, SW, LWs, LWa, Qh
real :: rho_soil, Cp_soil
real :: r 

integer :: istep,nsteps, istep_print
integer :: unit_out = 10

! Tham số lặp Picard
  integer :: iter
  integer, parameter :: max_iter = 100
  real, parameter :: tolerance = 1.0e-4
  real :: max_error
  real :: temp_guess(n)
  real :: temp_new(n)


! Khai báo hằng số
K = 3e-7
dt = 30.0
dz = 0.04
endtime = 72.0*3600.0

! Tính chất đất (giá trị lựa chọn đại diện)
rho_soil = 2000.0
Cp_soil = 1000.0


! Tính số bước và bước in
nsteps = int(endtime / dt)
istep_print = int(3600.0 / dt)


! Tính số Fourier (kiểm tra độ ổn định)
r = K*dt / (dz*dz)

! Gọi subroutine khởi tạo
call init(temp0, temp1, n)

! Ghi tiêu đề vào file TXT
      open(unit=unit_out, file='result_version4.txt', &
     status='replace',action='write')
      write(unit_out,'(A10,1X,A5,1X,A10,1X,A10,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A10,1X,A10)') &
'Time(h)','Day','Ts(C)','Ta(C)','SW(W/m2)','LWs(W/m2)', &
'LWa(W/m2)','Qh(W/m2)','E(W/m2)','T8cm(C)','T20cm(C)'

! In tiêu đề bảng Kết quả
    write(*,'(A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10)') &
'Time(h)','Day','Ts(do C)','Ta(do C)','SW(W/m2)','LWs(W/m2)', &
'LWa(W/m2)','Qh','E(W/m2)','T8cm(do C)','T20cm(do C)'

! Vòng lặp chính - 72 giờ
do istep = 0, nsteps
  time = real(istep) * dt
  
 call calculate_E_version4(temp0(1), time, SW, LWs, LWa, Qh, E)
   
! Ghi kết quả
 if(mod(istep, istep_print) == 0) then
    call write_results(unit_out, time, temp0, &
                         SW, LWs, LWa, Qh, E)
  end if

   if(istep >= nsteps) exit

   ! Cập nhật nhiệt độ bằng Crank-Nicolson
    ! LWs phụ thuộc vào Ts mới nên dùng Picard iteration
    call update_temperature_cn(temp0, temp1, K, dt, dz, n, &
                               rho_soil, Cp_soil, time, &
                               max_iter, tolerance)


end do

close(unit_out)
contains

subroutine init(temp0, temp1, n)
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

subroutine calculate_E_version4(Ts, time, SW, LWs, LWa, Qh, E)
 ! Tính E (E = SW - LWs + LWa + Qh)
implicit none
real, intent(in) :: Ts            ! nhiệt độ bề mặt đất, lấy từ temp0(1)
real, intent(in) :: time
real, intent(out) :: SW          ! bức xạ sóng ngắn
real, intent(out) :: LWa         ! bức xạ sóng dài khí quyển
real, intent(out) :: LWs         ! bức xạ sóng dài từ mặt đất
real, intent(out) :: Qh          ! thông lượng nhiệt cảm
real, intent(out) :: E           ! tổng năng lượng đưa vào phương trình nhiệt
real :: Ta                       ! nhiệt độ không khí
real :: Ta_K
real :: H                        ! thông lượng nhiệt cảm biến
real :: pi
real :: omega_day
real :: omega_solar
real :: time_day
real, parameter :: sigma = 5.67e-8
real, parameter :: epsilon_s = 0.95
real, parameter :: epsilon_a = 1.0
real, parameter :: SW_max = 1000.0
real, parameter :: V = 2.0              ! tốc độ gió
real, parameter :: rho_a = 1.0          ! khối lượng riêng không khí
real, parameter :: Cp_a = 1005.0        ! nhiệt dung riêng không khí
real, parameter :: Ch = 0.005           ! hệ số trao đổi nhiệt

! Lấy nhiệt độ bề mặt hiện tại 

pi = 4.0*atan(1.0)

omega_day = 2.0*pi/(24.0*3600.0)
omega_solar = pi/(12.0*3600.0)

time_day = modulo(time,24.0*3600.0)

! Tính SW
  if (time_day >= 0.0 .and. time_day <= 12.0*3600.0) then
    SW = SW_max*sin(omega_solar*time_day)

    if (SW < 0.0) SW = 0.0
else
    SW = 0.0
endif

! Tính LWa


    Ta = 23.0 + 5.0*sin(omega_day*time_day)

Ta_K = Ta + 273.15

    LWa = epsilon_a * sigma * Ta_K**4

! Tính LWs
    LWs = epsilon_s * sigma * (Ts**4.0)

! Tính H
    H = Ch * abs(V) * (Ta_K - Ts)

! Tính Qh
    Qh = rho_a * Cp_a * H


! Tính E
   E = SW - LWs + LWa + Qh

end subroutine calculate_E_version4

! CẬP NHẬT NHIỆT ĐỘ BẰNG CRANK-NICOLSON
 subroutine update_temperature_cn(temp0, temp1, K, dt, dz, n, &
                                   rho_soil, Cp_soil, time, &
                                   max_iter, tolerance)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: max_iter

    real, intent(in) :: K
    real, intent(in) :: dt
    real, intent(in) :: dz
    real, intent(in) :: rho_soil
    real, intent(in) :: Cp_soil
    real, intent(in) :: time
    real, intent(in) :: tolerance

    real, intent(inout) :: temp0(n)
    real, intent(inout) :: temp1(n)

    real :: a(n-1)
    real :: b(n-1)
    real :: c(n-1)
    real :: rhs(n-1)
    real :: solution(n-1)

    real :: temp_guess(n)
    real :: temp_new(n)

    real :: r
    real :: q
    real :: s
    real :: rho_cp_dz

    real :: SW_old
    real :: LWs_old
    real :: LWa_old
    real :: Qh_old
    real :: E_old

    real :: SW_new
    real :: LWs_guess
    real :: LWa_new
    real :: Qh_new
    real :: E_new_guess

    real :: max_error

    integer :: i
    integer :: iter

    ! Các hệ số
    r = K * dt / (dz * dz)
    q = 0.5 * r

    rho_cp_dz = rho_soil * Cp_soil * dz
    s = dt / (2.0 * rho_cp_dz)

    ! E tại thời điểm cũ
    call calculate_E_version4(temp0(1), time, &
                              SW_old, LWs_old, LWa_old, Qh_old, E_old)

    ! Nghiệm đoán ban đầu
    do i = 1, n
      temp_guess(i) = temp0(i)
    end do

    ! =======================================================
    ! PICARD ITERATION
    ! =======================================================

    do iter = 1, max_iter

      ! Tính E tại thời điểm mới dựa trên Ts đoán
      call calculate_E_version4(temp_guess(1), time + dt, &
                                SW_new, LWs_guess, LWa_new, Qh_new, &
                                E_new_guess)

      ! Khởi tạo hệ số
      do i = 1, n-1
        a(i) = 0.0
        b(i) = 0.0
        c(i) = 0.0
        rhs(i) = 0.0
      end do

      ! -----------------------------------------------------
      ! NODE 1 - BỀ MẶT
      !
      ! (1+q)T1_new - q*T2_new
      ! =
      ! T1_old + q(T2_old-T1_old)
      ! + s(E_old + E_new)
      ! -----------------------------------------------------

      a(1) = 0.0
      b(1) = 1.0 + q
      c(1) = -q

      rhs(1) = temp0(1) &
             + q * (temp0(2) - temp0(1)) &
             + s * (E_old + E_new_guess)

      ! -----------------------------------------------------
      ! NODE 2 ĐẾN NODE n-1
      ! -----------------------------------------------------

      do i = 2, n-1

        a(i) = -q
        b(i) = 1.0 + r
        c(i) = -q

        rhs(i) = q * temp0(i-1) &
               + (1.0 - r) * temp0(i) &
               + q * temp0(i+1)

      end do

      ! -----------------------------------------------------
      ! BIÊN SÂU: T(n) = 25°C
      ! -----------------------------------------------------

      rhs(n-1) = rhs(n-1) &
               + q * (25.0 + 273.15) &
               + q * (25.0 + 273.15)

      c(n-1) = 0.0

      ! Giải hệ ba đường chéo
      call thomas(n-1, a, b, c, rhs, solution)

      ! Lưu nghiệm mới
      do i = 1, n-1
        temp_new(i) = solution(i)
      end do

      temp_new(n) = 25.0 + 273.15

      ! Kiểm tra hội tụ
      max_error = 0.0

      do i = 1, n
        max_error = max(max_error, &
                        abs(temp_new(i) - temp_guess(i)))
      end do

      ! Cập nhật nghiệm đoán
      do i = 1, n
        temp_guess(i) = temp_new(i)
      end do

      if (max_error < tolerance) exit

    end do

    if (iter >= max_iter .and. max_error >= tolerance) then
      write(*,'(A,I5,A,E12.5)') &
           'Canh bao: Picard chua hoi tu sau ', &
           max_iter, ' vong lap. Sai so = ', max_error
    endif

    ! Cập nhật nhiệt độ
    do i = 1, n
      temp1(i) = temp_new(i)
      temp0(i) = temp_new(i)
    end do

  end subroutine update_temperature_cn


  ! =========================================================
  ! THOMAS ALGORITHM
  ! =========================================================

  subroutine thomas(nsys, a, b, c, d, x)
    implicit none

    integer, intent(in) :: nsys

    real, intent(in) :: a(nsys)
    real, intent(in) :: b(nsys)
    real, intent(in) :: c(nsys)
    real, intent(in) :: d(nsys)

    real, intent(out) :: x(nsys)

    real :: c_star(nsys)
    real :: d_star(nsys)
    real :: denominator

    integer :: i

    denominator = b(1)

    if (abs(denominator) < 1.0e-20) then
      write(*,*) 'LOI THOMAS: Pivot dau tien bang 0.'
      stop
    endif

    c_star(1) = c(1) / denominator
    d_star(1) = d(1) / denominator

    do i = 2, nsys

      denominator = b(i) - a(i) * c_star(i-1)

      if (abs(denominator) < 1.0e-20) then
        write(*,*) 'LOI THOMAS: Pivot bang 0 tai dong ', i
        stop
      endif

      c_star(i) = c(i) / denominator
      d_star(i) = (d(i) - a(i) * d_star(i-1)) / denominator

    end do

    x(nsys) = d_star(nsys)

    do i = nsys-1, 1, -1
      x(i) = d_star(i) - c_star(i) * x(i+1)
    end do

  end subroutine thomas


  ! =========================================================
  ! GHI KẾT QUẢ
  ! =========================================================

  subroutine write_results(unit_out, time, temp, SW, LWs, LWa, Qh, E)
    implicit none

    integer, intent(in) :: unit_out

    real, intent(in) :: time
    real, intent(in) :: temp(:)
    real, intent(in) :: SW
    real, intent(in) :: LWs
    real, intent(in) :: LWa
    real, intent(in) :: Qh
    real, intent(in) :: E

    integer :: day_number
    real :: Ta
    real :: pi
    real :: omega_day
    real :: time_day

    real, parameter :: day_seconds = 24.0 * 3600.0

    pi = 4.0 * atan(1.0)
    omega_day = 2.0 * pi / day_seconds
    time_day = modulo(time, day_seconds)

    Ta = 23.0 + 5.0 * sin(omega_day * time_day)

    day_number = int(time / day_seconds) + 1

    write(*,'(F7.2,A,I3,A,F8.3,A,F8.3,A,F10.3,A,F10.3,A,F10.3,A,F10.3,A,F10.3,A,F8.3,A,F8.3)') &
         time / 3600.0, ' | ', &
         day_number, ' | ', &
         temp(1) - 273.15, ' | ', &
         Ta, ' | ', &
         SW, ' | ', &
         LWs, ' | ', &
         LWa, ' | ', &
         Qh, ' | ', &
         E, ' | ', &
         temp(3) - 273.15, ' | ', &
         temp(6) - 273.15

    write(unit_out,'(F10.3,1X,I5,1X,F10.3,1X,F10.3,1X,F12.3,1X,F12.3,1X,F12.3,1X,F12.3,1X,F12.3,1X,F10.3,1X,F10.3)') &
  time/3600.0, day_number, temp(1)-273.15, Ta, SW, LWs, LWa, Qh, E, &
  temp(3)-273.15, temp(6)-273.15

  end subroutine write_results


end program test4