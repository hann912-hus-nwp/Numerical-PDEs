program test3
  implicit none

  integer, parameter :: n = 25
  integer, parameter :: m = n - 1
  integer, parameter :: unit_out = 10

  real :: temp0(n), temp1(n)

  real :: K
  real :: dt, dz, time, endtime
  real :: E, SW, LWs, LWa
  real :: rho_soil, Cp_soil
  real :: r

  real :: a(m), b(m), c(m), rhs(m)
  real :: solution(m)

  integer :: istep
  integer :: nsteps
  integer :: istep_print

  ! Tham số lặp Picard
  integer :: iter
  integer, parameter :: max_iter = 100
  real, parameter :: tolerance = 1.0e-4
  real :: max_error
  real :: temp_guess(n)
  real :: temp_new(n)

  ! =========================================================
  ! THÔNG SỐ MÔ PHỎNG
  ! =========================================================

  K       = 3.0e-7
  dt      = 30.0
  dz      = 0.04
  endtime = 72.0 * 3600.0

  ! Tính chất đất
  !
  ! rho_soil = 2000 kg/m3
  ! Cp_soil  = 1000 J/(kg K)
  ! rho*Cp   = 2.0e6 J/(m3 K)
  !
  rho_soil = 2000.0
  Cp_soil  = 1000.0

  ! Số Fourier
  r = K * dt / (dz * dz)

  ! Số bước thời gian
  nsteps = int(endtime / dt)

  ! In kết quả mỗi giờ
  istep_print = int(3600.0 / dt)

  ! Khởi tạo nhiệt độ
  call init(temp0, temp1, n)

  ! =========================================================
  ! MỞ FILE KẾT QUẢ
  ! =========================================================

  open(unit=unit_out, file='result_version3.txt', &
       status='replace', action='write')

  write(unit_out,'(A10,3X,A6,3X,A10,3X,A10,3X,A12,3X,A12,3X,A12,3X,A10,3X,A10)') &
       'Time(h)', 'Day', 'Ts(C)', 'Ta(C)', 'SW(W/m2)', &
       'LWs(W/m2)', 'LWa(W/m2)', 'E(W/m2)', 'T8cm(C)'

  ! =========================================================
  ! IN THÔNG TIN
  ! =========================================================

  write(*,*) ''
  write(*,*) '=========================================================='
  write(*,*) ' VERSION 3: E = SW - LWs + LWa'
  write(*,*) ' CRANK-NICOLSON + THOMAS ALGORITHM'
  write(*,*) '=========================================================='
  write(*,*) ''
  write(*,*) 'Quy uoc thoi gian:'
  write(*,*) '  time = 0 gio: 6AM'
  write(*,*) '  time = 6 gio: 12PM'
  write(*,*) '  time = 12 gio: 6PM'
  write(*,*) '  time = 24 gio: 6AM ngay tiep theo'
  write(*,*) ''
  write(*,*) 'Dieu kien ban dau:'
  write(*,*) '  T(z,0) = 25 do C'
  write(*,*) '  T(z=L) = 25 do C'
  write(*,*) ''
  write(*,*) 'Tinh chat dat:'
  write(*,*) '  rho_soil = 2000 kg/m3'
  write(*,*) '  Cp_soil  = 1000 J/(kg K)'
  write(*,*) '  rho*Cp   = 2.0e6 J/(m3 K)'
  write(*,*) ''
  write(*,*) 'Nhiet do khong khi:'
  write(*,*) '  Ta_mean = 23 do C'
  write(*,*) '  Ta_amp  = 5 do C'
  write(*,*) '  Ta dao dong chu ky 24 gio, bat dau tu 6AM'
  write(*,*) ''
  write(*,*) 'Thong so so:'
  write(*,'(A,E12.5,A)') '  K  = ', K, ' m2/s'
  write(*,'(A,F8.2,A)')  '  dt = ', dt, ' s'
  write(*,'(A,F8.4,A)')  '  dz = ', dz, ' m'
  write(*,'(A,F10.6)')   '  r  = ', r
  write(*,'(A,I8)')      '  So buoc = ', nsteps
  write(*,*) ''
  write(*,'(A)') &
       'Time(h) | Day | Ts(C) | Ta(C) | SW | LWs | LWa | E | T8cm(C)'
  write(*,'(A)') &
       '--------|-----|-------|-------|----|-----|-----|---|---------'

  ! =========================================================
  ! VÒNG LẶP THỜI GIAN
  ! =========================================================

  do istep = 0, nsteps

    time = real(istep) * dt

    ! Tính các thành phần năng lượng tại thời điểm hiện tại
    call calculate_E_version3(temp0(1), time, &
                              SW, LWs, LWa, E)

    ! In mỗi giờ
    if (mod(istep, istep_print) == 0) then
      call write_results(unit_out, time, temp0, &
                         SW, LWs, LWa, E)
    endif

    if (istep >= nsteps) exit

    ! Cập nhật nhiệt độ bằng Crank-Nicolson
    ! LWs phụ thuộc vào Ts mới nên dùng Picard iteration
    call update_temperature_cn(temp0, temp1, K, dt, dz, n, &
                               rho_soil, Cp_soil, time, &
                               max_iter, tolerance)

  end do

  close(unit_out)

  write(*,*) ''
  write(*,'(A,F8.2,A)') &
       'Ket thuc mo phong tai time = ', time / 3600.0, ' gio'
  write(*,*) 'Ket qua da luu vao file result_version3.txt'

contains

  ! =========================================================
  ! KHỞI TẠO NHIỆT ĐỘ
  ! =========================================================

  subroutine init(temp0, temp1, n)
    implicit none

    integer, intent(in) :: n
    real, intent(out) :: temp0(n)
    real, intent(out) :: temp1(n)

    integer :: i

    do i = 1, n
      temp0(i) = 25.0 + 273.15
      temp1(i) = 25.0 + 273.15
    end do

    write(*,*) 'Nhiet do ban dau: T(z,0) = 25 do C'

  end subroutine init


  ! =========================================================
  ! TÍNH SW, LW_s, LW_a VÀ E
  !
  ! Version 3:
  !
  ! E = SW - LWs + LWa
  !
  ! time = 0      : 6AM
  ! time = 21600  : 12PM
  ! time = 43200  : 6PM
  !
  ! SW:
  !   dao dong nua sin tu 6AM den 6PM
  !
  ! Ta:
  !   dao dong hinh sin chu ky 24 gio
  !   Ta = 23 + 5 sin(2*pi*time/24h)
  !
  ! LWs:
  !   epsilon_soil * sigma * Ts^4
  !
  ! LWa:
  !   epsilon_air * sigma * Ta^4
  ! =========================================================

  subroutine calculate_E_version3(Ts, time, SW, LWs, LWa, E)
    implicit none

    real, intent(in) :: Ts
    real, intent(in) :: time
    real, intent(out) :: SW
    real, intent(out) :: LWs
    real, intent(out) :: LWa
    real, intent(out) :: E

    real :: Ta
    real :: pi
    real :: omega_day    !tần số góc của Ta
    real :: omega_solar  !tần số góc của SW
    real :: time_day     !thời gian tính từ đầu ngày hiện tại
    real :: Ta_mean      !nhiệt độ không khí trung bình
    real :: Ta_amp       !biên độ

    real, parameter :: sigma = 5.67e-8
    real, parameter :: epsilon_soil = 0.95    !Lws
    real, parameter :: epsilon_air = 1.0      !Lwa
    real, parameter :: SW_max = 1000.0
    real, parameter :: day_seconds = 24.0 * 3600.0

    ! Nhiệt độ trung bình và biên độ không khí
    Ta_mean = 23.0
    Ta_amp  = 5.0

    pi = 4.0 * atan(1.0)

    ! Tần số dao động nhiệt độ không khí: chu kỳ 24 giờ
    omega_day = 2.0 * pi / day_seconds

    ! Tần số bức xạ mặt trời: ban ngày 12 giờ
    omega_solar = pi / (12.0 * 3600.0)

    ! Thời gian trong ngày
    time_day = modulo(time, day_seconds)

    ! =======================================================
    ! NHIỆT ĐỘ KHÔNG KHÍ Ta
    !
    ! time = 0 tương ứng với 6AM
    !
    ! Ta(6AM)  = 23°C
    ! Ta(12PM) = 28°C
    ! Ta(6PM)  = 23°C
    ! Ta(12AM) = 18°C
    ! =======================================================

    Ta = Ta_mean + Ta_amp * sin(omega_day * time_day)

    ! =======================================================
    ! BỨC XẠ SÓNG NGẮN SW
    !
    ! 6AM : 0
    ! 12PM: 1000 W/m2
    ! 6PM : 0
    ! Ban đêm: 0
    ! =======================================================

    if (time_day >= 0.0 .and. &
        time_day <= 12.0 * 3600.0) then

      SW = SW_max * sin(omega_solar * time_day)

      if (SW < 0.0) then
        SW = 0.0
      endif

    else
      SW = 0.0
    endif

    ! =======================================================
    ! BỨC XẠ SÓNG DÀI TỪ MẶT ĐẤT
    ! =======================================================

    LWs = epsilon_soil * sigma * Ts**4

    ! =======================================================
    ! BỨC XẠ SÓNG DÀI KHÍ QUYỂN
    ! =======================================================

    LWa = epsilon_air * sigma * (Ta + 273.15)**4

    ! =======================================================
    ! CÂN BẰNG NĂNG LƯỢNG
    ! =======================================================

    E = SW - LWs + LWa

  end subroutine calculate_E_version3


  ! =========================================================
  ! CẬP NHẬT NHIỆT ĐỘ BẰNG CRANK-NICOLSON
  !
  ! Node 1:
  !
  ! T1(n+1) - T1(n)
  ! ---------------- =
  !        dt
  !
  ! 1/2 [K(T2n-T1n)/dz2 + E(n)/(rho Cp dz)]
  !
  ! + 1/2 [K(T2n+1-T1n+1)/dz2
  !        + E(n+1)/(rho Cp dz)]
  !
  ! E(n+1) được tính trong Picard iteration.
  ! =========================================================

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
    real :: E_old

    real :: SW_new
    real :: LWs_guess
    real :: LWa_new
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
    call calculate_E_version3(temp0(1), time, &
                              SW_old, LWs_old, LWa_old, E_old)

    ! Nghiệm đoán ban đầu
    do i = 1, n
      temp_guess(i) = temp0(i)
    end do

    ! =======================================================
    ! PICARD ITERATION
    ! =======================================================

    do iter = 1, max_iter

      ! Tính E tại thời điểm mới dựa trên Ts đoán
      call calculate_E_version3(temp_guess(1), time + dt, &
                                SW_new, LWs_guess, LWa_new, &
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

  subroutine write_results(unit_out, time, temp, SW, LWs, LWa, E)
    implicit none

    integer, intent(in) :: unit_out

    real, intent(in) :: time
    real, intent(in) :: temp(:)
    real, intent(in) :: SW
    real, intent(in) :: LWs
    real, intent(in) :: LWa
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

    write(*,'(F7.2,A,I3,A,F8.3,A,F8.3,A,F10.3,A,F10.3,A,F10.3,A,F10.3,A,F8.3)') &
         time / 3600.0, ' | ', &
         day_number, ' | ', &
         temp(1) - 273.15, ' | ', &
         Ta, ' | ', &
         SW, ' | ', &
         LWs, ' | ', &
         LWa, ' | ', &
         E, ' | ', &
         temp(3) - 273.15

    write(unit_out, &
          '(F10.3,3X,I6,3X,F10.3,3X,F10.3,3X,F12.3,3X,F12.3,3X,F12.3,3X,F12.3,3X,F10.3)') &
         time / 3600.0, &
         day_number, &
         temp(1) - 273.15, &
         Ta, &
         SW, &
         LWs, &
         LWa, &
         E, &
         temp(3) - 273.15

  end subroutine write_results

end program test3