program version_4
  !====================================================================
  ! VERSION 4:  E = SW - LWs + LWa + Qh
  !   LWs = epsilon*sigma*Ts^4        (phi tuyen theo Ts)
  !   LWa = epsilon_air*sigma*Ta^4    (Ta dao dong sin: 23 + 5*sin(2*pi*t/24h))
  !   Qh  = rho_a*Cp_a*Ch*|V|*(Ta - Ts)   (nhiet cam bien, tuyen tinh theo Ts)
  !
  ! - So do an (theta = 1) va ban an Crank-Nicolson (theta = 0.5)
  ! - Picard: lap de xu ly E(Ts^{n+1}) (LWs phi tuyen, Qh phu thuoc Ts)
  ! - Thomas: giai he 3 duong cheo o moi vong lap Picard
  !====================================================================
  implicit none

  integer, parameter :: n = 25                 ! 25 node theo chieu sau
  integer, parameter :: nhours = 72            ! thoi gian mo phong [gio]
  integer, parameter :: nrows = nhours + 1     ! so dong ket qua (moi gio 1 dong)
  integer, parameter :: max_iter = 100         ! so vong lap Picard toi da
  real, parameter :: tolerance = 1.0e-4        ! sai so hoi tu Picard [K]

  real :: temp(n)
  real :: K, dt, dz, time, endtime, theta, r
  real :: E, SW, LWs, LWa, Qh, Ta
  real :: rho_soil, Cp_soil

  integer :: istep, nsteps, istep_print
  integer :: ischeme, unit_out, irow, hour_display
  integer :: hr(nrows)
  real :: res(nrows,5,2)    ! (dong theo gio, [Ts,Ta,E,T8cm,T96cm], so do 1/2)
  character(len=40) :: fname, label

  ! ---------------- THONG SO MO PHONG ----------------
  K       = 3.0e-7
  dt      = 30.0
  dz      = 0.04
  endtime = real(nhours) * 3600.0
  rho_soil = 2000.0
  Cp_soil  = 1000.0

  r = K * dt / (dz * dz)
  nsteps = nint(endtime / dt)
  istep_print = nint(3600.0 / dt)

  ! ischeme = 1 : so do an hoan toan (theta = 1)
  ! ischeme = 2 : so do ban an Crank-Nicolson (theta = 0.5)
  do ischeme = 1, 2

    if (ischeme == 1) then
      theta = 1.0
      fname = 'result_implicit_version4.txt'
      label = 'SO DO AN (Implicit)'
    else
      theta = 0.5
      fname = 'result_crank_nicolson_version4.txt'
      label = 'SO DO BAN AN CRANK-NICOLSON'
    end if

    call init(temp, n)
    time = 0.0
    irow = 0

    unit_out = 10
    open(unit=unit_out, file=trim(fname), status='replace', action='write')
    write(unit_out,'(A8,3X,A10,3X,A10,3X,A12,3X,A12,3X,A12)') &
         'Time(h)', 'Ts(C)', 'Ta(C)', 'E(W/m2)', 'T8cm(C)', 'T96cm(C)'

    write(*,*) ''
    write(*,*) '===== VERSION 4: ', trim(label), ' - GIAI BANG THOMAS + PICARD ====='
    write(*,*) 'r = K*dt/dz^2 = ', K*dt/(dz*dz)
    write(*,'(A)') 'Time(h) |  Ts(C)  |  Ta(C)  |     E      | T tai 8cm(do C) | T tai 96cm(do C)'
    write(*,'(A)') '========|=========|=========|============|=================|================='

    ! ---------------- VONG LAP THOI GIAN ----------------
    do istep = 0, nsteps

      time = real(istep) * dt

      ! E chan doan tai thoi diem hien tai: E = SW - LWs(Ts) + LWa + Qh
      call calculate_E_version4(temp(1), time, SW, LWs, LWa, Qh, Ta, E)

      if (mod(istep, istep_print) == 0) then
        hour_display = mod(6 + istep / istep_print, 24)

        write(*,'(I6,A,F7.2,A,F7.2,A,F11.1,A,F8.2,A,F8.2)') &
          hour_display, ' | ', temp(1)-273.15, ' | ', Ta, ' | ', E, ' | ', &
          temp(3)-273.15, ' | ', temp(n)-273.15

        write(unit_out,'(I6,3X,F10.3,3X,F10.3,3X,F12.3,3X,F12.3,3X,F12.3)') &
          hour_display, temp(1)-273.15, Ta, E, temp(3)-273.15, temp(n)-273.15

        ! Luu lai de ghi file chung cho 2 so do
        irow = irow + 1
        hr(irow) = hour_display
        res(irow,1,ischeme) = temp(1)-273.15
        res(irow,2,ischeme) = Ta
        res(irow,3,ischeme) = E
        res(irow,4,ischeme) = temp(3)-273.15
        res(irow,5,ischeme) = temp(n)-273.15
      end if

      if (istep >= nsteps) exit

      ! Cap nhat t^n -> t^(n+1): theta = 1 (an), theta = 0.5 (Crank-Nicolson)
      call update_temperature_theta(temp, K, dt, dz, n, rho_soil, Cp_soil, &
                                    time, theta, max_iter, tolerance)
    end do

    close(unit_out)
  end do

  ! =========================================================
  ! FILE CHUNG: ket qua 2 so do canh nhau (An = Implicit, CN = Crank-Nicolson)
  ! =========================================================
  open(unit=11, file='result_version4.txt', status='replace', action='write')
  write(11,'(A10,2X,A6,9(2X,A11))') 'Elapsed(h)', 'Hour', &
       'Ts_An', 'Ts_CN', 'Ta', 'E_An', 'E_CN', 'T8_An', 'T8_CN', 'T96_An', 'T96_CN'
  do irow = 1, nrows
    write(11,'(I10,2X,I6,9(2X,F11.3))') irow-1, hr(irow), &
         res(irow,1,1), res(irow,1,2), &
         res(irow,2,1), &
         res(irow,3,1), res(irow,3,2), &
         res(irow,4,1), res(irow,4,2), &
         res(irow,5,1), res(irow,5,2)
  end do
  close(11)

contains

  ! =========================================================
  subroutine init(temp, n)
    implicit none
    integer, intent(in) :: n
    real, intent(out) :: temp(n)
    temp = 25.0 + 273.15
  end subroutine init


  ! =========================================================
  ! SW, LWs, LWa, Qh, E = SW - LWs + LWa + Qh
  ! Qh = rho_a*Cp_a*H, H = Ch*|V|*(Ta - Ts)  (Ta, Ts tinh bang K)
  ! time = 0 : 6AM
  ! SW = SW_max*sin(omega*t) trong 12 gio ban ngay, 0 ban dem
  ! Ta = 23 + 5*sin(2*pi*t/24h)  [do C]  (Ta(6AM)=23, Ta(12PM)=28, Ta(12AM)=18)
  ! =========================================================
  subroutine calculate_E_version4(Ts, time, SW, LWs, LWa, Qh, Ta, E)
    implicit none
    real, intent(in)  :: Ts, time
    real, intent(out) :: SW, LWs, LWa, Qh, Ta, E   ! Ta [do C]

    real, parameter :: sigma = 5.67e-8
    real, parameter :: epsilon = 0.95        ! phat xa cua dat
    real, parameter :: epsilon_air = 1.0     ! phat xa cua khong khi
    real, parameter :: SW_max = 1000.0
    real, parameter :: Ta_mean = 23.0
    real, parameter :: Ta_amp  = 5.0
    real, parameter :: V     = 2.0           ! toc do gio [m/s]
    real, parameter :: rho_a = 1.0           ! khoi luong rieng khong khi [kg/m3]
    real, parameter :: Cp_a  = 1005.0        ! nhiet dung rieng khong khi [J/(kg K)]
    real, parameter :: Ch    = 0.005         ! he so trao doi nhiet
    real :: H
    real, parameter :: day_seconds = 24.0 * 3600.0
    real :: pi, omega, omega_day, time_day

    pi = 4.0 * atan(1.0)
    omega = pi / (12.0 * 3600.0)
    omega_day = 2.0 * pi / day_seconds
    time_day = modulo(time, day_seconds)

    ! Nhiet do khong khi [do C]
    Ta = Ta_mean + Ta_amp * sin(omega_day * time_day)

    if (time_day <= 12.0 * 3600.0) then
      SW = max(0.0, SW_max * sin(omega * time_day))
    else
      SW = 0.0
    end if

    LWs = epsilon * sigma * Ts**4
    LWa = epsilon_air * sigma * (Ta + 273.15)**4

    ! Nhiet cam bien: Ts (K) va Ta (K); Ta > Ts -> dat nhan nhiet
    H  = Ch * abs(V) * ((Ta + 273.15) - Ts)
    Qh = rho_a * Cp_a * H

    E = SW - LWs + LWa + Qh
  end subroutine calculate_E_version4


  ! =========================================================
  ! CAP NHAT NHIET DO - SO DO THETA (an: theta=1, CN: theta=0.5)
  !
  ! (T^{n+1}-T^n)/dt = theta*F(T^{n+1}) + (1-theta)*F(T^n)
  !
  ! Node 1 (be mat, an so):
  !   F = K*(T2-T1)/dz^2 + E/(rho*Cp*dz)
  !   (1+theta*r)T1 - theta*r*T2 = T1^n + (1-theta)*r*(T2^n-T1^n)
  !                                + s*(theta*E^{n+1} + (1-theta)*E^n)
  !   voi r = K*dt/dz^2, s = dt/(rho*Cp*dz)
  !
  ! Node 2..n-1:
  !   -theta*r*T_{i-1} + (1+2*theta*r)*T_i - theta*r*T_{i+1}
  !        = T_i^n + (1-theta)*r*(T_{i-1}^n - 2T_i^n + T_{i+1}^n)
  !
  ! Node n: T = 25 C (bien Dirichlet, da biet) -> chuyen ve phai
  !
  ! E^{n+1} = SW(t+dt) - LWs(T1^{n+1}) + LWa(t+dt) + Qh(T1^{n+1}) -> Picard
  ! (SW, LWa chi phu thuoc thoi gian; LWs, Qh phu thuoc T1):
  !   dung T1 cua vong lap truoc de tinh E^{n+1}, giai Thomas, lap den khi hoi tu.
  ! =========================================================
  subroutine update_temperature_theta(temp, K, dt, dz, n, rho_soil, Cp_soil, &
                                      time, theta, max_iter, tolerance)
    implicit none
    integer, intent(in) :: n, max_iter
    real, intent(inout) :: temp(n)
    real, intent(in) :: K, dt, dz, rho_soil, Cp_soil, time, theta, tolerance

    integer :: m, i, iter
    real :: a(n-1), b(n-1), c(n-1), rhs(n-1), solution(n-1)
    real :: temp_guess(n), temp_new(n)
    real :: r, s, rho_cp_dz, max_error
    real :: SW_old, LWs_old, LWa_old, Qh_old, Ta_old, E_old
    real :: SW_new, LWs_guess, LWa_new, Qh_new, Ta_new, E_new_guess
    real, parameter :: T_deep = 25.0 + 273.15

    m = n - 1                                   ! so an: node 1..n-1
    r = K * dt / (dz * dz)
    rho_cp_dz = rho_soil * Cp_soil * dz
    s = dt / rho_cp_dz

    ! E tai thoi diem cu
    call calculate_E_version4(temp(1), time, SW_old, LWs_old, LWa_old, Qh_old, Ta_old, E_old)

    temp_guess = temp
    temp_new   = temp

    do iter = 1, max_iter

      ! E tai thoi diem moi tu nghiem doan Picard
      call calculate_E_version4(temp_guess(1), time + dt, &
                                SW_new, LWs_guess, LWa_new, Qh_new, Ta_new, E_new_guess)

      ! ---- Node 1 ----
      a(1)   = 0.0
      b(1)   = 1.0 + theta * r
      c(1)   = -theta * r
      rhs(1) = temp(1) + (1.0 - theta) * r * (temp(2) - temp(1)) &
             + s * (theta * E_new_guess + (1.0 - theta) * E_old)

      ! ---- Node 2..n-1 ----
      do i = 2, n-1
        a(i)   = -theta * r
        b(i)   = 1.0 + 2.0 * theta * r
        c(i)   = -theta * r
        rhs(i) = temp(i) + (1.0 - theta) * r * &
                 (temp(i-1) - 2.0*temp(i) + temp(i+1))
      end do

      ! ---- Bien sau: T_n^{n+1} da biet -> chuyen ve phai o hang cuoi ----
      ! (phan T_n^n da nam trong rhs thong qua temp(i+1) o tren)
      rhs(m) = rhs(m) - c(m) * T_deep
      c(m)   = 0.0

      call thomas(m, a, b, c, rhs, solution)

      do i = 1, m
        temp_new(i) = solution(i)
      end do
      temp_new(n) = T_deep

      ! ---- Kiem tra hoi tu Picard ----
      max_error = 0.0
      do i = 1, n
        max_error = max(max_error, abs(temp_new(i) - temp_guess(i)))
      end do
      temp_guess = temp_new

      if (max_error < tolerance) exit
    end do

    if (iter > max_iter) then
      write(*,'(A,I5,A,E12.5)') &
           'Canh bao: Picard chua hoi tu sau ', max_iter, ' vong lap. Sai so = ', max_error
    end if

    temp = temp_new
  end subroutine update_temperature_theta


  ! =========================================================
  ! THOMAS: a(i)x(i-1) + b(i)x(i) + c(i)x(i+1) = d(i),  i = 1..nsys
  ! =========================================================
  subroutine thomas(nsys, a, b, c, d, x)
    implicit none
    integer, intent(in) :: nsys
    real, intent(in)  :: a(nsys), b(nsys), c(nsys), d(nsys)
    real, intent(out) :: x(nsys)

    real :: c_star(nsys), d_star(nsys), denominator
    integer :: i

    denominator = b(1)
    if (abs(denominator) < 1.0e-20) then
      write(*,*) 'LOI THOMAS: Pivot dau tien bang 0.'
      stop
    end if
    c_star(1) = c(1) / denominator
    d_star(1) = d(1) / denominator

    ! Khu xuoi
    do i = 2, nsys
      denominator = b(i) - a(i) * c_star(i-1)
      if (abs(denominator) < 1.0e-20) then
        write(*,*) 'LOI THOMAS: Pivot bang 0 tai dong ', i
        stop
      end if
      c_star(i) = c(i) / denominator
      d_star(i) = (d(i) - a(i) * d_star(i-1)) / denominator
    end do

    ! The nguoc
    x(nsys) = d_star(nsys)
    do i = nsys-1, 1, -1
      x(i) = d_star(i) - c_star(i) * x(i+1)
    end do
  end subroutine thomas

end program version_4