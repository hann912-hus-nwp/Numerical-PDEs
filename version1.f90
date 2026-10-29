program version_1
  implicit none

  integer, parameter :: n = 25
  real :: temp(n)
  real :: K, dt, dz, time, endtime, E, theta
  integer :: unit_out, hour_display, ischeme
  character(len=40) :: fname, label
  integer :: irow, hr(49)
  real :: res(49,4,2)   ! (hang theo gio, [Ts,E,T8cm,T96cm], so do 1/2)

  K       = 3e-7            ! He so khuech tan nhiet [m^2/s]
  dt      = 30.0            ! Buoc thoi gian [s]
  dz      = 0.04            ! Buoc luoi do sau [m]
  endtime = 48.0 * 3600.0   ! 48 gio [s]

  ! ischeme = 1 : so do an hoan toan (theta = 1)
  ! ischeme = 2 : so do ban an Crank-Nicolson (theta = 0.5)
  do ischeme = 1, 2

    if (ischeme == 1) then
      theta = 1.0
      fname = 'result_implicit_version1.txt'
      label = 'SO DO AN (Implicit)'
    else
      theta = 0.5
      fname = 'result_crank_nicolson_version1.txt'
      label = 'SO DO BAN AN CRANK-NICOLSON'
    end if

    call init(temp, n)
    time = 0.0
    irow = 0

    unit_out = 10
    open(unit=unit_out, file=trim(fname), status='replace', action='write')
    write(unit_out,'(A8,3X,A10,3X,A12,3X,A12,3X,A12)') &
         'Time(h)', 'Ts(C)', 'E(W/m2)', 'T8cm(C)', 'T96cm(C)'

    write(*,*) ''
    write(*,*) '===== ', trim(label), ' - GIAI BANG THOMAS ====='
    write(*,*) 'r = K*dt/dz^2 = ', K*dt/(dz*dz)
    write(*,'(A)') 'Time(h) |  Ts(C)  |     E      | T tai 8cm(do C) | T tai 96cm(do C)'
    write(*,'(A)') '========|=========|============|=================|================='

    do while (time <= endtime)

      temp(1) = tsfc(time)
      hour_display = mod(6 + int(time / 3600.0), 24)

      ! Tinh E nguoc tai thoi diem hien tai
      E = calculate_E_inverse(temp, K, dz, time)

      if (mod(int(time), 3600) == 0) then
        write(*,'(I6,A,F7.2,A,F11.1,A,F8.2,A,F8.2)') &
          hour_display, ' | ', temp(1)-273.15, ' | ', E, ' | ', &
          temp(3)-273.15, ' | ', temp(n)-273.15
        write(unit_out,'(I6,3X,F10.3,3X,F12.3,3X,F12.3,3X,F12.3)') &
          hour_display, temp(1)-273.15, E, temp(3)-273.15, temp(n)-273.15

        ! Luu lai de ghi file chung cho 2 so do
        irow = irow + 1
        hr(irow) = hour_display
        res(irow,1,ischeme) = temp(1)-273.15
        res(irow,2,ischeme) = E
        res(irow,3,ischeme) = temp(3)-273.15
        res(irow,4,ischeme) = temp(n)-273.15
      end if

      if (time >= endtime) exit

      ! Cap nhat nhiet do: theta = 1 (an), theta = 0.5 (Crank-Nicolson)
      call update_temperature_theta(temp, K, dt, dz, n, time, theta)

      time = time + dt
    end do

    close(unit_out)
  end do

  ! =========================================================
  ! FILE CHUNG: ket qua 2 so do canh nhau (An = Implicit, CN = Crank-Nicolson)
  ! =========================================================
  open(unit=11, file='result_version1.txt', status='replace', action='write')
  write(11,'(A8,2X,A7,2X,A9,2X,A11,2X,A11,2X,A9,2X,A9,2X,A10,2X,A10)') &
       'Elapsed(h)', 'Hour', 'Ts(C)', 'E_An', 'E_CN', 'T8_An', 'T8_CN', 'T96_An', 'T96_CN'
  do irow = 1, 49
    write(11,'(I8,2X,I7,2X,F9.3,2X,F11.3,2X,F11.3,2X,F9.3,2X,F9.3,2X,F10.3,2X,F10.3)') &
         irow-1, hr(irow), res(irow,1,1), &
         res(irow,2,1), res(irow,2,2), &
         res(irow,3,1), res(irow,3,2), &
         res(irow,4,1), res(irow,4,2)
  end do
  close(11)

contains

  subroutine init(temp, n)
    implicit none
    integer, intent(in) :: n
    real, intent(out) :: temp(n)
    temp = 25.0 + 273.15
  end subroutine init


  function tsfc(time) result(temp_surface)
    ! Dieu kien bien be mat: Ts(t) = 25 + 10*sin(omega*t) [do C] -> K
    implicit none
    real, intent(in) :: time
    real :: temp_surface
    real, parameter :: a0 = 10.0, T_mean = 25.0
    real, parameter :: pi = 4.0 * atan(1.0)
    real :: omega
    omega = (2.0 * pi) / (24.0 * 3600.0)
    temp_surface = a0 * sin(omega * time) + T_mean + 273.15
  end function tsfc


  function calculate_E_inverse(temp0, K, dz, time) result(E)
    ! E = rho*Cp*dz * [ dTs/dt - K*d2T/dz2 ]  (giu nguyen nhu code goc)
    implicit none
    real, intent(in) :: temp0(:), K, dz, time
    real :: E
    real :: dTdt, d2Tdz2, rho_cp_dz, a0, omega, pi

    rho_cp_dz = 2.0e6 * dz
    a0 = 10.0
    pi = 4.0 * atan(1.0)
    omega = (2.0 * pi) / (24.0 * 3600.0)

    dTdt   = a0 * omega * cos(omega * time)
    d2Tdz2 = (temp0(2) - temp0(1)) / (dz * dz)
    E = rho_cp_dz * (dTdt - K * d2Tdz2)
  end function calculate_E_inverse


  subroutine update_temperature_theta(temp, K, dt, dz, n, time, theta)
    !====================================================================
    ! Sơ đồ theta cho dT/dt = K d2T/dz2:
    !   (T^{n+1}-T^n)/dt = K*[ theta*L(T^{n+1}) + (1-theta)*L(T^n) ]
    !   theta = 1.0 : an hoan toan
    !   theta = 0.5 : Crank-Nicolson
    !
    ! Cho cac nut i = 2..N-1:
    !   a_i*T_{i-1}^{n+1} + b_i*T_i^{n+1} + c_i*T_{i+1}^{n+1} = d_i
    !   a_i = -theta*r ; b_i = 1 + 2*theta*r ; c_i = -theta*r
    !   d_i = T_i^n + (1-theta)*r*(T_{i-1}^n - 2T_i^n + T_{i+1}^n)
    !
    ! Bien da biet: T_1^{n+1} = Ts(t+dt), T_N^{n+1} = 25 C
    !   -> chuyen sang ve phai o hang dau (i=2) va hang cuoi (i=N-1)
    !====================================================================
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: temp(n)
    real, intent(in) :: K, dt, dz, time, theta

    integer :: i, m
    real :: r, T1_new, TN_new
    real :: a(n-2), b(n-2), c(n-2), d(n-2), x(n-2)

    r = K * dt / (dz * dz)
    m = n - 2                         ! so an: nut 2..N-1

    T1_new = tsfc(time + dt)          ! bien tren tai t+dt
    TN_new = 25.0 + 273.15            ! bien duoi

    do i = 2, n-1
      a(i-1) = -theta * r
      b(i-1) = 1.0 + 2.0 * theta * r
      c(i-1) = -theta * r
      d(i-1) = temp(i) + (1.0 - theta) * r * (temp(i-1) - 2.0*temp(i) + temp(i+1))
    end do

    ! Dua bien da biet sang ve phai
    d(1) = d(1) - a(1) * T1_new
    d(m) = d(m) - c(m) * TN_new

    call thomas(a, b, c, d, x, m)

    temp(1) = T1_new
    do i = 2, n-1
      temp(i) = x(i-1)
    end do
    temp(n) = TN_new
  end subroutine update_temperature_theta


  subroutine thomas(a, b, c, d, x, m)
    !====================================================================
    ! Thuat toan Thomas cho he 3 duong cheo:
    !   a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i),  i = 1..m
    ! (a(1) va c(m) khong dung den)
    !====================================================================
    implicit none
    integer, intent(in) :: m
    real, intent(in) :: a(m), b(m), c(m), d(m)
    real, intent(out) :: x(m)
    real :: cp(m), dp(m), den
    integer :: i

    ! Buoc xuoi (khu)
    cp(1) = c(1) / b(1)
    dp(1) = d(1) / b(1)
    do i = 2, m
      den   = b(i) - a(i) * cp(i-1)
      cp(i) = c(i) / den
      dp(i) = (d(i) - a(i) * dp(i-1)) / den
    end do

    ! Buoc nguoc (the)
    x(m) = dp(m)
    do i = m-1, 1, -1
      x(i) = dp(i) - cp(i) * x(i+1)
    end do
  end subroutine thomas

end program version_1