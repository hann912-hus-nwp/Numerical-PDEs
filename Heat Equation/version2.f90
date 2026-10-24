program test2
  ! Picard: Xử lý phi tuyến, lặp cho đến khi hội tụ
  !Vì bài toán có thành phần LWs= epsilon*sigma*T_s^4, làm phương trình trở thành phi tuyến. Picard được dùng để lặp và xử lý phần phi tuyến, mỗi vòng lặp Picard tạo ra một hệ tuyến tính được giải bằng thuật toán Thomas.
  ! Crank-Nicolson: Rời rạc hóa phương trình theo thời gian. Crank–Nicolson biến phương trình vi phân thành hệ phương trình đại số.
  ! Thomas: Giải hệ phương trình 3 đường chéo sau khi đã tuyến tính hóa

  implicit none

  integer, parameter :: n = 25         ! 25 node nhiệt độ theo chiều sâu
  integer, parameter :: unit_out = 10  ! số hiệu file = 10
  real :: temp0(n), temp1(n)           ! 2 mảng nhiệt độ tại thời điểm hiện tại và thời điểm mới

  real :: K                            ! hệ số khuếch tán nhiệt
  real :: dt, dz, time, endtime        
  real :: E, SW, LWs
  
  real :: rho_soil, Cp_soil
  real :: r
  
  integer :: istep         ! số bước thời gian hiện tại
  integer :: nsteps        ! tổng số bước thời gian
  integer :: istep_print   ! khoảng bao nhiêu bước thì in KQ/ lần

  ! Tham số lặp Picard
  
  integer, parameter :: max_iter = 100   ! số vòng lặp tối đa
  real, parameter :: tolerance = 1.0e-4  ! sai số cho phép
  
  ! =========================================================
  ! THÔNG SỐ MÔ PHỎNG
  ! =========================================================

  K       = 3.0e-7
  dt      = 30.0           ! bước thời gian 30s
  dz      = 0.04           ! khoảng cách giữa 2 node = 4cm
  endtime = 72.0 * 3600.0  ! thời gian kết thúc mô phỏng = 3 ngày

  ! Tính chất đất - giá trị đại diện
  !
  ! rho_soil = 2000 kg/m3       ! khối lượng riêng của đất
  ! Cp_soil  = 1000 J/(kg K)    ! nhiệt dung riêng
  ! rho*Cp   = 2.0e6 J/(m3 K)
  !
  rho_soil = 2000.0
  Cp_soil  = 1000.0

  ! Số Fourier (cho biết mức độ khuếch tán nhiệt trong một bước thời gian)
  r = K * dt / (dz * dz)

  ! Số bước thời gian
  nsteps = int(endtime / dt)

  ! In kết quả mỗi giờ
  istep_print = int(3600.0 / dt)

  ! =========================================================
  ! KHỞI TẠO NHIỆT ĐỘ
  ! =========================================================

  call init(temp0, temp1, n)

  ! =========================================================
  ! MỞ FILE KẾT QUẢ
  ! =========================================================

  open(unit=unit_out, file='result_version2.txt', &
       status='replace', action='write')

  write(unit_out,'(A10,3X,A6,3X,A10,3X,A12,3X,A12,3X,A12,3X,A10,3X,A10)') &
       'Time(h)', 'Day', 'Ts(C)', 'SW(W/m2)', 'LWs(W/m2)', &
       'E(W/m2)', 'T8cm(C)', 'T96cm(C)'

  ! =========================================================
  ! IN THÔNG TIN
  ! =========================================================

  write(*,'(A)') &
       'Time(h) | Day | Ts(C)   | SW       | LWs      | E        | T8cm    | T96cm'
  write(*,'(A)') &
       '--------|-----|---------|----------|----------|----------|---------|--------'

  ! =========================================================
  ! VÒNG LẶP THỜI GIAN
  ! =========================================================

  do istep = 0, nsteps

    time = real(istep) * dt     !real() : chuyển số nguyên thành số thực.

    ! Tính E tại trạng thái hiện tại để in kết quả
    call calculate_E_version2(temp0(1), time, SW, LWs, E)

    if (mod(istep, istep_print) == 0) then
      call write_results(unit_out, time, temp0, SW, LWs, E)
    endif

    if (istep >= nsteps) exit

    ! =======================================================
    ! CẬP NHẬT TỪ t^n ĐẾN t^(n+1)
    !
    ! Dùng Crank-Nicolson.
    ! LWs(t+dt) phụ thuộc vào Ts(t+dt), nên dùng Picard iteration.
    ! =======================================================

    call update_temperature_cn(temp0, temp1, K, dt, dz, n, &
                               rho_soil, Cp_soil, time, &
                               max_iter, tolerance)

  end do

  close(unit_out)

  write(*,*) ''
  write(*,'(A,F8.2,A)') &
       'Ket thuc mo phong tai time = ', time / 3600.0, ' gio'
  write(*,*) 'Ket qua da luu vao file result_version2.txt'

contains

  ! =========================================================
  ! KHỞI TẠO
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
  ! TÍNH SW, LW_s VÀ E
  !
  ! time = 0       : 6AM
  ! time = 21600   : 12PM
  ! time = 43200   : 6PM
  !
  ! SW lặp lại sau mỗi 24 giờ.
  !
  ! E = SW - LWs
  ! LWs = epsilon * sigma * Ts^4
  ! =========================================================

  subroutine calculate_E_version2(Ts, time, SW, LWs, E)
    implicit none

    real, intent(in) :: Ts
    real, intent(in) :: time
    real, intent(out) :: SW
    real, intent(out) :: LWs
    real, intent(out) :: E

    real :: pi
    real :: omega
    real :: time_day

    real, parameter :: sigma = 5.67e-8
    real, parameter :: epsilon = 0.95
    real, parameter :: SW_max = 1000.0
    real, parameter :: day_seconds = 24.0 * 3600.0

    pi = 4.0 * atan(1.0)

    ! Ban ngày dài 12 giờ
    omega = pi / (12.0 * 3600.0)

    ! Lấy thời gian trong ngày (Đây là phép lấy thời gian còn lại trong một chu kỳ 24 giờ)
    time_day = modulo(time, day_seconds)

    ! -------------------------------------------------------
    ! BỨC XẠ SÓNG NGẮN
    !
    ! 6AM : 0
    ! 12PM: SW_max
    ! 6PM : 0
    ! Ban đêm: 0
    ! -------------------------------------------------------

    if (time_day >= 0.0 .and. &
        time_day <= 12.0 * 3600.0) then

      SW = SW_max * sin(omega * time_day)

      if (SW < 0.0) then
        SW = 0.0
      endif

    else
      SW = 0.0
    endif

    ! -------------------------------------------------------
    ! BỨC XẠ SÓNG DÀI TỪ BỀ MẶT
    ! -------------------------------------------------------

    LWs = epsilon * sigma * Ts**4

    ! Cân bằng năng lượng
    E = SW - LWs

  end subroutine calculate_E_version2


  ! =========================================================
  ! CẬP NHẬT NHIỆT ĐỘ BẰNG CRANK-NICOLSON
  !
  ! Node 1:
  !
  ! T1(n+1) - T1(n)
  ! ---------------- = trung bình của:
  !        dt
  !
  ! K(T2-T1)/dz^2 + E/(rho Cp dz)
  !
  ! Các node 2 đến n-1:
  !
  ! dT/dt = K d2T/dz2
  !
  ! Node n:
  !
  ! Tn = 25°C
  !
  ! LWs(n+1) phụ thuộc vào T1(n+1).
  ! Dùng Picard iteration để xử lý phi tuyến.
  ! =========================================================

  subroutine update_temperature_cn(temp0, temp1, K, dt, dz, n, &
                                   rho_soil, Cp_soil, time, &
                                   max_iter, tolerance)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: max_iter  !số vòng lặp Picard tối đa

    real, intent(in) :: K
    real, intent(in) :: dt
    real, intent(in) :: dz
    real, intent(in) :: rho_soil
    real, intent(in) :: Cp_soil
    real, intent(in) :: time
    real, intent(in) :: tolerance   ! sai số hội tụ

    real, intent(inout) :: temp0(n)
    real, intent(inout) :: temp1(n)

    ! Các mảng của hệ phương trình THOMAS (dùng để lập và giải hệ phương trình ba đường chéo)
    real :: a(n-1)        !hệ số của node phía trên
    real :: b(n-1)        !hệ số của node đang xét
    real :: c(n-1)        !hệ số của node phía dưới
    real :: rhs(n-1)      !vế phải của hệ phương trình
    real :: solution(n-1) !kết quả sau khi Thomas giải hệ

    real :: temp_guess(n)   ! nghiệm đoán trong Picard
    real :: temp_new(n)     !nghiệm mới vừa được Thomas giải ra

    real :: r             !hệ số Fourier
    real :: q             !trung bình giữa thời điểm cũ và mới (Crank-Nicolson)
    real :: s             !s = dt / (2.0 * rho_cp_dz), dùng để đưa năng lượng E vào phương trình node 1
    real :: rho_cp_dz

    !Các biến năng lượng ở thời điểm cũ
    real :: SW_old
    real :: LWs_old
    real :: E_old
   
    !Các biến năng lượng ở thời điểm mới, sau đó được đưa vào rhs(1)
    real :: SW_new
    real :: LWs_guess
    real :: E_new_guess

    real :: max_error  !Dùng để kiểm tra Picard đã hội tụ chưa
    integer :: i       !chỉ số node
    integer :: iter    !số thứ tự vòng lặp Picard

    ! -------------------------------------------------------
    ! Các hệ số
    ! -------------------------------------------------------

    r = K * dt / (dz * dz)
    q = 0.5 * r             ! Do Crank–Nicolson lấy trung bình giữa

    rho_cp_dz = rho_soil * Cp_soil * dz

    ! Hệ số:
    !
    ! dt / (2 rho Cp dz)
    !
    s = dt / (2.0 * rho_cp_dz)

    ! E tại thời điểm cũ (thời điểm ngay trước bước thời gian đang tính)
    call calculate_E_version2(temp0(1), time, &
                              SW_old, LWs_old, E_old)

    ! Dự đoán ban đầu:
    ! lấy nhiệt độ ở thời điểm cũ làm nghiệm đoán đầu tiên cho thời điểm mới
    do i = 1, n
      temp_guess(i) = temp0(i)
    end do

    ! =======================================================
    ! PICARD ITERATION
    ! =======================================================

    do iter = 1, max_iter

      ! Tính SW tại thời điểm mới.
      ! SW chỉ phụ thuộc vào thời gian.
      call calculate_E_version2(temp_guess(1), time + dt, &
                                SW_new, LWs_guess, E_new_guess)

      ! -----------------------------------------------------
      ! Lập ma trận cho node 1 đến node n-1
      !
      ! Chỉ có n-1 ẩn số:
      ! T1, T2, ..., T(n-1)
      !
      ! T(n) là biên Dirichlet = 25°C.
      ! -----------------------------------------------------

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
      !
      ! E_new được tính từ nghiệm đoán Picard.
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
      ! Điều kiện biên sâu tại node n
      !
      ! Ở phương trình node n-1, đưa Tn sang vế phải.
      ! Vì Tn ở cả hai thời điểm bằng 25°C,
      ! đóng góp là q*Tn_old + q*Tn_new.
      ! -----------------------------------------------------

      rhs(n-1) = rhs(n-1) &
               + q * (25.0 + 273.15) &
               + q * (25.0 + 273.15)

      ! Không có hệ số phía trên ở dòng cuối
      c(n-1) = 0.0
      
      ! ==========================================================
! HỆ PHƯƠNG TRÌNH 3 ĐƯỜNG CHÉO
!
! a(i)*T(i-1) + b(i)*T(i) + c(i)*T(i+1) = rhs(i)
!
! Trong đó:
!   a(i)   : hệ số của nhiệt độ tại node phía trên i-1
!   b(i)   : hệ số của nhiệt độ tại chính node i
!   c(i)   : hệ số của nhiệt độ tại node phía dưới i+1
!   rhs(i) : vế phải của phương trình tại node i
!
!   T(i-1) : nhiệt độ tại node phía trên
!   T(i)   : nhiệt độ tại node đang xét
!   T(i+1) : nhiệt độ tại node phía dưới
!
! Các hệ số a, b, c và rhs được đưa vào thuật toán
! Thomas để giải ra T(i) ở bước thời gian mới.
! ==========================================================

      ! -----------------------------------------------------
      ! Giải hệ ba đường chéo bằng Thomas
      ! -----------------------------------------------------

      call thomas(n-1, a, b, c, rhs, solution)

      ! Gán nghiệm vào temp_new
      do i = 1, n-1
        temp_new(i) = solution(i)
      end do

      temp_new(n) = 25.0 + 273.15

      ! -----------------------------------------------------
      ! Kiểm tra hội tụ Picard
      ! -----------------------------------------------------

      max_error = 0.0
      ! tìm sai số lớn nhất trong 25 node
      do i = 1, n
        max_error = max(max_error, &
                        abs(temp_new(i) - temp_guess(i)))
      end do     ! max(a,b): lấy giá trị lớn hơn
                 ! abs() : dấu giá trị tuyệt đối
      ! Cập nhật nghiệm đoán
      do i = 1, n
        temp_guess(i) = temp_new(i)
      end do

      if (max_error < tolerance) exit

    end do

    ! Có thể bỏ phần này nếu không muốn cảnh báo
    if (iter >= max_iter .and. max_error >= tolerance) then
      write(*,'(A,I5,A,E12.5)') &
           'Canh bao: Picard chua hoi tu sau ', &
           max_iter, ' vong lap. Sai so = ', max_error
    endif

    ! Trả nghiệm mới
    do i = 1, n
      temp1(i) = temp_new(i)
      temp0(i) = temp_new(i)
    end do

  end subroutine update_temperature_cn


  ! =========================================================
  ! THOMAS ALGORITHM
  !
  ! Giải hệ:
  !
  ! a(i)x(i-1) + b(i)x(i) + c(i)x(i+1) = d(i)
  !
  ! Với ma trận ba đường chéo.
  ! =========================================================

  subroutine thomas(nsys, a, b, c, d, x)
    implicit none

    integer, intent(in) :: nsys    !số phương trình của hệ cần giải

    real, intent(in) :: a(nsys)    !hệ số nằm dưới đường chéo chính tương ứng với nhiệt độ của node phía trên: T(i-1)
    real, intent(in) :: b(nsys)    !hệ số trên đường chéo chính, là hệ số của ẩn x(i)
    real, intent(in) :: c(nsys)    !hệ số nằm trên đường chéo chính, tương ứng với nhiệt độ của node phía dưới: T(i+1)
    real, intent(in) :: d(nsys)    !vế phải của phương trình thứ i

    real, intent(out) :: x(nsys)   !nghiệm của hệ phương trình
! Sau khi Thomas giải xong:
! x(1) = nhiệt độ node 1 ở thời điểm mới
! x(2) = nhiệt độ node 2 ở thời điểm mới
! ...
! x(24) = nhiệt độ node 24 ở thời điểm mới

    ! mảng trung gian của thuật toán Thomas
    real :: c_star(nsys)  !hệ số c sau khi đã biến đổi trong bước khử xuôi, dùng để loại bỏ các phần tử bên dưới đường chéo chính.
    real :: d_star(nsys)  !vế phải d sau khi đã biến đổi trong bước khử xuôi, ! Sau khi khử xuôi, hệ phương trình được biến đổi thành dạng dễ tính nghiệm bằng thế ngược

    real :: denominator   !mẫu số trung gian trong quá trình khử xuôi

    integer :: i

    ! -------------------------------------------------------
    ! Khử xuôi, biến hệ phương trình ban đầu thành dạng đơn giản hơn
    ! -------------------------------------------------------

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

    ! -------------------------------------------------------
    ! Thế ngược, từ phương trình cuối cùng tính ngược lên để tìm toàn bộ nghiệm
    ! -------------------------------------------------------

    x(nsys) = d_star(nsys)

    do i = nsys-1, 1, -1
      x(i) = d_star(i) - c_star(i) * x(i+1)
    end do

  end subroutine thomas


  ! =========================================================
  ! GHI KẾT QUẢ
  ! =========================================================

  subroutine write_results(unit_out, time, temp, SW, LWs, E)
    implicit none

    integer, intent(in) :: unit_out

    real, intent(in) :: time
    real, intent(in) :: temp(:)
    real, intent(in) :: SW
    real, intent(in) :: LWs
    real, intent(in) :: E

    integer :: day_number

    day_number = int(time / (24.0 * 3600.0)) + 1

    write(*,'(F7.2,A,I3,A,F9.3,A,F10.3,A,F10.3,A,F10.3,A,F9.3,A,F9.3)') &
         time / 3600.0, ' | ', &
         day_number, ' | ', &
         temp(1) - 273.15, ' | ', &
         SW, ' | ', &
         LWs, ' | ', &
         E, ' | ', &
         temp(3) - 273.15, ' | ', &
         temp(25) - 273.15

    write(unit_out, &
          '(F10.3,3X,I6,3X,F10.3,3X,F12.3,3X,F12.3,3X,F12.3,3X,F10.3,3X,F10.3)') &
         time / 3600.0, &
         day_number, &
         temp(1) - 273.15, &
         SW, &
         LWs, &
         E, &
         temp(3) - 273.15, &
         temp(25) - 273.15

  end subroutine write_results

end program test2