program test1
  implicit none
  
  integer, parameter :: n = 25
  real :: temp0(n), temp1(n)
  real :: K              ! Hệ số khuếch tán nhiệt [m²/s]
  real :: dt, dz, time, endtime
  real :: E
  integer :: unit_out
  integer :: hour_display

  ! Khai báo hằng số
  K = 3e-7               ! Hệ số khuếch tán nhiệt [m²/s]
  dt = 30.0              ! Bước thời gian [s]
  dz = 0.04              ! Bước lưới độ sâu [m]
  time = 0.0             ! Thời gian ban đầu [s]
  endtime = 48.0 * 3600.0  ! 48 giờ [s]
  
  ! Gọi subroutine khởi tạo
  call init(temp0, temp1, n)

! =========================================================
! MỞ FILE ĐỂ LƯU KẾT QUẢ
! =========================================================
unit_out = 10

open(unit=unit_out, file='result_version1.txt', &
     status='replace', action='write')

write(unit_out,'(A8,3X,A10,3X,A12,3X,A12,3X,A12)') &
     'Time(h)', 'Ts(C)', 'E(W/m2)', 'T8cm(C)', 'T96cm(C)'

 ! In tiêu đề bảng kết quả
  write(*,*) ''
  write(*,*) '===== VERSION 1: TINH E NGUOC ====='
  write(*,*) 'LOGIC BAI TOAN NGUOC:'
  write(*,*) '  Ts(t) da biet -> tinh E -> cap nhat T cac lop sau'
  write(*,*) ''
  write(*,*) 'Du lieu dau vao:'
  write(*,*) '  - Ts(t) = 25 + 10*sin(omega(t - 6AM)) do C'
  write(*,*) '  - T(z=L) = 25 do C (hang so)'
  write(*,*) '  - K = 3e-7 m^2/s (he so khuech tan)'
  write(*,*) '  - r = K*dt/dz^2 = 0.005625 (on dinh < 0.5)'
  write(*,*) ''
  write(*,'(A)') 'Time(h) |  Ts(C)  |     E      | T tai 8cm(do C) | T tai 96cm(do C)'
  write(*,'(A)') '========|=========|============|=================|================='
  
  ! ========== VÒNG LẶP CHÍNH - 48 giờ ==========
  do while (time <= endtime)
    
    ! ========== BƯỚC 1: Gán điều kiện biên bề mặt tại thời điểm hiện tại ==========
    temp0(1) = tsfc(time)
    hour_display = mod(6 + int(time / 3600.0), 24)

    ! ========== BƯỚC 2: TÍNH E NGƯỢC tại thời điểm hiện tại ==========
    E = calculate_E_inverse(temp0, K, dz, time)
    
    ! ========== BƯỚC 3: IN KẾT QUẢ tại đúng thời điểm time ==========
    ! In TRƯỚC khi cập nhật, để nhãn thời gian chính xác
    ! Tại đây: temp0(1) = Ts(time), E = E(time), T = T(time)
    if (mod(int(time), 3600) == 0) then
      write(*,'(I6,A,F7.2,A,F11.1,A,F8.2,A,F8.2)') &
        hour_display, ' | ', temp0(1)-273.15, ' | ', E, ' | ', &
        temp0(3)-273.15, ' | ', temp0(25)-273.15
    endif

! Ghi kết quả vào file TXT
if (mod(int(time), 3600) == 0) then
  write(unit_out,'(I6,3X,F10.3,3X,F12.3,3X,F12.3,3X,F12.3)') &
    hour_display, &
    temp0(1)-273.15, &
    E, &
    temp0(3)-273.15, &
    temp0(25)-273.15
endif    
    ! ========== BƯỚC 4: Kiểm tra nếu đã đến cuối thời gian ==========
    if (time >= endtime) exit
    
    ! ========== BƯỚC 5: Cập nhật nhiệt độ các tầng sâu sang bước tiếp theo ==========
    ! Không dùng E trong cập nhật, chỉ khuếch tán
    call update_temperature(temp0, temp1, K, dt, dz, n, time)
    
    ! ========== BƯỚC 6: Tăng thời gian lên dt ==========
    time = time + dt
    
  end do
  
close(unit_out)
contains

  subroutine init(temp0, temp1, n)
    ! Khởi tạo nhiệt độ ban đầu
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: temp0(n), temp1(n)
    integer :: i
    
    do i = 1, n
      temp0(i) = 25.0 + 273.15   ! 25 độ C = 298.15 K
      temp1(i) = 25.0 + 273.15
    end do
    
    write(*,*) 'Nhiet do ban dau: ', temp0(1), ' K (25 do C)'
    
  end subroutine init


  function tsfc(time) result(temp_surface)
    !====================================================================
    ! Điều kiện biên bề mặt (BIẾT)
    ! Ts(t) = 25 + 10*sin(omega(t - 6AM)) (độ C)
    ! Bắt đầu từ 6AM với biên độ 10 độ C, trung bình 25 độ C
    !====================================================================
    
    implicit none
    real, intent(in) :: time
    real :: temp_surface
    real, parameter :: a0 = 10.0        ! Biên độ dao động = 10 độ C
    real, parameter :: T_mean = 25.0    ! Nhiệt độ trung bình = 25 độ C
    real, parameter :: pi = 4.0 * atan(1.0)
    real :: omega, time_from_6am
    
    ! Tần số góc: omega = 2pi/T (T = 48 giờ)
    omega = (2.0 * pi) / (24.0 * 3600.0)
    
    ! Thời gian tính từ 6AM
    time_from_6am = time
    
    ! Công thức dao động: Ts(t) = 25 + 10*sin(omega(t - 6AM)) (độ C)
    ! Chuyển sang Kelvin
    temp_surface = a0 * sin(omega * time_from_6am) + T_mean + 273.15
    
  end function tsfc


  function calculate_E_inverse(temp0, K, dz, time) result(E)
    !====================================================================
    ! TÍNH E NGƯỢC (Inverse Problem)
    !
    ! Version 1: E là kết quả "chẩn đoán" (diagnosed value)
    !
    ! Phương trình truyền nhiệt tại bề mặt:
    !   dT/dt|z=0 = c|z=0 + E/(ρ*Cp*dz)
    !
    ! Giải ngược để tính E [W/m²]:
    !   E = (ρ*Cp*dz) * [dT/dt|z=0 - K * d^2T/dz^2|z=0]
    !
    ! Ý nghĩa: Để bề mặt dao động theo Ts(t) = 25 + 10*sin(omega(t - 6AM)) (độ C)
    !          trong điều kiện khuếch tán K * d^2T/dz^2,
    !          cần có nguồn năng lượng bề mặt E bằng bao nhiêu.
    !====================================================================
    
    implicit none
    real, intent(in) :: temp0(:)  ! Mảng nhiệt độ hiện tại [K]
    real, intent(in) :: K         ! Hệ số khuếch tán [m²/s]
    real, intent(in) :: dz        ! Bước lưới độ sâu [m]
    real, intent(in) :: time      ! Thời gian hiện tại [s]
    real :: E
    
    real :: dTdt_at_surface       ! dT/dt tại z=0 [K/s]
    real :: d2Tdz2_at_surface     ! d^2T/dz^2 tại z=0 [K/m^2]
    real :: rho_cp_dz             ! (ρ*Cp*dz) [J/(m^2·K)]
    real :: a0, omega, pi, time_from_6am
    
    ! ========== HẰNG SỐ VẬT LÝ ==========
    ! ρ*Cp của đất: 2.0e6 J/(m³·K)
    ! (Giả định: đất khô ~ 1.8-2.2 MJ/(m³·K))
    ! Nhân với dz để chuyển sang [J/(m²·K)]
    rho_cp_dz = 2.0e6 * dz
    
    ! ========== TÍNH dT/dt TẠI BỀ MẶT (Đạo hàm giải tích) ==========
    ! Ts(t) = 25 + 10*sin(omega(t - 6AM)) (độ C)
    ! dTs/dt = 10*ω*cos(omega(t - 6AM)) [K/s]
    
    a0 = 10.0
    pi = 4.0 * atan(1.0)
    omega = (2.0 * pi) / (24.0 * 3600.0)
    time_from_6am = time
    
    dTdt_at_surface = a0 * omega * cos(omega * time_from_6am)
    
    ! ========== TÍNH d²T/dz² TẠI BỀ MẶT (Sai phân hữu hạn ở biên) ==========
    ! Công thức sai phân trung tâm đầy đủ: d²T/dz² = (T(2) - 2*T(1) + T(0))/dz²
    ! Tại biên z=0, sử dụng xấp xỉ:
    ! d²T/dz² ≈ [T(z=dz) - T(z=0)] / dz²
    ! 
    ! Đây là xấp xỉ biên (boundary approximation), độ chính xác O(dz)
    ! (không phải O(dz²) như sai phân trung tâm)
    
    d2Tdz2_at_surface = (temp0(2) - temp0(1)) / (dz * dz)
    
    ! ========== GIẢI NGƯỢC: TÍNH E ==========
    ! E = (ρ*Cp*dz) * [dT/dt - K * d²T/dz²]
    E = rho_cp_dz * (dTdt_at_surface - K * d2Tdz2_at_surface)
    
  end function calculate_E_inverse


  subroutine update_temperature(temp0, temp1, K, dt, dz, n, time)
    !====================================================================
    ! CẬP NHẬT NHIỆT ĐỘ BẰNG KHUẾCH TÁN (Không dùng E)
    !
    ! Version 1: E là kết quả "chẩn đoán", không dùng để dự báo T
    !
    ! Phương trình: dT/dt = K * d²T/dz² (chỉ khuếch tán)
    !
    ! Sơ đồ sai phân rõ ràng (Explicit Finite Difference):
    !   T(i, t+dt) = T(i, t) + r*[T(i-1, t) - 2*T(i, t) + T(i+1, t)]
    !
    ! Với r = K*dt/dz² = 0.005625 (số Fourier)
    ! Điều kiện ổn định: r ≤ 0.5 
    !
    ! Điều kiện biên:
    !   - Bề mặt (i=1): T(1, t+dt) = Ts(t+dt)
    !   - Tầng sâu (i=n): T(n) = 25°C (hằng số)
    !
    ! Node index → Độ sâu [cm]:
    !   - node 1 = 0 cm (bề mặt)
    !   - node 2 = 4 cm
    !   - node 3 = 8 cm
    !   - node 6 = 20 cm
    !====================================================================
    
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: temp0(n), temp1(n)
    real, intent(in) :: K        ! Hệ số khuếch tán [m²/s]
    real, intent(in) :: dt       ! Bước thời gian [s]
    real, intent(in) :: dz       ! Bước lưới độ sâu [m]
    real, intent(in) :: time     ! Thời gian hiện tại [s]
    
    integer :: i
    real :: r           ! Số Fourier = 0.005625
    
    ! ========== TÍNH HỆ SỐ FOURIER ==========
    ! r = K * dt / dz²
    r = K * dt / (dz * dz)
    
    ! ========== CẬP NHẬT BỀ MẶT (i=1) ==========
    ! Điều kiện biên Dirichlet: T được áp đặt theo Ts(t+dt)
    temp1(1) = tsfc(time + dt)
    
    ! ========== CẬP NHẬT CÁC LỚP TRUNG GIAN (i=2..n-1) ==========
    ! Công thức khuếch tán:
    ! T(i, t+dt) = T(i, t) + r*[T(i-1, t) - 2*T(i, t) + T(i+1, t)]
    
    do i = 2, n-1
      temp1(i) = temp0(i) + r * (temp0(i-1) - 2.0*temp0(i) + temp0(i+1))
    end do
    
    ! ========== CẬP NHẬT LỚP TẦNG SÂU (i=n) ==========
    ! Điều kiện biên Dirichlet: T = 25°C hằng số
    temp1(n) = 25.0 + 273.15
    
    ! ========== HOÁN ĐỔI MẢNG ==========
    ! Sao chép mảng temp1 vào temp0 cho bước tiếp theo
    do i = 1, n
      temp0(i) = temp1(i)
    end do
    
  end subroutine update_temperature

end program test1