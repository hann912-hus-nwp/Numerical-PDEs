program test2
  implicit none
  
  integer, parameter :: n = 25
  real :: temp0(n), temp1(n)
  real :: K              ! Hệ số khuếch tán nhiệt [m²/s]
  real :: dt, dz, time, endtime
  real :: E, SW, LWs
  real :: rho_soil, Cp_soil  ! Tính chất đất
  real :: r               ! Số Fourier
  integer :: unit_out
  
  ! Khai báo hằng số
  K = 3e-7               ! Hệ số khuếch tán nhiệt [m²/s]
  dt = 30.0              ! Bước thời gian [s]
  dz = 0.04              ! Bước lưới độ sâu [m]
  time = 0.0            
  endtime = 24.0 * 3600.0  ! 24 giờ [s]
  
  ! Tính chất đất (giá trị lựa chọn đại diện)
  ! rho_soil = 2000 kg/m³, Cp_soil = 1000 J/(kg·K)
  ! => rho*Cp = 2.0×10⁶ J/(m³·K)
  rho_soil = 2000.0
  Cp_soil  = 1000.0
  
  
  ! Tính số Fourier (kiểm tra ổn định)
  r = K * dt / (dz * dz)
  
  ! Gọi subroutine khởi tạo
  call init(temp0, temp1, n)

! =========================================================
! MỞ FILE ĐỂ LƯU KẾT QUẢ
! =========================================================
unit_out = 10

open(unit=unit_out, file='result_version2.txt', &
     status='replace', action='write')

! Ghi tên các cột vào file TXT
write(unit_out,'(A8,3X,A10,3X,A12,3X,A12,3X,A12,3X,A10,3X,A10)') &
     'Time(h)', 'Ts(C)', 'SW(W/m2)', 'LWs(W/m2)', &
     'E(W/m2)', 'T8cm(C)', 'T96cm(C)'

  ! In tiêu đề bảng kết quả
  write(*,*) ''
  write(*,*) '===== VERSION 2: E = SW - LWs (voi E tac dong vao be mat) ====='
  write(*,*) 'LOGIC BAI TOAN:'
  write(*,*) '  Tinh SW, LWs -> tinh E -> E anh huong den Ts'
  write(*,*) '  Ts khong duoc ap dat, ma tu thay doi theo nang luong'
  write(*,*) ''
  write(*,*) 'Dieu kien ban dau:'
  write(*,*) '  - T(z,0) = 25 do C cho toan bo lop dat'
  write(*,*) '  - T(z=L) = 25 do C (hang so, bien sau)'
  write(*,*) ''
  write(*,*) 'Tinh chat dat (gia dinh dai dien):'
  write(*,*) '  - p_soil = 2000 kg/m^3'
  write(*,*) '  - Cp_soil = 1000 J/(kg.K)'
  write(*,*) '  - p*Cp = 2.0×10⁶ J/(m^3.K)'
  write(*,*) ''
  write(*,*) 'Thong so mo phong:'
  write(*,'(A,E12.5,A)') '  - K = ', K, ' m^2/s'
  write(*,'(A,F6.2,A)') '  - delta t = ', dt, ' s'
  write(*,'(A,F6.4,A)') '  - delta z = ', dz, ' m'
  write(*,'(A,F8.6,A)') '  - r = K.dellta t/ delta z^2 = ', r, ' (on dinh < 0.5)'
  write(*,*) ''
  write(*,'(A)') 'Time(h) |  Ts(do C)  |    SW      |    LWs     |      E     | T8cm( do C) | T96cm( do C)'
  write(*,'(A)') '========|=========|============|============|============|===========|==========='
  
  ! ========== VÒNG LẶP CHÍNH - 24 giờ ==========
  do while (time <= endtime)
    
    ! Tính thời gian hiện tại
    temp0(1) = tsfc(time)
    
    ! ========== BƯỚC 1: TÍNH E từ cân bằng bức xạ ==========
    call calculate_E_version2(temp0, time, SW, LWs, E)
    
    ! ========== BƯỚC 2: IN KẾT QUẢ mỗi 1 giờ ==========
    if (mod(int(time), 3600) == 0) then
      write(*,'(F6.1,A,F8.2,A,F11.2,A,F11.2,A,F11.2,A,F9.2,A,F9.2)') &
        time/3600.0, ' | ', temp0(1)-273.15, ' | ', SW, ' | ', LWs, ' | ', E, ' | ', &
        temp0(3)-273.15, ' | ', temp0(25)-273.15
    endif

  ! Ghi kết quả vào file TXT
if (mod(int(time), 3600) == 0) then
  write(unit_out,'(F8.2,3X,F10.3,3X,F12.3,3X,F12.3,3X,F12.3,3X,F10.3,3X,F10.3)') &
    time/3600.0, &
    temp0(1)-273.15, &
    SW, &
    LWs, &
    E, &
    temp0(3)-273.15, &
    temp0(25)-273.15
endif
    
    ! ========== BƯỚC 3: Kiểm tra nếu đã đến cuối thời gian ==========
    if (time >= endtime) exit
    
    ! ========== BƯỚC 4: Cập nhật nhiệt độ các tầng ==========
    ! Node 1 (bề mặt): khuếch tán + E
    ! Node 2→n-1: chỉ khuếch tán
    ! Node n: T = 25°C (biên sâu)
    call update_temperature(temp0, temp1, K, dt, dz, n, E, rho_soil, Cp_soil)
    
 ! ========== BƯỚC 5: Tăng thời gian lên dt ==========
    time = time + dt

  end do
  
  write(*,'(A)') '========|=========|============|============|============|===========|==========='
  write(*,'(A,F6.1,A)') 'Ket thuc mo phong tai time = ', time/3600.0, ' gio'
  
  close(unit_out)
contains

  subroutine init(temp0, temp1, n)
    ! Khởi tạo nhiệt độ ban đầu: T(z,0) = 25°C cho toàn bộ lớp đất
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: temp0(n), temp1(n)
    integer :: i
    
    do i = 1, n
      temp0(i) = 25.0 + 273.15   ! 25°C = 298.15 K
      temp1(i) = 25.0 + 273.15
    end do
    
    write(*,*) 'Nhiệt độ ban đầu: T(z,0) = 25°C (298.15 K)'
    
  end subroutine init


  subroutine calculate_E_version2(temp0, time, SW, LWs, E)
    !====================================================================
    ! TÍNH E - VERSION 2: E = SW - LWs
    !
    ! Phương trình cân bằng năng lượng bề mặt:
    !   E = SW - LWs
    !
    ! Trong đó:
    !   SW: Bức xạ sóng ngắn [W/m²]
    !       Thay đổi theo thời gian trong ngày (6AM - 6PM)
    !       Cực đại ~1000 W/m² vào lúc trưa
    !       Bằng 0 vào đêm
    !
    !   LWs: Bức xạ sóng dài từ bề mặt (Stefan-Boltzmann)
    !       LWs = ε·σ·Ts⁴
    !       ε = 0.95 (hệ số phát xạ)
    !       σ = 5.67×10⁻⁸ W/(m²·K⁴) (hằng số Stefan-Boltzmann)
    !       Ts: Nhiệt độ bề mặt hiện tại [K]
    !
    ! E > 0: Năng lượng vào bề mặt (gia nhiệt)
    ! E < 0: Năng lượng thoát ra (mất nhiệt)
    !
    ! Version 2: E tác động vào phương trình truyền nhiệt tại bề mặt
    !====================================================================
    
    implicit none
    real, intent(in) :: temp0(:)  ! Mảng nhiệt độ [K]
    real, intent(in) :: time      ! Thời gian hiện tại [s]
    real, intent(out) :: SW       ! Bức xạ sóng ngắn [W/m²]
    real, intent(out) :: LWs      ! Bức xạ sóng dài từ bề mặt [W/m²]
    real, intent(out) :: E        ! Thông lượng năng lượng bề mặt [W/m²]
    
    real :: Ts                    ! Nhiệt độ bề mặt [K]
    real :: pi, omega, time_from_6am
    real, parameter :: sigma = 5.67e-8     ! Stefan-Boltzmann [W/(m²·K⁴)]
    real, parameter :: epsilon = 0.95      ! Hệ số phát xạ đất
    real, parameter :: SW_max = 1000.0     ! Bức xạ sóng ngắn cực đại [W/m²]
    
    ! ========== LẤY NHIỆT ĐỘ BỀ MẶT HIỆN TẠI ==========
    Ts = temp0(1)  ! [K]
    
    ! ========== TÍNH SW - BỨC XẠ SÓNG NGẮN ==========
    ! Mô hình đơn giản: SW = SW_max·sin(ω·(t - 6AM)) nếu 6AM ≤ t ≤ 6PM
    ! SW = 0 vào đêm
    ! Đạt cực đại tại lúc 12:00 (t = 6h từ 6AM)
    
    pi = 4.0 * atan(1.0)
    omega = pi / (12.0 * 3600.0)
    time_from_6am = time
    
    ! Kiểm tra thời gian có mặt trời (6AM - 6PM = 12 giờ)
    if (time >= 0.0 .and. time <= 12.0 * 3600.0) then
      SW = SW_max * sin(omega * time)
      if (SW < 0.0) SW = 0.0
    else
      SW = 0.0
    endif
    
    ! ========== TÍNH LWs - BỨC XẠ SÓNG DÀI ==========
    ! Phương trình Stefan-Boltzmann:
    ! LWs = ε·σ·Ts⁴
    LWs = epsilon * sigma * (Ts ** 4.0)
    
    ! ========== TÍNH E ==========
    ! E = SW - LWs
    E = SW - LWs
    
  end subroutine calculate_E_version2


  subroutine update_temperature(temp0, temp1, K, dt, dz, n, E, rho_soil, Cp_soil)
    !====================================================================
    ! CẬP NHẬT NHIỆT ĐỘ - VERSION 2
    !
    ! Phương trình truyền nhiệt:
    !   dT/dt = K·d²T/dz² + E/(ρ·Cp·Δz)  [tại z=0, node 1]
    !   dT/dt = K·d²T/dz²                  [tại z > 0, node 2 → n-1]
    !
    ! Sơ đồ sai phân rõ ràng (Explicit Finite Difference):
    !
    ! Node 1 (bề mặt): 
    !   T(1, t+Δt) = T(1, t) + Δt·[K·(T(2) - T(1))/Δz² + E/(ρ·Cp·Δz)]
    !   
    !   Giải thích: 
    ! - Thành phần khuếch tán phụ thuộc vào chênh lệch nhiệt độ giữa node 1 và node 2
    !   - Dòng năng lượng bề mặt E: được chia đều trên thể tích lớp
    !
    ! Node 2 → n-1 (trung gian):
    !   T(i, t+Δt) = T(i, t) + r·[T(i-1, t) - 2·T(i, t) + T(i+1, t)]
    !   với r = K·Δt/Δz² (số Fourier, phải < 0.5)
    !
    ! Node n (biên sâu):
    !   T(n) = 25°C = hằng số
    !
    ! Điểm khác biệt so với Version 1:
    ! - Không áp đặt Ts(t) = tsfc(time)
    ! - Node 1 được cập nhật bằng công thức với E
    ! - Ts tự thay đổi theo năng lượng (E > 0 → nóng lên, E < 0 → lạnh đi)
    !====================================================================
    
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: temp0(n), temp1(n)
    real, intent(in) :: K              ! Hệ số khuếch tán [m²/s]
    real, intent(in) :: dt             ! Bước thời gian [s]
    real, intent(in) :: dz             ! Bước lưới độ sâu [m]
    real, intent(in) :: E              ! Thông lượng năng lượng bề mặt [W/m²]
    real, intent(in) :: rho_soil       ! Mật độ đất [kg/m³]
    real, intent(in) :: Cp_soil        ! Nhiệt dung riêng [J/(kg·K)]
    
    integer :: i
    real :: r           ! Số Fourier = K·Δt/Δz²
    real :: rho_cp_dz   ! Tích ρ·Cp·Δz [J/(m²·K)]
    
    ! ========== TÍNH HỆ SỐ FOURIER ==========
    r = K * dt / (dz * dz)
    
    ! ========== TÍNH ρ·Cp·Δz ==========
    rho_cp_dz = rho_soil * Cp_soil * dz
    
    ! ========== CẬP NHẬT NODE 1 (BỀ MẶT) ==========
    ! Công thức: T1^(n+1) = T1^n + Δt·[K·(T2 - T1)/Δz² + E/(ρ·Cp·Δz)]
    ! 
    ! Diễn giải:
    ! - K·(T2 - T1)/Δz²: thành phần khuếch tán nhiệt phụ thuộc vào chênh lệch nhiệt độ giữa node 1 và node 2
    ! - E/(ρ·Cp·Δz): tác động của thông lượng năng lượng
    
        
    temp1(1) = temp0(1) + dt * ( &
        K * (temp0(2) - temp0(1)) / (dz * dz) + &
        E / rho_cp_dz )
    
    ! ========== CẬP NHẬT NODE 2 → n-1 (TRUNG GIAN) ==========
    ! Chỉ khuếch tán: T(i, t+Δt) = T(i, t) + r·[T(i-1, t) - 2·T(i, t) + T(i+1, t)]
    
    do i = 2, n-1
      temp1(i) = temp0(i) + r * (temp0(i-1) - 2.0*temp0(i) + temp0(i+1))
    end do
    
    ! ========== CẬP NHẬT NODE n (BIÊN SÂU) ==========
    ! Điều kiện biên Dirichlet: T = 25°C (hằng số)
    temp1(n) = 25.0 + 273.15
    
    ! ========== HOÁN ĐỔI MẢNG ==========
    ! Sao chép temp1 → temp0 cho bước tiếp theo
    do i = 1, n
      temp0(i) = temp1(i)
    end do
    
  end subroutine update_temperature

end program test2
