program shallow_water
  implicit none

  integer, parameter :: n = 40           
  integer, parameter :: na = n + 1       ! Kích thước mảng hai chiều của lưới A
  integer, parameter :: nb = n + 1       ! Kích thước mảng hai chiều của lưới C
  integer, parameter :: nt = 1800        ! 30 giờ / dt = 108000 s / 60 s = 1800 bước
  integer, parameter :: nout = 60        ! Ghi kết quả mỗi 60 bước = 1 giờ (để chọn vẽ tại t = 1, 6, 12, ..., 30 giờ)
  real(8), parameter :: km = 1.0d-4      ! Hệ số khuếch tán
  real(8), parameter :: dx=30000.0d0     ! Khoảng cách lưới theo phương x = 30 km
  real(8), parameter :: dy=30000.0d0     ! Khoảng cách lưới theo phương y = 30 km
  real(8), parameter :: dt=60.0d0        ! Bước thời gian
  real(8), parameter :: g=9.81d0
  real(8), parameter :: omega=7.2921150d-5
  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: phi0=35.0d0*pi/180.0d0     ! 35 độ
  real(8), parameter :: f=2.0d0*omega*sin(phi0)    ! Tham số Coriolis

  real(8), parameter :: tend=30.0d0*3600.0d0       ! Tổng thời gian mô phỏng = 30 giờ
  real(8), parameter :: z0=3000.0d0                ! Độ cao mặt ban đầu
  ! [MỚI] Thời gian "kéo" Z về Z mục tiêu trong số hạng nguồn S (s).
  ! Đặt rất lớn (ví dụ 1.0d30) nếu chỉ muốn S = dZ_mục_tiêu/dt (không kéo).
  real(8), parameter :: tau_s=3600.0d0

  ! [MỚI] z_mode = 1: Z ÁP ĐẶT Z*(x,y,t)  -> NaN của gió KHÔNG lan sang Z
  ! z_mode = 2: Z DỰ BÁO Z_t + (Zu)_x + (Zv)_y = S  -> NaN của gió sẽ lan sang Z

  integer, parameter :: z_mode = 1       !đổi sang biến dự báo thì đặt z_mode =2

  ! [MỚI] Tỉ lệ sai phân UPWIND trong số hạng đối lưu của LƯỚI A (0 = trung tâm thuần túy, 1 = upwind đầy đủ).
  ! Lưới C luôn dùng upwind đầy đủ. alpha_a nhỏ -> lưới A ít được khuếch tán số trị -> mất ổn định.
  ! (alpha_a >= ~0.3 thì lưới A cũng ổn định.)
  real(8), parameter :: alpha_a = 0.2d0

  real(8), parameter :: cs(4)=(/0.0d0,0.5d0,0.5d0,1.0d0/)   ! Hệ số RK4
  real(8), parameter :: ct(4)=(/0.0d0,0.5d0,0.5d0,1.0d0/)   ! Hệ số thời điểm của từng giai đoạn

  real(8) :: radius                ! Bán kính vùng trũng
  real(8) :: depth_max             ! Độ sâu cực đại của Z mục tiêu
  real(8) :: zmin                  ! Giá trị độ cao nhỏ nhất
  real(4) :: vmax                  ! Tốc độ gió lớn nhất 
  real(4) :: vrmin                 ! Vận tốc hướng tâm nhỏ nhất 
  real(8) :: time                  ! Thời gian mô phỏng
  integer :: it, is
  integer :: unit_a, unit_c        ! Tệp dữ liệu lưới A / lưới C
  integer :: unit_da, unit_dc      ! Tệp chuỗi thời gian chẩn đoán

 ! ---- Lưới A
real(4) :: ua(na,na), va(na,na)                   ! u, v
real(8) :: za(na,na)                              ! Z 
real(4) :: ut_a(na,na), vt_a(na,na)               ! Gió tạm thời trong RK4
real(8) :: zs_a(na,na)                            ! Z tạm thời
real(4) :: du_a(na,na), dv_a(na,na)               ! Xu thế gió
real(8) :: dz_a(na,na)                            ! Xu thế Z
real(8) :: sa(na,na)                              ! Số hạng nguồn S
real(4) :: kua(na,na,4), kva(na,na,4)             ! Hệ số RK4 của gió
real(8) :: kza(na,na,4)                           ! Hệ số RK4 của Z

 ! ---- Lưới C
real(4) :: uc(nb+1,nb), vc(nb,nb+1)
real(8) :: zc(nb,nb)
real(4) :: ut_c(nb+1,nb), vt_c(nb,nb+1)
real(8) :: zs_c(nb,nb)
real(4) :: du_c(nb+1,nb), dv_c(nb,nb+1)
real(8) :: dz_c(nb,nb)
real(8) :: sc(nb,nb)
real(4) :: kuc(nb+1,nb,4), kvc(nb,nb+1,4)
real(8) :: kzc(nb,nb,4)

 radius = 300000.0d0      ! Bán kính r = 300 km
 depth_max = 150.0d0      ! Biên độ giảm độ cao cực đại (Z mục tiêu) = 150 m

  write(*,'(a)') '=============================================='
  write(*,'(a)') ' SHALLOW-WATER MODEL: A-GRID / C-GRID'
  if z_mode = 2 then
    write(*,'(a)') ' Z la bien du bao'
  else
    write(*,'(a)') ' Z la bien ap dat'
  end if

  write(*,'(a)') '=============================================='

  write(*,'(a)') 'Grid          Zmin(m)       Vmax(m/s)     vrmin(m/s)'

  ! ==================================================================
  ! LƯỚI A
  ! ==================================================================
  open(newunit=unit_a,file='unit3_Agrid.txt',status='replace',action='write')
  write(unit_a,'(a)') 'time_h x_km y_km u v z'
  open(newunit=unit_da,file='diag_Agrid.txt',status='replace',action='write')
  write(unit_da,'(a)') 'time_h zmin vmax vrmin checker_z checker_u'


call init_a(ua,va,za,radius,depth_max)  !Khởi tạo trạng thái lưới A (t = 0)
call output_a(unit_a,0.0d0,ua,va,za)        ! ghi trạng thái ban đầu (t = 0) của lưới A (thời gian, tọa độ, u, v, Z)
call write_diag_a(unit_da,0.0d0,ua,va,za)   ! Ghi dòng chẩn đoán đầu tiên của lưới A vào diag_Agrid.txt

  do it=0,nt-1
    time=dble(it)*dt

   do is=1,4                                   ! Vòng lặp 4 giai đoạn của RK4
  ! Tính trạng thái tạm thời (u, v, Z) cho giai đoạn is: giai đoạn 1 lấy bằng giá trị hiện tại,
  ! giai đoạn 2-4 lấy bằng giá trị hiện tại + cs(is)*dt*(xu thế của giai đoạn trước)
  call forecast_a(ua,va,za,ut_a,vt_a,zs_a,kua,kva,kza,is,dt,cs)

  call init_tend_a(du_a,dv_a,dz_a)     ! Đặt các mảng xu thế du, dv, dz về 0 trước khi tính lại cho giai đoạn này

  ! Tính xu thế của u, v tại các điểm biên (và 4 góc) bằng công thức biên upwind,
  ! giá trị ngoài miền giả định bằng 0
  call boundary_a(ut_a,vt_a,du_a,dv_a)

  if (z_mode == 1) then
    ! Z áp đặt: ghi đè Z tạm thời bằng Z*(x,y,t) tại thời điểm của giai đoạn RK4
    ! (Z không phụ thuộc gió, nên NaN của gió không lan sang Z)
    call calc_z_a(zs_a,time+ct(is)*dt,radius,depth_max)
  else
    ! Z dự báo: tính số hạng nguồn S = dZ*/dt + div(Z* u) + (Z* - Z)/tau_s
    ! để tạo và giữ vùng trũng trong phương trình Z
    call calc_s_a(sa,zs_a,ut_a,vt_a,time+ct(is)*dt,radius,depth_max)
  end if

  call tend_u_a(ut_a,vt_a,zs_a,du_a)      ! Xu thế của u tại các điểm trong miền: đối lưu - g*dZ/dx + f*v + khuếch tán
  call tend_v_a(ut_a,vt_a,zs_a,dv_a)      ! Xu thế của v tại các điểm trong miền: đối lưu - g*dZ/dy - f*u + khuếch tán

  ! Xu thế của Z (chỉ khi Z là biến dự báo): Z_t = -(Zu)_x - (Zv)_y + S
  if (z_mode == 2) call tend_z_a(ut_a,vt_a,zs_a,sa,dz_a)

  kua(:,:,is)=du_a                          ! Lưu hệ số RK4 thứ is của u
  kva(:,:,is)=dv_a                          ! Lưu hệ số RK4 thứ is của v
  kza(:,:,is)=dz_a                          ! Lưu hệ số RK4 thứ is của Z (bằng 0 khi Z áp đặt)
end do

call update_a(ua,va,za,kua,kva,kza,dt)      ! Kết hợp 4 hệ số: x(n+1) = x(n) + dt*(k1 + 2*k2 + 2*k3 + k4)/6 cho u, v, Z

if (mod(it+1,nout)==0) then                 ! Cứ mỗi nout bước (= 1 giờ) thì ghi kết quả
  ! Khi Z áp đặt, za chưa được cập nhật nên phải tính lại Z*(t) trước khi ghi
  if (z_mode == 1) call calc_z_a(za,dble(it+1)*dt,radius,depth_max)

  call output_a(unit_a,dble(it+1)*dt,ua,va,za)     ! Ghi thời gian, tọa độ, u, v, Z của toàn lưới A 

  ! Ghi một dòng chẩn đoán (Zmin, Vmax, vrmin, chỉ số nhiễu bàn cờ) 
  call write_diag_a(unit_da,dble(it+1)*dt,ua,va,za)
end if
end do
close(unit_a)
close(unit_da)

! Chẩn đoán tại thời điểm cuối
if (z_mode == 1) call calc_z_a(za,tend,radius,depth_max)     ! Đặt lại Z = Z*(tend) nếu Z áp đặt

call diag_a(ua,va,za,zmin,vmax,vrmin)     ! Tính Zmin, Vmax và vận tốc hướng tâm nhỏ nhất vrmin

call output_profile_a(za)   ! Xuất mặt cắt Z dọc hàng y = 0 
call output_2d_a(za)       ! Xuất trường Z hai chiều tại thời điểm cuối 
  write(*,'(a8,3(2x,es13.5))') 'A-grid',zmin,vmax,vrmin

  ! ==================================================================
  ! LƯỚI C
  ! ==================================================================
  open(newunit=unit_c,file='unit3_Cgrid.txt',status='replace',action='write')
  write(unit_c,'(a)') 'time_h x_km y_km u v z'
  open(newunit=unit_dc,file='diag_Cgrid.txt',status='replace',action='write')
  write(unit_dc,'(a)') 'time_h zmin vmax vrmin checker_z checker_u'

  ! Lưới C: đặt u = v = 0 (u ở mặt tây/đông của ô, v ở mặt nam/bắc) và Z = 3000 m (tại tâm ô)
call init_c(uc,vc,zc,radius,depth_max)

call output_c(unit_c,0.0d0,uc,vc,zc)   ! Ghi trạng thái t = 0 của lưới C (u, v được nội suy về tâm ô khi ghi)

call write_diag_c(unit_dc,0.0d0,uc,vc,zc)   ! Ghi dòng chẩn đoán đầu tiên của lưới C 

do it=0,nt-1                                  ! Vòng lặp thời gian: nt bước, mỗi bước dt = 60 s
  time=dble(it)*dt

  do is=1,4                                   ! 4 giai đoạn của RK4
    ! Tính trạng thái tạm thời (u, v, Z) của giai đoạn is từ xu thế giai đoạn trước:
    ! giai đoạn 1 lấy giá trị hiện tại, giai đoạn 2-4 cộng thêm cs(is)*dt*xu thế trước
    call forecast_c(uc,vc,zc,ut_c,vt_c,zs_c,kuc,kvc,kzc,is,dt,cs)

    ! Đặt các mảng xu thế du, dv, dz về 0 trước khi tính lại cho giai đoạn này
    call init_tend_c(du_c,dv_c,dz_c)

    ! Tính xu thế của u, v tại các điểm biên (mặt biên và 4 góc) bằng công thức upwind,
    ! dùng thành phần vận tốc pháp tuyến đã nội suy về đúng vị trí so le của lưới C
    call boundary_c(ut_c,vt_c,du_c,dv_c)

    if (z_mode == 1) then
      ! Z áp đặt: ghi đè Z tạm thời bằng Z*(x,y,t) tại thời điểm của giai đoạn RK4
      call calc_z_c(zs_c,time+ct(is)*dt,radius,depth_max)
    else
      ! Z dự báo: tính số hạng nguồn S = dZ*/dt + div(Z* u) + (Z* - Z)/tau_s
      call calc_s_c(sc,zs_c,ut_c,vt_c,time+ct(is)*dt,radius,depth_max)
    end if

    ! Xu thế của u tại mặt tây của các ô trong miền: đối lưu upwind - g*dZ/dx + f*v + khuếch tán
    ! (v được nội suy từ 4 điểm lân cận về vị trí của u)
    call tend_u_c(ut_c,vt_c,zs_c,du_c)

    ! Xu thế của v tại mặt nam của các ô trong miền: đối lưu upwind - g*dZ/dy - f*u + khuếch tán
    ! (u được nội suy từ 4 điểm lân cận về vị trí của v)
    call tend_v_c(ut_c,vt_c,zs_c,dv_c)

    ! Xu thế của Z tại tâm ô (chỉ khi Z là biến dự báo): Z_t = -(Zu)_x - (Zv)_y + S,
    ! Z tại mặt ô lấy trung bình của hai ô kề nhau
    if (z_mode == 2) call tend_z_c(ut_c,vt_c,zs_c,sc,dz_c)

    kuc(:,:,is)=du_c                          ! Lưu hệ số RK4 thứ is của u
    kvc(:,:,is)=dv_c                          ! Lưu hệ số RK4 thứ is của v
    kzc(:,:,is)=dz_c                          ! Lưu hệ số RK4 thứ is của Z (bằng 0 khi Z áp đặt)
  end do

  ! Kết hợp 4 hệ số: x(n+1) = x(n) + dt*(k1 + 2*k2 + 2*k3 + k4)/6 cho u, v, Z
  call update_c(uc,vc,zc,kuc,kvc,kzc,dt)

  if (mod(it+1,nout)==0) then                 ! Cứ mỗi nout bước (= 1 giờ) thì ghi kết quả
    ! Khi Z áp đặt, zc chưa được cập nhật nên phải tính lại Z*(t) trước khi ghi
    if (z_mode == 1) call calc_z_c(zc,dble(it+1)*dt,radius,depth_max)

    ! Ghi thời gian, tọa độ, u, v (đã nội suy về tâm ô), Z 
    call output_c(unit_c,dble(it+1)*dt,uc,vc,zc)

    ! Ghi một dòng chẩn đoán (Zmin, Vmax, vrmin, chỉ số nhiễu bàn cờ) 
    call write_diag_c(unit_dc,dble(it+1)*dt,uc,vc,zc)
  end if
end do
close(unit_c)
close(unit_dc)

! ===== Chẩn đoán tại thời điểm cuối =====
! Đặt lại Z = Z*(tend) nếu Z áp đặt, để các chẩn đoán và file xuất dùng Z cuối
if (z_mode == 1) call calc_z_c(zc,tend,radius,depth_max)

call diag_c(uc,vc,zc,zmin,vmax,vrmin)   ! Tính Zmin, Vmax và vận tốc hướng tâm nhỏ nhất vrmin trên lưới C
call output_profile_c(zc) ! Xuất mặt cắt Z dọc hàng y = 0 
call output_2d_c(zc)     ! Xuất trường Z hai chiều tại thời điểm cuối 

  write(*,'(a8,3(2x,es13.5))') 'C-grid',zmin,vmax,vrmin

  write(*,'(a)') '=============================================='
  write(*,'(a)') 'Simulations completed.'

contains

  ! ==========================================================
  ! DỰ BÁO: LƯỚI A (u, v, Z)
  ! ==========================================================
 ! ====================================================================
  ! CÁC KÝ HIỆU DÙNG TRONG CHÚ THÍCH
  !   u, v : thành phần gió theo x, y            Z : độ cao mặt nước (3000 m lúc đầu)
  !   g = 9.81 m/s^2 ; f = 2*omega*sin(35 độ) : tham số Coriolis ; K = km : độ nhớt số trị
  !   dx = dy = 30 km ; dt = 60 s ; Z* : Z mục tiêu (vùng trũng hình nón)
  !   Hệ phương trình của đề (Z đã biết):
  !       u_t = -u u_x - v u_y - g Z_x + f v + d/dx(K du/dx)
  !       v_t = -u v_x - v v_y - g Z_y - f u + d/dy(K dv/dy)
  !   (khi Z là biến dự báo, thêm:  Z_t + (Zu)_x + (Zv)_y = S)
  ! ====================================================================


  ! ====================================================================
  ! FORECAST (lưới A): TÍNH TRẠNG THÁI TẠM THỜI CỦA TỪNG GIAI ĐOẠN RK4
  ! Công thức (c = cs = 0, 1/2, 1/2, 1):
  !     giai đoạn 1 :  x1 = x^n
  !     giai đoạn s :  xs = x^n + c(s)*dt*k(s-1)       (s = 2, 3, 4)
  !   tức là  u1 = u^n ,  u2 = u^n + (dt/2)*k1 ,  u3 = u^n + (dt/2)*k2 ,  u4 = u^n + dt*k3
  !   (đúng ghi chú R-K4 của đề: u^n -> u1, u2, u3, u^{n+1}; v^n -> v1, v2, v3, v^{n+1})
  ! Mục đích: RK4 cần tính xu thế (vế phải phương trình) tại 4 trạng thái trung gian;
  !   subroutine này tạo ra trạng thái đó (ut, vt, zt) cho giai đoạn is.
  ! ====================================================================
  subroutine forecast_a(u,v,z,ut,vt,zt,ku,kv,kz,is,dt,cs)
    implicit none
    integer, intent(in) :: is
    real(8), intent(in) :: dt, cs(4)
    real(8), intent(in) :: z(na,na)
    real(4), intent(in) :: u(na,na),v(na,na)
    real(8), intent(in) :: kz(na,na,4)
    real(4), intent(in) :: ku(na,na,4),kv(na,na,4)
    real(8), intent(out) :: zt(na,na)
    real(4), intent(out) :: ut(na,na),vt(na,na)

    if (is == 1) then
      ut = u
      vt = v
      zt = z
    else
      ut = u + cs(is)*dt*ku(:,:,is-1)
      vt = v + cs(is)*dt*kv(:,:,is-1)
      zt = z + cs(is)*dt*kz(:,:,is-1)
    end if
  end subroutine forecast_a

  ! ====================================================================
  ! FORECAST (lưới C): cùng công thức RK4 như lưới A, chỉ khác kích thước mảng:
  !   u có (nb+1) x nb điểm (ở mặt tây/đông của ô), v có nb x (nb+1) điểm (ở mặt nam/bắc),
  !   Z có nb x nb điểm (ở tâm ô).
  !     xs = x^n + c(s)*dt*k(s-1)
  ! [ĐÃ SỬA] Bản bạn dán thiếu 3 dòng khai báo: implicit none, is, dt và cs
  !   (do implicit none của chương trình chính áp dụng cả cho subroutine nội bộ, thiếu sẽ báo lỗi biên dịch).
  ! ====================================================================
  subroutine forecast_c(u,v,z,ut,vt,zt,ku,kv,kz,is,dt,cs)
    implicit none
    integer, intent(in) :: is
    real(8), intent(in) :: dt, cs(4)
    real(8), intent(in) :: z(nb,nb)
    real(4), intent(in) :: u(nb+1,nb),v(nb,nb+1)
    real(8), intent(in) :: kz(nb,nb,4)
    real(4), intent(in) :: ku(nb+1,nb,4),kv(nb,nb+1,4)
    real(8), intent(out) :: zt(nb,nb)
    real(4), intent(out) :: ut(nb+1,nb),vt(nb,nb+1)

    if (is == 1) then
      ut = u
      vt = v
      zt = z
    else
      ut = u + cs(is)*dt*ku(:,:,is-1)
      vt = v + cs(is)*dt*kv(:,:,is-1)
      zt = z + cs(is)*dt*kz(:,:,is-1)
    end if
  end subroutine forecast_c

  ! ==========================================================
  ! CẬP NHẬT (RK4)
  ! ==========================================================

  ! ====================================================================
  ! UPDATE (lưới A): GHÉP 4 HỆ SỐ RK4 THÀNH NGHIỆM Ở BƯỚC THỜI GIAN MỚI
  ! Công thức:
  !     x^{n+1} = x^n + (dt/6)*( k1 + 2*k2 + 2*k3 + k4 )      với x = u, v, Z
  ! Mục đích: tạo u^{n+1}, v^{n+1} (và Z^{n+1} khi Z là biến dự báo) có độ chính xác bậc 4
  !   theo thời gian, từ 4 xu thế k1..k4 đã tính ở 4 giai đoạn. Khi Z áp đặt thì kz = 0
  !   nên Z không đổi ở đây (Z được tính lại từ Z* khi ghi kết quả).
  ! ====================================================================
  subroutine update_a(u,v,z,ku,kv,kz,dt)
    implicit none
    real(8), intent(in) :: dt
    real(4), intent(in) :: ku(na,na,4),kv(na,na,4)
    real(8), intent(in) :: kz(na,na,4)
    real(4), intent(inout) :: u(na,na),v(na,na)
    real(8), intent(inout) :: z(na,na)

    u = u + dt*(ku(:,:,1) + 2.0d0*ku(:,:,2) + 2.0d0*ku(:,:,3) + ku(:,:,4))/6.0d0
    v = v + dt*(kv(:,:,1) + 2.0d0*kv(:,:,2) + 2.0d0*kv(:,:,3) + kv(:,:,4))/6.0d0
    z = z + dt*(kz(:,:,1) + 2.0d0*kz(:,:,2) + 2.0d0*kz(:,:,3) + kz(:,:,4))/6.0d0
  end subroutine update_a

  ! UPDATE (lưới C): cùng công thức  x^{n+1} = x^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
  !   áp dụng cho u (mặt tây/đông), v (mặt nam/bắc), Z (tâm ô).
  subroutine update_c(u,v,z,ku,kv,kz,dt)
    implicit none
    real(8), intent(in) :: dt
    real(4), intent(in) :: ku(nb+1,nb,4),kv(nb,nb+1,4)
    real(8), intent(in) :: kz(nb,nb,4)
    real(4), intent(inout) :: u(nb+1,nb),v(nb,nb+1)
    real(8), intent(inout) :: z(nb,nb)

    u = u + dt*(ku(:,:,1) + 2.0d0*ku(:,:,2) + 2.0d0*ku(:,:,3) + ku(:,:,4))/6.0d0
    v = v + dt*(kv(:,:,1) + 2.0d0*kv(:,:,2) + 2.0d0*kv(:,:,3) + kv(:,:,4))/6.0d0
    z = z + dt*(kz(:,:,1) + 2.0d0*kz(:,:,2) + 2.0d0*kz(:,:,3) + kz(:,:,4))/6.0d0
  end subroutine update_c

  ! ====================================================================
  ! CALC_Z (lưới A): Z MỤC TIÊU Z*(x,y,t) -- VÙNG TRŨNG HÌNH NÓN
  ! Công thức:
  !     r      = sqrt(x^2 + y^2)                    (khoảng cách tới tâm miền)
  !     bowl   = 1 - r/R   nếu r <= R ;  0 nếu r > R     (R = radius = 300 km)
  !     ramp   = t / T     (0 <= ramp <= 1, T = tend = 30 giờ, tăng TUYẾN TÍNH)
  !     Z*     = Z0 - H * ramp * bowl                (Z0 = 3000 m, H = depth_max = 150 m)
  !   => tâm vùng trũng hạ 5 m mỗi giờ, đạt 2850 m lúc t = 30 giờ; ngoài r > R vẫn là 3000 m.
  ! Mục đích trong bài:
  !   - z_mode = 1 : chính là trường Z ÁP ĐẶT trong phương trình động lượng (gradient Z gây ra gió);
  !   - z_mode = 2 : dùng làm điều kiện đầu và để tạo số hạng nguồn S của phương trình Z.
  ! ====================================================================
  subroutine calc_z_a(z,time,radius,depth_max)
    implicit none
    real(8), intent(out) :: z(na,na)
    real(8), intent(in)  :: time,radius,depth_max
    real(8) :: x,y,r,tau,ramp,bowl
    integer :: i,j

    tau  = min(1.0d0,max(0.0d0,time/tend))
    ramp = tau                                   ! [ĐÃ SỬA] ramp TUYẾN TÍNH: 150 m / 30 giờ = 5 m mỗi giờ

    do j=1,na
      y = dble(j-(na+1)/2)*dy
      do i=1,na
        x = dble(i-(na+1)/2)*dx
        r = sqrt(x*x+y*y)
        bowl = 0.0d0
        if (r <= radius) bowl = 1.0d0 - r/radius
        z(i,j) = z0-depth_max*ramp*bowl
      end do
    end do
  end subroutine calc_z_a

  ! CALC_Z (lưới C): cùng công thức  Z* = Z0 - H*(t/T)*bowl(r)  nhưng tính tại TÂM Ô của lưới C.
  subroutine calc_z_c(z,time,radius,depth_max)
    implicit none
    real(8), intent(out) :: z(nb,nb)
    real(8), intent(in)  :: time,radius,depth_max
    real(8) :: x,y,r,tau,ramp,bowl
    integer :: i,j

    tau  = min(1.0d0,max(0.0d0,time/tend))
    ramp = tau                                   ! [ĐÃ SỬA] ramp TUYẾN TÍNH: 150 m / 30 giờ = 5 m mỗi giờ

    do j=1,nb
      y = dble(j-(nb+1)/2)*dy
      do i=1,nb
        x = dble(i-(nb+1)/2)*dx
        r = sqrt(x*x+y*y)
        bowl = 0.0d0
        if (r <= radius) bowl = 1.0d0 - r/radius
        z(i,j) = z0-depth_max*ramp*bowl
      end do
    end do
  end subroutine calc_z_c

  ! ====================================================================
  ! CALC_S (lưới A): SỐ HẠNG NGUỒN S(x,y,t) CỦA PHƯƠNG TRÌNH Z  (CHỈ DÙNG KHI z_mode = 2)
  ! Công thức:
  !     S = dZ*/dt  +  div(Z* u)  +  (Z* - Z)/tau_s
  !     dZ*/dt = -H * bowl / T      (với 0 < t < T, bằng 0 ngoài khoảng này)
  !     div(Z* u) = (Z* u)_x + (Z* v)_y, tính bằng đúng toán tử sai phân của lưới (gọi tend_z_a với Z*, S=0)
  ! Mục đích trong bài: nếu u = v = 0 và Z = 3000 đều thì mọi gradient bằng 0 nên áp thấp
  !   KHÔNG tự xuất hiện (mục 1.2 của đề). S là cưỡng bức bên ngoài tạo vùng trũng:
  !   - dZ*/dt    : tốc độ hạ độ cao;
  !   - div(Z*u)  : bù phần hội tụ gió làm đầy vùng trũng, nên Z bám sát Z* (Zmin = 2850 m lúc t = 30 giờ);
  !   - (Z*-Z)/tau_s : kéo nhẹ Z về Z* để ổn định, tau_s = 3600 s (đặt 1.0d30 để bỏ số hạng này).
  ! ====================================================================
  subroutine calc_s_a(s,z,u,v,time,radius,depth_max)
    implicit none
    real(8), intent(out) :: s(na,na)
    real(8), intent(in)  :: z(na,na),time,radius,depth_max
    real(4), intent(in)  :: u(na,na),v(na,na)
    real(8) :: ztar(na,na),zer(na,na),dzf(na,na)
    real(8) :: x,y,r,tau,dramp,bowl
    integer :: i,j

    call calc_z_a(ztar,time,radius,depth_max)
    zer = 0.0d0
    call tend_z_a(u,v,ztar,zer,dzf)            ! dzf = -div(Z* u)
    tau = time/tend
    if (tau > 0.0d0 .and. tau < 1.0d0) then
      dramp = 1.0d0/tend                      ! d(ramp)/dt của ramp tuyến tính
    else
      dramp = 0.0d0
    end if

    do j=1,na
      y = dble(j-(na+1)/2)*dy
      do i=1,na
        x = dble(i-(na+1)/2)*dx
        r = sqrt(x*x+y*y)
        bowl = 0.0d0
        if (r <= radius) bowl = 1.0d0 - r/radius
        s(i,j) = -depth_max*dramp*bowl - dzf(i,j) + (ztar(i,j)-z(i,j))/tau_s
      end do
    end do
  end subroutine calc_s_a

  ! CALC_S (lưới C): cùng công thức  S = dZ*/dt + div(Z* u) + (Z*-Z)/tau_s,
  !   với div(Z*u) tính bằng toán tử thông lượng của lưới C (tend_z_c).
  subroutine calc_s_c(s,z,u,v,time,radius,depth_max)
    implicit none
    real(8), intent(out) :: s(nb,nb)
    real(8), intent(in)  :: z(nb,nb),time,radius,depth_max
    real(4), intent(in)  :: u(nb+1,nb),v(nb,nb+1)
    real(8) :: ztar(nb,nb),zer(nb,nb),dzf(nb,nb)
    real(8) :: x,y,r,tau,dramp,bowl
    integer :: i,j

    call calc_z_c(ztar,time,radius,depth_max)
    zer = 0.0d0
    call tend_z_c(u,v,ztar,zer,dzf)            ! dzf = -div(Z* u)
    tau = time/tend
    if (tau > 0.0d0 .and. tau < 1.0d0) then
      dramp = 1.0d0/tend                      ! d(ramp)/dt của ramp tuyến tính
    else
      dramp = 0.0d0
    end if

    do j=1,nb
      y = dble(j-(nb+1)/2)*dy
      do i=1,nb
        x = dble(i-(nb+1)/2)*dx
        r = sqrt(x*x+y*y)
        bowl = 0.0d0
        if (r <= radius) bowl = 1.0d0 - r/radius
        s(i,j) = -depth_max*dramp*bowl - dzf(i,j) + (ztar(i,j)-z(i,j))/tau_s
      end do
    end do
  end subroutine calc_s_c

  ! ====================================================================
  ! LƯỚI A: khởi tạo / biên / xu thế
  ! ====================================================================

  ! INIT_A: ĐIỀU KIỆN ĐẦU CỦA LƯỚI A
  !     u(x,y,0) = 0 ,  v(x,y,0) = 0 ,  Z(x,y,0) = Z*(x,y,0) = 3000 m
  ! Mục đích: bắt đầu từ trạng thái nghỉ (gió bằng 0, mặt nước phẳng). Gió chỉ xuất hiện khi
  !   vùng trũng hình thành và gradient Z kéo không khí vào.
  subroutine init_a(u,v,z,radius,depth_max)
    implicit none
    real(8), intent(in)  :: radius,depth_max
    real(8), intent(out) :: z(na,na)
    real(4), intent(out) :: u(na,na),v(na,na)
    u=0.0d0
    v=0.0d0
    call calc_z_a(z,0.0d0,radius,depth_max)    ! Z = 3000 tại t = 0
  end subroutine init_a

  ! INIT_TEND_A: ĐẶT LẠI CÁC XU THẾ  du = dv = dz = 0  trước mỗi giai đoạn RK4.
  ! Mục đích: các subroutine tính xu thế chỉ ghi vào điểm trong miền và điểm biên của chúng;
  !   khởi tạo bằng 0 để không còn sót giá trị của giai đoạn trước.
  subroutine init_tend_a(du,dv,dz)
    implicit none
    real(8), intent(out) :: dz(na,na)
    real(4), intent(out) :: du(na,na),dv(na,na)
    du=0.0d0
    dv=0.0d0
    dz=0.0d0
  end subroutine init_tend_a

  ! ====================================================================
  ! BOUNDARY_A: XU THẾ CỦA u, v TẠI CÁC ĐIỂM BIÊN (4 cạnh và 4 góc) CỦA LƯỚI A
  ! Với đại lượng q (= u hoặc v) và vận tốc pháp tuyến hướng RA NGOÀI un:
  !     dq/dt = - un^+ * (q_b - q_in)/dn  +  un^- * (q_b - q_ext)/dn
  !     un^+ = max(un,0) , un^- = min(un,0) , q_ext = 0 (giá trị giả định ngoài miền)
  !   - un > 0 (gió thổi ra ngoài): đối lưu upwind từ trong miền ra biên, q_b - q_in là chênh lệch với điểm kế bên trong;
  !   - un < 0 (gió thổi vào): giá trị ngoài miền bằng 0 được đưa vào, q_b bị kéo về q_ext.
  ! Mục đích: ở biên không có điểm lưới phía ngoài nên không dùng được sai phân của điểm trong;
  !   công thức này là điều kiện biên miền hở (theo đề) giữ cho biên không phản xạ/khuếch đại gió.
  ! Ở góc, cộng hai đóng góp (theo x và theo y).
  ! ====================================================================
  subroutine boundary_a(u,v,du,dv)
    implicit none
    real(4), intent(in)    :: u(na,na),v(na,na)
    real(4), intent(inout) :: du(na,na),dv(na,na)
    integer :: i,j

    do j=2,na-1
      du(1,j)=boundary_tendency(-u(1,j),u(1,j),u(2,j),dx)
      du(na,j)=boundary_tendency(u(na,j),u(na,j),u(na-1,j),dx)
      dv(1,j)=boundary_tendency(-u(1,j),v(1,j),v(2,j),dx)
      dv(na,j)=boundary_tendency(u(na,j),v(na,j),v(na-1,j),dx)
    end do

    do i=2,na-1
      du(i,1)=boundary_tendency(-v(i,1),u(i,1),u(i,2),dy)
      du(i,na)=boundary_tendency(v(i,na),u(i,na),u(i,na-1),dy)
      dv(i,1)=boundary_tendency(-v(i,1),v(i,1),v(i,2),dy)
      dv(i,na)=boundary_tendency(v(i,na),v(i,na),v(i,na-1),dy)
    end do

    du(1,1)=boundary_tendency(-u(1,1),u(1,1),u(2,1),dx)+ &
            boundary_tendency(-v(1,1),u(1,1),u(1,2),dy)
    du(na,1)=boundary_tendency(u(na,1),u(na,1),u(na-1,1),dx)+ &
             boundary_tendency(-v(na,1),u(na,1),u(na,2),dy)
    du(1,na)=boundary_tendency(-u(1,na),u(1,na),u(2,na),dx)+ &
             boundary_tendency(v(1,na),u(1,na),u(1,na-1),dy)
    du(na,na)=boundary_tendency(u(na,na),u(na,na),u(na-1,na),dx)+ &
              boundary_tendency(v(na,na),u(na,na),u(na,na-1),dy)

    dv(1,1)=boundary_tendency(-u(1,1),v(1,1),v(2,1),dx)+ &
            boundary_tendency(-v(1,1),v(1,1),v(1,2),dy)
    dv(na,1)=boundary_tendency(u(na,1),v(na,1),v(na-1,1),dx)+ &
             boundary_tendency(-v(na,1),v(na,1),v(na,2),dy)
    dv(1,na)=boundary_tendency(-u(1,na),v(1,na),v(2,na),dx)+ &
             boundary_tendency(v(1,na),v(1,na),v(1,na-1),dy)
    dv(na,na)=boundary_tendency(u(na,na),v(na,na),v(na-1,na),dx)+ &
              boundary_tendency(v(na,na),v(na,na),v(na,na-1),dy)
  end subroutine boundary_a

  ! ====================================================================
  ! TEND_U_A: XU THẾ CỦA u TẠI CÁC ĐIỂM TRONG MIỀN (LƯỚI A)
  ! Phương trình (1) của đề:
  !     u_t = -u u_x - v u_y  -  g Z_x  +  f v  +  d/dx( K du/dx )
  !   - đối lưu  -u u_x - v u_y : advect_a (pha trộn trung tâm/upwind, hệ số alpha_a);
  !   - lực gradient áp suất  -g Z_x : Z_x ~ (Z_{i+1,j} - Z_{i-1,j}) / (2 dx)   (sai phân trung tâm);
  !   - lực Coriolis  +f v ;
  !   - độ nhớt số trị (K không đổi nên = K u_xx) : u_xx ~ (u_{i+1,j} - 2 u_{i,j} + u_{i-1,j}) / dx^2 , K = km.
  ! Mục đích: Z thấp ở tâm tạo gradient Z, gradient này gia tốc gió vào tâm; Coriolis bẻ hướng
  !   thành xoáy thuận; đối lưu và độ nhớt là các hiệu ứng thứ cấp.
  ! ====================================================================
  subroutine tend_u_a(u,v,z,du)
    implicit none
    real(8), intent(in)    :: z(na,na)
    real(4), intent(in)    :: u(na,na),v(na,na)
    real(4), intent(inout) :: du(na,na)
    real(8) :: dzdx
    real(4) :: lapu,advu
    integer :: i,j

    do j=2,na-1
      do i=2,na-1
        advu=advect_a(u(i,j),u(i-1,j),u(i+1,j), &
          u(i,j-1),u(i,j+1),u(i,j),v(i,j))
        dzdx=(z(i+1,j)-z(i-1,j))/(2.0d0*dx)
        lapu=(u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx**2
        du(i,j)=advu-g*dzdx+f*v(i,j)+km*lapu
      end do
    end do
  end subroutine tend_u_a

  ! ====================================================================
  ! TEND_V_A: XU THẾ CỦA v TẠI CÁC ĐIỂM TRONG MIỀN (LƯỚI A)
  ! Phương trình (2) của đề:
  !     v_t = -u v_x - v v_y  -  g Z_y  -  f u  +  d/dy( K dv/dy )
  !   - Z_y ~ (Z_{i,j+1} - Z_{i,j-1}) / (2 dy) ;   v_yy ~ (v_{i,j+1} - 2 v_{i,j} + v_{i,j-1}) / dy^2
  !   - đối lưu qua advect_a (cùng pha trộn như với u).
  ! Mục đích: giống tend_u_a nhưng cho thành phần gió theo y; dấu -f u (đối nghịch với +f v)
  !   làm gió quay theo chiều xoáy thuận ở bán cầu bắc.
  ! ====================================================================
  subroutine tend_v_a(u,v,z,dv)
    implicit none
    real(8), intent(in)    :: z(na,na)
    real(4), intent(in)    :: u(na,na),v(na,na)
    real(4), intent(inout) :: dv(na,na)
    real(8) :: dzdy
    real(4) :: lapv,advv
    integer :: i,j

    do j=2,na-1
      do i=2,na-1
        advv=advect_a(v(i,j),v(i-1,j),v(i+1,j), &
          v(i,j-1),v(i,j+1),u(i,j),v(i,j))
        dzdy=(z(i,j+1)-z(i,j-1))/(2.0d0*dy)
        lapv=(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy**2
        dv(i,j)=advv-g*dzdy-f*u(i,j)+km*lapv
      end do
    end do
  end subroutine tend_v_a

  ! ====================================================================
  ! TEND_Z_A: XU THẾ CỦA Z (LƯỚI A)  -- CHỈ DÙNG KHI Z LÀ BIẾN DỰ BÁO (z_mode = 2)
  ! Phương trình bảo toàn khối lượng (6) của đề:
  !     Z_t = -(Z u)_x - (Z v)_y + S
  !     (Z u)_x ~ ( Z_{i+1,j} u_{i+1,j} - Z_{i-1,j} u_{i-1,j} ) / (2 dx)
  !     (Z v)_y ~ ( Z_{i,j+1} v_{i,j+1} - Z_{i,j-1} v_{i,j-1} ) / (2 dy)
  ! Mục đích: để Z phản ứng với dòng hội tụ/phân kỳ của gió (hội tụ làm Z tăng lên, phân kỳ làm Z giảm);
  !   S là nguồn cưỡng bức tạo vùng trũng. Ở điểm biên chỉ có S (xa tâm S ~ 0).
  ! Lưu ý: sai phân 2dx nối các điểm cách nhau 2 ô, nên mode bàn cờ không "thấy" nhau (nhược điểm lưới A).
  ! ====================================================================
  subroutine tend_z_a(u,v,z,s,dz)
    implicit none
    real(8), intent(in)  :: z(na,na),s(na,na)
    real(4), intent(in)  :: u(na,na),v(na,na)
    real(8), intent(out) :: dz(na,na)
    integer :: i,j

    dz = s                                   ! Biên: chỉ có nguồn S (S ~ 0 ở xa tâm)
    do j=2,na-1
      do i=2,na-1
        dz(i,j) = s(i,j) &
          - (z(i+1,j)*u(i+1,j)-z(i-1,j)*u(i-1,j))/(2.0d0*dx) &
          - (z(i,j+1)*v(i,j+1)-z(i,j-1)*v(i,j-1))/(2.0d0*dy)
      end do
    end do
  end subroutine tend_z_a

  ! ====================================================================
  ! LƯỚI C: khởi tạo / biên / xu thế
  !   Vị trí các biến: Z ở TÂM ô (i,j); u(i,j) ở MẶT TÂY của ô; v(i,j) ở MẶT NAM của ô.
  ! ====================================================================

  ! INIT_C: ĐIỀU KIỆN ĐẦU LƯỚI C:  u = 0 , v = 0 , Z = Z*(t=0) = 3000 m  (trạng thái nghỉ, như lưới A)
  subroutine init_c(u,v,z,radius,depth_max)
    implicit none
    real(8), intent(in)  :: radius,depth_max
    real(4), intent(out) :: u(nb+1,nb),v(nb,nb+1)
    real(8), intent(out) :: z(nb,nb)
    u=0.0d0
    v=0.0d0
    call calc_z_c(z,0.0d0,radius,depth_max)    ! Z = 3000 tại t = 0
  end subroutine init_c

  ! INIT_TEND_C: đặt lại  du = dv = dz = 0  trước mỗi giai đoạn RK4 (mục đích như init_tend_a).
  subroutine init_tend_c(du,dv,dz)
    implicit none
    real(8), intent(out) :: dz(nb,nb)
    real(4), intent(out) :: du(nb+1,nb),dv(nb,nb+1)
    du=0.0d0
    dv=0.0d0
    dz=0.0d0
  end subroutine init_tend_c

  ! ====================================================================
  ! BOUNDARY_C: XU THẾ CỦA u, v TẠI CÁC ĐIỂM BIÊN CỦA LƯỚI C
  ! Cùng công thức biên miền hở như lưới A (xem boundary_a):
  !     dq/dt = - un^+ (q_b - q_in)/dn + un^- (q_b - q_ext)/dn ,   q_ext = 0
  ! Khác biệt do lưới so le: vận tốc pháp tuyến phải lấy ĐÚNG VỊ TRÍ của biến đang xét.
  !   - u ở mặt tây/đông: pháp tuyến là chính u (nằm ngay trên biên);
  !   - u ở hàng nam/bắc: pháp tuyến là v, nội suy  vint_u = (v_{i-1,j}+v_{i,j}+v_{i-1,j+1}+v_{i,j+1})/4;
  !   - v ở mặt nam/bắc: pháp tuyến là chính v;
  !   - v ở cột tây/đông: pháp tuyến là u, nội suy  uint_v = (u_{i,j-1}+u_{i+1,j-1}+u_{i,j}+u_{i+1,j})/4.
  ! Mục đích: biên miền hở đúng với vị trí so le của lưới C (đề yêu cầu biên phù hợp lưới Arakawa C).
  ! ====================================================================
  subroutine boundary_c(u,v,du,dv)
    implicit none
    real(4), intent(in)    :: u(nb+1,nb),v(nb,nb+1)
    real(4), intent(inout) :: du(nb+1,nb),dv(nb,nb+1)
    real(4) :: vint_u,uint_v
    integer :: i,j

    ! ---- Điều kiện biên cho u: phía tây / phía đông
    do j=1,nb
      du(1,j)=boundary_tendency(-u(1,j),u(1,j),u(2,j),dx)
      du(nb+1,j)=boundary_tendency(u(nb+1,j),u(nb+1,j),u(nb,j),dx)
    end do

    ! ---- Điều kiện biên cho u: phía nam / phía bắc
    do i=2,nb
      vint_u=0.25d0*(v(i-1,1)+v(i,1)+v(i-1,2)+v(i,2))
      du(i,1)=boundary_tendency(-vint_u,u(i,1),u(i,2),dy)
      vint_u=0.25d0*(v(i-1,nb)+v(i,nb)+v(i-1,nb+1)+v(i,nb+1))
      du(i,nb)=boundary_tendency(vint_u,u(i,nb),u(i,nb-1),dy)
    end do

    ! ---- Bốn góc của u (đặt sau các cạnh)
    du(1,1)=boundary_tendency(-u(1,1),u(1,1),u(2,1),dx)+ &
            boundary_tendency(-v(1,1),u(1,1),u(1,2),dy)
    du(nb+1,1)=boundary_tendency(u(nb+1,1),u(nb+1,1),u(nb,1),dx)+ &
               boundary_tendency(-v(nb,1),u(nb+1,1),u(nb+1,2),dy)
    du(1,nb)=boundary_tendency(-u(1,nb),u(1,nb),u(2,nb),dx)+ &
             boundary_tendency(v(1,nb+1),u(1,nb),u(1,nb-1),dy)
    du(nb+1,nb)=boundary_tendency(u(nb+1,nb),u(nb+1,nb),u(nb,nb),dx)+ &
                boundary_tendency(v(nb,nb+1),u(nb+1,nb),u(nb+1,nb-1),dy)

    ! ---- Điều kiện biên cho v: phía nam / phía bắc
    do i=1,nb
      dv(i,1)=boundary_tendency(-v(i,1),v(i,1),v(i,2),dy)
      dv(i,nb+1)=boundary_tendency(v(i,nb+1),v(i,nb+1),v(i,nb),dy)
    end do

    ! ---- Điều kiện biên cho v: phía tây / phía đông
    do j=2,nb
      uint_v=0.25d0*(u(1,j-1)+u(2,j-1)+u(1,j)+u(2,j))
      dv(1,j)=boundary_tendency(-uint_v,v(1,j),v(2,j),dx)
      uint_v=0.25d0*(u(nb,j-1)+u(nb+1,j-1)+u(nb,j)+u(nb+1,j))
      dv(nb,j)=boundary_tendency(uint_v,v(nb,j),v(nb-1,j),dx)
    end do

    ! ---- Bốn góc của v (đặt sau các cạnh)
    dv(1,1)=boundary_tendency(-u(1,1),v(1,1),v(2,1),dx)+ &
            boundary_tendency(-v(1,1),v(1,1),v(1,2),dy)
    dv(nb,1)=boundary_tendency(u(nb+1,1),v(nb,1),v(nb-1,1),dx)+ &
             boundary_tendency(-v(nb,1),v(nb,1),v(nb,2),dy)
    dv(1,nb+1)=boundary_tendency(-u(1,nb),v(1,nb+1),v(2,nb+1),dx)+ &
               boundary_tendency(v(1,nb+1),v(1,nb+1),v(1,nb),dy)
    dv(nb,nb+1)=boundary_tendency(u(nb+1,nb),v(nb,nb+1),v(nb-1,nb+1),dx)+ &
                boundary_tendency(v(nb,nb+1),v(nb,nb+1),v(nb,nb),dy)
  end subroutine boundary_c

  ! ====================================================================
  ! TEND_U_C: XU THẾ CỦA u TẠI MẶT TÂY CỦA CÁC Ô TRONG MIỀN (LƯỚI C)
  !     u_t = -u u_x - v u_y - g Z_x + f v + K u_xx           (phương trình (1) của đề)
  !   - v tại vị trí của u (nội suy 4 điểm):  vint_u = (v_{i-1,j}+v_{i,j}+v_{i-1,j+1}+v_{i,j+1})/4
  !   - đối lưu UPWIND bậc 1 (a+ = max(a,0), a- = min(a,0), a = u, b = vint_u):
  !       -a+ (u_i - u_{i-1})/dx - a- (u_{i+1} - u_i)/dx - b+ (u_j - u_{j-1})/dy - b- (u_{j+1} - u_j)/dy
  !   - gradient áp suất:  Z_x ~ (Z_{i,j} - Z_{i-1,j}) / dx     (hai ô kề, đúng vị trí mặt tây — sai phân GỌN, độ rộng dx)
  !   - độ nhớt:  K (u_{i+1} - 2 u_i + u_{i-1}) / dx^2 ,  K = km
  ! Mục đích: cùng vật lý như tend_u_a nhưng trên lưới so le: upwind tiêu tán nên ổn định,
  !   và gradient Z gọn nối trực tiếp hai ô kề nên không có mode bàn cờ bị "vô hình".
  ! ====================================================================
  subroutine tend_u_c(u,v,z,du)
    implicit none
    real(8), intent(in)    :: z(nb,nb)
    real(4), intent(in)    :: u(nb+1,nb),v(nb,nb+1)
    real(4), intent(inout) :: du(nb+1,nb)
    real(4) :: vint_u,lapu
    real(8) :: dzdx_u
    real(4) :: ap,am,bp,bm,adv
    integer :: i,j

    do j=2,nb-1
      do i=2,nb
        vint_u=0.25d0*(v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1))
        ap=max(u(i,j),0.0)
        am=min(u(i,j),0.0)
        bp=max(vint_u,0.0)
        bm=min(vint_u,0.0)
        adv=-ap*(u(i,j)-u(i-1,j))/dx &
            -am*(u(i+1,j)-u(i,j))/dx &
            -bp*(u(i,j)-u(i,j-1))/dy &
            -bm*(u(i,j+1)-u(i,j))/dy
        dzdx_u=(z(i,j)-z(i-1,j))/dx
        lapu=(u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx**2
        du(i,j)=adv-g*dzdx_u+f*vint_u+km*lapu
      end do
    end do
  end subroutine tend_u_c

  ! ====================================================================
  ! TEND_V_C: XU THẾ CỦA v TẠI MẶT NAM CỦA CÁC Ô TRONG MIỀN (LƯỚI C)
  !     v_t = -u v_x - v v_y - g Z_y - f u + K v_yy           (phương trình (2) của đề)
  !   - u tại vị trí của v (nội suy 4 điểm):  uint_v = (u_{i,j-1}+u_{i+1,j-1}+u_{i,j}+u_{i+1,j})/4
  !   - đối lưu upwind bậc 1 với a = uint_v, b = v (cùng dạng như tend_u_c)
  !   - gradient áp suất:  Z_y ~ (Z_{i,j} - Z_{i,j-1}) / dy ;   độ nhớt:  K (v_{j+1} - 2 v_j + v_{j-1}) / dy^2
  ! Mục đích: giống tend_u_c nhưng cho thành phần gió theo y.
  ! ====================================================================
  subroutine tend_v_c(u,v,z,dv)
    implicit none
    real(8), intent(in)    :: z(nb,nb)
    real(4), intent(in)    :: u(nb+1,nb),v(nb,nb+1)
    real(4), intent(inout) :: dv(nb,nb+1)
    real(4) :: uint_v,lapv
    real(8) :: dzdy_v
    real(4) :: ap,am,bp,bm,adv
    integer :: i,j

    do j=2,nb
      do i=2,nb-1
        uint_v=0.25d0*(u(i,j-1)+u(i+1,j-1)+u(i,j)+u(i+1,j))
        ap=max(uint_v,0.0)
        am=min(uint_v,0.0)
        bp=max(v(i,j),0.0)
        bm=min(v(i,j),0.0)
        adv=-ap*(v(i,j)-v(i-1,j))/dx &
            -am*(v(i+1,j)-v(i,j))/dx &
            -bp*(v(i,j)-v(i,j-1))/dy &
            -bm*(v(i,j+1)-v(i,j))/dy
        dzdy_v=(z(i,j)-z(i,j-1))/dy
        lapv=(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy**2
        dv(i,j)=adv-g*dzdy_v-f*uint_v+km*lapv
      end do
    end do
  end subroutine tend_v_c

  ! ====================================================================
  ! TEND_Z_C: XU THẾ CỦA Z TẠI TÂM Ô (LƯỚI C) -- CHỈ KHI Z LÀ BIẾN DỰ BÁO (z_mode = 2)
  ! Dạng thông lượng, Z tại mặt ô lấy trung bình hai ô kề:
  !     Z_t = -[ Z_e u_{i+1,j} - Z_w u_{i,j} ]/dx - [ Z_n v_{i,j+1} - Z_s v_{i,j} ]/dy + S
  !     Z_w = (Z_{i-1,j}+Z_{i,j})/2 ,  Z_e = (Z_{i,j}+Z_{i+1,j})/2 ,  tương tự Z_s, Z_n
  !   (ở mép miền không có ô kề: lấy chính Z của ô đó)
  ! Mục đích: phương trình liên tục (6) của đề. Lượng Z đi vào/ra ô chỉ qua 4 mặt, hội tụ tại tâm
  !   thì Z tăng, nên bảo toàn khối lượng; S là nguồn cưỡng bức tạo vùng trũng.
  ! ====================================================================
  subroutine tend_z_c(u,v,z,s,dz)
    implicit none
    real(8), intent(in)  :: z(nb,nb),s(nb,nb)
    real(4), intent(in)  :: u(nb+1,nb),v(nb,nb+1)
    real(8), intent(out) :: dz(nb,nb)
    real(8) :: zw,ze,zsou,znor
    integer :: i,j

    do j=1,nb
      do i=1,nb
        if (i == 1) then
          zw = z(1,j)
        else
          zw = 0.5d0*(z(i-1,j)+z(i,j))
        end if
        if (i == nb) then
          ze = z(nb,j)
        else
          ze = 0.5d0*(z(i,j)+z(i+1,j))
        end if
        if (j == 1) then
          zsou = z(i,1)
        else
          zsou = 0.5d0*(z(i,j-1)+z(i,j))
        end if
        if (j == nb) then
          znor = z(i,nb)
        else
          znor = 0.5d0*(z(i,j)+z(i,j+1))
        end if
        dz(i,j) = s(i,j) &
          - (ze*u(i+1,j)-zw*u(i,j))/dx &
          - (znor*v(i,j+1)-zsou*v(i,j))/dy
      end do
    end do
  end subroutine tend_z_c

  ! ====================================================================
  ! Các hàm dùng chung
  ! ====================================================================

  ! ADVECT_A: SỐ HẠNG ĐỐI LƯU  -(u q_x + v q_y)  CHO LƯỚI A  (q = u hoặc v; a, b là gió u, v)
  ! Công thức (alpha = alpha_a):
  !     đối lưu = (1-alpha) * TRUNG TÂM + alpha * UPWIND
  !     trung tâm :  -a (q_e - q_w)/(2 dx)  -  b (q_n - q_s)/(2 dy)
  !     upwind    :  -a+ (q_c - q_w)/dx - a- (q_e - q_c)/dx - b+ (q_c - q_s)/dy - b- (q_n - q_c)/dy
  !                  với a+ = max(a,0), a- = min(a,0), b+ = max(b,0), b- = min(b,0)
  ! Mục đích: mô phỏng việc gió tự vận chuyển chính nó (số hạng phi tuyến của đề).
  !   Sai phân trung tâm không tiêu tán nên dễ sinh nhiễu 2dx; upwind tiêu tán nên ổn định.
  !   alpha_a điều chỉnh mức tiêu tán số trị của lưới A (alpha_a = 0: không tiêu tán;
  !   alpha_a >= ~0.3: ổn định). Lưới C luôn dùng upwind đầy đủ (alpha = 1).
  real(4) function advect_a(qc,qw,qe,qs,qn,a,b)
    implicit none
    real(4), intent(in) :: qc,qw,qe,qs,qn,a,b
    real(4) :: ap,am,bp,bm,cen,up
    ap=max(a,0.0)
    am=min(a,0.0)
    bp=max(b,0.0)
    bm=min(b,0.0)
    cen=-a*(qe-qw)/(2.0d0*dx)-b*(qn-qs)/(2.0d0*dy)
    up =-ap*(qc-qw)/dx-am*(qe-qc)/dx &
        -bp*(qc-qs)/dy-bm*(qn-qc)/dy
    advect_a=(1.0d0-alpha_a)*cen+alpha_a*up
  end function advect_a

  ! BOUNDARY_TENDENCY: XU THẾ MỘT CHIỀU CỦA ĐẠI LƯỢNG q TẠI ĐIỂM BIÊN (dùng cho boundary_a và boundary_c)
  ! Công thức:
  !     dq/dt = - un^+ (q_b - q_in)/dn  +  un^- (q_b - q_ext)/dn
  !     un^+ = max(un,0) , un^- = min(un,0) , q_ext = 0
  !   un : vận tốc pháp tuyến hướng ra ngoài ; q_b : q tại biên ; q_in : q tại điểm kế bên trong ; dn : dx hoặc dy
  ! Mục đích: điều kiện biên miền hở theo đề:
  !   gió thoát ra thì được tải ra ngoài (upwind), gió đi vào lấy giá trị ngoài miền bằng 0.
  real(4) function boundary_tendency(un,qb,qin,dn)
    implicit none
    real(4), intent(in) :: un,qb,qin
    real(8), intent(in) :: dn
    real(4) :: un_plus,un_minus
    real(4), parameter :: qext=0.0
    un_plus=max(un,0.0)
    un_minus=min(un,0.0)
    boundary_tendency=-un_plus*(qb-qin)/dn+un_minus*(qb-qext)/dn
  end function boundary_tendency

  ! [MỚI] CBMODE: ĐO BIÊN ĐỘ "MODE BÀN CỜ" 2dx CỦA TRƯỜNG Z (real(8))
  ! Công thức:
  !     C = | (1/N) * SUM_{i,j} (-1)^(i+j) * (q_{i,j} - q0) |      (N = n1*n2 điểm lưới)
  ! Mục đích: trường trơn cho C ~ 0, còn nhiễu xen kẽ (+,-,+,-) giữa các điểm kề nhau cho C lớn.
  !   Dùng để kiểm tra nhiễu 2dx có tăng dần theo thời gian không (dấu hiệu lưới sắp mất ổn định).
  real(8) function cbmode(q,n1,n2,q0)
    implicit none
    integer, intent(in) :: n1,n2
    real(8), intent(in) :: q(n1,n2),q0
    integer :: i,j
    real(8) :: s
    s=0.0d0
    do j=1,n2
      do i=1,n1
        s=s+(-1.0d0)**(i+j)*(q(i,j)-q0)
      end do
    end do
    cbmode=abs(s)/dble(n1*n2)
  end function cbmode

  ! CBMODE_U: giống cbmode nhưng cho mảng GIÓ (real(4)); q0 = 0 khi đo nhiễu bàn cờ của u.
  !     C = | (1/N) * SUM (-1)^(i+j) * (q_{i,j} - q0) |
  real(8) function cbmode_u(q,n1,n2,q0)
    implicit none
    integer, intent(in) :: n1,n2
    real(4), intent(in) :: q(n1,n2),q0
    integer :: i,j
    real(8) :: s
    s=0.0d0
    do j=1,n2
      do i=1,n1
        s=s+(-1.0d0)**(i+j)*dble(q(i,j)-q0)
      end do
    end do
    cbmode_u=abs(s)/dble(n1*n2)
  end function cbmode_u

  ! ====================================================================
  ! Ghi và xuất số liệu
  ! ====================================================================
  subroutine output_a(unit_out,time,u,v,z)
    implicit none
    integer, intent(in) :: unit_out
    real(8), intent(in) :: time,z(na,na)
    real(4), intent(in) :: u(na,na),v(na,na)
    integer :: i,j
    real(8) :: x,y

    do j=1,na
      y=dble(j-(na+1)/2)*dy/1000.0d0
      do i=1,na
        x=dble(i-(na+1)/2)*dx/1000.0d0
        write(unit_out,'(f7.2,1x,f9.2,1x,f9.2,1x,es14.6,1x,es14.6,1x,f10.3)') &
          time/3600.0d0,x,y,u(i,j),v(i,j),z(i,j)
      end do
    end do
    write(unit_out,*)
  end subroutine output_a

  subroutine output_c(unit_out,time,u,v,z)
    implicit none
    integer, intent(in) :: unit_out
    real(8), intent(in) :: time,z(nb,nb)
    real(4), intent(in) :: u(nb+1,nb),v(nb,nb+1)
    integer :: i,j
    real(8) :: x,y
    real(4) :: uc,vc

    do j=1,nb
      y=dble(j-(nb+1)/2)*dy/1000.0d0
      do i=1,nb
        x=dble(i-(nb+1)/2)*dx/1000.0d0
        uc=0.5d0*(u(i,j)+u(i+1,j))
        vc=0.5d0*(v(i,j)+v(i,j+1))
        write(unit_out,'(f7.2,1x,f9.2,1x,f9.2,1x,es14.6,1x,es14.6,1x,f10.3)') &
          time/3600.0d0,x,y,uc,vc,z(i,j)
      end do
    end do
    write(unit_out,*)
  end subroutine output_c

  subroutine write_diag_a(unit_out,time,u,v,z)
    implicit none
    integer, intent(in) :: unit_out
    real(8), intent(in) :: time,z(na,na)
    real(4), intent(in) :: u(na,na),v(na,na)
    real(8) :: zm
    real(4) :: vm,vr
    call diag_a(u,v,z,zm,vm,vr)
    write(unit_out,'(f7.2,3(1x,es14.6),2(1x,es14.6))') time/3600.0d0,zm,vm,vr, &
      cbmode(z,na,na,z0),cbmode_u(u,na,na,0.0)
  end subroutine write_diag_a

  subroutine write_diag_c(unit_out,time,u,v,z)
    implicit none
    integer, intent(in) :: unit_out
    real(8), intent(in) :: time,z(nb,nb)
    real(4), intent(in) :: u(nb+1,nb),v(nb,nb+1)
    real(8) :: zm
    real(4) :: vm,vr
    call diag_c(u,v,z,zm,vm,vr)
    write(unit_out,'(f7.2,3(1x,es14.6),2(1x,es14.6))') time/3600.0d0,zm,vm,vr, &
      cbmode(z,nb,nb,z0),cbmode_u(u(1:nb,1:nb),nb,nb,0.0)
  end subroutine write_diag_c

  subroutine output_profile_a(z)
    implicit none
    real(8), intent(in) :: z(na,na)
    integer :: i,j0,un
    real(8) :: x

    open(newunit=un,file='Z_profile_Agrid.txt',status='replace',action='write')
    write(un,'(a)') 'x_km z_m'
    j0=(na+1)/2
    do i=1,na
      x=dble(i-(na+1)/2)*dx/1000.0d0
      write(un,'(f12.3,1x,f12.5)') x,z(i,j0)
    end do
    close(un)
  end subroutine output_profile_a

  subroutine output_profile_c(z)
    implicit none
    real(8), intent(in) :: z(nb,nb)
    integer :: i,j0,un
    real(8) :: x

    open(newunit=un,file='Z_profile_Cgrid.txt',status='replace',action='write')
    write(un,'(a)') 'x_km z_m'
    j0=(nb+1)/2
    do i=1,nb
      x=dble(i-(nb+1)/2)*dx/1000.0d0
      write(un,'(f12.3,1x,f12.5)') x,z(i,j0)
    end do
    close(un)
  end subroutine output_profile_c

  subroutine output_2d_a(z)
    implicit none
    real(8), intent(in) :: z(na,na)
    integer :: i,j,un
    real(8) :: x,y

    open(newunit=un,file='Z_final_Agrid.txt',status='replace',action='write')
    write(un,'(a)') 'x_km y_km z_m'
    do j=1,na
      y=dble(j-(na+1)/2)*dy/1000.0d0
      do i=1,na
        x=dble(i-(na+1)/2)*dx/1000.0d0
        write(un,'(f12.3,1x,f12.3,1x,f12.5)') x,y,z(i,j)
      end do
    end do
    close(un)
  end subroutine output_2d_a

  subroutine output_2d_c(z)
    implicit none
    real(8), intent(in) :: z(nb,nb)
    integer :: i,j,un
    real(8) :: x,y

    open(newunit=un,file='Z_final_Cgrid.txt',status='replace',action='write')
    write(un,'(a)') 'x_km y_km z_m'
    do j=1,nb
      y=dble(j-(nb+1)/2)*dy/1000.0d0
      do i=1,nb
        x=dble(i-(nb+1)/2)*dx/1000.0d0
        write(un,'(f12.3,1x,f12.3,1x,f12.5)') x,y,z(i,j)
      end do
    end do
    close(un)
  end subroutine output_2d_c

  subroutine diag_a(u,v,z,zmin,vmax,vrmin)
    implicit none
    real(8), intent(in)  :: z(na,na)
    real(4), intent(in)  :: u(na,na),v(na,na)
    real(8), intent(out) :: zmin
    real(4), intent(out) :: vmax,vrmin
    integer :: i,j
    real(8) :: x,y,r
    real(4) :: speed,vr,snan
    logical :: hasnan

    zmin=minval(z)
    vmax=0.0
    vrmin=huge(1.0)
    hasnan=.false.
    do j=1,na
      y=dble(j-(na+1)/2)*dy
      do i=1,na
        x=dble(i-(na+1)/2)*dx
        speed=sqrt(u(i,j)**2+v(i,j)**2)
        if (speed /= speed) then            ! Giá trị NaN
          hasnan=.true.
          snan=speed
        end if
        vmax=max(vmax,speed)
        r=sqrt(x*x+y*y)
        if (r>0.0d0) then
          vr=(x*u(i,j)+y*v(i,j))/r
          vrmin=min(vrmin,vr)
        end if
      end do
    end do
    if (hasnan) vmax=snan                   ! Báo NaN nếu có tại bất kỳ điểm nào
  end subroutine diag_a

  subroutine diag_c(u,v,z,zmin,vmax,vrmin)
    implicit none
    real(8), intent(in)  :: z(nb,nb)
    real(4), intent(in)  :: u(nb+1,nb),v(nb,nb+1)
    real(8), intent(out) :: zmin
    real(4), intent(out) :: vmax,vrmin
    integer :: i,j
    real(8) :: x,y,r
    real(4) :: uc,vc,speed,vr,snan
    logical :: hasnan

    zmin=minval(z)
    vmax=0.0
    vrmin=huge(1.0)
    hasnan=.false.
    do j=1,nb
      y=dble(j-(nb+1)/2)*dy
      do i=1,nb
        x=dble(i-(nb+1)/2)*dx
        uc=0.5d0*(u(i,j)+u(i+1,j))
        vc=0.5d0*(v(i,j)+v(i,j+1))
        speed=sqrt(uc**2+vc**2)
        if (speed /= speed) then            ! Giá trị NaN
          hasnan=.true.
          snan=speed
        end if
        vmax=max(vmax,speed)
        r=sqrt(x*x+y*y)
        if (r>0.0d0) then
          vr=(x*uc+y*vc)/r
          vrmin=min(vrmin,vr)
        end if
      end do
    end do
    if (hasnan) vmax=snan
  end subroutine diag_c

end program shallow_water
