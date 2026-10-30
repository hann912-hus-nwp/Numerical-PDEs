import os
folder = r"D:\fortran"
print(os.listdir(folder))


#Vẽ mặt cắt Z của lưới A và lưới C
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. Khai báo file dữ liệu
# ============================================================

folder = r"D:\fortran"

files = {
    "Lưới A (A-grid)": folder + r"\unit3_Agrid.txt",
    "Lưới C (C-grid)": folder + r"\unit3_Cgrid.txt",
}

# ============================================================
# 2. Thông số mô hình
# ============================================================

times = [6, 18, 30]
colors = ["blue", "green", "red"]

T_END = 30.0
Z0 = 3000.0
DEPTH = 150.0
R = 300.0
XLIM = 600

# ============================================================
# 3. Hàm tính độ cao mục tiêu
# ============================================================

def z_target(x_km, t_h):
    tau = np.clip(t_h / T_END, 0, 1)

    bowl = np.where(
        np.abs(x_km) <= R,
        1 - np.abs(x_km) / R,
        0.0
    )

    return Z0 - DEPTH * tau * bowl

# ============================================================
# 4. Hàm đọc dữ liệu
# ============================================================

def load(path):
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=0,
        usecols=range(6),
        names=["time_h", "x_km", "y_km", "u", "v", "z"]
    )

    df = df.apply(pd.to_numeric, errors="coerce")

    return df.dropna(
        subset=["time_h", "x_km", "y_km", "z"]
    )

# ============================================================
# 5. Lấy mặt cắt Z dọc theo y = 0
# ============================================================

def profile(df, t):

    d = df[np.isclose(df["time_h"], t)]

    if d.empty:
        return None

    # Chọn hàng có y gần 0 nhất
    y0 = d.loc[d["y_km"].abs().idxmin(), "y_km"]

    return (
        d[np.isclose(d["y_km"], y0)]
        .sort_values("x_km")
    )

# ============================================================
# 6. Đọc dữ liệu hai lưới
# ============================================================

data = {
    title: load(path)
    for title, path in files.items()
}

# ============================================================
# 7. Vẽ biểu đồ
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(15, 6),
    constrained_layout=True
)

xs = np.linspace(-XLIM, XLIM, 601)

for ax, (title, df) in zip(axes, data.items()):

    # Vẽ ba thời điểm
    for t, color in zip(times, colors):

        p = profile(df, t)

        if p is None:
            print(f"Không có dữ liệu t = {t} giờ ({title})")
            continue

        ax.plot(
            p["x_km"],
            p["z"],
            color=color,
            lw=2,
            label=f"t = {t} giờ"
        )

    # Thông tin biểu đồ
    ax.set_title(
        title,
        fontsize=14,
        fontweight="bold"
    )

    ax.set_xlabel(
        "Khoảng cách x (km)",
        fontsize=12
    )

    ax.set_xlim(-XLIM, XLIM)

    ax.set_ylim(
        Z0 - DEPTH - 10,
        Z0 + 5
    )

    ax.grid(
        True,
        ls="--",
        alpha=0.5
    )

    ax.legend(fontsize=10)

# Nhãn trục tung
axes[0].set_ylabel(
    "Độ cao Z (m)",
    fontsize=12
)

# Tiêu đề chung
fig.suptitle(
    "Mặt cắt độ cao Z dọc y = 0 tại các thời điểm 6, 18 và 30 giờ",
    fontsize=15,
    fontweight="bold"
)

plt.savefig(
    r"D:\fortran\Mat_cat_Z_Agrid_Cgrid.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


#Vẽ hình phân bố Z và xu thế gió các thời điểm của lưới A và lưới C
import numpy as np, pandas as pd, matplotlib.pyplot as plt

folder = r"D:\fortran"
def load(p):
    df = pd.read_csv(p, sep=r"\s+", header=0, usecols=range(6),
                     names=["time_h", "x_km", "y_km", "u", "v", "z"])
    return df.apply(pd.to_numeric, errors="coerce").dropna(subset=["time_h", "x_km", "y_km", "z"])

def plot_times(path, title, times, ncols=5, km_per_ms=25):
    df = load(path)
    nrows = int(np.ceil(len(times) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6*ncols, 4.3*nrows), constrained_layout=True)
    for ax in np.ravel(axes)[len(times):]:
        ax.axis("off")
    for ax, t in zip(np.ravel(axes), times):
        d = df[np.isclose(df["time_h"], t)]
        if d.empty:
            ax.set_title(f"t = {t} h: không có dữ liệu"); continue
        P = lambda c: (d.pivot(index="y_km", columns="x_km", values=c)
                         .sort_index(axis=0).sort_index(axis=1))
        Z, U, V = P("z"), P("u"), P("v")
        X, Y = U.columns.to_numpy() + 600, U.index.to_numpy() + 600     # dời về 0..1200 km
        XX, YY = np.meshgrid(X, Y)
        u, v = U.to_numpy(), V.to_numpy()
        fin = np.isfinite(u) & np.isfinite(v)
        vmax = np.hypot(u, v).max() if fin.all() else np.nan
        u, v = np.where(fin, u, np.nan), np.where(fin, v, np.nan)
        cf = ax.contourf(XX, YY, Z.to_numpy(), levels=20, cmap="viridis")
        fig.colorbar(cf, ax=ax, label="Z (m)")
        ax.quiver(XX, YY, u, v, color="white", angles="xy", scale_units="xy",
                  scale=1/km_per_ms, width=0.0025)
        ax.set_title(f"t = {t:g} h  (|V|max = {vmax:.1f} m/s)", fontsize=10)
        ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)"); ax.set_aspect("equal")
    fig.suptitle(title, fontsize=14, fontweight="bold")
    plt.show()

times = [1, 6, 12, 13, 14, 15, 16, 17, 18, 30]
plot_times(folder + r"\unit3_Agrid.txt", "Lưới A (A-grid)", times)
plot_times(folder + r"\unit3_Cgrid.txt", "Lưới C (C-grid)", times)


#Vẽ hình trường gió của lưới A và lưới C
import numpy as np, pandas as pd, matplotlib.pyplot as plt

folder = r"D:\fortran"
files = {"Lưới A (A-grid)": folder + r"\unit3_Agrid.txt",
         "Lưới C (C-grid)": folder + r"\unit3_Cgrid.txt"}
times = [1, 6, 12, 13, 14, 15, 16, 17, 18, 30]

STEP  = 2      # vẽ mũi tên mỗi 2 ô lưới (60 km); đặt 1 nếu muốn dày hơn
L_MAX = 55     # km: độ dài mũi tên tối đa
L_MIN = 10     # km: độ dài tối thiểu
V_REF = 30     # m/s: từ mức này trở lên mũi tên đạt L_MAX
V_MIN = 0.05   # m/s: gió nhỏ hơn mức này không vẽ

def load(p):
    df = pd.read_csv(p, sep=r"\s+", header=0, usecols=range(6),
                     names=["time_h", "x_km", "y_km", "u", "v", "z"])
    return df.apply(pd.to_numeric, errors="coerce").dropna(subset=["time_h", "x_km", "y_km", "z"])

DATA = {name: load(p) for name, p in files.items()}

def fields(df, t):
    d = df[np.isclose(df["time_h"], t)]
    if d.empty:
        return None
    P = lambda c: (d.pivot(index="y_km", columns="x_km", values=c)
                     .sort_index(axis=0).sort_index(axis=1))
    U, V, Z = P("u"), P("v"), P("z")
    X = U.columns.to_numpy() + 600          # dời về 0..1200 km
    Y = U.index.to_numpy() + 600
    XX, YY = np.meshgrid(X, Y)
    return XX, YY, Z.to_numpy(), U.to_numpy(), V.to_numpy()

for t in times:
    F = {name: fields(df, t) for name, df in DATA.items()}
    zs = [f[2] for f in F.values() if f is not None]
    if not zs:
        print(f"Không có dữ liệu t = {t} giờ"); continue
    zlo, zhi = min(np.nanmin(z) for z in zs), max(np.nanmax(z) for z in zs)
    levels = np.linspace(zlo, zhi, 21) if zhi - zlo > 1e-6 else 21

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.8), constrained_layout=True)
    for ax, (name, f) in zip(axes, F.items()):
        if f is None:
            ax.set_title(f"{name}, t = {t:g} h: không có dữ liệu"); continue
        XX, YY, Z, u, v = f
        cf = ax.contourf(XX, YY, Z, levels=levels, cmap="turbo")

        fin = np.isfinite(u) & np.isfinite(v)
        spd = np.where(fin, np.hypot(u, v), np.nan)
        vmax = np.nanmax(spd) if fin.all() else np.nan

        # chọn điểm vẽ đối xứng quanh tâm, bỏ gió quá nhỏ hoặc NaN
        ic = np.argmin(np.abs(XX[0] - 600)); jc = np.argmin(np.abs(YY[:, 0] - 600))
        I, J = np.meshgrid(np.arange(XX.shape[1]), np.arange(XX.shape[0]))
        sel = fin & (((I - ic) % STEP) == 0) & (((J - jc) % STEP) == 0) & (spd >= V_MIN)

        # độ dài mũi tên (km) nén theo căn bậc hai, có chặn trên và dưới
        L = L_MIN + (L_MAX - L_MIN) * np.clip(spd / V_REF, 0, 1) ** 0.5
        with np.errstate(invalid="ignore", divide="ignore"):
            uu, vv = u / spd * L, v / spd * L
        ax.quiver(XX[sel], YY[sel], uu[sel], vv[sel], angles="xy", scale_units="xy",
                  scale=1, pivot="mid", color="black", edgecolor="white",
                  linewidth=0.4, width=0.004, headwidth=4, headlength=5)

        if not fin.all():
            ax.text(0.5, 0.04, "gió = NaN (mô hình đã vỡ)", transform=ax.transAxes,
                    ha="center", color="red", fontsize=11, fontweight="bold")
        ax.set_title(f"{name}, t = {t:g} h\n|V|max = {vmax:.1f} m/s", fontsize=12, fontweight="bold")
        ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
        ax.set_xlim(0, 1200); ax.set_ylim(0, 1200); ax.set_aspect("equal")
        ax.grid(alpha=0.25, ls="--")

    fig.colorbar(cf, ax=axes, shrink=0.85, label="Z (m)")
    
       # Lưu từng hình theo thời điểm
    plt.savefig(
        folder + rf"\Truong_gio_Z_t{t:g}h.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
