import os
folder = r"D:\fortran"
print(os.listdir(folder))


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cases = [
    (r"D:\fortran\Z_final_Agrid.txt", "Lưới A (A-grid)", 2850),
    (r"D:\fortran\Z_final_Cgrid.txt", "Lưới C (C-grid)", 2850),
]
fig, axes = plt.subplots(
    1, 2, figsize=(13, 5.5), constrained_layout=True
)

for ax, (filename, title, z_center) in zip(axes, cases):

    # Đọc dữ liệu độ cao
    df = pd.read_csv(filename, sep=r"\s+", comment="#")

    # Chuẩn hóa tên ba cột dữ liệu
    df = df.iloc[:, :3]
    df.columns = ["x_km", "y_km", "z_m"]

    # Lấy vùng quanh tâm từ -600 đến 600 km
    df = df[
        df["x_km"].between(-600, 600)
        & df["y_km"].between(-600, 600)
    ]

    # Chuyển dữ liệu thành lưới hai chiều
    Z = df.pivot(
        index="y_km",
        columns="x_km",
        values="z_m"
    )
    Z = Z.sort_index(axis=0).sort_index(axis=1)

    X = Z.columns.to_numpy()
    Y = Z.index.to_numpy()
    XX, YY = np.meshgrid(X, Y)
    ZZ = Z.to_numpy()

    # Vẽ trường độ cao và đường đồng mức
    cf = ax.contourf(
        XX, YY, ZZ,
        levels=20,
        cmap="viridis"
    )

    cs = ax.contour(
        XX, YY, ZZ,
        levels=10,
        colors="white",
        linewidths=0.7
    )

    ax.clabel(cs, inline=True, fontsize=8, fmt="%.0f")

    # Đánh dấu tâm miền
    ax.plot(
        0, 0, "r+",
        markersize=10,
        markeredgewidth=2,
        label="Tâm miền"
    )

    ax.set_title(
        f"{title}\nĐộ cao nhỏ nhất dự kiến: {z_center} m"
    )
    ax.set_xlabel("Khoảng cách theo hướng X (km)")
    ax.set_ylabel("Khoảng cách theo hướng Y (km)")
    ax.set_xlim(-600, 600)
    ax.set_ylim(-600, 600)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right")

    fig.colorbar(cf, ax=ax, label="Độ cao Z (m)")

fig.suptitle(
    "Phân bố độ cao mặt nước ở thời điểm kết thúc mô phỏng",
    fontsize=14
)

plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. File dữ liệu
folder = r"D:\fortran"

files = {
    "Lưới A (A-grid)": folder + r"\unit3_Agrid.txt",
    "Lưới C (C-grid)": folder + r"\unit3_Cgrid.txt",
}

times = [0, 6, 18, 30]
colors = ["black", "blue", "green", "red"]

# 2. Thiết lập đồ thị
fig, axes = plt.subplots(
    1, 2,
    figsize=(14, 5.5),
    sharey=True,
    constrained_layout=True
)

for ax, (title, file_path) in zip(axes, files.items()):

    df = pd.read_csv(file_path, sep=r"\s+", header=0)
    df = df.iloc[:, :6]
    df.columns = ["time_h", "x_km", "y_km", "u", "v", "z"]

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(
        subset=["time_h", "x_km", "y_km", "z"]
    )

    for t, color in zip(times, colors):

        # Chọn thời điểm cần vẽ
        d = df[np.isclose(df["time_h"], t)].copy()

        if d.empty:
            print(f"Không có dữ liệu t={t} giờ trong {file_path}")
            continue

        # Chọn hàng có y gần 0 nhất
        y_nearest = d.loc[d["y_km"].abs().idxmin(), "y_km"]

        d = d[
            np.isclose(d["y_km"], y_nearest)
            & d["x_km"].between(-800, 800)
        ].sort_values("x_km")

        if d.empty:
            continue

        ax.plot(
            d["x_km"],
            d["z"],
            color=color,
            linewidth=2,
            label=f"t = {t} giờ"
        )

    # Tiêu đề và nhãn trục
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Khoảng cách x (km)", fontsize=11)
    ax.set_xlim(-800, 800)

    # Giới hạn trục Z
    ax.set_ylim(2850, 3000)
    ax.set_yticks(np.arange(2850, 3001, 25))

    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()

# Tiêu đề chung
fig.suptitle(
    "So sánh mặt cắt độ cao Z giữa lưới A và lưới C",
    fontsize=15,
    fontweight="bold"
)

axes[0].set_ylabel("Độ cao Z (m)", fontsize=11)

plt.show()

