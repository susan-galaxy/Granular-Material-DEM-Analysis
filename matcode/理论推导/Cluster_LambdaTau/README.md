# Cluster λ-τ 分析工具箱

本目录给出一套自洽的 MATLAB 代码，用来从你的 DEM 数据中提取
**特征长度 λ** 与 **特征时间 τ**，并按理论框架（弹道-热化双流体机制）
生成下面三张论文级别图：

1. **Fig 1** &nbsp; 特征长度 λ(φ) 随体积分数变化：理论曲线 vs DEM 测量值。
2. **Fig 2** &nbsp; 相图（两块面板）：
   - (A) (φ, e_n) 平面，含理论临界曲线 λ = L/2，DEM 散点；
   - (B) 无量纲 (Λ* = λ/(L/2), τ* = τ_H/τ_mfp) 平面，理论轨迹 + DEM 轨迹 + 状态分区。
3. **Fig 3** &nbsp; 状态演化预测：
   - (A) T_center / T_wall 随 φ 的预测曲线（含线性 cosh 近似 + 非线性 ODE 精确解，
        后者给出与 PPT 中“骤减”一致的尖锐相变）；
   - (B) 多个 φ 下的 T(z)/T_w 空间剖面（理论 + DEM）。

---

## 文件清单

| 文件 | 作用 |
|------|------|
| `chi_CS.m`              | Carnahan-Starling pair-correlation χ(φ) |
| `theory_lambda.m`       | 理论特征长度 λ(φ, e_n, d)，含可调前因子 `Cscale` |
| `theory_tau.m`          | 理论 Haff 冷却时间 τ_H 与平均自由时间 τ_mfp |
| `solve_T_profile.m`     | 解一维稳态能量平衡 ODE 得到 T(z)/T_w（首积分法 + cosh fallback） |
| `phi_from_name.m`       | 从文件夹/文件名解析 `phi`（如 `phi0.02`） |
| `dem_extract_run.m`     | 对单个 .mat 文件提取：n, T, T_wall, T_center, T(r_wall) 剖面、σ²_n、λ_DEM、τ_H, τ_mfp |
| `dem_load_all.m`        | 批量扫描目录 + 调 `dem_extract_run` |
| `Fig1_LengthVsPhi.m`    | 生成 Fig 1 |
| `Fig2_PhaseDiagram.m`   | 生成 Fig 2 |
| `Fig3_Prediction.m`     | 生成 Fig 3 |
| `Main_run.m`            | 一键运行：加载 → 提取 → 出图 |

---

## 使用方法

1. 用 `Pos_Vel_Extrace_v3.m` 生成 `Data_*.mat`（你已经做了）。
2. 打开 `Main_run.m`，修改 **USER SETTINGS** 区块：
   - `DATA_DIR`：上一步保存的 .mat 文件夹路径（Linux/Windows 都行）。
   - `L`, `d`, `rho`, `en`：与你 DEM 设置一致即可。
3. 在 MATLAB 里运行 `Main_run`。三张图会保存到 `figures/` 子目录下（.png + .fig）。

> 如果 `DATA_DIR` 不存在，脚本会回退到一组人造数据，给你看图的视觉效果。

---

## 物理含义速览

理论核心（与你的 PPT/Word 一致）：

- 1D 稳态能量平衡：&nbsp; `d(κ dT/dz)/dz − Γ + S = 0`
- κ = κ₀ · n · d · √T ; &nbsp; Γ = γ₀ · n² d² χ(φ) (1−e_n²) T^{3/2}
- 微重力下 P = nT = const, 代入后得到 χ(φ) 主导的耗散；
- 量纲分析给出 **能量穿透深度**：

  λ²  =  κ₀ / [ γ₀ (1−e_n²) n d χ(φ) ]  ⇒
  **λ(φ) = C · d / √[(1−e_n²) φ χ(φ)]**

- 临界判据 λ ≈ L/2 → 中心温度坍缩 → 团簇形成；
- 自然时间尺度：Haff 冷却时间 τ_H = 1 / [γ₀(1−e²) n d² χ(φ) √(T/m)]。

`Fig1_LengthVsPhi` 会自动从 DEM 数据拟合前因子 `Cscale`，
并把它回传给 Fig 2 / Fig 3 保证一致性。

---

## 预测怎么用

把 `Fig3_Prediction` 输出的 `phi_c^{predict}` 与 `Fig1` 的 `phi_c^{th}`
配合，就可以在 (φ, e_n) 相图上读出：
**“给定 e_n 与系统尺寸 L，预测临界 φ_c”**。

如果做新参数（不同 L、不同 e_n、不同 d），只需：

```matlab
phi_grid = logspace(-3, log10(0.3), 200);
lam      = theory_lambda(phi_grid, en_new, d_new, 'Cscale', C_fitted);
phi_c    = interp1(lam, phi_grid, L_new/2);
```

即可立刻给出该工况下的预测临界体积分数。
