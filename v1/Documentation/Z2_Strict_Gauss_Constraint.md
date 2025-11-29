# 严格 Gauss 约束实现文档

## 概述

本文档描述了在 ALF 框架中实现的 Z₂ 规范场耦合费米子模型的**严格 Gauss 约束**投影，对标 **PRX 10, 041057 (2020) Appendix A**。

---

## 理论基础

### 1. Gauss 算符定义

在 orthogonal-fermion/slave-spin 构造中：

$$G_r = Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b$$

**注意**：这里**没有** $(-1)^{n_f}$ 项（已被吸收到 τ 结构中）。

### 2. λ 场（PRX A5-A6）

λ 是 **τ-independent** 的空间场：
$$\lambda_i = \pm 1, \quad i \in \text{sites}$$

### 3. Gauss 投影权重（PRX A6）

$$W_i(\lambda_i; \tau^z_{i,0}, \tau^z_{i,M-1}) \propto e^{\gamma \cdot \tau^z_{i,0} \cdot \lambda_i \cdot \tau^z_{i,M-1}}$$

其中：
$$\gamma = -\frac{1}{2}\ln[\tanh(\epsilon \cdot h)]$$

### 4. 费米子行列式修正

$$\det M = \det(1 + P[\lambda] \cdot \mathcal{B})$$

其中对角矩阵 $P_{ij}[\lambda] = \lambda_i \cdot \delta_{ij}$

---

## 使用方法

### 参数设置

```fortran
&VAR_Z2_Matter
UseStrictGauss = .T.
GaussSector    = "even"    ! "even", "odd", "staggered"
/
```

### GaussSector 定义

| GaussSector | Q_r 定义 | 适用场景 |
|-------------|----------|----------|
| `"even"` | $Q_i = +1$ 所有 site | 标准物理 sector |
| `"odd"` | $Q_i = -1$ 所有 site | 全局奇 sector |
| `"staggered"` | $Q_{x,y} = (-1)^{x+y}$ | A/B 子格交替 |

---

## 实现状态

### ✅ 已完成

| 功能 | 文件 | 说明 |
|------|------|------|
| λ 场初始化 | `Hamiltonian_Z2_Matter_smod.F90` | `lambda_field(site)` 初始化为 +1 |
| Q_background 初始化 | `Hamiltonian_Z2_Matter_smod.F90` | 根据 GaussSector 设置 |
| Gamma_Gauss 计算 | `Hamiltonian_Z2_Matter_smod.F90` | $\gamma = -\frac{1}{2}\ln[\tanh(\epsilon h)]$ |
| λ 翻转玻色权重 | `Compute_Gauss_Weight_Ratio_Lambda_PRX` | $R = e^{-2\gamma \tau^z_0 \tau^z_{M-1} \lambda_{\text{old}}}$ |
| Sweep_Lambda 循环 | `main.F90` | 遍历 site（不遍历 τ），Metropolis 接受 |
| Gauss 观测量 | `Obser` | 测量 $\langle G_r \rangle$ 和 GaussViol |
| GaussViol 诊断 | `Measure_GaussViolation_Diagnostic` | 实时检查 Gauss 约束 |

### ⚠️ 已禁用（数值不稳定）

| 功能 | 文件 | 问题 |
|------|------|------|
| P[λ] 应用到 B 矩阵 | `Apply_P_Lambda_To_B` | 函数为空，不修改 B |
| 费米子行列式比率 | `Lambda_Ferm_Ratio_site` | 直接返回 1 |
| Sherman-Morrison 更新 | `Lambda_Update_Green_site` | 调用被注释掉 |

**当前替代方案**：在 `Sweep_Lambda` 后调用 CGR 重建 Green 函数。

### 🔴 核心问题

1. **σ 更新没有 Gauss 约束**
   - σ 翻转时 $X_r = \prod \sigma^x$ 改变，可能导致 $G_r$ 从 +1 变为 -1
   - 当前实现中，σ 更新完全没有 Gauss 权重检查
   - 这导致 $\langle G_r \rangle \approx 0$ 而不是 +1

2. **P[λ] 修改与 ALF wrap 机制不兼容**
   - ALF 的 UDV 分解和稳定化方案与 P[λ] 修改冲突
   - 需要更深层次的集成才能正确工作

---

## 文件修改列表

### `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`

**新增变量**：
- `lambda_field(:)` - λ 场（一维，site-only）
- `Q_background(:)` - 背景电荷
- `Gamma_Gauss` - PRX A6 耦合常数
- `B_lambda_slice(:,:)` - 保存的 B 矩阵（目前未使用）

**新增函数**：
- `Setup_Gauss_constraint()` - 初始化
- `Get_Tau_Z_At_Time_0(I)` / `Get_Tau_Z_At_Time_M1(I)` - 获取边界 τ^z
- `Compute_Gauss_Action_PRX(I)` - 计算单点 Gauss 作用量
- `Compute_Gauss_Weight_Ratio_Lambda_PRX(I)` - λ 翻转玻色权重
- `Compute_Star_Product_X(I, nt)` - 计算 star product
- `Compute_Gauss_Operator_Int(I, nt)` - 计算 Gauss 算符
- `Apply_P_Lambda_To_B(B, nf)` - P[λ] 应用（**目前禁用**）
- `Lambda_Ferm_Ratio_site(i, G, R)` - 费米子比率（**返回 1**）
- `Lambda_Update_Green_site(i, G, R)` - SM 更新（**未调用**）
- `Sweep_Lambda(G, Phase)` - λ sweep 主循环
- `Measure_GaussViolation_Diagnostic(sweep)` - 诊断输出

### `Prog/wrapur_mod.F90`

- 在 `nt == Ltrot` 时调用 `ham%Apply_P_Lambda_To_B`（**目前该函数为空**）

### `Prog/main.F90`

- CGR 后调用 `ham%Sweep_Lambda(GR, Phase)`
- λ sweep 后调用 CGR 重建 G

---

## 待解决问题

### 高优先级

1. **实现 σ 更新的 Gauss 约束**
   - 在 `S0` 函数中检查 σ 翻转是否违反 Gauss 约束
   - 如果 $G_r$ 从 +1 变为 -1，拒绝该更新

2. **正确实现 P[λ] 修改**
   - 需要与 ALF 的 wrap 机制兼容
   - 可能需要修改 CGR/WRAPUR 的核心逻辑

### 低优先级

- 时空 plaquette 项（如需要 3D gauge action）
- GaussSector odd/staggered 测试

---

## 数值稳定性

### γ 参数

当 $h \to 0$ 时，$\gamma \to \infty$。实现中使用：
- 小 $\epsilon h$ 渐近展开
- 最大值截断 `Gamma_max = 100`

### λ 翻转权重

$$R_{\text{bose}} = e^{-2\gamma \cdot \tau^z_0 \cdot \tau^z_{M-1} \cdot \lambda_{\text{old}}}$$

**注意负号**：当前配置"好"时 $R < 1$（保持），"坏"时 $R > 1$（翻转）。

---

## 诊断输出示例

```
============================================================
 GAUSS CONSTRAINT DIAGNOSTIC - Sweep      100
============================================================
   <G_r>         (should be ~1): 0.10000000E+01
   GaussViol     (should be ~0): 0.12345678E-11
   Lambda_BC_sum (PRX A6 check): 0.50000000E+00
   Gamma_Gauss:                    1.234567
------------------------------------------------------------
```

如果 `<G_r> ≈ 0` 而不是 `+1`，说明 σ 更新没有被 Gauss 约束限制。

---

## 参考文献

- PRX 10, 041057 (2020) Appendix A

---

## 作者

ALF Collaboration
