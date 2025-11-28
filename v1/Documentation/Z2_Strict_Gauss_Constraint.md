# 严格 Gauss 约束实现文档

## 概述

本文档描述了在 ALF 框架中实现的 Z₂ 规范场耦合费米子模型的**严格 Gauss 约束**投影，严格对标 **PRX 10, 041057 (2020) Appendix A**。

---

## 模块 0：路径积分离散化（PRX A1–A5）

### 0.1 Trotter 分解

$$e^{-\beta H} = \left(e^{-\epsilon H}\right)^M, \quad \epsilon = \beta/M = \Delta\tau$$

### 0.2 时间片完备基插入

在每个时间片 $\tau$ 插入 $\tau^z$ 的完备基：
$$\sum_{\{\tau^z(\tau)\}} |\tau^z(\tau)\rangle\langle\tau^z(\tau)|$$

---

## 模块 1：λ 场的引入方式（PRX A5–A6 核心）

### 1.1 Gauss 算符定义

在 orthogonal-fermion/slave-spin 构造中，费米子奇偶性 $(-1)^{n_f}$ 被吸收到 τ 自旋结构中。Gauss 算符为：

$$G_r = Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b$$

**注意**：这里**没有** $(-1)^{n_f}$ 项！

### 1.2 λ 是 τ-independent 的空间场

🚨 **关键点**：λ 只有空间索引，**没有时间索引**！

$$\lambda_i = \pm 1, \quad i \in \text{sites}$$

**不是** `lambda_field(site, tau)`，**而是** `lambda_field(site)`。

### 1.3 Gauss 投影权重（PRX A6）

从 Gauss projector $\hat{P}_i = \frac{1}{2}(1+G_i)$ 出发，经过路径积分推导（PRX A5），得到权重：

$$W_i(\lambda_i; \tau^z_{i,0}, \tau^z_{i,M-1}) \propto e^{\gamma \cdot \tau^z_{i,0} \cdot \lambda_i \cdot \tau^z_{i,M-1}} \tag{A6}$$

其中：
$$\gamma = -\frac{1}{2}\ln[\tanh(\epsilon \cdot h)]$$

这里：
- $\tau^z_{i,0}$：格点 $i$ 在 $\tau=0$ 的 τ 自旋
- $\tau^z_{i,M-1}$：格点 $i$ 在 $\tau=M-1$ 的 τ 自旋
- $h$：横场强度 (Ham_h)
- $\epsilon = \Delta\tau$

### 1.4 物理含义

- $\lambda_i = +1$ → **周期边界条件 (PBC)**：$\tau^z_{i,0}$ 和 $\tau^z_{i,M-1}$ 同号有利
- $\lambda_i = -1$ → **反周期边界条件 (APBC)**：$\tau^z_{i,0}$ 和 $\tau^z_{i,M-1}$ 异号有利

---

## 模块 2：时空 Plaquette 项

PRX Appendix A 明确指出，Gauss 约束在 σ 场上产生额外的**时空 plaquette** 作用量：

$$S_{\text{plaq}} = -K_{\text{plaq}} \sum_{\Box_{i,\tau}} \sigma^z_{\Box}$$

其中 $\sigma^z_{\Box} = \prod_{b \in \Box} \sigma^z_b$ 是时空 plaquette 上的 gauge 场乘积。

---

## 模块 3：费米子行列式的修正（PRX A6 后段）

### 3.1 传播子结构

整条虚时间传播子：
$$\mathcal{B} = B(M) B(M-1) \cdots B(1)$$

其中每个 B-slice：
$$B(\tau) = e^{-\Delta\tau K} \cdot e^{-\Delta\tau V(\sigma(\tau), \tau(\tau))}$$

### 3.2 λ 修正方式（关键！）

🚨 **关键点**：λ **只在时间闭合处**修改费米子行列式，**不是**逐 τ 乘 P(τ)！

❌ 错误写法：$B(\tau) = P(\tau) B_0(\tau)$

✅ 正确写法：
$$\det M = \det(1 + P[\lambda] \cdot \mathcal{B})$$

其中对角矩阵：
$$P_{ij}[\lambda] = \lambda_i \cdot \delta_{ij}$$

### 3.3 物理解释

- $\lambda_i = +1$：费米子在格点 $i$ 满足 PBC
- $\lambda_i = -1$：费米子在格点 $i$ 满足 APBC

λ 通过修改**时间边界条件**影响费米子行列式。

---

## 模块 4：玻色作用量

### 4.1 总玻色作用量

$$S_{\text{total}} = S_{\text{Z2-gauge}} + S_{\tau\text{-path}} + S_{\text{plaq-time}} + S_{\text{Gauss-}\lambda}$$

### 4.2 Gauss λ 作用量

$$S_{\text{Gauss-}\lambda} = -\sum_i \gamma \cdot \tau^z_{i,0} \cdot \lambda_i \cdot \tau^z_{i,M-1}$$

### 4.3 对应权重

$$W_{\text{Gauss}} = \prod_i e^{\gamma \cdot \tau^z_{i,0} \cdot \lambda_i \cdot \tau^z_{i,M-1}}$$

---

## 模块 5：蒙特卡洛更新

### 5.1 更新 λ(i)

翻转 $\lambda_i \to -\lambda_i$：

**玻色权重比率**：
$$R_{\text{bose}}^{(\lambda)} = \exp\left(2\gamma \cdot \tau^z_{i,0} \cdot \tau^z_{i,M-1} \cdot \lambda_i^{\text{old}}\right)$$

**费米子行列式比率**：
$$R_{\text{ferm}}^{(\lambda)} = \frac{\det(1 + P[\lambda^{\text{new}}] \mathcal{B})}{\det(1 + P[\lambda^{\text{old}}] \mathcal{B})}$$

**总比率**：
$$R^{(\lambda)} = R_{\text{bose}}^{(\lambda)} \cdot R_{\text{ferm}}^{(\lambda)}$$

### 5.2 更新 τ 自旋

τ 自旋更新可能改变 $\tau^z_{i,0}$ 或 $\tau^z_{i,M-1}$，从而改变 Gauss 作用量：

$$\Delta S_{\text{Gauss}} = \gamma \left[\tau^z_{i,0}^{\text{new}} \lambda_i \tau^z_{i,M-1}^{\text{new}} - \tau^z_{i,0}^{\text{old}} \lambda_i \tau^z_{i,M-1}^{\text{old}}\right]$$

**玻色权重比率**：
$$R_{\text{bose}}^{(\tau)} = e^{-\Delta S_{\text{Gauss}}}$$

### 5.3 更新 σ 自旋

σ 更新影响 star product，但通常不直接改变 $\tau^z_{i,0}$ 或 $\tau^z_{i,M-1}$（除非通过耦合）。

如果有时空 plaquette 项，需要计算：
$$\Delta S_{\text{plaq}} = -K_{\text{plaq}} \left[\sigma^z_{\Box}^{\text{new}} - \sigma^z_{\Box}^{\text{old}}\right]$$

---

## 模块 6：观测量

### 6.1 Gauss 算符期望值

$$\langle G_r \rangle = \left\langle Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b \right\rangle$$

应接近 $+1$（或 $Q_r$）。

### 6.2 Gauss 约束违反度

$$\langle (G_r - Q_r)^2 \rangle \approx 0$$

---

## 使用方法

### 参数设置

```
UseStrictGauss = .true.
GaussSector = "even"    ! "even", "odd", "staggered"
```

### 示例参数文件

```
Model = Z2_Matter
Lattice_type = Square
L1 = 6
L2 = 6
ham_T = 1.0
ham_TZ2 = 1.0
Ham_K = 1.0
Ham_h = 1.0
Ham_g = 1.0
Beta = 10.0
Dtau = 0.1
UseStrictGauss = .true.
GaussSector = "even"
```

---

## 实现细节

### 新增场变量

```fortran
! Lambda 场：τ-independent，只有空间索引
Integer, allocatable :: lambda_field(:)  ! lambda_field(site) = +1 或 -1

! tau^z 场（已存在，需要访问首尾）
! tau_z(site, tau=0) 和 tau_z(site, tau=M-1)

! 背景电荷数组
Integer, allocatable :: Q_background(:)

! Gauss 耦合常数
Real (Kind=Kind(0.d0)) :: Gamma_Gauss  ! γ = -0.5 * ln(tanh(ε*h))
```

### 核心函数

| 函数名 | 功能 | 公式 |
|--------|------|------|
| `Setup_Gauss_constraint()` | 初始化 λ 场和计算 γ | $\gamma = -\frac{1}{2}\ln[\tanh(\epsilon h)]$ |
| `Get_Tau_Z_At_Time_0(I)` | 获取 τ=0 处的 τ^z | $\tau^z_{i,0}$ |
| `Get_Tau_Z_At_Time_M1(I)` | 获取 τ=M-1 处的 τ^z | $\tau^z_{i,M-1}$ |
| `Compute_Gauss_Action_PRX(I)` | 计算单点 Gauss 作用量 | $S_i = -\gamma \tau^z_{i,0} \lambda_i \tau^z_{i,M-1}$ |
| `Compute_Gauss_Weight_Ratio_Lambda_PRX(I)` | λ 翻转的权重比 | $R = e^{2\gamma \tau^z_{i,0} \tau^z_{i,M-1} \lambda_i^{\text{old}}}$ |
| `Compute_Delta_S_Gauss_Tau_Update(...)` | τ 更新的 ΔS | $\Delta S = S^{\text{new}} - S^{\text{old}}$ |
| `Compute_Gauss_Operator(I, nt, GRC)` | 计算 Gauss 算符（观测量） | $G_r = Q_r \tau_r^x X_r$（无 $(-1)^{n_f}$） |

### 费米子行列式修正

```fortran
! 计算完整传播子
Btotal = B(M) * B(M-1) * ... * B(1)

! 构造 P[λ] 对角矩阵
P_lambda(i,i) = lambda_field(i)

! 修正后的 Green 函数逆
Ginv = I + P_lambda * Btotal

! 行列式
detM = det(Ginv)
```

---

## 与文献的对应关系

| 本文档内容 | 对应 PRX 公式 |
|------------|---------------|
| λ 是 τ-independent | Appendix A 整体结构 |
| $W_i = e^{\gamma \tau^z_0 \lambda \tau^z_{M-1}}$ | (A6) |
| $\gamma = -\frac{1}{2}\ln[\tanh(\epsilon h)]$ | (A6) |
| $\det(1 + P[\lambda]\mathcal{B})$ | A6 后段 |
| 时空 plaquette | A6 后 "spatiotemporal plaquette" |

---

## 文件修改列表

主要修改的文件：
- `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`
  - **关键修改**：`lambda_field(:)` 改为 τ-independent（一维数组）
  - 添加 `Gamma_Gauss` 参数，计算公式 $\gamma = -\frac{1}{2}\ln[\tanh(\epsilon h)]$
  - 添加 `Get_Tau_Z_At_Time_0(I)` 和 `Get_Tau_Z_At_Time_M1(I)` 函数
  - 添加 `Compute_Gauss_Action_PRX(I)` 函数
  - 添加 `Compute_Gauss_Weight_Ratio_Lambda_PRX(I)` 函数
  - 添加 `Compute_Delta_S_Gauss_Tau_Update(...)` 函数
  - 修改 `Compute_Gauss_Operator` 去除 $(-1)^{n_f}$（PRX orthogonal-fermion 构造）
  - 修改 `Setup_Gauss_constraint` 初始化 τ-independent λ 场
  - 修改 `S0` 函数中 λ 更新使用 PRX A6 公式
  - 修改 `Hamiltonian_set_nsigma` 中 λ 初始化

**待完成**（需要 ALF 核心框架支持）：
- 费米子行列式修正：$\det(1 + P[\lambda]\mathcal{B})$ 而非逐 τ 乘 P(τ)
- 时空 plaquette 项 $S_{\text{plaq}}$

---

## 注意事项

1. **λ 不是逐 τ 的**：这是最关键的点。λ 只有空间索引。

2. **费米子边界条件**：λ 通过修改时间边界条件（PBC/APBC）影响费米子行列式，不是逐 τ 乘对角矩阵。

3. **γ 的计算**：需要 $h > 0$ 才能定义 γ。当 $h \to 0$ 时，$\gamma \to \infty$。

4. **初始化**：初始配置应满足 Gauss 约束。

5. **时空 plaquette**：如果模型包含 gauge 场动力学，需要添加时空 plaquette 项。

---

## 作者

ALF Collaboration
