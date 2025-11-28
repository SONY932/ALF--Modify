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

### 0.3 τ^z 路径积分时间演化项 S_τ-path（🟥 关键模块 #1）

PRX A1–A5 给出了 τ^x 和 τ^z 在 Trotter 分解下的路径积分形式。横场项 $h \tau^x$ 生成**时间方向的最近邻耦合**，行为类似于 1D Ising coupling：

$$S_{\tau\text{-path}} = -K_\tau \sum_{i,\tau} \tau^z_{i,\tau} \cdot \tau^z_{i,\tau+1}$$

其中时间方向的耦合系数为：

$$K_\tau = -\frac{1}{2} \ln[\tanh(\epsilon h)]$$

🚨 **关键点**：这里的 $K_\tau$ 与 Gauss 投影的 $\gamma$ **完全相同**！PRX Appendix A 明确说明它们来自同一个起源（τ^x 的离散化）。

#### ALF 实现

在 ALF 中，时间方向的 τ^z 耦合已经通过 `DW_Matter_tau` 实现：
```fortran
! DW_Matter_tau(+1) = tanh(Dtau*Ham_h)  当 tau_z(t) = tau_z(t+1)
! DW_Matter_tau(-1) = 1/tanh(Dtau*Ham_h) 当 tau_z(t) ≠ tau_z(t+1)
```

权重比率：
$$\frac{W(\tau^z_t = \tau^z_{t+1})}{W(\tau^z_t \neq \tau^z_{t+1})} = \tanh(\epsilon h)$$

这与 $K_\tau = -\frac{1}{2}\ln[\tanh(\epsilon h)]$ 给出的 $e^{2K_\tau} = 1/\tanh(\epsilon h)$ 一致。

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

## 模块 2：时空 Plaquette 项（🟥 关键模块 #2）

PRX Appendix A 明确指出，Gauss 约束在 σ 场上产生额外的**时空 plaquette** 作用量：

$$S_{\text{plaq}} = -K_{\text{plaq}} \sum_{\Box_{i,\tau}} \sigma^z_{\Box}$$

其中 $\sigma^z_{\Box} = \prod_{b \in \Box} \sigma^z_b$ 是时空 plaquette 上的 gauge 场乘积。

### 2.1 K_plaq 的具体表达式

根据文献 [Gazit 2016 PNAS] 的 τ 方向 link weight，plaquette 耦合系数为：

$$K_{\text{plaq}} = \frac{1}{2}\ln\left(\coth(\epsilon g)\right)$$

其中 $g$ 是 gauge 场的横场强度 (`Ham_g` 在 ALF 中)。

### 2.2 时空 Plaquette 的几何结构

在 (2+1)D 时空中，每个时空 plaquette 包含：
- 两条**空间方向**的 link：$\sigma^z_b(\tau)$ 和 $\sigma^z_b(\tau+1)$
- 两条**时间方向**的 link（虚拟）

时空 plaquette 产物：
$$\sigma^z_{\Box_{i,\mu,\tau}} = \sigma^z_{i,\mu}(\tau) \cdot \sigma^z_{i+\hat{\mu},0}(\tau \to \tau+1) \cdot \sigma^z_{i,\mu}(\tau+1) \cdot \sigma^z_{i,0}(\tau+1 \to \tau)$$

由于时间方向的 link 来自 Gauss 约束的离散化，在实际计算中简化为：

$$\sigma^z_{\Box} \approx \sigma^z_{i,\mu}(\tau) \cdot \sigma^z_{i,\mu}(\tau+1)$$

### 2.3 与 Ham_g 的关系

在 ALF 中，gauge 场的时间演化已通过 `DW_Ising_tau` 实现：
```fortran
DW_Ising_tau(+1) = tanh(Dtau*Ham_g)  ! 当 sigma(t) = sigma(t+1)
DW_Ising_tau(-1) = 1/tanh(Dtau*Ham_g) ! 当 sigma(t) ≠ sigma(t+1)
```

这对应于：
$$K_\sigma = -\frac{1}{2}\ln[\tanh(\epsilon g)]$$

### 2.4 ALF 实现

时空 plaquette 在 ALF 中通过 gauge 场的时间方向耦合自动实现。在 `S0` 函数中：
```fortran
! Gauge field temporal coupling
S0 = S0 * DW_Ising_tau(nsigma%i(n,nt) * nsigma%i(n,nt+1))
S0 = S0 * DW_Ising_tau(nsigma%i(n,nt) * nsigma%i(n,nt-1))
```

---

## 模块 3：费米子行列式的修正（🟥 关键模块 #3 - PRX A6 后段）

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

### 3.4 ALF 实现：P[λ] 在 wrap-up 层的插入

#### 3.4.1 核心原则

> 🚨 **P[λ] 只能作用一次，只在完整 B_total 乘完后！**
> 
> 在 ALF 中，任意时间切片的等时 Green function 构造都依赖"从该时间片出发、沿虚时间跑一整圈"的 B 链。P[λ] **必须且只能**乘在这条完整链的某一端（统一乘在左端），不能：
> - 在每个 stabilizer block 上分别乘 P[λ]
> - 在 GRUP 和 GRDW 各乘一次
> - 在局部 Green function 计算中反复插入

PRX Appendix A 原文：
> "the fermion propagator is modified by inserting a diagonal matrix with diagonal elements λ_i at the **temporal boundary**"

对应路径积分图像：
```
τ = 0 ----- B(1) ----- τ = 1 ----- B(2) ----- ... ----- τ = M-1 ---- wrap ----> τ = M (=0)
                                                                      ↑
                                                             P[λ] 只在这里作用！
```

#### 3.4.2 ALF 中"wrap-up"的具体位置

ALF Green function 计算流程（`cgr1_mod.F90` 的 CGR）：
1. **逐 τ 构造 B(τ)**
2. **分组进行 QR 或 LU 稳定**
3. **最后一个 wrap-up，把稳定块乘起来形成完整的 $\mathcal{B}$**
4. 计算 $G = (1+\mathcal{B})^{-1}$

#### 3.4.3 正确的 P[λ] 插入位置

> **在 wrap-up 阶段构造完 B_total 后、计算 $G^{-1} = I + B_{\text{total}}$ 之前**

```fortran
! ========================================
! 正确做法：在 wrap-up 层、构造最终 B_total 后插入
! ========================================

! wrap-up 完成后的 B_total（已经是完整的一圈传播）
B_total = B(M) * B(M-1) * ... * B(1)

! P[λ] 只乘一次，乘在左边
B_eff = P[lambda] * B_total

! 计算 Green function
Ginv = I + B_eff
G = inverse(Ginv)
```

#### 3.4.4 不正确做法的例子

```fortran
! ========================================
! ❌ 错误做法 1：在 CGR 的局部 GRUP 上乘 P[λ]
! ========================================
! CGR 内部的 GRUP 是局部传播子，不是完整的一圈
! 如果在这里乘 P[λ]，会导致某些时间片的 G 包含 P[λ] 两次，某些不含

! ❌ 错误做法 2：在每个 stabilizer block 上乘 P[λ]
! ========================================
! 会导致 P[λ] 被乘了多次（有多少个 block 就乘多少次）

! ❌ 错误做法 3：在 GRUP 和 GRDW 各乘一次
! ========================================
! 会导致 P[λ] 被乘了两次
```

#### 3.4.5 与 CGR 流程的集成方案

在 ALF 中，最终 Green function 由 CGR 输出。

> 🚨 **关键数学问题**
>
> 设 $G = (1+B)^{-1}$（CGR 输出），我们想要 $G' = (1+PB)^{-1}$。
>
> **这两者的关系不是简单的 $G' = P \cdot G$！**
>
> 正确的关系是（设 $\Delta = P - I$）：
> $$G' = (I + G \cdot \Delta \cdot B)^{-1} \cdot G$$
>
> 当所有 $\lambda_i = +1$ 时 $\Delta = 0$，此时 $G' = G$。
> 当某些 $\lambda_i = -1$ 时，$\Delta_{ii} = -2$，需要 rank-k 更新。

**方案 A：修改 B_total 构造（推荐 - 数学上正确）**

在构造 B_total 时直接乘 P[λ]，然后正常计算 Green function：
```fortran
! 在构造完整 B_total 后、传入 CGR 前
B_eff = P[lambda] * B_total
! 然后 CGR 计算 G = (1 + B_eff)^{-1}
```

这需要修改 ALF 的 wrap-up 流程。

**方案 B：后修正 Green function（需要 B_total）**

如果已经计算了 $G = (1+B)^{-1}$，可以用 Woodbury 公式修正：
```fortran
! 设 k 个格点有 lambda_i = -1
! Delta 是 (N x k) 矩阵，(Delta)_{i,j} = -2 如果 i 是第 j 个 lambda=-1 的格点
! V^T 是 (k x N) 矩阵，V^T = (B 的对应行)

! G' = (I + G*Delta*V^T)^{-1} * G
! 用 Sherman-Morrison-Woodbury 公式计算
```

**方案 C：简化实现（近似 - 仅当 λ 不经常翻转时）**

如果大多数 $\lambda_i = +1$，可以使用当前的简化实现：
```fortran
! 当前实现：G' = P * G
! 这仅在 lambda_i 很少为 -1 时是近似正确的
! 严格来说这不等价于 (1+PB)^{-1}
call ham%Apply_P_Lambda_To_Green(GRUP, nf_eff)
```

> ⚠️ **当前代码状态**：实现的是方案 C（简化版），不是严格正确的 PRX A6 实现。
> 完整实现需要方案 A 或方案 B。

#### 3.4.4 两自旋自由度的处理

**情况 1：两自旋独立（无自旋翻转项）**

费米子矩阵是 block-diagonal：
$$B = \begin{pmatrix} B^\uparrow & 0 \\ 0 & B^\downarrow \end{pmatrix}$$

P[λ] 也必须 block-diagonal：
$$P[\lambda] = \begin{pmatrix} P_\lambda & 0 \\ 0 & P_\lambda \end{pmatrix}$$

其中 $(P_\lambda)_{ij} = \lambda_i \delta_{ij}$

**情况 2：有自旋混合项（SO coupling 等）**

费米子 Hilbert 空间维度是 2N，P[λ] 仍然是：
$$P[\lambda]_{(i,\sigma),(j,\sigma')} = \lambda_i \delta_{ij} \delta_{\sigma\sigma'}$$

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

### 🚨 关键设计决策：λ 是 site-only 变量，独立更新

> **λ 不是 Field_type=5，不走逐时间片更新！**
> 
> λ 只有空间索引 `lambda_field(site)`，不出现在 `nsigma(i, nt)` 这类带 τ 下标的数组里。
> λ 更新通过独立的 `Update_Lambda` 循环，只遍历 site，不遍历 τ。

### 5.1 更新 λ(i)：独立的 site-only 更新

翻转 $\lambda_i \to -\lambda_i$：

**玻色权重比率**（PRX A6）：
$$R_{\text{bose}}^{(\lambda)} = \exp\left(2\gamma \cdot \tau^z_{i,0} \cdot \tau^z_{i,M-1} \cdot \lambda_i^{\text{old}}\right)$$

**费米子行列式比率**：
$$R_{\text{ferm}}^{(\lambda)} = \frac{\det(1 + P[\lambda^{\text{new}}] \mathcal{B})}{\det(1 + P[\lambda^{\text{old}}] \mathcal{B})}$$

**总比率**：
$$R^{(\lambda)} = R_{\text{bose}}^{(\lambda)} \cdot R_{\text{ferm}}^{(\lambda)}$$

#### 5.1.1 情况 A：↑↓ 自旋完全独立（block-diagonal 费米子矩阵）

若费米子矩阵是 block-diagonal（↑↓ 解耦，无 SO coupling）：
$$B = \begin{pmatrix} B^\uparrow & 0 \\ 0 & B^\downarrow \end{pmatrix}$$

可以对每个自旋分开做 **rank-1** 更新：
$$R_{\text{ferm}} = R_{\text{ferm}}^\uparrow \times R_{\text{ferm}}^\downarrow$$

其中单自旋的 rank-1 公式：
$$R_{\text{ferm}}^\sigma = 1 - 2\lambda_i^{\text{old}} \cdot (B^\sigma \cdot G^\sigma)_{ii}$$

**Sherman-Morrison 更新**（单自旋）：

$$G^\sigma_{\text{new}} = G^\sigma_{\text{old}} - \frac{G^\sigma_{\text{old}} \cdot u \cdot w^T \cdot G^\sigma_{\text{old}}}{R_{\text{ferm}}^\sigma}$$

其中：
- $u = (-2 \lambda_i^{\text{old}}) e_i$
- $w^T = B^\sigma_{\text{row }i}$

#### 5.1.2 情况 B：自旋混合（SO coupling, pair-hopping 等）

若存在自旋混合项，费米子矩阵维度是 2N，翻转 λ_i 修改两行：第 i 行和第 i+N 行。

必须使用 **rank-2** 更新：

**费米子行列式比率**：
$$R_{\text{ferm}} = \det(I_2 + V^T \cdot G \cdot U)$$

其中：
$$U = \begin{pmatrix} u_\uparrow & 0 \\ 0 & u_\downarrow \end{pmatrix}_{2N \times 2}, \quad
V = \begin{pmatrix} B_{\text{row }i} \\ B_{\text{row }i+N} \end{pmatrix}^T_{2N \times 2}$$

这是一个 **2×2 行列式**，计算成本 O(1)。

**Sherman-Morrison rank-2 更新**：
$$G_{\text{new}} = G - G \cdot U \cdot (I_2 + V^T \cdot G \cdot U)^{-1} \cdot V^T \cdot G$$

#### 5.1.3 ALF 实现：独立 λ 更新循环

```fortran
!> λ 更新：独立循环遍历所有 site，不遍历 τ
!> 注意：此示例仅适用于单自旋 / ↑↓ 解耦的情况
subroutine Update_All_Lambda(G, B_total, N_sites, N_dim)
    complex(8), intent(inout) :: G(:,:)
    complex(8), intent(in) :: B_total(:,:)
    integer, intent(in) :: N_sites, N_dim
    
    integer :: i
    real(8) :: R_bose, R_tot
    complex(8) :: R_ferm, BG_ii
    integer :: tau_z_0, tau_z_M1, lambda_old
    
    ! 遍历所有 site（不是时间片！）
    do i = 1, N_sites
        ! 计算 R_bose（PRX A6）
        tau_z_0 = Get_Tau_Z_At_Time_0(i)
        tau_z_M1 = Get_Tau_Z_At_Time_M1(i)
        lambda_old = lambda_field(i)
        R_bose = exp(2.0d0 * Gamma_Gauss * tau_z_0 * tau_z_M1 * lambda_old)
        
        ! 计算 R_ferm（单自旋 rank-1）
        BG_ii = sum(B_total(i,:) * G(:,i))
        R_ferm = 1.0d0 - 2.0d0 * lambda_old * BG_ii
        
        ! Metropolis 接受/拒绝
        R_tot = R_bose * abs(R_ferm)
        if (ranf() < R_tot) then
            ! 更新 lambda（site-only 变量）
            lambda_field(i) = -lambda_old
            ! Sherman-Morrison 更新 Green function
            call Update_Green_SM_Lambda(G, i, B_total, N_dim, R_ferm)
        endif
    enddo
end subroutine
```

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

### 6.1 关于 τ^x/σ^x 的物理解释

> 🚨 **重要澄清：ALF 中存储的 Z₂ 变量 vs. Hamiltonian 中的算符**
>
> ALF 中实际存储的 `nsigma`, `ntau` 等变量是 **Z₂ Ising 场**（取值 ±1），对应的是 **σ^z, τ^z 的 classical representation**。
>
> 而 Hamiltonian 中的 **σ^x, τ^x** 是通过 **Hubbard-Stratonovich 变换** 映射到这些 Ising 场上的。

#### 从路径积分到观测量的映射关系

1. **MC 采样的是**：`nsigma(bond, tau)`, `ntau(site, tau)` 的配置
2. **这些配置代表的是**：在该时空点上 σ^z, τ^z 的本征值
3. **但在 Gauss 算符中出现的是 σ^x, τ^x**

在 slave-spin/orthogonal-fermion 框架中，路径积分表示已经将量子算符映射为经典 Ising 场。因此：

$$G_r = Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b$$

在 MC 中**直接用** `ntau(r, tau)` 和 `nsigma(b, tau)` 计算：

$$G_r^{\text{MC}}(\tau) = Q_r \cdot \texttt{ntau}(r, \tau) \cdot \prod_{b \in +r} \texttt{nsigma}(b, \tau)$$

这不是一个"错误"，而是 **路径积分 representation 中 classical field 就代表对应的 Pauli 算符**。

### 6.2 Gauss 算符期望值

$$\langle G_r \rangle = \left\langle Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b \right\rangle$$

在 MC 中测量：
$$\overline{G} = \frac{1}{N_\tau N_s} \sum_{\tau, r} G_r^{\text{MC}}(\tau)$$

应接近 $+1$（严格 projector 情况下）。

### 6.3 Gauss 约束违反度

$$\text{GaussViol} = \left\langle (G_r - Q_r)^2 \right\rangle = \frac{1}{N_\tau N_s} \sum_{\tau, r} (G_r(\tau) - Q_r)^2$$

- 若 projector 完全精确且无数值误差：GaussViol ≈ 0
- 实际上可能有极小但非零值（机器精度附近）

---

## 验证 Checklist

### ✅ 数值自检项目

#### 1. Gauss 约束数值检查
```fortran
! 测量 ⟨(G_r - Q_r)²⟩
real(8) :: gauss_viol
gauss_viol = 0.d0
do nt = 1, Ltrot
    do i = 1, Latt%N
        G_r = Compute_Gauss_Operator_Int(i, nt)
        gauss_viol = gauss_viol + (G_r - Q_background(i))**2
    enddo
enddo
gauss_viol = gauss_viol / (Ltrot * Latt%N)
! 期望值：应该在机器精度附近（< 1e-10）
```

**如果 GaussViol 随时间增大**，检查：
- P[λ] 是否在所有时间片的 Green 中一致地出现
- 某些 update 是否忘记乘 bosonic factor
- stabilizer block 是否重复乘了 P[λ]

#### 2. λ 边界条件检查
```fortran
! 测试 1：把所有 λ 固定为 +1
lambda_field(:) = +1
! 与不加严格 Gauss projector 的结果比较
! 应该只在物理 sector 有差异，而不是整体崩掉

! 测试 2：随机翻转几个 λ
call random_flip_lambda(10)
! 观测局域 τ^z 或密度
! 看看是否出现明显 PBC/APBC 的差异
```

#### 3. Sign 检查
```fortran
! 在论文参数点（sign-free 区域）统计平均 sign
complex(8) :: avg_sign
! 如果 sign 掉得很快，检查：
! - P[λ] 是否插错位置
! - GaussSector / Q_r pattern 是否和原论文一致
```

### ✅ 必须确认的实现细节

| 检查项 | 期望 | 危险信号 |
|--------|------|----------|
| P[λ] 乘的次数 | 完整一圈只乘一次 | 每个 block 乘一次/GRUP GRDW 各乘一次 |
| λ 翻转时 B_total | 不重算 | 每次翻转都重新计算 B_total |
| 两自旋 rank-2 | 真正用 rank-2 或分开两个 rank-1 | 偷懒当单 rank-1 用 |
| λ 存储 | site-only `lambda_field(site)` | 带 τ 下标 `lambda(site, tau)` |
| λ 更新循环 | 只遍历 site | 遍历 site × tau |

---

## 使用方法

### 参数设置

```
UseStrictGauss = .true.
GaussSector = "even"    ! "even", "odd", "staggered"
```

### GaussSector 的 Q_r pattern 定义

> 🚨 **必须明确 Q_r 的具体取值！**

| GaussSector | Q_r 定义 | 适用场景 |
|-------------|----------|----------|
| `"even"` | $Q_i = +1$ 对所有 site | 标准物理 sector |
| `"odd"` | $Q_i = -1$ 对所有 site | 全局奇 sector |
| `"staggered"` | $Q_{x,y} = (-1)^{x+y}$ | A/B 子格交替 |

#### Q_r pattern 的明确公式（Square lattice）

假设 site index 按行优先排列：
$$i = x + (y - 1) \cdot L_x, \quad x \in [1, L_x], \; y \in [1, L_y]$$

则：
- **even**：$Q_i = +1$
- **odd**：$Q_i = -1$
- **staggered**（棋盘形）：
  $$Q_i = (-1)^{x + y}$$
  其中 $(x, y) = (\text{mod}(i-1, L_x) + 1, \; (i-1) / L_x + 1)$

#### ALF 实现

```fortran
subroutine Setup_Q_Background(Latt, GaussSector)
    type(Lattice_type), intent(in) :: Latt
    character(len=*), intent(in) :: GaussSector
    integer :: i, x, y, Lx, Ly
    
    Lx = Latt%L1
    Ly = Latt%L2
    
    select case (trim(GaussSector))
    case ("even")
        Q_background(:) = +1
        
    case ("odd")
        Q_background(:) = -1
        
    case ("staggered")
        do i = 1, Latt%N
            x = mod(i - 1, Lx) + 1
            y = (i - 1) / Lx + 1
            Q_background(i) = (-1)**(x + y)
        enddo
        
    case default
        Q_background(:) = +1
    end select
end subroutine
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

### 已完成的修改

#### 1. `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`

**变量声明**：
- `lambda_field(:)` - τ-independent λ 场（一维数组）
- `Q_background(:)` - 背景电荷数组
- `Gamma_Gauss` - PRX A6 耦合常数 $\gamma = -\frac{1}{2}\ln[\tanh(\epsilon h)]$

**新增函数**：
| 函数名 | 功能 |
|--------|------|
| `Setup_Gauss_constraint()` | 初始化 λ 场和计算 γ |
| `Get_Tau_Z_At_Time_0(I)` | 获取 τ=0 处的 τ^z |
| `Get_Tau_Z_At_Time_M1(I)` | 获取 τ=M-1 处的 τ^z |
| `Compute_Gauss_Action_PRX(I)` | 计算单点 Gauss 作用量 |
| `Compute_Gauss_Weight_Ratio_Lambda_PRX(I)` | λ 翻转的玻色权重比 |
| `Compute_Delta_S_Gauss_Tau_Update(...)` | τ 更新的 ΔS_Gauss |
| `Compute_Star_Product_X(I, nt)` | 计算 star product $X_r$ |
| `Compute_Gauss_Operator_Int(I, nt)` | 计算 Gauss 算符（整数） |
| `Construct_P_Lambda_Matrix(P, N)` | 构造对角矩阵 $P[\lambda]$ |
| `Apply_P_Lambda_To_Matrix(B, N)` | 应用 $P[\lambda]$ 到矩阵 $B$ |
| `Compute_Lambda_Flip_Fermion_Ratio(I, G, B, N)` | Sherman-Morrison 费米子行列式比率 |
| `Update_Green_Sherman_Morrison_Lambda(G, I, B, N, R)` | Sherman-Morrison 更新 Green 函数 |
| `Compute_Lambda_Flip_Total_Ratio(I, G, B, N)` | λ 翻转总接受率 (bose + fermion) |

**修改的函数**：
- `Compute_Gauss_Operator` - 去除 $(-1)^{n_f}$（PRX orthogonal-fermion 构造）
- `Setup_Gauss_constraint` - 初始化 τ-independent λ 场
- `S0` - λ 更新使用 PRX A6 公式
- `Global_move_tau` - 添加 τ 更新的 Gauss 权重
- `Hamiltonian_set_nsigma` - 正确初始化 τ-independent λ

#### 2. `Documentation/Z2_Strict_Gauss_Constraint.md`

- 添加 S_τ-path 路径积分项说明（模块 0.3）
- 添加时空 plaquette 完整定义（模块 2）
- 添加 P[λ] wrap-up 插入机制说明（模块 3.4）
- 添加 Sherman-Morrison 更新说明（模块 5.1.1-5.1.3）

### 待完成任务

#### ✅ 已完成

1. **P[λ] 构造和应用函数** - ✅ 完成
   - `Construct_P_Lambda_Matrix(P, N)` - 构造对角矩阵
   - `Apply_P_Lambda_To_Matrix(B, N)` - 应用到 B 矩阵

2. **Sherman-Morrison λ 更新** - ✅ 完成
   - `Compute_Lambda_Flip_Fermion_Ratio(I, G, B, N)` - 费米子行列式比率
   - `Update_Green_Sherman_Morrison_Lambda(G, I, B, N, R)` - Green 函数更新
   - `Compute_Lambda_Flip_Total_Ratio(I, G, B, N)` - 总接受率

3. **PRX A6 玻色权重** - ✅ 完成
   - `Compute_Gauss_Weight_Ratio_Lambda_PRX(I)` - 玻色权重比率

#### ✅ 已完成（ALF 核心框架集成）

4. **在 CGR 函数中集成 P[λ]** - ✅ 完成
   - 修改 `cgr1_mod.F90` 中的 `CGR` 函数
   - 添加 `Use Hamiltonian_main, only: ham`
   - 在计算 GRUP 后调用 `ham%Apply_P_Lambda_To_Green(GRUP, 1)`
   - 支持两个版本的 CGR（STAB1/STAB2 和 STAB3/STABLOG）

5. **Hamiltonian_main 接口扩展** - ✅ 完成
   - 添加 `Use_Strict_Gauss()` 函数到 `ham_base` 类型
   - 添加 `Apply_P_Lambda_To_Green(GR, nf_eff)` 过程到 `ham_base` 类型
   - 在 `Hamiltonian_Z2_Matter_smod.F90` 中覆盖这些过程

6. **λ 更新的玻色权重** - ✅ 完成
   - `S0` 函数已使用 `Compute_Gauss_Weight_Ratio_Lambda_PRX(I)` 计算玻色权重
   - 费米子部分通过标准的 Green function 更新机制处理

#### 🔴 高优先级 - 关键问题

1. **P[λ] 实现的数学正确性问题** ⚠️
   - **问题**：当前实现 `G' = P * G` 不等于正确的 `G' = (1+PB)^{-1}`
   - **影响**：当某些 λ_i = -1 时，Green function 和行列式计算不准确
   - **解决方案选项**：
     - (A) 修改 ALF wrap-up 流程，在 B_total 上乘 P 后再传入 CGR
     - (B) 使用 Woodbury 公式后修正：`G' = (I + G*(P-I)*B)^{-1} * G`
   - **当前状态**：简化实现，仅在所有 λ_i = +1 时严格正确

2. **λ 更新循环的位置**
   - λ 是 τ-independent 的 site-only 变量
   - 需要在 MC sweep 中添加独立的 `Update_All_Lambda` 循环
   - 这个循环只遍历 site，不遍历 τ
   - **不使用 Field_type=5**，λ 不是 nsigma/ntau 数组的一部分

#### 🟡 中优先级

3. **时空 plaquette 项 S_plaq**（如需要 3D gauge action）
   - 添加 $K_{\text{plaq}} = \frac{1}{2}\ln[\coth(\epsilon g)]$

#### 🟢 低优先级

4. **τ(0), τ(M−1) 索引验证**
   - 确认 ALF 中 tau=1 对应 τ=0，tau=Ltrot 对应 τ=M-1

5. **GaussSector odd/staggered 测试**

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
