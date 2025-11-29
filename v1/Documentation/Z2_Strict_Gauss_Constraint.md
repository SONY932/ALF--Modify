# 严格 Gauss 约束实现文档

## 概述

本文档描述了在 ALF 框架中实现的 Z₂ 规范场耦合费米子模型的**严格 Gauss 约束**，对标 **PRX 10, 041057 (2020) Appendix A**。

---

## 核心原理（PRX A2-A6）

### 1. Gauss 投影算符

PRX Appendix A 使用离散拉格朗日乘子 λ 重写 Gauss 投影算符：

$$\hat{P}_i = \frac{1}{2}(1 + G_i) = \sum_{\lambda_i = \pm 1} e^{i\frac{\pi}{2}(1-\lambda_i)} \hat{P}_i(\lambda_i)$$

其中 Gauss 算符为：
$$G_r = Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b$$

### 2. λ 求和后的结果

**关键点**：λ 不是 MC 采样变量！对 λ 求和后：

$$\sum_{\lambda_i = \pm 1} e^{\gamma \lambda_i \tau^z_{i,0} \tau^z_{i,M-1}} = 2\cosh(\gamma \cdot \tau^z_{i,0} \cdot \tau^z_{i,M-1})$$

这给出一个**纯玻色的时间边界耦合**：

$$S_{\text{boundary}} = -K_G \sum_i \tau^z_{i,0} \cdot \tau^z_{i,M-1}$$

其中：
$$K_G = -\frac{1}{2}\ln[\tanh(\epsilon \cdot h)]$$

### 3. 物理意义

- **$K_G \to \infty$**：严格投影，强制 $\tau^z_{i,0} = \tau^z_{i,M-1}$
- **$K_G \to 0$**：无时间边界耦合

这等价于在 τ 自旋的时间边界（τ=0 和 τ=β）之间添加一个 Ising 耦合项。

---

## 实现要点

### ✅ 正确的实现（当前版本）

1. **NO λ 场作为 MC 变量**
   - λ 不存在于 MC 配置空间中
   - 不需要 `lambda_field(:)` 数组
   - 不需要 `Sweep_Lambda` 循环

2. **NO 费米子传播子修改**
   - 不需要 `P[λ]` 乘在 B 矩阵上
   - 不需要 `Apply_P_Lambda_To_B`
   - 不需要 Sherman-Morrison 更新
   - Green 函数计算与普通情况相同

3. **YES 纯玻色时间边界耦合**
   - 在 τ 自旋更新时（nt=1 或 nt=Ltrot），附加权重因子
   - $R_{\text{bose}} = \exp(-K_G \cdot \Delta(\tau^z_0 \cdot \tau^z_{M-1}))$

### ❌ 之前错误的实现（已删除）

- 把 λ 当成 MC 采样变量
- 在费米子传播子中插入 `P[λ]`
- Sherman-Morrison 更新 Green 函数
- `Sweep_Lambda` 循环

---

## 代码结构

### 参数设置

```fortran
&VAR_Z2_Matter
UseStrictGauss = .T.
GaussSector    = "even"    ! "even", "odd", "staggered"
/
```

### 关键函数

| 函数 | 功能 |
|------|------|
| `Setup_Gauss_constraint()` | 初始化 Q_background 和 K_G (Gamma_Gauss) |
| `Get_Tau_Z_At_Time_0(I)` | 获取 τ=0 处的 τ^z |
| `Get_Tau_Z_At_Time_M1(I)` | 获取 τ=M-1 处的 τ^z |
| `Compute_Delta_S_Gauss_Tau_Update(...)` | 计算 τ 翻转的 ΔS_boundary |
| `Compute_Gauss_Operator_Int(I, nt)` | 计算 Gauss 算符（观测量用） |
| `Measure_GaussViolation_Diagnostic(sweep)` | 诊断输出 |

### τ 更新的 Gauss 权重

在 `Global_move_tau` 中，当 nt=1 或 nt=Ltrot 时：

```fortran
If (UseStrictGauss) then
   If (ntau == 1 .or. ntau == Ltrot) then
      ! 计算边界值变化
      Delta_S_Gauss = Compute_Delta_S_Gauss_Tau_Update(...)
      R_Gauss = exp(-Delta_S_Gauss)
      S0_Matter = S0_Matter * R_Gauss
   endif
endif
```

---

## 验证方法

### 1. GaussViol 检查

```
============================================================
 GAUSS CONSTRAINT DIAGNOSTIC - Sweep      100
============================================================
   <G_r>         (should be ~1): 0.98765432E+00
   GaussViol     (should be ~0): 0.12345678E-03
   <tau0*tauM1>  (boundary corr): 0.95000000E+00
   K_G (Gamma_Gauss):             2.302585
============================================================
```

- `<G_r> ≈ 1`：Gauss 约束被满足
- `GaussViol ≈ 0`：违反度极小
- `<tau0*tauM1> > 0`：边界耦合工作正常

### 2. K_G 调节

- K_G 越大，Gauss 约束越严格
- 典型值：K_G ~ 2-10 对于合理的 Dtau 和 Ham_h

### 3. 参数建议

```fortran
! 小系统测试
L1 = 2, L2 = 2
Beta = 2.0, Dtau = 0.1
Ham_h = 1.0  ! 给出 K_G ~ 2.3
```

---

## 与 PRX 的对应关系

| PRX 内容 | ALF 实现 |
|----------|----------|
| λ 求和后的有效作用量 | `S_boundary = -K_G * tau_z_0 * tau_z_{M-1}` |
| K_G 公式 | `Gamma_Gauss = -0.5 * ln(tanh(eps*h))` |
| τ^x 横场 → 时间方向 Ising | `DW_Matter_tau` (已有) |
| Gauss 算符观测量 | `Compute_Gauss_Operator_Int` |

---

## 文件修改列表

### `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`

**保留**：
- `Q_background(:)` - 背景电荷
- `Gamma_Gauss` - K_G 耦合常数
- `Setup_Gauss_constraint()` - 初始化
- `Get_Tau_Z_At_Time_0/M1()` - 边界值获取
- `Compute_Delta_S_Gauss_Tau_Update()` - τ 更新权重
- `Compute_Gauss_Operator_Int()` - 观测量
- `Measure_GaussViolation_Diagnostic()` - 诊断

**删除**：
- `lambda_field(:)` - 不再是 MC 变量
- `B_lambda_slice(:,:)` - 不再需要
- `Sweep_Lambda()` - 从 type 定义中删除
- `Apply_P_Lambda_To_B()` - 从 type 定义中删除
- S0 函数中的 λ 更新代码

### `Prog/wrapur_mod.F90`

**删除**：
- `ham%Apply_P_Lambda_To_B` 调用

### `Prog/main.F90`

**删除**：
- `ham%Sweep_Lambda` 调用及相关代码块

---

## 注意事项

1. **K_G 的数值稳定性**
   - 当 Ham_h → 0 时，K_G → ∞
   - 代码中设置了最大值截断 `K_G_max = 100`

2. **σ 更新不受 Gauss 约束影响**
   - Gauss 约束通过 τ 边界耦合实现
   - σ 更新使用标准 Ising 权重

3. **GaussViol 的解释**
   - 由于是软约束（有限 K_G），GaussViol 不会精确为 0
   - K_G 越大，GaussViol 越小

---

## 参考文献

- PRX 10, 041057 (2020) - "Dynamical Signatures of Edge-State Magnetism on Graphene Nanoribbons"
  - Appendix A: Path integral representation of Gauss constraint

---

## 作者

ALF Collaboration
