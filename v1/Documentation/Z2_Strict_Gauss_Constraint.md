# 严格 Gauss 约束实现文档

## 概述

本文档描述了在 ALF 框架中实现的 Z₂ 规范场耦合费米子模型的**严格 Gauss 约束**，对标 **PRX 10, 041057 (2020) Appendix A** 和 **PNAS 115, E6987 (2018)**。

---

## 核心原理

### 1. Gauss 算符定义

对每个格点 r：
$$G_r = Q_r \cdot \tau_r^x \cdot \prod_{b \in +r} \sigma^x_b$$

严格 Gauss 约束要求：
$$G_r = +1 \quad \forall r, \forall \tau$$

其中 $Q_r = \pm 1$ 是背景电荷。

### 2. 完整的严格 Gauss 约束 = 两个部分

**严格 Gauss 约束需要同时实现两个作用量项：**

#### Part 1: τ 时间边界耦合 (PRX A5-A6)

对 λ 求和后得到：
$$S_{\tau}^{\rm boundary} = -K_G \sum_r \tau^z_{r,0} \cdot \tau^z_{r,M-1}$$

其中：
$$K_G = -\frac{1}{2}\ln[\tanh(\epsilon \cdot h_\tau)]$$

- $h_\tau$ = `Ham_h` (τ 自旋的横场)
- $\epsilon$ = `Dtau` (虚时间步长)

#### Part 2: σ 星乘积时间耦合 (电场约束)

$$S_{\rm time}^{\rm plaq} = -K_{\rm time} \sum_{r,n} X_r(n) \cdot X_r(n+1)$$

其中：
- $X_r(n) = \prod_{b \in +r} \sigma^z_{b,n}$ 是时间片 n 上格点 r 处的星乘积
- $K_{\rm time} = -\frac{1}{2}\ln[\tanh(\epsilon \cdot h_\sigma)] = \frac{1}{2}\ln[\coth(\epsilon \cdot g)]$
- $h_\sigma$ = `Ham_g` (σ 自旋的横场)

### 3. 物理意义

**τ 边界耦合**：强制 $\tau^z_{r,0} = \tau^z_{r,M-1}$
- 保证 τ 在虚时间方向的周期性/反周期性

**σ 星乘积耦合**：强制 $X_r(n) = X_r(n+1)$ 对所有 n
- 保证星乘积（电场项）在虚时间方向一致
- 这是 Gauss 算符中 $\prod \sigma^x_b$ 部分的约束

**两者结合**：完整的 Gauss 约束
$$G_r(\tau) = Q_r \quad \forall r, \forall \tau$$

---

## 实现要点

### ✅ 当前完整实现

1. **τ 时间边界耦合** (Part 1)
   - 在 `Global_move_tau` 中
   - 当 nt=1 或 nt=Ltrot 时，附加权重 `R = exp(-Delta_S_tau)`
   - `Delta_S_tau = -K_G * (tau0_new * tauM1_new - tau0_old * tauM1_old)`

2. **σ 星乘积时间耦合** (Part 2)
   - 在 `S0` 函数中，当 field_type=1 (gauge field) 时
   - 附加权重 `R = exp(-Compute_Delta_S_Star_Time(n, nt))`
   - 影响星乘积两端点（link 两端格点）的时间耦合

3. **NO λ 场作为 MC 变量**
   - λ 已被求和消除
   - 不需要 `Sweep_Lambda`

4. **NO 费米子传播子修改**
   - 不需要 `P[λ]` 乘在 B 矩阵上
   - Green 函数计算与普通情况相同

---

## 代码结构

### 参数设置

```fortran
&VAR_Z2_Matter
UseStrictGauss = .T.
GaussSector    = "even"    ! "even", "odd", "staggered"
/
```

### 关键变量

| 变量 | 含义 |
|------|------|
| `Gamma_Gauss` | K_G (τ 边界耦合强度) |
| `Gamma_Gauss_Sigma` | K_time (σ 星乘积耦合强度) |
| `Q_background(:)` | 背景电荷 Q_r |

### 关键函数

| 函数 | 功能 |
|------|------|
| `Setup_Gauss_constraint()` | 初始化 Q_r, K_G, K_time |
| `Compute_Delta_S_Gauss_Tau_Update(...)` | 计算 τ 翻转的 ΔS_tau |
| `Compute_Delta_S_Star_Time(n, nt)` | 计算 σ 翻转的 ΔS_star |
| `Compute_Star_Product_X(I, nt)` | 计算星乘积 X_r |
| `Compute_Gauss_Operator_Int(I, nt)` | 计算 Gauss 算符（观测量） |
| `Measure_GaussViolation_Diagnostic(sweep)` | 诊断输出 |

### τ 更新的 Gauss 权重

在 `Global_move_tau` 中，当 nt=1 或 nt=Ltrot 时：

```fortran
If (UseStrictGauss) then
   If (ntau == 1 .or. ntau == Ltrot) then
      Delta_S_Gauss = Compute_Delta_S_Gauss_Tau_Update(...)
      R_Gauss = exp(-Delta_S_Gauss)
      S0_Matter = S0_Matter * R_Gauss
   endif
endif
```

### σ 更新的 Gauss 权重

在 `S0` 函数中，当 field_type=1 (gauge field) 时：

```fortran
If (UseStrictGauss) then
   R_Gauss = exp(-Compute_Delta_S_Star_Time(n, nt))
   S0 = S0 * R_Gauss
endif
```

---

## 验证方法

### 1. GaussViol 检查

```
============================================================
 GAUSS CONSTRAINT DIAGNOSTIC - Sweep      100
============================================================
   <G_r>         (should be ~1): 0.99999999E+00
   GaussViol     (should be ~0): 0.12345678E-08
   <tau0*tauM1>  (boundary corr): 0.99999000E+00
   K_G (Gamma_Gauss):             2.302585
============================================================
```

- `<G_r> ≈ 1`：Gauss 约束被严格满足
- `GaussViol ≈ 0`：违反度极小（< 10⁻⁶ 表示成功）
- `<tau0*tauM1> ≈ 1`：τ 边界耦合工作正常

### 2. 参数建议

```fortran
! 小系统测试
L1 = 2, L2 = 2
Beta = 2.0, Dtau = 0.1
Ham_h = 1.0  ! 给出 K_G ~ 2.3
Ham_g = 1.0  ! 给出 K_time ~ 2.3
```

### 3. 严格约束的标志

当 K_G 和 K_time 都足够大时（> 2），应该观察到：
- `GaussViol < 10⁻⁶`
- `<G_r> ≈ 1.0` (误差 < 10⁻⁴)
- 配置在物理 Hilbert 子空间内

---

## 与 PRX/PNAS 的对应关系

| 文献内容 | ALF 实现 |
|----------|----------|
| λ 求和后的 τ 边界耦合 | `Gamma_Gauss`, `Compute_Delta_S_Gauss_Tau_Update` |
| Gauss 投影 → σ 时间一致性 | `Gamma_Gauss_Sigma`, `Compute_Delta_S_Star_Time` |
| τ^x 横场 → 时间方向 Ising | `DW_Matter_tau` (已有) |
| σ^x 横场 → 时间方向 Ising | `DW_Ising_tau` (已有) |
| Gauss 算符观测量 | `Compute_Gauss_Operator_Int` |
| 完整作用量 | S = S_tau + S_star |

---

## 文件修改列表

### `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`

**保留/添加**：
- `Q_background(:)` - 背景电荷
- `Gamma_Gauss` - K_G (τ 边界耦合)
- `Gamma_Gauss_Sigma` - K_time (σ 星乘积耦合) **[新增]**
- `Setup_Gauss_constraint()` - 初始化（计算 K_G 和 K_time）
- `Compute_Delta_S_Gauss_Tau_Update()` - τ 更新权重
- `Compute_Delta_S_Star_Time()` - σ 更新权重 **[新增]**
- `Compute_Star_Product_X()` - 星乘积计算
- `Compute_Gauss_Operator_Int()` - 观测量
- `Measure_GaussViolation_Diagnostic()` - 诊断

**删除**：
- `lambda_field(:)` - 不再是 MC 变量
- `Sweep_Lambda()` - 不再需要
- `Apply_P_Lambda_To_B()` - 不再需要

### `Prog/wrapur_mod.F90`

**删除**：
- `ham%Apply_P_Lambda_To_B` 调用

### `Prog/main.F90`

**删除**：
- `ham%Sweep_Lambda` 调用

---

## 注意事项

1. **K_G 和 K_time 的数值稳定性**
   - 当 Ham_h → 0 或 Ham_g → 0 时，K → ∞
   - 代码中设置了最大值截断 `K_max = 100`

2. **两个耦合都是必须的**
   - 只有 K_G（τ 边界）：不是严格 Gauss，σ 不受约束
   - 只有 K_time（σ 星乘积）：不完整
   - 两者都有：完整的严格 Gauss 约束

3. **性能考虑**
   - `Compute_Delta_S_Star_Time` 在每次 σ 更新时被调用
   - 它计算 6 个星乘积值，有一定开销
   - 对于大系统可能需要优化

---

## 参考文献

- PRX 10, 041057 (2020) - "Dynamical Signatures of Edge-State Magnetism on Graphene Nanoribbons"
  - Appendix A: Path integral representation of Gauss constraint
- PNAS 115, E6987 (2018) - "Monte Carlo studies of the Z₂ gauge-Higgs model"
  - Gauss law enforcement methods

---

## 作者

ALF Collaboration
