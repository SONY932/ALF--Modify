# 严格 Gauss 约束实现文档

## 概述

本文档描述了在 ALF (Algorithms for Lattice Fermions) 框架中实现的 Z₂ 规范场耦合费米子模型的**严格 Gauss 约束**投影。实现对标 PRX 10, 041057 和 Gazit et al. (2016) 的方法。

---

## 模块 1：Gauss 算符的定义与 Sector 选择

### 1.1 完整的 Gauss 算符形式

在 orthogonal-fermion/Z₂-Matter 模型中，Gauss 算符的完整定义为（参考 PRX 附录 A）：

$$G_r(\tau) = Q_r \cdot (-1)^{n^f_r(\tau)} \cdot \tau_r^x(\tau) \cdot \prod_{b\in +r} \sigma^x_b(\tau)$$

其中：
- $n^f_r(\tau)$：格点 $r$ 在虚时间 $\tau$ 处的费米子占据数
- $\tau_r^x(\tau)$：格点 $r$ 上的 $\tau$ 自旋的 $x$ 分量（Ising 变量，取值 $\pm 1$）
- $\sigma^x_b(\tau)$：连接格点 $r$ 的 bond $b$ 上的规范场变量
- $\prod_{b\in +r}$：对以 $r$ 为中心的 star（所有相邻 bond）求积
- $Q_r = \pm 1$：背景电荷，决定 Gauss sector

### 1.2 Gauss Sector 选择

物理态满足局域 Gauss 约束：
$$G_r(\tau) = +1, \quad \forall r, \tau$$

投影算符定义为：
$$\hat{P}_r = \frac{1}{2}(1 + G_r), \qquad \hat{P} = \prod_r \hat{P}_r$$

**Sector 参数**：
- `GaussSector = "even"`：对应 $Q_r = +1$（所有格点）
- `GaussSector = "odd"`：对应 $Q_r = -1$（所有格点）
- `GaussSector = "staggered"`：对应交错 pattern（如子格 A 取 $+1$，子格 B 取 $-1$）

### 1.3 ALF 实现要点

1. **Gauss 算符计算**需要包含所有组成部分：
   - 费米子占据数奇偶性 $(-1)^{n^f_r}$
   - τ 自旋 $\tau_r^x$
   - Star product $\prod_{b\in +r} \sigma^x_b$
   - 背景电荷 $Q_r$

2. **参数层接口**：
   ```fortran
   ! Sector 选择参数
   Character(len=64) :: GaussSector = "even"  ! "even", "odd", "staggered"
   
   ! 背景电荷数组
   Integer, allocatable :: Q_background(:)  ! Q_r for each site
   ```

---

## 模块 2：投影算符在路径积分中的插入方式

### 2.1 起点：带投影的配分函数

配分函数从投影到物理 Hilbert 空间出发：
$$Z = \mathrm{Tr}\Big[\hat{P} \, e^{-\beta H}\Big]$$

其中 $\hat{P} = \prod_r \hat{P}_r$，$\hat{P}_r = \frac{1}{2}(1 + G_r)$。

### 2.2 Trotter 分解

设 $\beta = M \Delta\tau$，则：
$$e^{-\beta H} \approx \Big(e^{-\Delta\tau K} e^{-\Delta\tau V}\Big)^M$$

### 2.3 投影算符的均匀分布

将投影算符均匀分散到每个时间片（PRX 附录中明确采用此做法）：
$$Z \approx \mathrm{Tr}\left[\prod_{\tau=1}^{M}\left(\prod_r \hat{P}_r(\tau) \cdot e^{-\Delta\tau K} e^{-\Delta\tau V(\tau)}\right)\right]$$

### 2.4 路径积分形式

引入 $\lambda$ 场后，配分函数写成：
$$Z \approx \sum_{\{\sigma, \tau, \lambda\}} W_{\text{Bose}}[\sigma, \tau, \lambda] \cdot \det M[\sigma, \tau, \lambda]$$

其中玻色权重为：
$$W_{\text{Bose}}[\sigma, \tau, \lambda] = \prod_{\tau, r} W_r\bigl(\lambda_r(\tau), G_r(\tau)\bigr)$$

**关键点**：Gauss 投影贡献的是**局域权重因子** $W_r$，而不是简单地乘一个相位因子。$\lambda$ 场同时进入玻色 action 和费米子行列式。

---

## 模块 3：投影算符的 λ 展开

### 3.1 严格等价的 λ 展开公式

投影算符 $\hat{P}_r = \frac{1}{2}(1 + G_r)$ 可以用 Z₂ 辅助场 $\lambda_r = \pm 1$ 展开：

$$\frac{1 + G_r}{2} = \frac{1}{4} \sum_{\lambda_r = \pm 1} (1 + \lambda_r)(1 + \lambda_r G_r)$$

**验证**：
- 若 $G_r = +1$：
  - $\lambda_r = +1$：$(1+1)(1+1)/4 = 1$
  - $\lambda_r = -1$：$(1-1)(1-1)/4 = 0$
  - 求和 = 1 ✓
- 若 $G_r = -1$：
  - $\lambda_r = +1$：$(1+1)(1-1)/4 = 0$
  - $\lambda_r = -1$：$(1-1)(1+1)/4 = 0$
  - 求和 = 0 ✓

### 3.2 局域 Gauss 权重

定义局域权重函数：
$$W_r\bigl(\lambda_r(\tau), G_r(\tau)\bigr) \equiv \frac{1}{4}\bigl(1 + \lambda_r(\tau)\bigr)\bigl(1 + \lambda_r(\tau) G_r(\tau)\bigr)$$

### 3.3 整体投影因子

整个虚时间轨迹的投影因子：
$$W_{\text{Gauss}}[\lambda, G] = \prod_{\tau, r} W_r\bigl(\lambda_r(\tau), G_r(\tau)\bigr)$$

### 3.4 作用量形式

对应的 Gauss 作用量：
$$S_{\text{Gauss}}[\lambda, G] = -\sum_{\tau, r} \ln W_r\bigl(\lambda_r(\tau), G_r(\tau)\bigr)$$

### 3.5 关键性质

| $G_r$ | $\lambda_r$ | $W_r$ | 物理含义 |
|-------|-------------|-------|----------|
| $+1$ | $+1$ | $1$ | 满足 Gauss 约束，有限权重 |
| $+1$ | $-1$ | $0$ | 满足约束但 $\lambda$ 选择导致零权重 |
| $-1$ | $+1$ | $0$ | 违反约束，被投影杀掉 |
| $-1$ | $-1$ | $0$ | 违反约束，被投影杀掉 |

**核心要点**：
- 对 $G_r = +1$ 的配置，只有 $\lambda_r = +1$ 有非零权重
- 对 $G_r = -1$ 的配置，无论 $\lambda_r$ 取何值，权重均为 0
- 这实现了对违反 Gauss law 轨迹的**严格排除**

---

## 模块 4：费米子传播子中 λ 的作用方式

### 4.1 对角矩阵 P[λ]

根据 PRX 附录 A，$\lambda$ 场通过一个对角矩阵修改费米子 Green 函数。在单粒子 Hilbert 空间中定义：
$$P_{ij}(\tau) = \lambda_i(\tau) \cdot \delta_{ij}$$

即 $P(\tau)$ 是一个对角矩阵，对角元为各格点的 $\lambda_r(\tau)$。

### 4.2 修改后的 B-slice

原始的单时间片传播子：
$$B_0(\tau) = e^{-\Delta\tau K} \cdot e^{-\Delta\tau V(\sigma(\tau), \tau(\tau))}$$

加入 $\lambda$ 场后：
$$B(\tau) = P(\tau) \cdot B_0(\tau)$$

或等价地（因 $P$ 是对角矩阵，与某些项可交换）：
$$B(\tau) = P(\tau) \cdot e^{-\Delta\tau K} \cdot e^{-\Delta\tau V(\sigma(\tau), \tau(\tau))}$$

### 4.3 完整传播子

整条虚时间传播子：
$$\mathcal{B} = B(M) B(M-1) \cdots B(1)$$

费米子行列式：
$$\det M[\sigma, \tau, \lambda] = \det(1 + \mathcal{B})$$

$\lambda$ 场通过 $P(\tau)$ 矩阵进入行列式，影响费米子权重。

### 4.4 ALF 实现要点

```fortran
! 在 Op_V 或相关函数中，构造 P(tau) 对角矩阵
! P_diag(i) = lambda_field(i, nt)

! 修改 B-slice 的计算
! B(tau) = P(tau) * exp(-dtau*K) * exp(-dtau*V)
```

---

## 模块 5：蒙特卡洛更新的接受率

### 5.1 总接受率结构

任意场更新的接受率为：
$$R_{\text{tot}} = R_{\text{Gauss}} \cdot R_{\text{fermion}}$$

其中：
- $R_{\text{Gauss}}$：来自局域 Gauss 权重 $W_r$ 的贡献
- $R_{\text{fermion}}$：来自费米子行列式的贡献

### 5.2 翻转 λ 场：$\lambda_r(\tau) \to -\lambda_r(\tau)$

**受影响区域**：仅单点 $(r, \tau)$

**Gauss 权重比率**：
$$R_{\text{Gauss}}^{(\lambda)} = \frac{W_r(\lambda_r^{\text{new}}, G_r)}{W_r(\lambda_r^{\text{old}}, G_r)}$$

注意：
- 若 $G_r = -1$，则 $W_r^{\text{new}} = W_r^{\text{old}} = 0$，对应投影完全杀掉此配置
- 若 $G_r = +1$ 且 $\lambda_r^{\text{old}} = +1$，则翻转后 $\lambda_r^{\text{new}} = -1$，$W_r^{\text{new}} = 0$，更新被拒绝

**费米子行列式比率**：
$$R_{\text{fermion}}^{(\lambda)} = \frac{\det(1 + \mathcal{B}^{\text{new}})}{\det(1 + \mathcal{B}^{\text{old}})}$$

差别来自 $P(\tau)$ 在第 $r$ 个对角元的 $\pm 1$ 翻转。

**总比率**：
$$R_{\text{tot}}^{(\lambda)} = R_{\text{Gauss}}^{(\lambda)} \cdot R_{\text{fermion}}^{(\lambda)}$$

### 5.3 翻转 gauge link：$\sigma_b^x(\tau) \to -\sigma_b^x(\tau)$

**受影响区域**：link $b$ 连接的两个格点 $r_1, r_2$ 的 Gauss 算符

**Gauss 权重比率**：
$$R_{\text{Gauss}}^{(\sigma)} = \prod_{r \in \{r_1, r_2\}} \frac{W_r(\lambda_r, G_r^{\text{new}})}{W_r(\lambda_r, G_r^{\text{old}})}$$

关键点：
- 翻转 $\sigma_b^x$ 会改变 $G_{r_1}$ 和 $G_{r_2}$ 的 star product 部分
- 若翻转后某个 $G_r^{\text{new}} = -1$，则 $W_r = 0$，整个更新被拒绝

**费米子行列式比率**：
$$R_{\text{fermion}}^{(\sigma)} = \frac{\det(1 + \mathcal{B}^{\text{new}})}{\det(1 + \mathcal{B}^{\text{old}})}$$

这与标准 DQMC 相同，但使用的是 $\lambda$ 修正后的 B-slice。

### 5.4 翻转 τ 自旋：$\tau_r^x(\tau) \to -\tau_r^x(\tau)$

**受影响区域**：单点 $(r, \tau)$ 的 Gauss 算符

**Gauss 算符变化**：
$$G_r^{\text{new}} = -G_r^{\text{old}}$$

因为 $\tau_r^x$ 是 Gauss 算符的直接组成部分。

**Gauss 权重比率**：
$$R_{\text{Gauss}}^{(\tau)} = \frac{W_r(\lambda_r, G_r^{\text{new}})}{W_r(\lambda_r, G_r^{\text{old}})}$$

关键点：
- 若 $G_r^{\text{old}} = +1$，翻转后 $G_r^{\text{new}} = -1$，则 $W_r^{\text{new}} = 0$，更新被**严格禁止**
- 这保证了满足 Gauss 约束的配置不会演化到违反约束的配置

**费米子行列式比率**：
$$R_{\text{fermion}}^{(\tau)} = \frac{\det(1 + \mathcal{B}^{\text{new}})}{\det(1 + \mathcal{B}^{\text{old}})}$$

### 5.5 更新策略总结

| 更新类型 | 受影响的 $G_r$ | $R_{\text{Gauss}}$ | 约束强制机制 |
|----------|----------------|---------------------|--------------|
| $\lambda_r(\tau)$ 翻转 | 无（$G_r$ 不变） | $W_r^{\text{new}}/W_r^{\text{old}}$ | 若 $G_r = -1$ 则权重为 0 |
| $\sigma_b^x(\tau)$ 翻转 | $G_{r_1}, G_{r_2}$ | $\prod_{r} W_r^{\text{new}}/W_r^{\text{old}}$ | 若任一 $G_r^{\text{new}} = -1$ 则拒绝 |
| $\tau_r^x(\tau)$ 翻转 | $G_r$ 变号 | $W_r^{\text{new}}/W_r^{\text{old}}$ | 若 $G_r^{\text{new}} = -1$ 则严格拒绝 |

---

## 模块 6：观测量定义

### 6.1 局域 Gauss 期望值

定义时空平均：
$$\overline{G} \equiv \frac{1}{N_\tau N_s} \sum_{\tau, r} G_r(\tau)$$

对于 even sector（$Q_r = +1$），蒙特卡洛平均应满足：
$$\langle \overline{G} \rangle \approx +1$$

对于非平凡 sector，定义：
$$\tilde{G}_r = Q_r \cdot G_r$$

则仍期望 $\langle \tilde{G}_r \rangle \simeq 1$。

### 6.2 Gauss 约束违反度

**标准定义**（even sector）：
$$\overline{\Delta G^2} \equiv \frac{1}{N_\tau N_s} \sum_{\tau, r} (G_r(\tau) - 1)^2$$

**一般 sector**：
$$\overline{\Delta G^2} \equiv \frac{1}{N_\tau N_s} \sum_{\tau, r} (G_r(\tau) - Q_r)^2$$

**期望值**：
- 若投影完全精确且无数值误差：$\langle \overline{\Delta G^2} \rangle = 0$
- 实际模拟中可能有极小非零值，用于数值自检

### 6.3 ALF 输出文件

| 文件名 | 内容 | 预期值 |
|--------|------|--------|
| `Gauss_scal.dat` | $\langle \overline{G} \rangle$ | $\approx Q_r$（通常为 $+1$） |
| `GaussViol_scal.dat` | $\langle \overline{\Delta G^2} \rangle$ | $\approx 0$ |

---

## 使用方法

### 参数设置

在参数文件中：
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
! Lambda 场数组：lambda_field(site, tau) = +1 或 -1
Integer, allocatable :: lambda_field(:,:)

! 背景电荷数组
Integer, allocatable :: Q_background(:)

! Star product 缓存
Integer, allocatable :: star_product_cache(:,:)
```

### 场类型

Lambda 场为第 5 种场类型（Field_type = 5）：
- Field_type = 1: Z₂ 规范场 $\sigma$
- Field_type = 2: Bond matter 场
- Field_type = 3: Hubbard HS 场
- Field_type = 4: Site matter 场 $\tau$
- **Field_type = 5: Gauss lambda 场** (新增)

### 核心函数

| 函数名 | 功能 | 公式 |
|--------|------|------|
| `Compute_Star_Product_X(I, nt)` | 计算 star product | $X_r = \prod_{b \in +r} \sigma_b^x$ |
| `Compute_Gauss_Operator_Int(I, nt)` | 计算整数 Gauss 算符（用于 MC 更新） | $G_r^{\text{bose}} = Q_r \cdot \tau_r^x \cdot X_r$ |
| `Compute_Gauss_Operator(I, nt, GRC)` | 计算完整 Gauss 算符（用于观测量） | $G_r = Q_r \cdot (-1)^{n_r^f} \cdot \tau_r^x \cdot X_r$ |
| `Compute_Gauss_Weight(lambda, G_r)` | 计算局域权重 | $W_r = \frac{1}{4}(1+\lambda)(1+\lambda G_r)$ |
| `Compute_Gauss_Weight_Ratio(...)` | 计算权重比率 | $R = W_r^{\text{new}} / W_r^{\text{old}}$ |

**注**：P[λ] 对角矩阵通过 `Ham_V` 中 Field_type=5 的 `Op_V` 设置实现，详见代码注释。

---

## 与文献的对应关系

| 本文档内容 | 对应文献 |
|------------|----------|
| Gauss 算符定义 | PRX 10, 041057 附录 A |
| λ 展开公式 | Gazit et al. (2016) |
| P[λ] 对角矩阵 | PRX 附录 A 最后一段 |
| 接受率公式 | 标准 DQMC + 投影权重 |

---

## 验证方法

### 1. 检查 Gauss 约束满足程度

运行模拟后检查：
- `Gauss_scal.dat`：应接近 $Q_r$（通常为 $+1$）
- `GaussViol_scal.dat`：应接近 $0$（典型值 $< 10^{-10}$）

### 2. 配置有效性检查

在每次更新后可选地验证：
$$\forall (r, \tau): \quad G_r(\tau) = Q_r$$

若任何配置违反此条件，说明实现有误。

### 3. 对比测试

比较以下两种情况：
1. `UseStrictGauss = .false.`（原有实现）
2. `UseStrictGauss = .true.`（严格 Gauss 约束）

检查相图结构、临界指数、费米面结构的差异。

---

## 注意事项

1. **权重为零的处理**：当 $W_r = 0$ 时，$\ln W_r$ 发散。实际实现中直接拒绝使权重为零的更新，无需计算对数。

2. **数值稳定性**：对角矩阵 $P(\tau)$ 的元素为 $\pm 1$，不会引入额外的数值不稳定性。

3. **符号问题**：在论文参数下应是 sign-free，但在新参数区域需要验证。

4. **初始化**：初始配置必须满足 $G_r = Q_r$，否则权重为零。建议从满足约束的配置开始（如 $\lambda_r = +1$，所有 $G_r = +1$）。

5. **性能**：$\lambda$ 场更新与其他 Ising 场更新开销相当。

---

## 文件修改列表

主要修改的文件：
- `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`
  - 添加 `UseStrictGauss`、`GaussSector` 参数
  - 添加 `lambda_field`、`Q_background`、`DW_Gauss_weight` 存储
  - 添加 `Compute_Gauss_Operator_Int`、`Compute_Gauss_Weight`、`Compute_Gauss_Weight_Ratio` 函数
  - 修改 `Compute_Gauss_Operator` 加入 τ 自旋和 Q_r
  - 修改 `Setup_Gauss_constraint` 初始化背景电荷和权重表
  - 修改 `Ham_V` 中 Field_type=5 处理 P[λ] 对角算符
  - 修改 `S0` 函数加入 λ 翻转和 σ 翻转的 Gauss 权重比率
  - 修改 `Global_move_tau` 加入 τ 翻转的 Gauss 权重比率
  - 修改 `Obser` 测量 Gauss 和 GaussViol 观测量

---

## 作者

ALF Collaboration
