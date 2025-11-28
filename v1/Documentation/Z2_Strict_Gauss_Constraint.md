# 严格 Gauss 约束实现文档

## 概述

本文档描述了在 ALF (Algorithms for Lattice Fermions) 框架中实现的 Z₂ 规范场耦合费米子模型的**严格 Gauss 约束**投影。

### 物理背景

Gauss 约束的核心要求是：

$$G_r(\tau) \equiv (-1)^{n_r(\tau)} \prod_{b\in +r} \sigma^x_b(\tau) = 1,\quad \forall r,\tau$$

通过在虚时间方向引入 Z₂ 拉格朗日乘子场 $\lambda_r(\tau)=\pm1$，可以实现对物理 Hilbert 空间的显式投影。

投影算符可以写成：
$$P_r = \frac{1}{2}(1 + G_r) \propto \sum_{\lambda_r=\pm1} \exp\left[ i\lambda_r(\tau)\frac{\pi}{2}\,G_r(\tau)\right]$$

## 使用方法

### 参数设置

在参数文件中设置以下参数来启用严格 Gauss 约束：

```
UseStrictGauss = .true.
```

当 `UseStrictGauss = .false.`（默认值）时，代码恢复到原有的 Z₂_Matter 模型行为。

### 新增观测量

当启用严格 Gauss 约束时，以下额外观测量会被自动添加：

1. **Gauss**: Gauss 算符的期望值 $\langle G_r \rangle$
   - 如果 Gauss 约束被严格满足，此值应该接近 1

2. **GaussViol**: Gauss 约束违反 $\langle (G_r - 1)^2 \rangle$
   - 如果 Gauss 约束被严格满足，此值应该接近 0

## 实现细节

### 1. 新增场变量

在 `Hamiltonian_Z2_Matter_smod.F90` 中添加了以下存储：

```fortran
! Lambda 场数组：lambda_field(site, tau) = +1 或 -1
Integer, allocatable :: lambda_field(:,:)

! Star product 缓存
Integer, allocatable :: star_product_cache(:,:)
```

### 2. 场类型

Lambda 场被添加为第 5 种场类型（Field_type = 5）：
- Field_type = 1: Z₂ 规范场
- Field_type = 2: Bond matter 场
- Field_type = 3: Hubbard HS 场
- Field_type = 4: Site matter 场
- **Field_type = 5: Gauss lambda 场** (新增)

### 3. 核心函数

#### `Compute_Star_Product_X(I, nt)`
计算 star product $X_r(\tau) = \prod_{b\in +r} \sigma^x_b(\tau)$

#### `Compute_Gauss_Operator(I, nt, GRC)`
计算 Gauss 算符 $G_r = (-1)^{n_r} \cdot X_r$

#### `Compute_Gauss_Phase(I, nt)`
计算 Gauss 相位因子 $\exp(i \lambda_r \pi/2 \cdot X_r)$

### 4. 蒙特卡洛更新

Lambda 场通过标准的顺序更新（sequential updates）进行采样，与其他 Ising 场相同。更新范围在 `Overide_global_tau_sampling_parameters` 中设置。

### 5. 传播子修改

Gauss 相位因子作为对角的一体算符插入到传播子 $B(\tau)$ 中：

$$B(\tau) = e^{-\Delta\tau K}\, e^{-\Delta\tau V(\sigma(\tau))}\, B^{(\text{Gauss})}(\tau)$$

其中：
$$B^{(\text{Gauss})}(\tau) = \prod_r \exp\left[ i\lambda_r(\tau)\frac{\pi}{2} X_r(\tau)\right]$$

## 与文献的对应关系

本实现与以下文献中的方法保持一致：
- Nat. Phys. 2017: Gazit et al.
- PNAS 2018: 相关工作
- PRX 2020: Assaad 等

## 验证方法

### 1. 检查 Gauss 约束满足程度

运行模拟后，检查输出文件：
- `Gauss_scal.dat`: 应该接近 1
- `GaussViol_scal.dat`: 应该接近 0

### 2. 对比测试

建议进行以下对比：
1. `UseStrictGauss = .false.` (原有实现)
2. `UseStrictGauss = .true.` (严格 Gauss 约束)

比较两种情况下的：
- 相图结构
- 临界指数
- 费米面结构

## 注意事项

1. **数值稳定性**: 新的相位因子可能在某些参数区域影响数值稳定性，建议从较小的系统开始测试。

2. **符号问题**: 在论文参数下应该是 sign-free，但建议在新参数区域进行验证。

3. **性能**: Lambda 场的更新会增加计算量，但与原有 Ising 场更新的开销相当。

## 示例参数文件

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
```

## 文件修改列表

主要修改的文件：
- `Prog/Hamiltonians/Hamiltonian_Z2_Matter_smod.F90`
  - 添加 `UseStrictGauss` 参数
  - 添加 `lambda_field` 存储
  - 添加辅助函数
  - 修改 `Ham_V`、`S0`、`Hamiltonian_set_nsigma` 等

## 作者

ALF Collaboration
