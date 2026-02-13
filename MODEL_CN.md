# 薄壁容器内压计算模型（平衡弯矩 + 稳定性）

## 1. 工程假设
- 容器为薄壁圆柱壳，满足 \(t/r \ll 1\)。
- 材料线弹性，小变形。
- 主要考虑：外部弯矩 \(M\)、轴向载荷 \(N\)、内压 \(p\)。
- 忽略局部应力集中、焊缝缺陷、开孔补强、热应力等二阶效应（用于初步方案估算）。

## 2. 力学模型

### 2.1 膜应力（内压）
- 环向应力：\(\sigma_\theta = p r / t\)
- 轴向应力：\(\sigma_{z,p} = p r / (2t)\)

### 2.2 弯矩引起的轴向应力
薄壁圆环截面二次矩：
\[
I = \pi r^3 t
\]
弯曲最大应力：
\[
\sigma_b = \frac{M r}{I}
\]

### 2.3 轴向载荷引起应力
薄壁截面积近似 \(A \approx 2\pi r t\)：
\[
\sigma_N = N/A
\]

## 3. 设计准则

### 3.1 弯矩平衡目标
定义抵消比例 \(\eta\in[0,1]\)：
\[
\sigma_{z,p} = \eta \sigma_b
\]
由此得到目标内压：
\[
p_{balance} = \frac{2t}{r}\eta\sigma_b
\]

### 3.2 强度/稳定性约束（简化）
采用平面应力 von Mises：
\[
\sigma_{vm} = \sqrt{\sigma_z^2 - \sigma_z\sigma_\theta + \sigma_\theta^2}
\]
其中最不利纤维取：
\[
\sigma_z = \frac{pr}{2t} + \sigma_N + \sigma_b
\]
要求 \(\sigma_{vm} \le \sigma_{allow}\)。

## 4. 计算策略
1. 由弯矩平衡计算 \(p_{balance}\)
2. 由许用应力反算内压上限 \(p_{strength,max}\)
3. 在工艺边界 \([p_{min}, p_{max}]\) 内取推荐值：
\[
p_{rec} = \min\big(\max(p_{min}, p_{balance}),\; p_{max},\; p_{strength,max}\big)
\]

> 说明：该模型中的“稳定性”由简化强度裕度体现；若涉及真实壳体屈曲（外压/轴压/几何缺陷敏感），应补充 ASME/EN 规范屈曲校核或有限元特征屈曲与非线性分析。

## 5. 代码
实现见 `thin_wall_pressure.py`，包含：
- 输入/输出数据结构
- 弯矩平衡内压计算
- 强度上限反算（二分搜索）
- 综合推荐内压计算
- 示例运行
