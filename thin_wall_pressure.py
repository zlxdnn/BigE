"""薄壁容器内压估算（考虑弯矩平衡与稳定性约束）.

模型说明（工程假设）：
1) 容器近似为薄壁圆柱壳体，满足 t / r << 1。
2) 载荷主要来源：
   - 外部弯矩 M（例如支撑偏心、外部管线载荷等）
   - 轴向附加载荷 N（可选）
3) 使用线弹性、小变形近似，材料在弹性范围内。
4) 薄壁圆柱在内压 p 作用下：
   - 环向膜应力 sigma_hoop = p * r / t
   - 轴向膜应力 sigma_axial_p = p * r / (2t)
5) 外部弯矩引起的最大轴向弯曲应力（薄壁环形截面）
   sigma_b = M * r / I, 其中 I = pi * r^3 * t。
6) 设计准则：
   a) 弯矩平衡目标（可选参数 balance_ratio）:
      令内压产生的轴向拉应力抵消一定比例的弯曲应力，
      即 sigma_axial_p = balance_ratio * sigma_b。
   b) 稳定性/强度约束：
      使用简化等效应力（von Mises）不超过许用应力 sigma_allow。

注意：
- 该模型用于方案估算，不替代完整壳体有限元、屈曲校核和规范设计。
- 若工况涉及高温、循环疲劳、局部屈曲、开孔补强等，应进行更详细分析。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from typing import Optional


@dataclass
class VesselInput:
    """输入参数.

    Attributes:
        radius: 容器中径半径 r (m)
        thickness: 壳体厚度 t (m)
        bending_moment: 外部弯矩 M (N·m)
        allowable_stress: 许用应力 sigma_allow (Pa)
        axial_force: 轴向附加载荷 N (N), 拉为正, 压为负
        balance_ratio: 弯矩抵消比例 [0, 1], 1 表示完全抵消
        pressure_min: 工艺给定最小内压 (Pa)
        pressure_max: 工艺给定最大内压 (Pa)
    """

    radius: float
    thickness: float
    bending_moment: float
    allowable_stress: float
    axial_force: float = 0.0
    balance_ratio: float = 1.0
    pressure_min: float = 0.0
    pressure_max: Optional[float] = None


@dataclass
class VesselResult:
    """计算结果."""

    pressure_from_balance: float
    pressure_from_strength_limit: float
    recommended_pressure: float
    hoop_stress: float
    axial_stress_total: float
    von_mises_stress: float
    is_within_allowable: bool


def _validate_inputs(inp: VesselInput) -> None:
    if inp.radius <= 0:
        raise ValueError("radius 必须 > 0")
    if inp.thickness <= 0:
        raise ValueError("thickness 必须 > 0")
    if inp.allowable_stress <= 0:
        raise ValueError("allowable_stress 必须 > 0")
    if not (0.0 <= inp.balance_ratio <= 1.0):
        raise ValueError("balance_ratio 必须在 [0, 1] 内")
    if inp.pressure_min < 0:
        raise ValueError("pressure_min 不能为负")
    if inp.pressure_max is not None and inp.pressure_max < inp.pressure_min:
        raise ValueError("pressure_max 不能小于 pressure_min")


def section_second_moment(radius: float, thickness: float) -> float:
    """薄壁圆环截面对中性轴的二次矩 I = pi * r^3 * t."""
    return pi * radius**3 * thickness


def bending_stress_from_moment(moment: float, radius: float, thickness: float) -> float:
    """外弯矩导致的最大轴向弯曲应力（绝对值）."""
    ixx = section_second_moment(radius, thickness)
    return abs(moment) * radius / ixx


def axial_stress_from_force(axial_force: float, radius: float, thickness: float) -> float:
    """轴向附加载荷导致的平均轴向应力 N/A, A≈2*pi*r*t."""
    area = 2.0 * pi * radius * thickness
    return axial_force / area


def pressure_for_moment_balance(
    moment: float,
    radius: float,
    thickness: float,
    balance_ratio: float = 1.0,
) -> float:
    """按弯矩平衡目标计算所需内压.

    由: p*r/(2t) = balance_ratio * sigma_b
    得: p = 2t/r * balance_ratio * sigma_b
    """
    sigma_b = bending_stress_from_moment(moment, radius, thickness)
    return (2.0 * thickness / radius) * balance_ratio * sigma_b


def von_mises_plane_stress(sigma_x: float, sigma_y: float, tau_xy: float = 0.0) -> float:
    """平面应力下 von Mises 等效应力."""
    return sqrt(sigma_x**2 - sigma_x * sigma_y + sigma_y**2 + 3.0 * tau_xy**2)


def pressure_strength_upper_bound(inp: VesselInput) -> float:
    """由强度约束反推可接受内压上限（简化搜索）.

    约束: sigma_vm(p) <= sigma_allow
    其中:
      sigma_hoop = p*r/t
      sigma_axial = p*r/(2t) + sigma_N +/- sigma_b
    对最不利纤维取 "+ sigma_b".
    """
    sigma_b = bending_stress_from_moment(inp.bending_moment, inp.radius, inp.thickness)
    sigma_n = axial_stress_from_force(inp.axial_force, inp.radius, inp.thickness)

    def vm_at_pressure(p: float) -> float:
        sigma_hoop = p * inp.radius / inp.thickness
        sigma_axial = p * inp.radius / (2.0 * inp.thickness) + sigma_n + sigma_b
        return von_mises_plane_stress(sigma_axial, sigma_hoop)

    # 先指数扩展上界，再二分。
    lo, hi = 0.0, max(inp.allowable_stress * inp.thickness / inp.radius, 1.0)
    while vm_at_pressure(hi) < inp.allowable_stress:
        hi *= 2.0
        if hi > 1e10:  # 防止异常工况无限扩张
            break

    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if vm_at_pressure(mid) <= inp.allowable_stress:
            lo = mid
        else:
            hi = mid

    return lo


def evaluate_pressure(inp: VesselInput) -> VesselResult:
    """综合弯矩平衡与稳定性约束，给出推荐内压。"""
    _validate_inputs(inp)

    p_balance = pressure_for_moment_balance(
        moment=inp.bending_moment,
        radius=inp.radius,
        thickness=inp.thickness,
        balance_ratio=inp.balance_ratio,
    )

    p_strength_max = pressure_strength_upper_bound(inp)

    p_rec = max(inp.pressure_min, p_balance)
    if inp.pressure_max is not None:
        p_rec = min(p_rec, inp.pressure_max)
    p_rec = min(p_rec, p_strength_max)

    sigma_b = bending_stress_from_moment(inp.bending_moment, inp.radius, inp.thickness)
    sigma_n = axial_stress_from_force(inp.axial_force, inp.radius, inp.thickness)
    sigma_hoop = p_rec * inp.radius / inp.thickness
    sigma_axial = p_rec * inp.radius / (2.0 * inp.thickness) + sigma_n + sigma_b
    sigma_vm = von_mises_plane_stress(sigma_axial, sigma_hoop)

    return VesselResult(
        pressure_from_balance=p_balance,
        pressure_from_strength_limit=p_strength_max,
        recommended_pressure=p_rec,
        hoop_stress=sigma_hoop,
        axial_stress_total=sigma_axial,
        von_mises_stress=sigma_vm,
        is_within_allowable=sigma_vm <= inp.allowable_stress,
    )


def _demo() -> None:
    """示例运行。"""
    inp = VesselInput(
        radius=0.8,
        thickness=0.012,
        bending_moment=1.8e5,
        allowable_stress=1.5e8,
        axial_force=5.0e4,
        balance_ratio=0.9,
        pressure_min=2.0e5,
        pressure_max=2.5e6,
    )

    res = evaluate_pressure(inp)
    print("=== 薄壁容器内压估算结果 ===")
    print(f"弯矩平衡需求内压: {res.pressure_from_balance:,.1f} Pa")
    print(f"强度约束内压上限: {res.pressure_from_strength_limit:,.1f} Pa")
    print(f"推荐内压: {res.recommended_pressure:,.1f} Pa")
    print(f"环向应力: {res.hoop_stress:,.1f} Pa")
    print(f"轴向总应力(最不利): {res.axial_stress_total:,.1f} Pa")
    print(f"von Mises 等效应力: {res.von_mises_stress:,.1f} Pa")
    print(f"是否满足许用应力: {res.is_within_allowable}")


if __name__ == "__main__":
    _demo()
