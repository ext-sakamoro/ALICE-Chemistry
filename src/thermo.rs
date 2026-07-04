//! 熱力学 + pH

use crate::constants::GAS_CONSTANT;

/// 熱力学状態
#[derive(Debug, Clone, Copy)]
pub struct ThermodynamicState {
    /// 温度 (K)
    pub temperature: f64,
    /// 圧力 (Pa)
    pub pressure: f64,
    /// モル数
    pub moles: f64,
}

/// ギブズ自由エネルギー
#[must_use]
pub fn gibbs_free_energy(enthalpy: f64, temperature: f64, entropy: f64) -> f64 {
    temperature.mul_add(-entropy, enthalpy)
}

/// 反応が自発的か判定
#[must_use]
pub fn is_spontaneous(delta_g: f64) -> bool {
    delta_g < 0.0
}

/// 理想気体の状態方程式 PV = nRT
#[must_use]
pub fn ideal_gas_volume(moles: f64, temperature: f64, pressure: f64) -> f64 {
    moles * GAS_CONSTANT * temperature / pressure
}

/// 理想気体の圧力
#[must_use]
pub fn ideal_gas_pressure(moles: f64, temperature: f64, volume: f64) -> f64 {
    moles * GAS_CONSTANT * temperature / volume
}

/// エントロピー変化（可逆過程）
#[must_use]
pub fn entropy_change_reversible(heat: f64, temperature: f64) -> f64 {
    heat / temperature
}

/// 等温膨張でのエントロピー変化
#[must_use]
pub fn entropy_isothermal_expansion(moles: f64, v1: f64, v2: f64) -> f64 {
    moles * GAS_CONSTANT * (v2 / v1).ln()
}

/// 定圧熱容量による温度変化のエントロピー
#[must_use]
pub fn entropy_temperature_change(moles: f64, cp: f64, t1: f64, t2: f64) -> f64 {
    moles * cp * (t2 / t1).ln()
}

/// 反応エンタルピー（生成エンタルピーの差）
#[must_use]
pub fn reaction_enthalpy(
    product_enthalpies: &[(f64, f64)],
    reactant_enthalpies: &[(f64, f64)],
) -> f64 {
    let prod: f64 = product_enthalpies.iter().map(|(h, n)| h * n).sum();
    let react: f64 = reactant_enthalpies.iter().map(|(h, n)| h * n).sum();
    prod - react
}

/// ヘスの法則: 中間反応のエンタルピーを加算
#[must_use]
pub fn hess_law(enthalpies: &[f64]) -> f64 {
    enthalpies.iter().sum()
}

/// 温度依存エンタルピー変化（キルヒホッフの式）
#[must_use]
pub fn kirchhoff_enthalpy(delta_h_t1: f64, delta_cp: f64, t1: f64, t2: f64) -> f64 {
    delta_cp.mul_add(t2 - t1, delta_h_t1)
}

/// カルノーサイクルの効率
#[must_use]
pub fn carnot_efficiency(t_hot: f64, t_cold: f64) -> f64 {
    1.0 - t_cold / t_hot
}

/// 内部エネルギー変化（第一法則）
#[must_use]
pub fn internal_energy_change(heat: f64, work: f64) -> f64 {
    heat + work
}

/// pH = -log10([H+])
#[must_use]
pub fn ph(h_concentration: f64) -> f64 {
    -h_concentration.log10()
}

/// pOH = -log10([OH-])
#[must_use]
pub fn poh(oh_concentration: f64) -> f64 {
    -oh_concentration.log10()
}

/// [H+] from pH
#[must_use]
pub fn h_concentration_from_ph(ph_val: f64) -> f64 {
    10.0_f64.powf(-ph_val)
}

/// Henderson-Hasselbalch equation
#[must_use]
pub fn henderson_hasselbalch(pka: f64, conjugate_base: f64, weak_acid: f64) -> f64 {
    pka + (conjugate_base / weak_acid).log10()
}
