//! 反応速度論 (Arrhenius) + 化学平衡 + 化学量論

use crate::constants::{AVOGADRO, GAS_CONSTANT};

/// アレニウスの式
#[derive(Debug, Clone, Copy)]
pub struct ArrheniusParams {
    /// 頻度因子 A (1/s)
    pub pre_exponential: f64,
    /// 活性化エネルギー Ea (J/mol)
    pub activation_energy: f64,
}

impl ArrheniusParams {
    #[must_use]
    pub const fn new(pre_exponential: f64, activation_energy: f64) -> Self {
        Self {
            pre_exponential,
            activation_energy,
        }
    }

    /// 速度定数 k(T)
    #[must_use]
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        self.pre_exponential * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp()
    }

    /// 2つの温度での速度定数の比 k2/k1
    #[must_use]
    pub fn rate_ratio(&self, t1: f64, t2: f64) -> f64 {
        (self.activation_energy / GAS_CONSTANT * (1.0 / t1 - 1.0 / t2)).exp()
    }

    /// 2つの温度での速度定数から活性化エネルギーを逆算
    #[must_use]
    pub fn activation_energy_from_rates(k1: f64, k2: f64, t1: f64, t2: f64) -> f64 {
        GAS_CONSTANT * (k2 / k1).ln() / (1.0 / t1 - 1.0 / t2)
    }
}

/// 反応次数ごとの濃度変化
#[derive(Debug, Clone, Copy)]
pub enum ReactionOrder {
    /// 零次
    Zero,
    /// 一次
    First,
    /// 二次
    Second,
}

impl ReactionOrder {
    /// 時刻tでの濃度
    #[must_use]
    pub fn concentration(self, initial: f64, k: f64, t: f64) -> f64 {
        match self {
            Self::Zero => k.mul_add(-t, initial).max(0.0),
            Self::First => initial * (-k * t).exp(),
            Self::Second => {
                let denom = k.mul_add(t, 1.0 / initial);
                if denom > 0.0 {
                    1.0 / denom
                } else {
                    0.0
                }
            }
        }
    }

    /// 半減期
    #[must_use]
    pub fn half_life(self, initial: f64, k: f64) -> f64 {
        match self {
            Self::Zero => initial / (2.0 * k),
            Self::First => core::f64::consts::LN_2 / k,
            Self::Second => 1.0 / (k * initial),
        }
    }
}

/// 化学平衡定数から自由エネルギー変化を算出
#[must_use]
pub fn gibbs_from_equilibrium(keq: f64, temperature: f64) -> f64 {
    -GAS_CONSTANT * temperature * keq.ln()
}

/// 自由エネルギー変化から平衡定数を算出
#[must_use]
pub fn equilibrium_constant(delta_g: f64, temperature: f64) -> f64 {
    (-delta_g / (GAS_CONSTANT * temperature)).exp()
}

/// ヴァント・ホッフの式
#[must_use]
pub fn vant_hoff(k1: f64, delta_h: f64, t1: f64, t2: f64) -> f64 {
    k1 * (-delta_h / GAS_CONSTANT * (1.0 / t2 - 1.0 / t1)).exp()
}

/// 反応商 Q の計算
#[must_use]
pub fn reaction_quotient(products: &[(f64, f64)], reactants: &[(f64, f64)]) -> f64 {
    let numerator: f64 = products
        .iter()
        .map(|(conc, coeff)| conc.powf(*coeff))
        .product();
    let denominator: f64 = reactants
        .iter()
        .map(|(conc, coeff)| conc.powf(*coeff))
        .product();
    if denominator.abs() < 1e-30 {
        return f64::INFINITY;
    }
    numerator / denominator
}

/// Q と K の比較に基づく反応方向
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionDirection {
    Forward,
    Reverse,
    Equilibrium,
}

/// Q と K を比較して反応方向を判定
#[must_use]
pub fn predict_direction(q: f64, keq: f64) -> ReactionDirection {
    let ratio = q / keq;
    if ratio < 0.999 {
        ReactionDirection::Forward
    } else if ratio > 1.001 {
        ReactionDirection::Reverse
    } else {
        ReactionDirection::Equilibrium
    }
}

/// 化学量論: モル数から質量 (g)
#[must_use]
pub fn moles_to_grams(moles: f64, molar_mass: f64) -> f64 {
    moles * molar_mass
}

/// 化学量論: 質量からモル数
#[must_use]
pub fn grams_to_moles(grams: f64, molar_mass: f64) -> f64 {
    grams / molar_mass
}

/// 化学量論: モル数から分子数
#[must_use]
pub fn moles_to_molecules(moles: f64) -> f64 {
    moles * AVOGADRO
}

/// 化学量論: 分子数からモル数
#[must_use]
pub fn molecules_to_moles(molecules: f64) -> f64 {
    molecules / AVOGADRO
}

/// モル濃度 (mol/L)
#[must_use]
pub fn molarity(moles: f64, volume_liters: f64) -> f64 {
    moles / volume_liters
}

/// 希釈の式 M1*V1 = M2*V2
#[must_use]
pub fn dilution_volume(m1: f64, v1: f64, m2: f64) -> f64 {
    m1 * v1 / m2
}

/// 制限試薬の特定
#[must_use]
pub fn limiting_reagent(reagents: &[(f64, f64)]) -> usize {
    reagents
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let ra = a.0 / a.1;
            let rb = b.0 / b.1;
            ra.partial_cmp(&rb).unwrap_or(core::cmp::Ordering::Equal)
        })
        .map_or(0, |(i, _)| i)
}

/// 理論収量 (mol)
#[must_use]
pub fn theoretical_yield(limiting_moles: f64, limiting_coeff: f64, product_coeff: f64) -> f64 {
    limiting_moles / limiting_coeff * product_coeff
}

/// 収率 (%)
#[must_use]
pub fn percent_yield(actual: f64, theoretical: f64) -> f64 {
    actual / theoretical * 100.0
}

/// 経験式から分子式の倍数を求める
#[must_use]
pub fn molecular_formula_multiplier(empirical_mass: f64, molecular_mass: f64) -> u32 {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = (molecular_mass / empirical_mass).round() as u32;
    n.max(1)
}
