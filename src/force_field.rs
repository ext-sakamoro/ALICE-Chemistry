//! 力場: Lennard-Jones + Coulomb + 結合エネルギー

use crate::constants::COULOMB_CONSTANT;
use crate::vec3::Vec3;

/// Lennard-Jonesポテンシャルのパラメータ
#[derive(Debug, Clone, Copy)]
pub struct LennardJonesParams {
    /// 井戸の深さ epsilon (J)
    pub epsilon: f64,
    /// 衝突直径 sigma (m)
    pub sigma: f64,
}

impl LennardJonesParams {
    #[must_use]
    pub const fn new(epsilon: f64, sigma: f64) -> Self {
        Self { epsilon, sigma }
    }

    /// Lennard-Jonesポテンシャルエネルギー
    #[must_use]
    pub fn potential(&self, r: f64) -> f64 {
        let sr = self.sigma / r;
        let sr6 = sr * sr * sr * sr * sr * sr;
        let sr12 = sr6 * sr6;
        4.0 * self.epsilon * (sr12 - sr6)
    }

    /// Lennard-Jones力の大きさ
    #[must_use]
    pub fn force_magnitude(&self, r: f64) -> f64 {
        let sr = self.sigma / r;
        let sr6 = sr * sr * sr * sr * sr * sr;
        let sr12 = sr6 * sr6;
        24.0 * self.epsilon / r * 2.0f64.mul_add(sr12, -sr6)
    }

    /// 平衡距離
    #[must_use]
    pub fn equilibrium_distance(&self) -> f64 {
        self.sigma * (1.0_f64 / 6.0).exp2()
    }

    /// 2粒子間の力ベクトル
    #[must_use]
    pub fn force_vector(&self, pos_a: Vec3, pos_b: Vec3) -> Vec3 {
        let diff = pos_b - pos_a;
        let r = diff.length();
        if r < 1e-15 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        let f_mag = self.force_magnitude(r);
        diff.scale(f_mag / r)
    }
}

/// クーロン力
pub struct CoulombForce;

impl CoulombForce {
    /// クーロンポテンシャル
    #[must_use]
    pub fn potential(q1: f64, q2: f64, r: f64) -> f64 {
        COULOMB_CONSTANT * q1 * q2 / r
    }

    /// クーロン力の大きさ
    #[must_use]
    pub fn force_magnitude(q1: f64, q2: f64, r: f64) -> f64 {
        COULOMB_CONSTANT * (q1 * q2).abs() / (r * r)
    }

    /// クーロン力ベクトル
    #[must_use]
    pub fn force_vector(q1: f64, q2: f64, pos_a: Vec3, pos_b: Vec3) -> Vec3 {
        let diff = pos_a - pos_b;
        let r = diff.length();
        if r < 1e-15 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        let sign = if q1 * q2 > 0.0 { 1.0 } else { -1.0 };
        let f_mag = Self::force_magnitude(q1, q2, r);
        diff.scale(sign * f_mag / r)
    }

    /// 電場の大きさ
    #[must_use]
    pub fn electric_field_magnitude(q: f64, r: f64) -> f64 {
        COULOMB_CONSTANT * q.abs() / (r * r)
    }
}

/// 結合の種類
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}

/// 代表的な結合エネルギー (kJ/mol)
#[must_use]
pub fn bond_energy(atom1: &str, atom2: &str, bond_type: BondType) -> Option<f64> {
    match (atom1, atom2, bond_type) {
        ("H", "H", BondType::Single) => Some(436.0),
        ("C", "H", BondType::Single) | ("H", "C", BondType::Single) => Some(413.0),
        ("C", "C", BondType::Single) => Some(348.0),
        ("C", "C", BondType::Double) => Some(614.0),
        ("C", "C", BondType::Triple) => Some(839.0),
        ("C", "C", BondType::Aromatic) => Some(518.0),
        ("C", "O", BondType::Single) | ("O", "C", BondType::Single) => Some(360.0),
        ("C", "O", BondType::Double) | ("O", "C", BondType::Double) => Some(743.0),
        ("C", "N", BondType::Single) | ("N", "C", BondType::Single) => Some(305.0),
        ("C", "N", BondType::Double) | ("N", "C", BondType::Double) => Some(615.0),
        ("C", "N", BondType::Triple) | ("N", "C", BondType::Triple) => Some(891.0),
        ("O", "H", BondType::Single) | ("H", "O", BondType::Single) => Some(463.0),
        ("O", "O", BondType::Single) => Some(146.0),
        ("O", "O", BondType::Double) => Some(497.0),
        ("N", "H", BondType::Single) | ("H", "N", BondType::Single) => Some(391.0),
        ("N", "N", BondType::Single) => Some(163.0),
        ("N", "N", BondType::Double) => Some(418.0),
        ("N", "N", BondType::Triple) => Some(941.0),
        ("H", "F", BondType::Single) | ("F", "H", BondType::Single) => Some(567.0),
        ("H", "Cl", BondType::Single) | ("Cl", "H", BondType::Single) => Some(431.0),
        ("H", "Br", BondType::Single) | ("Br", "H", BondType::Single) => Some(366.0),
        ("C", "Cl", BondType::Single) | ("Cl", "C", BondType::Single) => Some(339.0),
        ("C", "F", BondType::Single) | ("F", "C", BondType::Single) => Some(485.0),
        ("S", "H", BondType::Single) | ("H", "S", BondType::Single) => Some(363.0),
        ("C", "S", BondType::Single) | ("S", "C", BondType::Single) => Some(272.0),
        ("C", "S", BondType::Double) | ("S", "C", BondType::Double) => Some(573.0),
        _ => None,
    }
}

/// 結合エネルギーから反応エンタルピーを概算
#[must_use]
pub fn reaction_enthalpy_from_bonds(reactant_bonds: f64, product_bonds: f64) -> f64 {
    reactant_bonds - product_bonds
}
