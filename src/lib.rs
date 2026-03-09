#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::similar_names)]
#![allow(clippy::doc_markdown)]
#![allow(clippy::cast_precision_loss)]

//! ALICE-Chemistry: 分子動力学・化学反応シミュレーションライブラリ
//!
//! - 周期表データ（元素情報）
//! - 力場（Lennard-Jones, Coulomb）
//! - 反応速度論（Arrhenius式）
//! - 化学平衡
//! - 分子構造・結合エネルギー
//! - 熱力学（エンタルピー、エントロピー、ギブズ自由エネルギー）
//! - 化学量論

use core::f64::consts::PI;
use core::ops::{Add, Sub};

// ============================================================
// 定数
// ============================================================

/// ボルツマン定数 (J/K)
pub const BOLTZMANN: f64 = 1.380_649e-23;

/// アボガドロ数 (1/mol)
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// 気体定数 R (J/(mol K))
pub const GAS_CONSTANT: f64 = 8.314_462_618;

/// クーロン定数 (N m^2 / C^2)
pub const COULOMB_CONSTANT: f64 = 8.987_551_792e9;

/// 素電荷 (C)
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// 真空の誘電率 (F/m)
pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_817e-12;

// ============================================================
// 周期表
// ============================================================

/// 元素データ
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Element {
    /// 原子番号
    pub atomic_number: u32,
    /// 元素記号
    pub symbol: &'static str,
    /// 元素名
    pub name: &'static str,
    /// 原子量 (u)
    pub atomic_mass: f64,
    /// 電気陰性度 (Pauling)
    pub electronegativity: Option<f64>,
}

/// 最初の36元素（H-Kr）のデータ
const ELEMENTS: [Element; 36] = [
    Element {
        atomic_number: 1,
        symbol: "H",
        name: "Hydrogen",
        atomic_mass: 1.008,
        electronegativity: Some(2.20),
    },
    Element {
        atomic_number: 2,
        symbol: "He",
        name: "Helium",
        atomic_mass: 4.003,
        electronegativity: None,
    },
    Element {
        atomic_number: 3,
        symbol: "Li",
        name: "Lithium",
        atomic_mass: 6.941,
        electronegativity: Some(0.98),
    },
    Element {
        atomic_number: 4,
        symbol: "Be",
        name: "Beryllium",
        atomic_mass: 9.012,
        electronegativity: Some(1.57),
    },
    Element {
        atomic_number: 5,
        symbol: "B",
        name: "Boron",
        atomic_mass: 10.81,
        electronegativity: Some(2.04),
    },
    Element {
        atomic_number: 6,
        symbol: "C",
        name: "Carbon",
        atomic_mass: 12.011,
        electronegativity: Some(2.55),
    },
    Element {
        atomic_number: 7,
        symbol: "N",
        name: "Nitrogen",
        atomic_mass: 14.007,
        electronegativity: Some(3.04),
    },
    Element {
        atomic_number: 8,
        symbol: "O",
        name: "Oxygen",
        atomic_mass: 15.999,
        electronegativity: Some(3.44),
    },
    Element {
        atomic_number: 9,
        symbol: "F",
        name: "Fluorine",
        atomic_mass: 18.998,
        electronegativity: Some(3.98),
    },
    Element {
        atomic_number: 10,
        symbol: "Ne",
        name: "Neon",
        atomic_mass: 20.180,
        electronegativity: None,
    },
    Element {
        atomic_number: 11,
        symbol: "Na",
        name: "Sodium",
        atomic_mass: 22.990,
        electronegativity: Some(0.93),
    },
    Element {
        atomic_number: 12,
        symbol: "Mg",
        name: "Magnesium",
        atomic_mass: 24.305,
        electronegativity: Some(1.31),
    },
    Element {
        atomic_number: 13,
        symbol: "Al",
        name: "Aluminium",
        atomic_mass: 26.982,
        electronegativity: Some(1.61),
    },
    Element {
        atomic_number: 14,
        symbol: "Si",
        name: "Silicon",
        atomic_mass: 28.086,
        electronegativity: Some(1.90),
    },
    Element {
        atomic_number: 15,
        symbol: "P",
        name: "Phosphorus",
        atomic_mass: 30.974,
        electronegativity: Some(2.19),
    },
    Element {
        atomic_number: 16,
        symbol: "S",
        name: "Sulfur",
        atomic_mass: 32.06,
        electronegativity: Some(2.58),
    },
    Element {
        atomic_number: 17,
        symbol: "Cl",
        name: "Chlorine",
        atomic_mass: 35.45,
        electronegativity: Some(3.16),
    },
    Element {
        atomic_number: 18,
        symbol: "Ar",
        name: "Argon",
        atomic_mass: 39.948,
        electronegativity: None,
    },
    Element {
        atomic_number: 19,
        symbol: "K",
        name: "Potassium",
        atomic_mass: 39.098,
        electronegativity: Some(0.82),
    },
    Element {
        atomic_number: 20,
        symbol: "Ca",
        name: "Calcium",
        atomic_mass: 40.078,
        electronegativity: Some(1.00),
    },
    Element {
        atomic_number: 21,
        symbol: "Sc",
        name: "Scandium",
        atomic_mass: 44.956,
        electronegativity: Some(1.36),
    },
    Element {
        atomic_number: 22,
        symbol: "Ti",
        name: "Titanium",
        atomic_mass: 47.867,
        electronegativity: Some(1.54),
    },
    Element {
        atomic_number: 23,
        symbol: "V",
        name: "Vanadium",
        atomic_mass: 50.942,
        electronegativity: Some(1.63),
    },
    Element {
        atomic_number: 24,
        symbol: "Cr",
        name: "Chromium",
        atomic_mass: 51.996,
        electronegativity: Some(1.66),
    },
    Element {
        atomic_number: 25,
        symbol: "Mn",
        name: "Manganese",
        atomic_mass: 54.938,
        electronegativity: Some(1.55),
    },
    Element {
        atomic_number: 26,
        symbol: "Fe",
        name: "Iron",
        atomic_mass: 55.845,
        electronegativity: Some(1.83),
    },
    Element {
        atomic_number: 27,
        symbol: "Co",
        name: "Cobalt",
        atomic_mass: 58.933,
        electronegativity: Some(1.88),
    },
    Element {
        atomic_number: 28,
        symbol: "Ni",
        name: "Nickel",
        atomic_mass: 58.693,
        electronegativity: Some(1.91),
    },
    Element {
        atomic_number: 29,
        symbol: "Cu",
        name: "Copper",
        atomic_mass: 63.546,
        electronegativity: Some(1.90),
    },
    Element {
        atomic_number: 30,
        symbol: "Zn",
        name: "Zinc",
        atomic_mass: 65.38,
        electronegativity: Some(1.65),
    },
    Element {
        atomic_number: 31,
        symbol: "Ga",
        name: "Gallium",
        atomic_mass: 69.723,
        electronegativity: Some(1.81),
    },
    Element {
        atomic_number: 32,
        symbol: "Ge",
        name: "Germanium",
        atomic_mass: 72.630,
        electronegativity: Some(2.01),
    },
    Element {
        atomic_number: 33,
        symbol: "As",
        name: "Arsenic",
        atomic_mass: 74.922,
        electronegativity: Some(2.18),
    },
    Element {
        atomic_number: 34,
        symbol: "Se",
        name: "Selenium",
        atomic_mass: 78.971,
        electronegativity: Some(2.55),
    },
    Element {
        atomic_number: 35,
        symbol: "Br",
        name: "Bromine",
        atomic_mass: 79.904,
        electronegativity: Some(2.96),
    },
    Element {
        atomic_number: 36,
        symbol: "Kr",
        name: "Krypton",
        atomic_mass: 83.798,
        electronegativity: Some(3.00),
    },
];

/// 原子番号から元素を取得（1-indexed）
#[must_use]
pub const fn element_by_number(atomic_number: u32) -> Option<&'static Element> {
    if atomic_number >= 1 && atomic_number <= 36 {
        Some(&ELEMENTS[(atomic_number - 1) as usize])
    } else {
        None
    }
}

/// 元素記号から元素を取得
#[must_use]
pub fn element_by_symbol(symbol: &str) -> Option<&'static Element> {
    ELEMENTS.iter().find(|e| e.symbol == symbol)
}

// ============================================================
// 3Dベクトル
// ============================================================

/// 3次元ベクトル
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    #[must_use]
    pub fn length_squared(self) -> f64 {
        self.x
            .mul_add(self.x, self.y.mul_add(self.y, self.z * self.z))
    }

    #[must_use]
    pub fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    #[must_use]
    pub fn scale(self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.x
            .mul_add(other.x, self.y.mul_add(other.y, self.z * other.z))
    }

    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        (self - other).length()
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// ============================================================
// 力場: Lennard-Jones
// ============================================================

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

// ============================================================
// 力場: Coulomb
// ============================================================

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

// ============================================================
// 反応速度論: Arrhenius
// ============================================================

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

// ============================================================
// 化学平衡
// ============================================================

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

// ============================================================
// 結合エネルギー
// ============================================================

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

// ============================================================
// 分子構造
// ============================================================

/// 原子
#[derive(Debug, Clone)]
pub struct Atom {
    pub element: &'static str,
    pub position: Vec3,
    pub charge: f64,
    pub mass: f64,
}

/// 結合
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub atom_a: usize,
    pub atom_b: usize,
    pub bond_type: BondType,
}

/// 分子
#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl Molecule {
    #[must_use]
    pub const fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }

    pub fn add_atom(
        &mut self,
        element: &'static str,
        position: Vec3,
        charge: f64,
        mass: f64,
    ) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(Atom {
            element,
            position,
            charge,
            mass,
        });
        idx
    }

    pub fn add_bond(&mut self, atom_a: usize, atom_b: usize, bond_type: BondType) {
        self.bonds.push(Bond {
            atom_a,
            atom_b,
            bond_type,
        });
    }

    /// 分子量
    #[must_use]
    pub fn molecular_weight(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass).sum()
    }

    /// 重心
    #[must_use]
    pub fn center_of_mass(&self) -> Vec3 {
        let total_mass: f64 = self.atoms.iter().map(|a| a.mass).sum();
        if total_mass < 1e-30 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        let wx: f64 = self.atoms.iter().map(|a| a.mass * a.position.x).sum();
        let wy: f64 = self.atoms.iter().map(|a| a.mass * a.position.y).sum();
        let wz: f64 = self.atoms.iter().map(|a| a.mass * a.position.z).sum();
        Vec3::new(wx / total_mass, wy / total_mass, wz / total_mass)
    }

    /// 総結合エネルギー (kJ/mol)
    #[must_use]
    pub fn total_bond_energy(&self) -> f64 {
        self.bonds
            .iter()
            .filter_map(|b| {
                let a1 = &self.atoms[b.atom_a];
                let a2 = &self.atoms[b.atom_b];
                bond_energy(a1.element, a2.element, b.bond_type)
            })
            .sum()
    }

    /// 結合距離
    #[must_use]
    pub fn bond_length(&self, bond_index: usize) -> Option<f64> {
        self.bonds.get(bond_index).map(|b| {
            self.atoms[b.atom_a]
                .position
                .distance(self.atoms[b.atom_b].position)
        })
    }

    /// 結合角（3原子 i-j-k の角度をラジアンで返す）
    #[must_use]
    pub fn bond_angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let v1 = self.atoms[i].position - self.atoms[j].position;
        let v2 = self.atoms[k].position - self.atoms[j].position;
        let cos_angle = v1.dot(v2) / (v1.length() * v2.length());
        cos_angle.clamp(-1.0, 1.0).acos()
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================
// 熱力学
// ============================================================

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

// ============================================================
// 化学量論
// ============================================================

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

// ============================================================
// 分子動力学シミュレーション
// ============================================================

/// 粒子
#[derive(Debug, Clone)]
pub struct Particle {
    pub position: Vec3,
    pub velocity: Vec3,
    pub force: Vec3,
    pub mass: f64,
    pub charge: f64,
}

impl Particle {
    #[must_use]
    pub const fn new(position: Vec3, velocity: Vec3, mass: f64, charge: f64) -> Self {
        Self {
            position,
            velocity,
            force: Vec3::new(0.0, 0.0, 0.0),
            mass,
            charge,
        }
    }
}

/// Velocity Verlet積分器による1ステップ
pub fn velocity_verlet_step(particles: &mut [Particle], dt: f64, lj: &LennardJonesParams) {
    let n = particles.len();

    // 位置の更新
    for p in particles.iter_mut() {
        let ax = p.force.x / p.mass;
        let ay = p.force.y / p.mass;
        let az = p.force.z / p.mass;
        p.position.x += p.velocity.x.mul_add(dt, 0.5 * ax * dt * dt);
        p.position.y += p.velocity.y.mul_add(dt, 0.5 * ay * dt * dt);
        p.position.z += p.velocity.z.mul_add(dt, 0.5 * az * dt * dt);
    }

    // 古い力を保存
    let old_forces: Vec<Vec3> = particles.iter().map(|p| p.force).collect();

    // 力の再計算
    for p in particles.iter_mut() {
        p.force = Vec3::new(0.0, 0.0, 0.0);
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let f = lj.force_vector(particles[i].position, particles[j].position);
            particles[i].force = particles[i].force + f;
            particles[j].force = particles[j].force - f;
        }
    }

    // 速度の更新
    for (i, p) in particles.iter_mut().enumerate() {
        let old_accel_x = old_forces[i].x / p.mass;
        let old_accel_y = old_forces[i].y / p.mass;
        let old_accel_z = old_forces[i].z / p.mass;
        let new_accel_x = p.force.x / p.mass;
        let new_accel_y = p.force.y / p.mass;
        let new_accel_z = p.force.z / p.mass;
        p.velocity.x += 0.5 * (old_accel_x + new_accel_x) * dt;
        p.velocity.y += 0.5 * (old_accel_y + new_accel_y) * dt;
        p.velocity.z += 0.5 * (old_accel_z + new_accel_z) * dt;
    }
}

/// 系の運動エネルギー
#[must_use]
pub fn kinetic_energy(particles: &[Particle]) -> f64 {
    particles
        .iter()
        .map(|p| 0.5 * p.mass * p.velocity.length_squared())
        .sum()
}

/// 系の温度（等分配定理）
#[must_use]
pub fn system_temperature(particles: &[Particle]) -> f64 {
    let ek = kinetic_energy(particles);
    let n = particles.len() as f64;
    if n < 1.0 {
        return 0.0;
    }
    2.0 * ek / (3.0 * n * BOLTZMANN)
}

/// マクスウェル-ボルツマン分布の最確速度
#[must_use]
pub fn most_probable_speed(mass: f64, temperature: f64) -> f64 {
    (2.0 * BOLTZMANN * temperature / mass).sqrt()
}

/// マクスウェル-ボルツマン分布の平均速度
#[must_use]
pub fn mean_speed(mass: f64, temperature: f64) -> f64 {
    (8.0 * BOLTZMANN * temperature / (PI * mass)).sqrt()
}

/// マクスウェル-ボルツマン分布のRMS速度
#[must_use]
pub fn rms_speed(mass: f64, temperature: f64) -> f64 {
    (3.0 * BOLTZMANN * temperature / mass).sqrt()
}

// ============================================================
// pH
// ============================================================

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

// ============================================================
// テスト
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-6;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn assert_approx(a: f64, b: f64, tol: f64) {
        assert!(
            approx_eq(a, b, tol),
            "assertion failed: {a} != {b} (tol={tol})"
        );
    }

    // --- 周期表 ---
    #[test]
    fn test_element_hydrogen() {
        let h = element_by_number(1).unwrap();
        assert_eq!(h.symbol, "H");
        assert_eq!(h.name, "Hydrogen");
        assert_approx(h.atomic_mass, 1.008, 0.01);
    }

    #[test]
    fn test_element_carbon() {
        let c = element_by_number(6).unwrap();
        assert_eq!(c.symbol, "C");
        assert_approx(c.atomic_mass, 12.011, 0.01);
    }

    #[test]
    fn test_element_by_symbol() {
        let o = element_by_symbol("O").unwrap();
        assert_eq!(o.atomic_number, 8);
        assert_approx(o.atomic_mass, 15.999, 0.01);
    }

    #[test]
    fn test_element_not_found() {
        assert!(element_by_number(0).is_none());
        assert!(element_by_number(37).is_none());
        assert!(element_by_symbol("Xx").is_none());
    }

    #[test]
    fn test_noble_gas_no_electronegativity() {
        let he = element_by_symbol("He").unwrap();
        assert!(he.electronegativity.is_none());
    }

    #[test]
    fn test_element_iron() {
        let fe = element_by_symbol("Fe").unwrap();
        assert_eq!(fe.atomic_number, 26);
        assert_approx(fe.atomic_mass, 55.845, 0.01);
    }

    #[test]
    fn test_element_krypton() {
        let kr = element_by_number(36).unwrap();
        assert_eq!(kr.symbol, "Kr");
    }

    // --- Vec3 ---
    #[test]
    fn test_vec3_length() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert_approx(v.length(), 5.0, EPSILON);
    }

    #[test]
    fn test_vec3_distance() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        assert_approx(a.distance(b), 1.0, EPSILON);
    }

    #[test]
    fn test_vec3_dot() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_approx(a.dot(b), 32.0, EPSILON);
    }

    #[test]
    fn test_vec3_scale() {
        let v = Vec3::new(1.0, 2.0, 3.0).scale(2.0);
        assert_approx(v.x, 2.0, EPSILON);
        assert_approx(v.y, 4.0, EPSILON);
        assert_approx(v.z, 6.0, EPSILON);
    }

    #[test]
    fn test_vec3_add() {
        let v = Vec3::new(1.0, 2.0, 3.0) + Vec3::new(4.0, 5.0, 6.0);
        assert_approx(v.x, 5.0, EPSILON);
    }

    // --- Lennard-Jones ---
    #[test]
    fn test_lj_equilibrium() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        let r_eq = lj.equilibrium_distance();
        assert_approx(r_eq, 2.0_f64.powf(1.0 / 6.0), EPSILON);
    }

    #[test]
    fn test_lj_potential_at_sigma() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        assert_approx(lj.potential(1.0), 0.0, EPSILON);
    }

    #[test]
    fn test_lj_potential_at_equilibrium() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        let r_eq = lj.equilibrium_distance();
        assert_approx(lj.potential(r_eq), -1.0, EPSILON);
    }

    #[test]
    fn test_lj_force_zero_at_equilibrium() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        let r_eq = lj.equilibrium_distance();
        assert_approx(lj.force_magnitude(r_eq), 0.0, 1e-10);
    }

    #[test]
    fn test_lj_repulsive_close() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        assert!(lj.potential(0.5) > 0.0);
    }

    #[test]
    fn test_lj_force_vector() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(2.0, 0.0, 0.0);
        let f = lj.force_vector(a, b);
        assert!(f.x < 0.0);
    }

    #[test]
    fn test_lj_force_vector_zero_distance() {
        let lj = LennardJonesParams::new(1.0, 1.0);
        let a = Vec3::new(0.0, 0.0, 0.0);
        let f = lj.force_vector(a, a);
        assert_approx(f.length(), 0.0, EPSILON);
    }

    // --- Coulomb ---
    #[test]
    fn test_coulomb_potential_like_charges() {
        let v = CoulombForce::potential(1.0, 1.0, 1.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_coulomb_potential_unlike_charges() {
        let v = CoulombForce::potential(1.0, -1.0, 1.0);
        assert!(v < 0.0);
    }

    #[test]
    fn test_coulomb_force_magnitude() {
        let f = CoulombForce::force_magnitude(1e-6, 1e-6, 1.0);
        assert_approx(f, COULOMB_CONSTANT * 1e-12, 1.0);
    }

    #[test]
    fn test_coulomb_inverse_square() {
        let f1 = CoulombForce::force_magnitude(1.0, 1.0, 1.0);
        let f2 = CoulombForce::force_magnitude(1.0, 1.0, 2.0);
        assert_approx(f1 / f2, 4.0, EPSILON);
    }

    #[test]
    fn test_coulomb_force_vector_repulsion() {
        let f = CoulombForce::force_vector(
            1.0,
            1.0,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
        );
        assert!(f.x > 0.0);
    }

    #[test]
    fn test_coulomb_force_vector_attraction() {
        let f = CoulombForce::force_vector(
            1.0,
            -1.0,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
        );
        assert!(f.x < 0.0);
    }

    #[test]
    fn test_electric_field() {
        let e = CoulombForce::electric_field_magnitude(1.0, 1.0);
        assert_approx(e, COULOMB_CONSTANT, 1.0);
    }

    // --- Arrhenius ---
    #[test]
    fn test_arrhenius_rate_constant() {
        let arr = ArrheniusParams::new(1e13, 75000.0);
        let k = arr.rate_constant(300.0);
        let expected = 1e13 * (-75000.0 / (GAS_CONSTANT * 300.0)).exp();
        assert_approx(k, expected, k * 1e-10);
    }

    #[test]
    fn test_arrhenius_higher_temp_faster() {
        let arr = ArrheniusParams::new(1e13, 75000.0);
        let k1 = arr.rate_constant(300.0);
        let k2 = arr.rate_constant(400.0);
        assert!(k2 > k1);
    }

    #[test]
    fn test_arrhenius_rate_ratio() {
        let arr = ArrheniusParams::new(1e13, 50000.0);
        let ratio = arr.rate_ratio(300.0, 310.0);
        let k1 = arr.rate_constant(300.0);
        let k2 = arr.rate_constant(310.0);
        assert_approx(ratio, k2 / k1, ratio * 1e-8);
    }

    #[test]
    fn test_activation_energy_from_rates() {
        let ea_original = 50000.0;
        let arr = ArrheniusParams::new(1e13, ea_original);
        let k1 = arr.rate_constant(300.0);
        let k2 = arr.rate_constant(350.0);
        let ea_calc = ArrheniusParams::activation_energy_from_rates(k1, k2, 300.0, 350.0);
        assert_approx(ea_calc, ea_original, 0.1);
    }

    #[test]
    fn test_arrhenius_zero_ea() {
        let arr = ArrheniusParams::new(1e10, 0.0);
        assert_approx(arr.rate_constant(300.0), 1e10, 1.0);
    }

    // --- 反応次数 ---
    #[test]
    fn test_zero_order_concentration() {
        let c = ReactionOrder::Zero.concentration(1.0, 0.1, 5.0);
        assert_approx(c, 0.5, EPSILON);
    }

    #[test]
    fn test_zero_order_no_negative() {
        let c = ReactionOrder::Zero.concentration(1.0, 0.1, 20.0);
        assert_approx(c, 0.0, EPSILON);
    }

    #[test]
    fn test_first_order_concentration() {
        let c = ReactionOrder::First.concentration(1.0, 0.1, 10.0);
        let expected = (-1.0_f64).exp();
        assert_approx(c, expected, EPSILON);
    }

    #[test]
    fn test_second_order_concentration() {
        let c = ReactionOrder::Second.concentration(1.0, 0.5, 2.0);
        assert_approx(c, 0.5, EPSILON);
    }

    #[test]
    fn test_first_order_half_life() {
        let t = ReactionOrder::First.half_life(1.0, 0.1);
        assert_approx(t, core::f64::consts::LN_2 / 0.1, EPSILON);
    }

    #[test]
    fn test_zero_order_half_life() {
        let t = ReactionOrder::Zero.half_life(2.0, 0.5);
        assert_approx(t, 2.0, EPSILON);
    }

    #[test]
    fn test_second_order_half_life() {
        let t = ReactionOrder::Second.half_life(1.0, 0.5);
        assert_approx(t, 2.0, EPSILON);
    }

    // --- 化学平衡 ---
    #[test]
    fn test_gibbs_from_equilibrium() {
        let dg = gibbs_from_equilibrium(1.0, 298.15);
        assert_approx(dg, 0.0, EPSILON);
    }

    #[test]
    fn test_equilibrium_constant_zero_dg() {
        let k = equilibrium_constant(0.0, 298.15);
        assert_approx(k, 1.0, EPSILON);
    }

    #[test]
    fn test_gibbs_equilibrium_roundtrip() {
        let keq = 100.0;
        let dg = gibbs_from_equilibrium(keq, 298.15);
        let k_back = equilibrium_constant(dg, 298.15);
        assert_approx(k_back, keq, 1e-8);
    }

    #[test]
    fn test_vant_hoff_endothermic() {
        let k2 = vant_hoff(1.0, 50000.0, 300.0, 350.0);
        assert!(k2 > 1.0);
    }

    #[test]
    fn test_vant_hoff_exothermic() {
        let k2 = vant_hoff(1.0, -50000.0, 300.0, 350.0);
        assert!(k2 < 1.0);
    }

    #[test]
    fn test_reaction_quotient() {
        let q = reaction_quotient(&[(2.0, 2.0)], &[(1.0, 1.0)]);
        assert_approx(q, 4.0, EPSILON);
    }

    #[test]
    fn test_predict_direction_forward() {
        assert_eq!(predict_direction(0.1, 1.0), ReactionDirection::Forward);
    }

    #[test]
    fn test_predict_direction_reverse() {
        assert_eq!(predict_direction(10.0, 1.0), ReactionDirection::Reverse);
    }

    #[test]
    fn test_predict_direction_equilibrium() {
        assert_eq!(predict_direction(1.0, 1.0), ReactionDirection::Equilibrium);
    }

    // --- 結合エネルギー ---
    #[test]
    fn test_bond_energy_ch() {
        assert_approx(
            bond_energy("C", "H", BondType::Single).unwrap(),
            413.0,
            EPSILON,
        );
    }

    #[test]
    fn test_bond_energy_cc_double() {
        assert_approx(
            bond_energy("C", "C", BondType::Double).unwrap(),
            614.0,
            EPSILON,
        );
    }

    #[test]
    fn test_bond_energy_nn_triple() {
        assert_approx(
            bond_energy("N", "N", BondType::Triple).unwrap(),
            941.0,
            EPSILON,
        );
    }

    #[test]
    fn test_bond_energy_unknown() {
        assert!(bond_energy("X", "Y", BondType::Single).is_none());
    }

    #[test]
    fn test_bond_energy_symmetric() {
        let e1 = bond_energy("C", "H", BondType::Single);
        let e2 = bond_energy("H", "C", BondType::Single);
        assert_eq!(e1, e2);
    }

    #[test]
    fn test_reaction_enthalpy_from_bonds() {
        let dh = reaction_enthalpy_from_bonds(678.0, 862.0);
        assert_approx(dh, -184.0, EPSILON);
    }

    // --- 分子構造 ---
    #[test]
    fn test_molecule_weight() {
        let mut mol = Molecule::new();
        mol.add_atom("O", Vec3::new(0.0, 0.0, 0.0), -0.82, 15.999);
        mol.add_atom("H", Vec3::new(0.96, 0.0, 0.0), 0.41, 1.008);
        mol.add_atom("H", Vec3::new(-0.24, 0.93, 0.0), 0.41, 1.008);
        assert_approx(mol.molecular_weight(), 18.015, 0.01);
    }

    #[test]
    fn test_molecule_center_of_mass() {
        let mut mol = Molecule::new();
        mol.add_atom("C", Vec3::new(0.0, 0.0, 0.0), 0.0, 12.0);
        mol.add_atom("C", Vec3::new(2.0, 0.0, 0.0), 0.0, 12.0);
        let com = mol.center_of_mass();
        assert_approx(com.x, 1.0, EPSILON);
    }

    #[test]
    fn test_molecule_bond_length() {
        let mut mol = Molecule::new();
        mol.add_atom("H", Vec3::new(0.0, 0.0, 0.0), 0.0, 1.008);
        mol.add_atom("H", Vec3::new(0.74, 0.0, 0.0), 0.0, 1.008);
        mol.add_bond(0, 1, BondType::Single);
        assert_approx(mol.bond_length(0).unwrap(), 0.74, EPSILON);
    }

    #[test]
    fn test_molecule_bond_angle() {
        let mut mol = Molecule::new();
        mol.add_atom("H", Vec3::new(1.0, 0.0, 0.0), 0.0, 1.008);
        mol.add_atom("O", Vec3::new(0.0, 0.0, 0.0), 0.0, 15.999);
        mol.add_atom("H", Vec3::new(0.0, 1.0, 0.0), 0.0, 1.008);
        let angle = mol.bond_angle(0, 1, 2);
        assert_approx(angle, PI / 2.0, EPSILON);
    }

    #[test]
    fn test_molecule_total_bond_energy() {
        let mut mol = Molecule::new();
        mol.add_atom("H", Vec3::new(0.0, 0.0, 0.0), 0.0, 1.008);
        mol.add_atom("H", Vec3::new(0.74, 0.0, 0.0), 0.0, 1.008);
        mol.add_bond(0, 1, BondType::Single);
        assert_approx(mol.total_bond_energy(), 436.0, EPSILON);
    }

    #[test]
    fn test_molecule_default() {
        let mol = Molecule::default();
        assert!(mol.atoms.is_empty());
        assert!(mol.bonds.is_empty());
    }

    // --- 熱力学 ---
    #[test]
    fn test_gibbs_free_energy() {
        let dg = gibbs_free_energy(-100.0, 298.15, 0.2);
        assert_approx(dg, -100.0 - 298.15 * 0.2, EPSILON);
    }

    #[test]
    fn test_is_spontaneous() {
        assert!(is_spontaneous(-10.0));
        assert!(!is_spontaneous(10.0));
        assert!(!is_spontaneous(0.0));
    }

    #[test]
    fn test_ideal_gas_volume() {
        let v = ideal_gas_volume(1.0, 273.15, 101325.0);
        assert_approx(v, 0.02241, 0.001);
    }

    #[test]
    fn test_ideal_gas_pressure() {
        let p = ideal_gas_pressure(1.0, 273.15, 0.02241);
        assert_approx(p, 101325.0, 500.0);
    }

    #[test]
    fn test_entropy_change_reversible() {
        let ds = entropy_change_reversible(1000.0, 500.0);
        assert_approx(ds, 2.0, EPSILON);
    }

    #[test]
    fn test_entropy_isothermal_expansion() {
        let ds = entropy_isothermal_expansion(1.0, 1.0, 2.0);
        assert_approx(ds, GAS_CONSTANT * 2.0_f64.ln(), EPSILON);
    }

    #[test]
    fn test_entropy_temperature_change() {
        let ds = entropy_temperature_change(1.0, 29.1, 300.0, 600.0);
        assert_approx(ds, 29.1 * (2.0_f64).ln(), 0.01);
    }

    #[test]
    fn test_reaction_enthalpy() {
        let dh = reaction_enthalpy(&[(0.0, 2.0), (0.0, 1.0)], &[(-285.8, 2.0)]);
        assert_approx(dh, 571.6, EPSILON);
    }

    #[test]
    fn test_hess_law() {
        let total = hess_law(&[-100.0, 50.0, -30.0]);
        assert_approx(total, -80.0, EPSILON);
    }

    #[test]
    fn test_kirchhoff() {
        let dh2 = kirchhoff_enthalpy(-100.0, 10.0, 298.0, 398.0);
        assert_approx(dh2, -100.0 + 10.0 * 100.0, EPSILON);
    }

    #[test]
    fn test_carnot_efficiency() {
        let eff = carnot_efficiency(500.0, 300.0);
        assert_approx(eff, 0.4, EPSILON);
    }

    #[test]
    fn test_internal_energy_change() {
        assert_approx(internal_energy_change(100.0, -30.0), 70.0, EPSILON);
    }

    // --- 化学量論 ---
    #[test]
    fn test_moles_to_grams() {
        assert_approx(moles_to_grams(2.0, 18.015), 36.03, 0.01);
    }

    #[test]
    fn test_grams_to_moles() {
        assert_approx(grams_to_moles(36.03, 18.015), 2.0, 0.001);
    }

    #[test]
    fn test_moles_to_molecules() {
        let n = moles_to_molecules(1.0);
        assert_approx(n, AVOGADRO, 1e16);
    }

    #[test]
    fn test_molecules_to_moles() {
        assert_approx(molecules_to_moles(AVOGADRO), 1.0, EPSILON);
    }

    #[test]
    fn test_molarity() {
        assert_approx(molarity(0.5, 0.25), 2.0, EPSILON);
    }

    #[test]
    fn test_dilution() {
        let v2 = dilution_volume(2.0, 0.5, 0.5);
        assert_approx(v2, 2.0, EPSILON);
    }

    #[test]
    fn test_limiting_reagent() {
        let idx = limiting_reagent(&[(3.0, 2.0), (2.0, 1.0)]);
        assert_eq!(idx, 0);
    }

    #[test]
    fn test_limiting_reagent_second() {
        let idx = limiting_reagent(&[(10.0, 1.0), (1.0, 2.0)]);
        assert_eq!(idx, 1);
    }

    #[test]
    fn test_theoretical_yield() {
        assert_approx(theoretical_yield(3.0, 2.0, 2.0), 3.0, EPSILON);
    }

    #[test]
    fn test_percent_yield() {
        assert_approx(percent_yield(2.5, 3.0), 83.333_333, 0.001);
    }

    #[test]
    fn test_molecular_formula_multiplier() {
        assert_eq!(molecular_formula_multiplier(30.0, 180.0), 6);
    }

    #[test]
    fn test_molecular_formula_multiplier_one() {
        assert_eq!(molecular_formula_multiplier(44.0, 44.0), 1);
    }

    // --- MD シミュレーション ---
    #[test]
    fn test_kinetic_energy() {
        let p = Particle::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), 2.0, 0.0);
        assert_approx(kinetic_energy(&[p]), 1.0, EPSILON);
    }

    #[test]
    fn test_system_temperature_single() {
        let p = Particle::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(100.0, 0.0, 0.0),
            1e-26,
            0.0,
        );
        let t = system_temperature(&[p]);
        assert!(t > 0.0);
    }

    #[test]
    fn test_velocity_verlet_conservation() {
        let lj = LennardJonesParams::new(1e-21, 3.4e-10);
        let mut particles = vec![
            Particle::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(10.0, 0.0, 0.0),
                6.63e-26,
                0.0,
            ),
            Particle::new(
                Vec3::new(1e-9, 0.0, 0.0),
                Vec3::new(-10.0, 0.0, 0.0),
                6.63e-26,
                0.0,
            ),
        ];
        let dt = 1e-15;
        for _ in 0..10 {
            velocity_verlet_step(&mut particles, dt, &lj);
        }
        assert!((particles[0].position.x - 0.0).abs() > 1e-20);
    }

    #[test]
    fn test_most_probable_speed() {
        let v = most_probable_speed(6.63e-26, 300.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_mean_speed() {
        let v = mean_speed(6.63e-26, 300.0);
        assert!(v > most_probable_speed(6.63e-26, 300.0));
    }

    #[test]
    fn test_rms_speed() {
        let v = rms_speed(6.63e-26, 300.0);
        assert!(v > mean_speed(6.63e-26, 300.0));
    }

    #[test]
    fn test_speed_ordering() {
        let m = 6.63e-26;
        let t = 300.0;
        assert!(most_probable_speed(m, t) < mean_speed(m, t));
        assert!(mean_speed(m, t) < rms_speed(m, t));
    }

    // --- pH ---
    #[test]
    fn test_ph_neutral() {
        assert_approx(ph(1e-7), 7.0, EPSILON);
    }

    #[test]
    fn test_ph_acidic() {
        assert_approx(ph(1e-3), 3.0, EPSILON);
    }

    #[test]
    fn test_poh() {
        assert_approx(poh(1e-7), 7.0, EPSILON);
    }

    #[test]
    fn test_h_from_ph() {
        assert_approx(h_concentration_from_ph(7.0), 1e-7, 1e-12);
    }

    #[test]
    fn test_henderson_hasselbalch() {
        let result = henderson_hasselbalch(4.75, 1.0, 1.0);
        assert_approx(result, 4.75, EPSILON);
    }

    #[test]
    fn test_henderson_hasselbalch_excess_base() {
        let result = henderson_hasselbalch(4.75, 10.0, 1.0);
        assert_approx(result, 5.75, EPSILON);
    }

    // --- 定数の検証 ---
    #[test]
    fn test_gas_constant() {
        assert_approx(GAS_CONSTANT, BOLTZMANN * AVOGADRO, 0.01);
    }

    #[test]
    fn test_constants_positive() {
        assert!(BOLTZMANN > 0.0);
        assert!(AVOGADRO > 0.0);
        assert!(COULOMB_CONSTANT > 0.0);
        assert!(ELEMENTARY_CHARGE > 0.0);
        assert!(VACUUM_PERMITTIVITY > 0.0);
    }

    // --- 追加テスト ---
    #[test]
    fn test_lj_symmetric_potential() {
        let lj = LennardJonesParams::new(0.01, 3.4);
        let v1 = lj.potential(4.0);
        let v2 = lj.potential(4.0);
        assert_approx(v1, v2, EPSILON);
    }

    #[test]
    fn test_coulomb_symmetry() {
        let v1 = CoulombForce::potential(1.0, -2.0, 3.0);
        let v2 = CoulombForce::potential(-2.0, 1.0, 3.0);
        assert_approx(v1, v2, EPSILON);
    }

    #[test]
    fn test_equilibrium_large_k() {
        let dg = gibbs_from_equilibrium(1e10, 298.15);
        assert!(dg < 0.0);
    }

    #[test]
    fn test_equilibrium_small_k() {
        let dg = gibbs_from_equilibrium(1e-10, 298.15);
        assert!(dg > 0.0);
    }

    #[test]
    fn test_vant_hoff_same_temp() {
        let k2 = vant_hoff(5.0, 50000.0, 300.0, 300.0);
        assert_approx(k2, 5.0, EPSILON);
    }

    #[test]
    fn test_bond_energy_hf() {
        assert_approx(
            bond_energy("H", "F", BondType::Single).unwrap(),
            567.0,
            EPSILON,
        );
    }

    #[test]
    fn test_bond_energy_cs_double() {
        assert_approx(
            bond_energy("C", "S", BondType::Double).unwrap(),
            573.0,
            EPSILON,
        );
    }

    #[test]
    fn test_ph_poh_sum() {
        let ph_val = 4.0;
        let h = h_concentration_from_ph(ph_val);
        let kw = 1e-14;
        let oh = kw / h;
        let poh_val = poh(oh);
        assert_approx(ph_val + poh_val, 14.0, 0.01);
    }

    #[test]
    fn test_carnot_same_temp() {
        assert_approx(carnot_efficiency(300.0, 300.0), 0.0, EPSILON);
    }

    #[test]
    fn test_ideal_gas_double_moles() {
        let v1 = ideal_gas_volume(1.0, 300.0, 101325.0);
        let v2 = ideal_gas_volume(2.0, 300.0, 101325.0);
        assert_approx(v2 / v1, 2.0, EPSILON);
    }

    #[test]
    fn test_zero_order_half_exhaustion() {
        let t_half = ReactionOrder::Zero.half_life(1.0, 0.1);
        let c = ReactionOrder::Zero.concentration(1.0, 0.1, 2.0 * t_half);
        assert_approx(c, 0.0, EPSILON);
    }

    #[test]
    fn test_first_order_independent_of_initial() {
        let t1 = ReactionOrder::First.half_life(1.0, 0.5);
        let t2 = ReactionOrder::First.half_life(100.0, 0.5);
        assert_approx(t1, t2, EPSILON);
    }

    #[test]
    fn test_molecule_multi_bond() {
        let mut mol = Molecule::new();
        let c = mol.add_atom("C", Vec3::new(0.0, 0.0, 0.0), 0.0, 12.011);
        let h1 = mol.add_atom("H", Vec3::new(1.0, 0.0, 0.0), 0.0, 1.008);
        let h2 = mol.add_atom("H", Vec3::new(0.0, 1.0, 0.0), 0.0, 1.008);
        let h3 = mol.add_atom("H", Vec3::new(0.0, 0.0, 1.0), 0.0, 1.008);
        let h4 = mol.add_atom("H", Vec3::new(-1.0, 0.0, 0.0), 0.0, 1.008);
        mol.add_bond(c, h1, BondType::Single);
        mol.add_bond(c, h2, BondType::Single);
        mol.add_bond(c, h3, BondType::Single);
        mol.add_bond(c, h4, BondType::Single);
        assert_approx(mol.total_bond_energy(), 4.0 * 413.0, EPSILON);
        assert_approx(mol.molecular_weight(), 16.043, 0.01);
    }

    #[test]
    fn test_particle_creation() {
        let p = Particle::new(Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.0, 0.0, 0.0), 1.0, 0.5);
        assert_approx(p.position.x, 1.0, EPSILON);
        assert_approx(p.charge, 0.5, EPSILON);
    }

    #[test]
    fn test_system_temperature_empty() {
        let t = system_temperature(&[]);
        assert_approx(t, 0.0, EPSILON);
    }

    #[test]
    fn test_reaction_quotient_zero_reactant() {
        let q = reaction_quotient(&[(1.0, 1.0)], &[(0.0, 1.0)]);
        assert!(q.is_infinite());
    }
}
