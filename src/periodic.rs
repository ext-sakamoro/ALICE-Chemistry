//! 周期表 (元素データ)

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
