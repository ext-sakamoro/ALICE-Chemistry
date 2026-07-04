//! 分子構造 (`Molecule` / `Bond` 等)

use crate::force_field::{bond_energy, BondType};
use crate::vec3::Vec3;

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
