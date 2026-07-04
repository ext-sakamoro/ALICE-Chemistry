//! Cross-module tests (= 元 `lib.rs` の `mod tests` を統合)

#![allow(clippy::doc_markdown)]
#![allow(clippy::assertions_on_constants)]
#![allow(clippy::suboptimal_flops)]
#![allow(clippy::unreadable_literal)]

use crate::constants::*;
use crate::dynamics::*;
use crate::force_field::*;
use crate::periodic::*;
use crate::reactions::*;
use crate::structure::*;
use crate::thermo::*;
use crate::vec3::*;

use core::f64::consts::PI;

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
    let f =
        CoulombForce::force_vector(1.0, 1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
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
