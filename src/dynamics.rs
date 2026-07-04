//! 分子動力学シミュレーション

use crate::constants::BOLTZMANN;
use crate::force_field::LennardJonesParams;
use crate::vec3::Vec3;
use core::f64::consts::PI;

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
