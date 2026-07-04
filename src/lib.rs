//! ALICE-Chemistry: 分子動力学・化学反応シミュレーションライブラリ
//!
//! - 周期表データ (元素情報)
//! - 力場 (Lennard-Jones, Coulomb, 結合エネルギー)
//! - 反応速度論 (Arrhenius 式) + 化学平衡 + 化学量論
//! - 分子構造
//! - 熱力学 (エンタルピー / エントロピー / ギブズ自由エネルギー) + pH
//! - 分子動力学シミュレーション
//!
//! # Module 構成
//!
//! | Module | 内容 |
//! |--------|------|
//! | [`constants`] | 物理・化学定数 |
//! | [`periodic`] | 周期表 (元素データ) |
//! | [`vec3`] | 3D ベクトル primitive |
//! | [`force_field`] | Lennard-Jones + Coulomb + 結合エネルギー |
//! | [`reactions`] | Arrhenius + 化学平衡 + 化学量論 |
//! | [`structure`] | 分子構造 |
//! | [`thermo`] | 熱力学 + pH |
//! | [`dynamics`] | 分子動力学 |
//! | [`prelude`] | 主要 API 一括 re-export |
//!
//! # Backward compatibility
//!
//! v1.0.0 まで crate ルート直下に定義していた項目は module 移動後も
//! `pub use` re-export で crate root から到達可能

#![warn(clippy::all, clippy::pedantic, clippy::nursery)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::similar_names)]
#![allow(clippy::doc_markdown)]
#![allow(clippy::cast_precision_loss)]

pub mod constants;
pub mod dynamics;
pub mod force_field;
pub mod periodic;
pub mod prelude;
pub mod reactions;
pub mod structure;
pub mod thermo;
pub mod vec3;

#[cfg(test)]
mod integration_tests;

// Backward-compatible re-exports
pub use crate::constants::*;
pub use crate::dynamics::*;
pub use crate::force_field::*;
pub use crate::periodic::*;
pub use crate::reactions::*;
pub use crate::structure::*;
pub use crate::thermo::*;
pub use crate::vec3::*;
