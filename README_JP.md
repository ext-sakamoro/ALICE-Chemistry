[English](README.md) | **日本語**

# ALICE-Chemistry

[Project A.L.I.C.E.](https://github.com/anthropics/alice) の分子動力学・化学シミュレーションライブラリ

## 概要

`alice-chemistry` は純Rustによる化学シミュレーションツールキットです。分子動力学、反応速度論、熱力学、周期表データを提供します。

## 機能

- **周期表** — H〜Krの元素データ（原子番号、質量、電気陰性度）
- **力場** — レナード・ジョーンズポテンシャルとクーロンポテンシャル
- **反応速度論** — アレニウス式による速度定数計算
- **化学平衡** — 平衡定数と濃度の計算
- **分子構造** — 結合エネルギーと分子幾何学
- **熱力学** — エンタルピー、エントロピー、ギブズ自由エネルギー計算
- **化学量論** — 平衡反応係数の取り扱い
- **物理定数** — ボルツマン定数、アボガドロ数、気体定数、クーロン定数など

## クイックスタート

```rust
use alice_chemistry::{lookup_element, BOLTZMANN, GAS_CONSTANT};

let hydrogen = lookup_element(1).unwrap();
assert_eq!(hydrogen.symbol, "H");
assert!((hydrogen.atomic_mass - 1.008).abs() < 0.001);

// アレニウス速度定数: k = A * exp(-Ea / (R * T))
let k = 1e13_f64 * (-75000.0 / (GAS_CONSTANT * 298.15)).exp();
```

## アーキテクチャ

```
alice-chemistry
├── Element           — 元素データ（原子番号、質量、電気陰性度）
├── ELEMENTS          — 周期表データ（H-Kr）
├── lennard_jones     — LJポテンシャルと力の計算
├── coulomb           — 静電ポテンシャル
├── arrhenius         — 反応速度定数
├── equilibrium       — 化学平衡ソルバー
├── thermodynamics    — エンタルピー、エントロピー、ギブズエネルギー
└── stoichiometry     — 反応の平衡化とモル比
```

## ライセンス

MIT OR Apache-2.0
