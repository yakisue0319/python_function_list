# MetPy 関数・モジュール一覧

## 🔹 基本情報
| 名前 | 概要 |
|------|------|
| `units` | 単位の扱い（Pintベース） |
| `units.Quantity` | 単位付きのデータオブジェクト |
| `units.units()` | 単位の指定（例：`units('degC')`） |

## 🔹 熱力学・大気安定性
| 名前 | 概要 |
|------|------|
| `dewpoint_from_relative_humidity()` | 相対湿度から露点を計算 |
| `saturation_vapor_pressure()` | 飽和水蒸気圧 |
| `relative_humidity_from_dewpoint()` | 露点から相対湿度を計算 |
| `potential_temperature()` | ポテンシャル温位 |
| `equivalent_potential_temperature()` | 相当温位 |
| `virtual_temperature()` | 仮温度 |

## 🔹 降水・凝結関連
| 名前 | 概要 |
|------|------|
| `lifted_index()` | LI（持ち上げ指数）の計算 |
| `parcel_profile()` | 気塊上昇時の温度プロファイル |
| `lcl_pressure()`, `lcl_temperature()` | LCL（凝結高度）の気圧・気温 |
| `lcl()` | LCL気圧と気温を同時に計算 |

## 🔹 風・運動量関連
| 名前 | 概要 |
|------|------|
| `wind_components()` | 風速と風向からu,v成分を計算 |
| `wind_speed()` | u,v成分から風速を計算 |
| `wind_direction()` | u,v成分から風向を計算 |
| `shear()` | 風の鉛直シア計算 |
| `bulk_richardson_number()` | バルクリチャードソン数（対流判定） |

## 🔹 成層・安定度解析
| 名前 | 概要 |
|------|------|
| `brunt_vaisala_frequency()` | ブリュン・ヴァイサラ周波数 |
| `static_stability()` | 静的安定度の計算 |
| `moist_static_energy()` | 湿静的エネルギー |
| `cape_cin()` | CAPE・CINの同時計算 |

## 🔹 その他
| 名前 | 概要 |
|------|------|
| `interpolate_1d()` | 高度や気圧に沿った補間 |
| `smooth_n_point()`, `smooth_window()` | データの平滑化 |
| `geostrophic_wind()` | 地衡風の計算 |
| `absolute_vorticity()` | 絶対渦度 |

