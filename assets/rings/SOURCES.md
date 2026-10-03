# 行星环数据来源

## `saturn_rss.ring`：土星环径向光学深度剖面

| 项目 | 内容 |
|---|---|
| 来源 | NASA Planetary Data System, Ring-Moon Systems Node，卷 CORSS_8001（Cassini Radio Science 环掩星，E. Marouf 等），`RSS_2005_123_X43_E_TAU_10KM.TAB` |
| 下载地址 | <https://pds-rings.seti.org/holdings/volumes/CORSS_8xxx/CORSS_8001/data/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_TAU_10KM.TAB>（5 MB） |
| 观测 | 第 7 圈出射段，X 波段（3.6 cm），DSN 43 站，2005-05-03；环开口角约 −23.6°；径向分辨率 10 km，采样 2.5 km，72,770–144,988 km |
| 许可 | NASA PDS 数据，公有领域 |
| 处理 | `tools/rings/make_saturn_rings.py`：半径加上两项改正（改进的极轴、时间偏移）；信号丢失（功率 ≤ 0）或光学深度达到探测上限（约 5.06，B 环核心）处取上限值；负值噪声取 0；平均到 8192 个等宽区间（每个 8.82 km） |

文件格式（小端）：`char[8] "AXRING1\0"`，`float64` 内、外半径（km），`uint32` 样本数，`uint32` 保留，随后是 `float32` 法向光学深度。

处理后的几个特征值：C 环（80,000 km）0.175，B 环核心（110,000 km）3.87，卡西尼缝（119,000 km）0.16，A 环（130,000 km）0.63，恩克缝（133,580 km）0.001（`tests/sim_tests.cpp` 会检查）。

局限：
- 3.6 cm 无线电波测的是厘米到米级颗粒，可见光下还有更多细尘，光学深度会略有差别
- 不含 D 环（< 72,770 km，光学深度约 10⁻³）以及 E、G 等外侧尘埃环
- 这是一次掩星的单一经度剖面，A 环的方位不对称、F 环的扭结、B 环外缘随经度的变化都没有

## 木星环

没有数据文件，各环带直接写在 `assets/scenes/jupiter.toml` 里（来源：PDS Ring-Moon Systems Node 的木星环与内卫星统计表），加载时栅格化成同样的剖面（4096 个样本）。
