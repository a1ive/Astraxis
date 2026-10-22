# 地磁数据来源

| 文件 | 来源 | 许可 |
|---|---|---|
| `kp.bin` | GFZ Helmholtz Centre for Geosciences，行星际地磁指数 Kp（每 3 小时一个值，1932-01-01 起）。Matzka, J., Stolle, C., Yamazaki, Y., Bronkalla, O. and Morschhauser, A., 2021, The geomagnetic Kp index and derived indices of geomagnetic activity, Space Weather 19, e2020SW002641；数据集 Matzka et al. 2021, Geomagnetic Kp index, V. 1.0, GFZ Data Services, <https://doi.org/10.5880/Kp.0001> | CC BY 4.0（请注明 GFZ Helmholtz Centre for Geosciences） |

## 处理方式

`tools/geomag/kp_bake.py` 读 GFZ 的年度 WDC 格式文件（`kp1932.wdc` .. `kp2026.wdc`），只保留 Kp（以 1/3 为单位，0–27，每个三小时区间一个字节），从第一年的 1 月 1 日 00:00 UT 起连续排列，到最后一个文件的最后一天为止。WDC 文件里的太阳黑子数（SILSO，CC BY-NC 4.0）和 10.7 cm 射电流量不保留。

```
python -I tools/geomag/kp_bake.py <放 kpYYYY.wdc 的目录> assets/geomag/kp.bin
```

下载：GFZ 自己的服务器（`kp.gfz.de`、`ftp.gfz-potsdam.de`）在本机 TLS 握手失败，用 GWDG 的镜像 <https://ftp.gwdg.de/pub/geophys/kp-ap/wdc/yearly/>（格式说明在同目录的 `wdc_fmt.txt`）。当前文件到 2026-09-30 为止（2026-10-22 下载，没有缺值），更新时重新下载全部年份并重跑工具：GFZ 会修订最近的值。

## 程序里的其他地磁、极光常数

都写在 `src/scene/aurora.cpp`，注明了出处：

- IGRF-14 的偶极项 g10、g11、h11（1900–2025 每 5 年，加 2025–30 的长期变化）：IAGA Working Group V-MOD，<https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt>
- 极光卵边界：Starkov 1994 的系数，取自 Sigernes et al. 2011, Opt. Pura Apl. 44, 599, Table 1（开放获取，<https://folk.uib.no/st07350/papers/sigernes2011OPA599.pdf>）。表里 poleward 边界 A2 的 b0 印成“-4-44”，按 −4.44 用
- 各谱线的颜色（每 kR 的线性 RGB）：`tools/geomag/aurora_colors.py` 用 CIE 1931 2° 配色函数（CVRL，<http://cvrl.ucl.ac.uk/database/data/cmfs/ciexyz31_1.csv>，不入库）算出，以 5772 K 黑体代表的阳光为白色
