# 3DZH-Inversion

基于环境噪声 Rayleigh 波 **ZH 振幅比** 的三维 S 波速度（Vs）反演程序。  
项目使用 Fortran 实现，支持真实数据反演与合成数据检验（checkerboard/恢复测试）。


## 主要功能

- 支持两种模式：
  - `synthetic_flag=0`：真实观测反演
  - `synthetic_flag=1`：合成数据测试（由 `MOD.true` 生成）
- 自适应平滑半径（按台站空间分布计算）
- Tikhonov 正则化（0/1/2 阶，支持方向控制）
- LSMR 稀疏反演求解
- 每次迭代输出模型与日志（RMS、方程规模、非零元数量等）

## 项目结构

```text
.
├── src/                 # Fortran 源码
├── makefile             # 编译入口
├── bin/                 # 可执行文件（编译后生成 3DZHTomo）
├── obj/                 # 中间编译产物
├── datasets_Yunnan/     # 示例数据与示例参数
│   ├── 3DZHTomo.in
│   ├── MOD
│   ├── MOD.true
│   └── ZHdata/
└── tools/               # 绘图/分析 notebook
