# scTenifoldKnk 龙芯 SDAA 适配与验证报告

## 1. 结论摘要

scTenifoldKnk 2.1.0 已在 LoongArch CPU 环境跑通官方最小推理示例，固定源码 commit 为
`08118337f6fc22cca5a22c443ccb78ef920dfa3b`。官方 R/C++ 算法源码未修改。

当前验证结论：

- 官方实现使用 R、RcppEigen、RSpectra 与 BLAS/LAPACK 执行 CPU 计算；源码没有 CUDA、Torch 或 SDAA 调用。
- 本文将已有 x86 CPU 与 LoongArch CPU 对照单独列出。
- CPU 跨架构输出最大绝对误差低于 `1e-6`，显著基因集合相同。

## 2. 软件环境

| 环境 | 架构与系统 | R / scTenifoldKnk | BLAS/LAPACK | 设备与实际后端 |
| --- | --- | --- | --- | --- |
| x86 基线 | x86_64，Ubuntu 22.04 | R 4.4.3 / 2.1.0 | OpenBLAS 0.3.34 / LAPACK 3.12.0 | 容器未申请 CUDA 设备；实际运行 CPU 路径 |
| 龙芯基线 | loongarch64，Loongnix 23.1 | R 4.4.1 / 2.1.0 | FlexiBLAS/OpenBLAS / LAPACK 3.11.0 | 容器未映射 SDAA 设备；实际运行 CPU 路径 |

两端使用相同的 R 包版本：Matrix 1.7-6、MASS 7.3-66、Rcpp 1.0.14、RcppEigen 0.3.4.0.2、RSpectra 0.16-2。
官方示例为 `100 × 3000` 输入、目标基因 `G100`、随机种子 1。

## 3. 源码适配点

### 3.1 计算后端

官方实现没有 GPU 后端，也没有模型权重。LoongArch 适配是在目标系统编译 R/C++ 依赖，并运行原有 CPU 算法；没有添加伪造的 CUDA/SDAA 设备选择，也没有修改算法。

### 3.2 第三方依赖与 wheel

包的运行依赖在 `DESCRIPTION` 中声明：Matrix、RSpectra、MASS、Rcpp、RcppEigen，以及 R 自带的 stats、parallel、methods、utils。它们是 R 包，不是 Python 包，不能用 `.whl` 替代。对应固定版本源码包放在 `vendor/`，并由 `vendor/checksums.sha256` 校验。

独立 Python venv 中没有该算法需要的第三方 Python 包，因此没有生成 wheel。若后续增加 Python/SDAA 桥接依赖，再针对 Python 3.12 与 LoongArch 构建 wheel 并锁定其哈希。

### 3.3 镜像与环境锁

基础镜像只有本机 tag，没有 Registry RepoDigest；构建脚本因此固定并检查本机镜像 ID。R 依赖、系统 RPM、Python 与最终镜像信息记录在 `source.lock` 和 `locks/` 中。

## 4. 适配过程中遇到的问题

| 问题/报错 | 原因 | 处理方法 |
| --- | --- | --- |
| LoongArch 缺少可直接使用的 R 二进制扩展 | 目标架构的预编译包有限 | 在目标镜像中编译固定版本的 R 包源码，并记录 SHA-256 |
| `R-core-devel` 依赖 `redhat-rpm-config`，DNF 无法解析 | Loongnix 仓库提供 `system-rpm-config`，没有该 capability | DNF 安装实际依赖后，校验固定 RPM 哈希，再对该 RPM 使用 `rpm --nodeps` |
| EPLL RPM 缺少可用的 GPG 校验链，仓库配置 `gpgcheck=0` | RPM 签名密钥不在目标系统信任库中 | 使用官方 HTTPS 来源并固定 RPM SHA-256，来源记录在 `source.lock` |
## 5. 精度对比


### 5.1 x86 CPU 与 LoongArch CPU

两端使用相同源码 commit、R 包版本、官方输入、随机种子与算法参数。误差表保存在 `../work/loong64-smoke-test/x86-comparison.csv`。

| 输出字段 | 最大绝对误差 |
| --- | ---: |
| WT 网络矩阵 | 0 |
| KO 网络矩阵 | 0 |
| 流形对齐矩阵 | `1.77e-14` |
| distance | `4.41e-16` |
| Z | `1.29e-9` |
| FC | `1.10e-7` |
| p.value | `3.67e-11` |
| p.adj | `8.23e-12` |

显著基因集合相同；FC 和 p.adj 排名 Spearman 相关系数均为 1。最大误差低于 `1e-6`。
## 6. 性能对比



### 6.1 x86 CPU 与 LoongArch CPU

计时窗口只覆盖 `scTenifoldKnk()`，不包括容器启动与 R 进程启动；两端均固定单线程。此前每端先预热一次，再记录 10 次正式运行：

| 平台/后端 | 10 次平均 | 中位数 | 范围 |
| --- | ---: | ---: | ---: |
| x86_64 Xeon Platinum 8468 / CPU | `0.2275 s` | `0.2290 s` | `0.214–0.243 s` |
| Loongson-3C6000/D / CPU | `0.6224 s` | `0.5965 s` | `0.568–0.689 s` |

该批次龙芯 CPU 平均耗时约为 x86 CPU 的 `2.74x`。之后的一次重跑分别记录为 x86 CPU `0.204 s`、LoongArch CPU `0.568 s`；目录迁移后的 LoongArch CPU 验证耗时为 `0.698 s`。单次时间受服务器负载影响。
## 7. 复现命令

从 `scTenifoldKnk` 仓库根目录执行：

```bash
./sdaa/scripts/build.sh
./sdaa/scripts/export-locks.sh
./sdaa/scripts/run-smoke-test.sh
```

源码与环境锁位于 `source.lock`、`vendor/` 和 `locks/`；示例输出位于仓库上层的 `../work/loong64-smoke-test/`。构建、运行细节见同目录的 `README.md`。
