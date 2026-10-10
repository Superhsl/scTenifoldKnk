# 龙芯 LoongArch / SDAA 复现教程

本文从已上传的 scTenifoldKnk 仓库开始，在龙芯服务器构建容器并直接运行官方最小示例。当前 scTenifoldKnk 使用 R/C++ 与 BLAS/LAPACK 执行 CPU 计算，没有 CUDA 或 SDAA 后端；服务器装有 SDAA 卡不代表本程序会调用它。

## 1. 前置条件

- 服务器架构为 `loongarch64`，已安装 Docker 并可访问 Docker daemon。
- 本机已加载基础镜像 `lasermpnn-sdaa321-loong:torch2.12.0-sdaa3.2.1-py312`，镜像 ID 为 `sha256:3e425da47a74a03b589f8e1c18a67dcb90c8bcbde977cc1df09fa01209b7a524`。
- Docker build 阶段需要 Loongnix 软件源安装 R 和系统编译依赖；固定版本 R 源码包及 RPM 已随仓库提供。
- 上传仓库时保留官方源码 commit `08118337f6fc22cca5a22c443ccb78ef920dfa3b` 的 Git 历史。构建脚本会确认该版本存在，并验证其源码未在 `sdaa/` 之外被修改。

本流程不需要模型权重、CUDA/Torch-SDAA 或 SDAA 设备映射。R 包依赖不是 Python 包，当前没有需要制作的 wheel。

## 2. 获取仓库并检查环境

把仓库地址替换为上传后的地址：

```bash
REPOSITORY_URL='<已上传仓库地址>'
git clone --branch adapt/sdaa "$REPOSITORY_URL" scTenifoldKnk
cd scTenifoldKnk

test "$(uname -m)" = loongarch64
git cat-file -e '08118337f6fc22cca5a22c443ccb78ef920dfa3b^{commit}'
test "$(docker image inspect lasermpnn-sdaa321-loong:torch2.12.0-sdaa3.2.1-py312 --format '{{.Id}}')" = \
  sha256:3e425da47a74a03b589f8e1c18a67dcb90c8bcbde977cc1df09fa01209b7a524
```

若官方 commit 不存在，检查克隆是否为浅克隆或上传时是否丢失历史；需要补全历史，不能 checkout 回官方 commit 来代替，因为那会切掉新提交的 `sdaa/` 文件。

## 3. 校验锁文件并构建

在仓库根目录执行：

```bash
(cd sdaa/vendor && sha256sum --check checksums.sha256)
./sdaa/scripts/build.sh
```

构建脚本会检查架构、官方源码基线、基础镜像 ID 和 vendor 哈希，然后生成 `sctenifoldknk-loong64:r2.1.0-cpu-adapt-sdaa`。镜像安装 R 4.4.1、固定版本的 Matrix/MASS/Rcpp/RcppEigen/RSpectra，并创建独立 Python 3.12 venv。

## 4. 运行官方最小示例

在仓库根目录运行：

```bash
./sdaa/scripts/run-smoke-test.sh
```

脚本使用固定的新镜像标签 `sctenifoldknk-loong64:r2.1.0-cpu-adapt-sdaa` 和独立容器
`sctenifoldknk-loong64-sdaa-adapt-smoke`。它把源码所在项目目录只读挂载到容器，单独把输出目录可写挂载；网络关闭、不映射 CUDA 或 SDAA 设备，并使用当前宿主机 UID/GID。输出目录非空时推理脚本会拒绝覆盖。

输出位于仓库上层的 `../work/loong64-smoke-test-adapt-sdaa/`，包括 `result.rds`、
`diff_regulation.csv`、`session-info.txt`、`summary.txt`、校验和及设备审计文件。本流程只运行推理，不集成 x86 对比。

## 5. 检查输出

```bash
OUTPUT_DIR=../work/loong64-smoke-test-adapt-sdaa
ls -lh "$OUTPUT_DIR"
cat "$OUTPUT_DIR/summary.txt"
sha256sum --check "$OUTPUT_DIR/output-checksums.sha256"
```
