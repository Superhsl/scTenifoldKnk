# 固定的 R 源码依赖

这些文件来自 CRAN 官方站点，由 `checksums.sha256` 固定：

- `Matrix_1.7-6.tar.gz`: `https://cran.r-project.org/src/contrib/Matrix_1.7-6.tar.gz`
- `MASS_7.3-66.tar.gz`: `https://cran.r-project.org/src/contrib/MASS_7.3-66.tar.gz`
- `Rcpp_1.0.14.tar.gz`: `https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.14.tar.gz`
- `RcppEigen_0.3.4.0.2.tar.gz`: `https://cran.r-project.org/src/contrib/RcppEigen_0.3.4.0.2.tar.gz`
- `RSpectra_0.16-2.tar.gz`: `https://cran.r-project.org/src/contrib/RSpectra_0.16-2.tar.gz`

scTenifoldKnk 不是预训练模型，不存在模型权重下载。

`R-core-devel-4.4.1-6.lns23.loongarch64.rpm` 来自 Loongnix 官方 EPLL：

`https://pkg.loongnix.cn/loongnix-server/23.1/epll/loongarch64/Packages/R-core-devel-4.4.1-6.lns23.loongarch64.rpm`

Loongnix EPLL 的 `R-core-devel` 错误依赖了不存在的 `redhat-rpm-config` capability，
而仓库自身配置为 `gpgcheck=0`。构建先通过 DNF 安装全部真实依赖，再校验该 RPM 的
固定 SHA-256，最后只对这一 RPM 使用 `rpm --nodeps`。
