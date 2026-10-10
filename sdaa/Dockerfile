ARG BASE_IMAGE=lasermpnn-sdaa321-loong:torch2.12.0-sdaa3.2.1-py312
FROM ${BASE_IMAGE}

ARG BASE_IMAGE_ID=sha256:3e425da47a74a03b589f8e1c18a67dcb90c8bcbde977cc1df09fa01209b7a524
ARG HOST_UID=1006
ARG HOST_GID=1006

USER 0:0
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

LABEL org.opencontainers.image.title="scTenifoldKnk LoongArch CPU baseline" \
      org.opencontainers.image.source="https://github.com/Zaoqu-Liu/scTenifoldKnk" \
      org.opencontainers.image.revision="08118337f6fc22cca5a22c443ccb78ef920dfa3b" \
      org.opencontainers.image.version="2.1.0" \
      org.sctenifoldknk.base.image-id="${BASE_IMAGE_ID}" \
      org.sctenifoldknk.compute-backend="CPU_R_RcppEigen_RSpectra" \
      org.sctenifoldknk.gpu-required="false"

COPY sdaa/vendor/R-core-devel-4.4.1-6.lns23.loongarch64.rpm /tmp/R-core-devel.rpm
RUN test "$(uname -m)" = "loongarch64" \
    && dnf --assumeyes --setopt=install_weak_deps=False install \
         R-core-4.4.1-6.lns23 \
         system-rpm-config \
         bzip2-devel \
         flexiblas-devel \
         gcc-c++ \
         gcc-gfortran \
         libRmath-devel \
         libX11-devel \
         libdeflate-devel \
         libicu-devel \
         libtirpc-devel \
         make \
         pcre2-devel \
         pkgconf-pkg-config \
         tcl-devel \
         tk-devel \
         tre-devel \
         xz-devel \
         zlib-devel \
    && printf '%s  %s\n' \
         'f3d38155ccfbfe7b7083e52f040c2329c72ebe8e7d111b1245ae5a53dc1afbcb' \
         '/tmp/R-core-devel.rpm' | sha256sum --check - \
    && rpm -Uvh --nodeps /tmp/R-core-devel.rpm \
    && rm -f /tmp/R-core-devel.rpm \
    && dnf clean all \
    && rm -rf /var/cache/dnf

ENV SCTENIFOLDKNK_PREFIX=/opt/scTenifoldKnk \
    VIRTUAL_ENV=/opt/venvs/scTenifoldKnk \
    R_LIBS_SITE=/opt/scTenifoldKnk/R/library \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    TZ=UTC \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN /opt/venvs/lasermpnn-sdaa321-py312/bin/python3 -m venv --copies "${VIRTUAL_ENV}" \
    && test "$("${VIRTUAL_ENV}/bin/python3" -c 'import sys; print(sys.version_info[:2])')" = "(3, 12)" \
    && "${VIRTUAL_ENV}/bin/python3" -c 'import sys; assert sys.prefix == "/opt/venvs/scTenifoldKnk"'

ENV PATH=/opt/venvs/scTenifoldKnk/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

COPY sdaa/vendor/ /tmp/scTenifoldKnk-cran/
RUN cd /tmp/scTenifoldKnk-cran \
    && sha256sum --check checksums.sha256 \
    && install -d -m 0755 "${R_LIBS_SITE}" \
    && MAKEFLAGS=-j8 R CMD INSTALL --no-multiarch --library="${R_LIBS_SITE}" Rcpp_1.0.14.tar.gz \
    && MAKEFLAGS=-j8 R CMD INSTALL --no-multiarch --library="${R_LIBS_SITE}" Matrix_1.7-6.tar.gz \
    && MAKEFLAGS=-j8 R CMD INSTALL --no-multiarch --library="${R_LIBS_SITE}" MASS_7.3-66.tar.gz \
    && MAKEFLAGS=-j8 R CMD INSTALL --no-multiarch --library="${R_LIBS_SITE}" RcppEigen_0.3.4.0.2.tar.gz \
    && MAKEFLAGS=-j8 R CMD INSTALL --no-multiarch --library="${R_LIBS_SITE}" RSpectra_0.16-2.tar.gz

COPY . /tmp/scTenifoldKnk
RUN MAKEFLAGS=-j8 R CMD INSTALL --no-multiarch --with-keep.source \
      --library="${R_LIBS_SITE}" /tmp/scTenifoldKnk \
    && Rscript -e 'expected <- c(scTenifoldKnk="2.1.0", Matrix="1.7.6", MASS="7.3.66", Rcpp="1.0.14", RcppEigen="0.3.4.0.2", RSpectra="0.16.2"); actual <- vapply(names(expected), function(x) as.character(packageVersion(x)), character(1)); print(actual); stopifnot(identical(unname(actual), unname(expected)))' \
    && rm -rf /tmp/scTenifoldKnk /tmp/scTenifoldKnk-cran \
    && install -d -m 0755 -o "${HOST_UID}" -g "${HOST_GID}" \
         /workspace/repro /workspace/output

WORKDIR /workspace/repro
USER ${HOST_UID}:${HOST_GID}
CMD ["sleep", "infinity"]
