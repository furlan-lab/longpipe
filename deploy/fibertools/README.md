# Rhino GPU Fibertools Build

This directory catalogs the Furlan Lab's GPU-enabled Fibertools deployment.
The binary and LibTorch libraries are deployment artifacts and are not tracked
in this repository.

## Validated Build

- Fibertools release: `v0.12.0`
- Git commit: `c881786`
- Rust: `rustc 1.90.0`
- Rhino modules:
  - `GCC/13.3.0`
  - `CMake/3.29.3-GCCcore-13.3.0`
  - `CUDA/12.4.0`
- LibTorch: `2.6.0+cu124`
- Source checkout: `/home/sfurlan/develop/fibertools-rs`
- Cargo target: `target-rhino-v0.12.0-cu124`
- Versioned executable:
  `/home/sfurlan/.local/opt/fibertools-rs/0.12.0-cu124/bin/ft`
- SHA-256:
  `89514bac7b94f6f1e9c09b8718d8d97ac6b444d56dbf9c279e8aba60c927c744`

The executable has this embedded `RUNPATH`:

```text
/home/sfurlan/lib/libtorch-2.6.0-cu124/lib
/app/software/GCCcore/13.3.0/lib64
/app/software/CUDA/12.4.0/lib
```

No existing Fibertools binary, LibTorch installation, Conda environment, or
module-managed software was replaced.

## Build Recipe

The source checkout was placed at the official `v0.12.0` release tag. LibTorch
was downloaded from the official PyTorch distribution:

```text
https://download.pytorch.org/libtorch/cu124/libtorch-shared-with-deps-2.6.0%2Bcu124.zip
```

The release build used:

```bash
module purge
module load GCC/13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load CUDA/12.4.0

export LIBTORCH=/home/sfurlan/lib/libtorch-2.6.0-cu124
export LIBTORCH_CXX11_ABI=0
export LD_LIBRARY_PATH="$LIBTORCH/lib:${LD_LIBRARY_PATH:-}"
export CARGO_TARGET_DIR=/home/sfurlan/develop/fibertools-rs/target-rhino-v0.12.0-cu124
export RUSTFLAGS="-C link-arg=-Wl,-rpath,$LIBTORCH/lib \
-C link-arg=-Wl,-rpath,$EBROOTGCC/lib64 \
-C link-arg=-Wl,-rpath,$EBROOTCUDA/lib"

cd /home/sfurlan/develop/fibertools-rs
cargo build --release --all-features
```

## Validation

The installed binary was tested with an empty environment. `ft --version` and
`ft predict-m6a --help` succeeded, and `ldd` reported no unresolved direct
libraries.

Slurm job `62646509` then ran GPU prediction on 100 reads from the SNF_4N
kinetics BAM. It completed successfully in 28 seconds using one GPU, retained
all 100 reads, and emitted `MM`/`ML` modification tags plus `nuc` and `msp`
annotations in the current `MA` encoding. Longpipe's input detector classified
the result as `use_existing` with 99% m6A and annotation coverage.

Configure a run to use this deployment with:

```yaml
ft_predict_backend: external
ft_predict_executable: /home/sfurlan/.local/opt/fibertools-rs/0.12.0-cu124/bin/ft
ft_predict_expected_version: 0.12.0
```
