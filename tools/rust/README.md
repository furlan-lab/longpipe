# Longpipe Rust tools

This crate contains native helpers used only by the Longpipe workflow:

- `copy_pbWDL` materializes the file outputs from a MiniWDL JSON result;
- `shrink_fire_bed` converts a FIRE pileup to the four-column bedGraph used
  for bigWig generation, merging adjacent sentinel intervals.

Snakemake builds both binaries from these sources into the retained scratch
run. They are not user-configured external dependencies.

For development, run:

```bash
cargo test --locked --manifest-path tools/rust/Cargo.toml
```
