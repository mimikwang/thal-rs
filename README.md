# thal-rs

A rust rewrite of the thermal alignment module in [primer3](https://github.com/primer3-org/primer3). This is largely a learning exercise.

## Build

Build the tool by running:

```
cargo build --release
```

This should create `thal-rs` in the `target/release` directory.

## Run

### Dimers

```
./target/release/thal-rs dimer CCCCCATCCGATCAGGGGG GGGGGACTAGCCTACCCCC
```

You should see the following:
```
dg=-9865.89800000003 dh=-101200.0 ds=-294.4836433983556
```

### Monomers

```
./target/release/thal-rs monomer CCCCCATCCGATCAGGGGG
```

You should see the following:
```
dg=-3796.279999999999 dh=-36300.0 ds=-104.80000000000001
```

## Python bindings

Python bindings live in `python/` as a separate crate (`thal-py`) that depends on
this one, built with [maturin](https://www.maturin.rs/) and managed with
[uv](https://docs.astral.sh/uv/). The main crate has no Python-related
dependencies; `cargo build` here is unaffected.

To build and install the `thal_rs` module into a venv at `python/.venv`:

```
cd python
uv run maturin develop
```

Then use it:

```
uv run python3 -c "
import thal_rs
print(thal_rs.monomer('CCCCCATCCGATCAGGGGG'))
print(thal_rs.dimer('CCCCCATCCGATCAGGGGG', 'GGGGGACTAGCCTACCCCC'))
"
```

You should see the same `(dg, dh, ds)` values as the CLI:
```
(-3796.279999999999, -36300.0, -104.80000000000001)
(-9865.89800000003, -101200.0, -294.4836433983556)
```

To use `thal_rs` from another uv project, add it as an editable path dependency
pointing at `python/`, e.g. in `pyproject.toml`:

```toml
[tool.uv.sources]
thal-rs = { path = "../python", editable = true }
```
