# thal-rs

A rust rewrite of the thermal alignment module in [primer3](https://github.com/primer3-org/primer3). This is largely a learning exercise.

## Build

Build the tool by running:

```
cargo build --release
```

This should create `thal-rs` in the `target/release` directory.

## Run

Right now, only dimer thermal alignment is available. Run the tool by running

```
./target/release/thal-rs dimer CCCCCATCCGATCAGGGGG GGGGGACTAGCCTACCCCC
```

You should see the following:

```
dg=-9865.89800000003 dh=-101200.0 ds=-294.4836433983556
```
