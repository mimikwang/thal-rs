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
