use std::{
    env, fs,
    io::Write,
    path::{Path, PathBuf},
};

const COMMENT_PREFIX: &str = "#";
const THERMO_DIR: &str = "thermo";
const STACK_DH: &str = "stack.dh";
const STACK_DS: &str = "stack.ds";
const STACKMM_DH: &str = "stackmm.dh";
const STACKMM_DS: &str = "stackmm.ds";
const TSTACK_DH: &str = "tstack.dh";
const TSTACK_DS: &str = "tstack.ds";
const DANGLE_DH: &str = "dangle.dh";
const DANGLE_DS: &str = "dangle.ds";
const LOOPS_DH: &str = "loops.dh";
const LOOPS_DS: &str = "loops.ds";
const TRILOOP_DH: &str = "triloop.dh";
const TRILOOP_DS: &str = "triloop.ds";
const TETRALOOP_DH: &str = "tetraloop.dh";
const TETRALOOP_DS: &str = "tetraloop.ds";

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let thermo_dir = Path::new(&manifest_dir).join(THERMO_DIR);

    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("gen_params.rs");
    let mut file = fs::File::create(&dest_path).unwrap();

    write_stack(&mut file, &thermo_dir);
    write_stack_mm(&mut file, &thermo_dir);
    write_tstack(&mut file, &thermo_dir);
    write_dangle(&mut file, &thermo_dir);
    write_triloop(&mut file, &thermo_dir);
    write_tetraloop(&mut file, &thermo_dir);
    write_loops(&mut file, &thermo_dir);
}

fn write_stack(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(STACK_DH));
    let ds_data = read_file(thermo_dir.join(STACK_DS));

    let mut values = [(f64::INFINITY, f64::INFINITY); 625];

    for line in dh_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[pair_key_index(&key)].1 = value;
    }

    for line in ds_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[pair_key_index(&key)].0 = value;
    }

    let mut code = String::from("pub static STACK_ARRAY: [Thermo; 625] = [\n");
    build_thermo_lines(&values, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn write_stack_mm(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(STACKMM_DH));
    let ds_data = read_file(thermo_dir.join(STACKMM_DS));

    let mut values = [(f64::INFINITY, f64::INFINITY); 625];

    for line in dh_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[pair_key_index(&key)].1 = value;
    }

    for line in ds_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[pair_key_index(&key)].0 = value;
    }

    let mut code = String::from("pub static STACK_MM_ARRAY: [Thermo; 625] = [\n");
    build_thermo_lines(&values, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn write_tstack(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(TSTACK_DH));
    let ds_data = read_file(thermo_dir.join(TSTACK_DS));

    let mut values = [(f64::INFINITY, f64::INFINITY); 625];

    for line in dh_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[pair_key_index(&key)].1 = value;
    }

    for line in ds_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[pair_key_index(&key)].0 = value;
    }

    let mut code = String::from("pub static TSTACK_ARRAY: [Thermo; 625] = [\n");
    build_thermo_lines(&values, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn write_dangle(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(DANGLE_DH));
    let ds_data = read_file(thermo_dir.join(DANGLE_DS));

    let mut dangle3_values = [(f64::INFINITY, f64::INFINITY); 125];
    let mut dangle5_values = [(f64::INFINITY, f64::INFINITY); 125];

    for line in dh_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        set_dangle_value(&key, value, true, &mut dangle3_values, &mut dangle5_values);
    }

    for line in ds_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        set_dangle_value(&key, value, false, &mut dangle3_values, &mut dangle5_values);
    }

    let mut code = String::from("pub static DANGLE3_ARRAY: [Thermo; 125] = [\n");
    build_thermo_lines(&dangle3_values, &mut code);
    code.push_str("];\n\n");

    code.push_str("pub static DANGLE5_ARRAY: [Thermo; 125] = [\n");
    build_thermo_lines(&dangle5_values, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn set_dangle_value(
    key: &str,
    value: f64,
    is_dh: bool,
    dangle3_values: &mut [(f64, f64); 125],
    dangle5_values: &mut [(f64, f64); 125],
) {
    let (b1, b2) = key.split_once('_').expect("invalid dangle key");
    let b1 = b1.as_bytes();
    let b2 = b2.as_bytes();

    let (values, index) = if b1.len() == 2 && b2.len() == 1 {
        let index = triple_index(base_to_num(b1[0]), base_to_num(b1[1]), base_to_num(b2[0]));
        (dangle3_values, index)
    } else if b1.len() == 1 && b2.len() == 2 {
        let index = triple_index(base_to_num(b1[0]), base_to_num(b2[0]), base_to_num(b2[1]));
        (dangle5_values, index)
    } else {
        panic!("invalid dangle key length: {key}");
    };

    if is_dh {
        values[index].1 = value;
    } else {
        values[index].0 = value;
    }
}

fn write_triloop(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(TRILOOP_DH));
    let ds_data = read_file(thermo_dir.join(TRILOOP_DS));

    let mut values = [(f64::INFINITY, f64::INFINITY); 3125];

    for line in dh_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[triloop_key_index(&key)].1 = value;
    }

    for line in ds_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[triloop_key_index(&key)].0 = value;
    }

    let mut code = String::from("pub static TRILOOP_ARRAY: [Thermo; 3125] = [\n");
    build_thermo_lines(&values, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn write_tetraloop(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(TETRALOOP_DH));
    let ds_data = read_file(thermo_dir.join(TETRALOOP_DS));

    let mut values = [(f64::INFINITY, f64::INFINITY); 15625];

    for line in dh_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[tetraloop_key_index(&key)].1 = value;
    }

    for line in ds_data.lines() {
        if should_skip_line(line) {
            continue;
        }
        let (key, value) = parse_line(line);
        values[tetraloop_key_index(&key)].0 = value;
    }

    let mut code = String::from("pub static TETRALOOP_ARRAY: [Thermo; 15625] = [\n");
    build_thermo_lines(&values, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn write_loops(file: &mut fs::File, thermo_dir: &PathBuf) {
    let dh_data = read_file(thermo_dir.join(LOOPS_DH));
    let ds_data = read_file(thermo_dir.join(LOOPS_DS));

    let mut internal = [(f64::INFINITY, f64::INFINITY); 30];
    let mut bulge = [(f64::INFINITY, f64::INFINITY); 30];
    let mut hairpin = [(f64::INFINITY, f64::INFINITY); 30];

    fill_loop_values(&dh_data, true, &mut internal, &mut bulge, &mut hairpin);
    fill_loop_values(&ds_data, false, &mut internal, &mut bulge, &mut hairpin);

    let mut code = String::from("pub static INTERNAL_ARRAY: [Thermo; 30] = [\n");
    build_thermo_lines(&internal, &mut code);
    code.push_str("];\n\n");

    code.push_str("pub static BULGE_ARRAY: [Thermo; 30] = [\n");
    build_thermo_lines(&bulge, &mut code);
    code.push_str("];\n\n");

    code.push_str("pub static HAIRPIN_ARRAY: [Thermo; 30] = [\n");
    build_thermo_lines(&hairpin, &mut code);
    code.push_str("];\n\n");

    file.write_all(code.as_bytes()).unwrap();
}

fn fill_loop_values(
    data: &str,
    is_dh: bool,
    internal: &mut [(f64, f64); 30],
    bulge: &mut [(f64, f64); 30],
    hairpin: &mut [(f64, f64); 30],
) {
    let mut index = 0;
    for line in data.lines() {
        if should_skip_line(line) {
            continue;
        }

        let mut parts = line.split_whitespace();
        parts.next().expect("missing size column");
        let internal_value = parse_value(parts.next().expect("missing internal value"));
        let bulge_value = parse_value(parts.next().expect("missing bulge value"));
        let hairpin_value = parse_value(parts.next().expect("missing hairpin value"));

        if is_dh {
            internal[index].1 = internal_value;
            bulge[index].1 = bulge_value;
            hairpin[index].1 = hairpin_value;
        } else {
            internal[index].0 = internal_value;
            bulge[index].0 = bulge_value;
            hairpin[index].0 = hairpin_value;
        }
        index += 1;
    }
}

fn triloop_key_index(key: &str) -> usize {
    let bases = key.as_bytes();
    assert_eq!(bases.len(), 5, "invalid triloop key length: {key}");
    triloop_index(
        base_to_num(bases[0]),
        base_to_num(bases[1]),
        base_to_num(bases[2]),
        base_to_num(bases[3]),
        base_to_num(bases[4]),
    )
}

fn tetraloop_key_index(key: &str) -> usize {
    let bases = key.as_bytes();
    assert_eq!(bases.len(), 6, "invalid tetraloop key length: {key}");
    tetraloop_index(
        base_to_num(bases[0]),
        base_to_num(bases[1]),
        base_to_num(bases[2]),
        base_to_num(bases[3]),
        base_to_num(bases[4]),
        base_to_num(bases[5]),
    )
}

fn triple_index(b0: usize, b1: usize, b2: usize) -> usize {
    (b0 * 5 + b1) * 5 + b2
}

fn quad_index(b0: usize, b1: usize, b2: usize, b3: usize) -> usize {
    triple_index(b0, b1, b2) * 5 + b3
}

fn triloop_index(b0: usize, b1: usize, b2: usize, b3: usize, b4: usize) -> usize {
    quad_index(b0, b1, b2, b3) * 5 + b4
}

fn tetraloop_index(b0: usize, b1: usize, b2: usize, b3: usize, b4: usize, b5: usize) -> usize {
    triloop_index(b0, b1, b2, b3, b4) * 5 + b5
}

fn thermo_line(ds: String, dh: String) -> String {
    format!("    Thermo {{ ds: {ds}, dh: {dh} }},\n")
}

fn build_thermo_lines(values: &[(f64, f64)], code: &mut String) {
    for (ds, dh) in values {
        let ds = format_value(*ds);
        let dh = format_value(*dh);
        code.push_str(&thermo_line(ds, dh));
    }
}

fn parse_line(line: &str) -> (String, f64) {
    let mut parts = line.split_whitespace();
    let key = parts.next().expect("missing key").to_owned();
    let value = parse_value(parts.next().expect("missing value"));
    (key, value)
}

fn parse_value(val: &str) -> f64 {
    if val == "inf" {
        f64::INFINITY
    } else {
        val.parse().expect("value should be a number")
    }
}

fn base_to_num(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'G' => 1,
        b'C' => 2,
        b'T' => 3,
        _ => 4,
    }
}

fn pair_key_index(key: &str) -> usize {
    let (b1, b2) = key.split_once('_').expect("invalid key");
    let b1 = b1.as_bytes();
    let b2 = b2.as_bytes();
    assert!(b1.len() == 2 && b2.len() == 2, "invalid key length: {key}");

    let (b0, b1) = (base_to_num(b1[0]), base_to_num(b1[1]));
    let (b2, b3) = (base_to_num(b2[0]), base_to_num(b2[1]));
    ((b0 * 5 + b1) * 5 + b2) * 5 + b3
}

fn format_value(val: f64) -> String {
    if val == f64::INFINITY {
        "f64::INFINITY".to_owned()
    } else if val == f64::NEG_INFINITY {
        "f64::NEG_INFINITY".to_owned()
    } else {
        format!("{val:?}")
    }
}

fn read_file(path: PathBuf) -> String {
    println!("cargo::rerun-if-changed={}", path.display());
    fs::read_to_string(&path).unwrap()
}

fn should_skip_line(line: &str) -> bool {
    line.starts_with(COMMENT_PREFIX) || line.is_empty()
}
