use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join(name)
}

fn run(args: &[&str]) -> std::process::Output {
    let output = Command::new(env!("CARGO_BIN_EXE_refmat"))
        .args(args)
        .output()
        .expect("refmat should run");
    assert!(
        output.status.success(),
        "refmat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn assert_reference(input: &str, expected: &str, extra: &[&str]) {
    let output = std::env::temp_dir().join(format!("refmat-{}-{input}.tsv", std::process::id()));
    let input_path = fixture(input);
    let mut args = vec![
        "build",
        input_path.to_str().unwrap(),
        "--column",
        "cell_type",
    ];
    args.extend_from_slice(extra);
    args.extend(["--output", output.to_str().unwrap()]);
    run(&args);

    let actual = fs::read_to_string(&output).unwrap();
    let expected = fs::read_to_string(fixture(expected)).unwrap();
    let actual_rows = numeric_rows(&actual);
    let expected_rows = numeric_rows(&expected);
    assert_eq!(actual_rows.len(), expected_rows.len());
    for ((actual_name, actual_values), (expected_name, expected_values)) in
        actual_rows.iter().zip(&expected_rows)
    {
        assert_eq!(actual_name, expected_name);
        assert_eq!(actual_values.len(), expected_values.len());
        for (actual, expected) in actual_values.iter().zip(expected_values) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }
    fs::remove_file(output).unwrap();
}

fn numeric_rows(contents: &str) -> Vec<(String, Vec<f64>)> {
    contents
        .lines()
        .skip(1)
        .map(|line| {
            let mut fields = line.split('\t');
            let name = fields.next().unwrap().to_owned();
            let values = fields.map(|value| value.parse().unwrap()).collect();
            (name, values)
        })
        .collect()
}

#[test]
fn detects_and_prints_sce_and_h5ad_metadata() {
    for input in ["sce.rds", "anndata.h5ad"] {
        let path = fixture(input);
        let output = run(&["head", path.to_str().unwrap(), "-n", "2"]);
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.contains("| cell  | cell_type | batch | score |"));
        assert!(stdout.contains("| Cell1 | T cell    | one   | 0.1   |"));
        assert!(stdout.lines().all(|line| !line.contains('\t')));
    }
}

#[test]
fn counts_cells_by_metadata_column() {
    for input in ["sce.rds", "anndata.h5ad"] {
        let path = fixture(input);
        let output = run(&["col", path.to_str().unwrap(), "cell_type"]);
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.contains("| cell_type | cells |"));
        assert!(stdout.contains("| T cell    | 2     |"));
        assert!(stdout.contains("| B cell    | 2     |"));
    }
}

#[test]
fn builds_sce_sparse_and_dense_assays() {
    assert_reference("sce.rds", "sce.expected.tsv", &[]);
    assert_reference(
        "sce.rds",
        "sce.expected.tsv",
        &["--assay", "dense_logcounts"],
    );
}

#[test]
fn builds_h5ad_dense_x_and_sparse_counts_layer() {
    assert_reference("anndata.h5ad", "anndata.expected.tsv", &[]);
    assert_reference(
        "anndata.h5ad",
        "anndata.expected.tsv",
        &["--layer", "counts"],
    );
}
