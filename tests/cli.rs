use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

static OUTPUT_SEQUENCE: AtomicUsize = AtomicUsize::new(0);

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
    assert_reference_with_column(input, expected, &["--column", "cell_type"], extra);
}

fn assert_reference_with_column(input: &str, expected: &str, column_args: &[&str], extra: &[&str]) {
    let sequence = OUTPUT_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let output = std::env::temp_dir().join(format!(
        "refmat-{}-{sequence}-{input}.tsv",
        std::process::id()
    ));
    let input_path = fixture(input);
    let mut args = vec!["build", input_path.to_str().unwrap()];
    args.extend_from_slice(column_args);
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
        for column_args in [
            &["cell_type"][..],
            &["-c", "cell_type"][..],
            &["--column", "cell_type"][..],
        ] {
            let mut args = vec!["col", path.to_str().unwrap()];
            args.extend_from_slice(column_args);
            let output = run(&args);
            let stdout = String::from_utf8(output.stdout).unwrap();
            assert!(stdout.contains("| cell_type | cells |"));
            assert!(stdout.contains("| T cell    | 2     |"));
            assert!(stdout.contains("| B cell    | 2     |"));
        }
    }
}

#[test]
fn build_accepts_each_column_syntax() {
    assert_reference_with_column("anndata.h5ad", "anndata.expected.tsv", &["cell_type"], &[]);
    assert_reference_with_column(
        "anndata.h5ad",
        "anndata.expected.tsv",
        &["-c", "cell_type"],
        &[],
    );
}

#[test]
fn rejects_multiple_column_arguments() {
    let path = fixture("anndata.h5ad");
    let output = Command::new(env!("CARGO_BIN_EXE_refmat"))
        .args([
            "col",
            path.to_str().unwrap(),
            "cell_type",
            "--column",
            "batch",
        ])
        .output()
        .expect("refmat should run");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot be used with"));
}

#[test]
fn detects_and_prints_legacy_seurat_metadata() {
    for (input, version) in [("seurat-v3.rds", "3.2.3"), ("seurat-v4.rds", "4.4.0")] {
        let path = fixture(input);

        let inspect = run(&["inspect", path.to_str().unwrap(), "--depth", "1"]);
        let inspect_stdout = String::from_utf8(inspect.stdout).unwrap();
        assert!(inspect_stdout.contains(&format!("detected: Seurat {version}")));
        assert!(inspect_stdout.contains("Assay (v3/v4 layout)"));

        let head = run(&["head", path.to_str().unwrap(), "-n", "2"]);
        let head_stdout = String::from_utf8(head.stdout).unwrap();
        let head_stderr = String::from_utf8(head.stderr).unwrap();
        assert!(head_stderr.contains(&format!("detected: Seurat {version}")));
        assert!(head_stderr.contains("Assay (v3/v4 layout)"));
        assert!(head_stdout.contains("| cell  | cell_type | batch | score |"));
        assert!(head_stdout.contains("| Cell1 | T cell    | one   | 0.1   |"));

        let col = run(&["col", path.to_str().unwrap(), "cell_type"]);
        let col_stdout = String::from_utf8(col.stdout).unwrap();
        assert!(col_stdout.contains("| B cell    | 2     |"));
        assert!(col_stdout.contains("| T cell    | 2     |"));
    }
}

#[test]
fn builds_legacy_seurat_sparse_and_dense_layers() {
    for input in ["seurat-v3.rds", "seurat-v4.rds"] {
        assert_reference(input, "seurat-legacy.expected.tsv", &[]);
        assert_reference(input, "seurat-legacy.expected.tsv", &["--layer", "counts"]);
        assert_reference(
            input,
            "seurat-legacy.expected.tsv",
            &["--layer", "scale.data"],
        );
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
