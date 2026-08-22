use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail, ensure};
use hdf5::types::{TypeDescriptor, VarLenAscii, VarLenUnicode};
use hdf5::{Dataset, File, Group, LocationType};

use crate::{
    InputScale, default_output, escape_tsv, finish_group_means, to_linear, write_reference,
};

pub(crate) fn inspect(path: &Path, _depth: usize) -> Result<()> {
    let file = open(path)?;
    validate_root(&file)?;
    let obs = file.group("obs").context("H5AD has no /obs dataframe")?;
    let var = file.group("var").context("H5AD has no /var dataframe")?;
    let columns = dataframe_columns(&obs)?;
    let obs_names = dataframe_index(&obs)?;
    let var_names = dataframe_index(&var)?;
    println!("root: AnnData (H5AD)");
    println!(
        "  obs: data.frame[{} rows x {} columns]",
        obs_names.len(),
        columns.len()
    );
    for name in columns {
        println!("    {name}: {}", column_description(&obs, &name)?);
    }
    println!("  var: data.frame[{} rows]", var_names.len());
    if file.link_exists("X") {
        println!("  X: {}", matrix_description(&file, "X")?);
    }
    if file.link_exists("layers") {
        let layers = file.group("layers")?;
        println!("  layers: [{}]", layers.member_names()?.join(", "));
    }
    Ok(())
}

pub(crate) fn head(path: &Path, rows: usize) -> Result<()> {
    let file = open(path)?;
    validate_root(&file)?;
    let obs = file.group("obs").context("H5AD has no /obs dataframe")?;
    let names = dataframe_index(&obs)?;
    let take = rows.min(names.len());
    let columns = dataframe_columns(&obs)?;
    let rendered = columns
        .iter()
        .map(|name| {
            read_column_strings(&obs, name, take)
                .with_context(|| format!("failed to read obs column '{name}'"))
        })
        .collect::<Result<Vec<_>>>()?;
    print!("cell");
    for name in &columns {
        print!("\t{name}");
    }
    println!();
    for row in 0..take {
        print!("{}", escape_tsv(&names[row]));
        for values in &rendered {
            print!("\t{}", escape_tsv(&values[row]));
        }
        println!();
    }
    Ok(())
}

pub(crate) fn build(
    path: &Path,
    column: &str,
    layer: Option<&str>,
    scale: InputScale,
    output: Option<&Path>,
) -> Result<()> {
    let file = open(path)?;
    validate_root(&file)?;
    let obs = file.group("obs").context("H5AD has no /obs dataframe")?;
    let (group_names, cell_groups) = read_groups(&obs, column)?;
    ensure!(
        !group_names.is_empty(),
        "obs column '{column}' has no non-missing groups"
    );
    let var = file.group("var").context("H5AD has no /var dataframe")?;
    let feature_names = dataframe_index(&var)?;
    let matrix_path = layer.map_or_else(|| "X".to_owned(), |name| format!("layers/{name}"));
    ensure!(
        file.link_exists(&matrix_path),
        "matrix '/{matrix_path}' does not exist"
    );
    let matrix_name = layer.unwrap_or("X");
    let sums = aggregate_matrix(
        &file,
        &matrix_path,
        matrix_name,
        &cell_groups,
        group_names.len(),
        feature_names.len(),
        scale,
    )?;
    let output_path: PathBuf = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| default_output(path));
    write_reference(&output_path, &feature_names, &group_names, &sums)?;
    eprintln!(
        "wrote {} features x {} groups to {}",
        feature_names.len(),
        group_names.len(),
        output_path.display()
    );
    Ok(())
}

fn open(path: &Path) -> Result<File> {
    hdf5::silence_errors(true);
    File::open(path).with_context(|| format!("failed to open H5AD file {}", path.display()))
}

fn validate_root(file: &File) -> Result<()> {
    let encoding = attr_string(file, "encoding-type").unwrap_or_default();
    ensure!(
        encoding == "anndata" || (file.link_exists("obs") && file.link_exists("var")),
        "HDF5 file is not an AnnData object"
    );
    Ok(())
}

fn dataframe_columns(group: &Group) -> Result<Vec<String>> {
    if let Ok(attr) = group.attr("column-order") {
        if let Ok(values) = attr.read_raw::<VarLenUnicode>() {
            return Ok(values.iter().map(ToString::to_string).collect());
        }
        if let Ok(values) = attr.read_raw::<VarLenAscii>() {
            return Ok(values.iter().map(ToString::to_string).collect());
        }
    }
    let index = attr_string(group, "_index").unwrap_or_else(|_| "_index".into());
    Ok(group
        .member_names()?
        .into_iter()
        .filter(|name| name != &index)
        .collect())
}

fn dataframe_index(group: &Group) -> Result<Vec<String>> {
    let index = attr_string(group, "_index").unwrap_or_else(|_| "_index".into());
    read_string_dataset(
        &group
            .dataset(&index)
            .with_context(|| format!("dataframe index dataset '{index}' is missing"))?,
    )
}

fn attr_string(location: &hdf5::Location, name: &str) -> Result<String> {
    let attr = location.attr(name)?;
    if let Ok(value) = attr.read_scalar::<VarLenUnicode>() {
        return Ok(value.to_string());
    }
    if let Ok(value) = attr.read_scalar::<VarLenAscii>() {
        return Ok(value.to_string());
    }
    bail!("attribute '{name}' is not a scalar string")
}

fn read_string_dataset(dataset: &Dataset) -> Result<Vec<String>> {
    if let Ok(values) = dataset.read_raw::<VarLenUnicode>() {
        return Ok(values.iter().map(ToString::to_string).collect());
    }
    if let Ok(values) = dataset.read_raw::<VarLenAscii>() {
        return Ok(values.iter().map(ToString::to_string).collect());
    }
    bail!(
        "dataset '{}' is not a supported string array",
        dataset.name()
    )
}

fn read_column_strings(group: &Group, name: &str, take: usize) -> Result<Vec<String>> {
    match group.loc_type_by_name(name)? {
        LocationType::Group => {
            let column = group.group(name)?;
            match attr_string(&column, "encoding-type")
                .unwrap_or_default()
                .as_str()
            {
                "categorical" => {
                    let categories = read_string_dataset(&column.dataset("categories")?)?;
                    let codes = read_integer_dataset(&column.dataset("codes")?)?;
                    Ok(codes
                        .into_iter()
                        .take(take)
                        .map(|code| {
                            if code >= 0 {
                                categories
                                    .get(code as usize)
                                    .cloned()
                                    .unwrap_or_else(|| "NA".into())
                            } else {
                                "NA".into()
                            }
                        })
                        .collect())
                }
                "nullable-integer" | "nullable-boolean" => {
                    let values = read_numeric_dataset(&column.dataset("values")?)?;
                    let mask = read_bool_dataset(&column.dataset("mask")?)?;
                    Ok(values
                        .into_iter()
                        .zip(mask)
                        .take(take)
                        .map(|(value, missing)| {
                            if missing {
                                "NA".into()
                            } else {
                                format_number(value)
                            }
                        })
                        .collect())
                }
                encoding => bail!("unsupported AnnData column encoding '{encoding}' for '{name}'"),
            }
        }
        LocationType::Dataset => {
            let dataset = group.dataset(name)?;
            if is_string(&dataset)? {
                Ok(read_string_dataset(&dataset)?
                    .into_iter()
                    .take(take)
                    .collect())
            } else if dataset.dtype()?.is::<bool>() {
                Ok(read_bool_dataset(&dataset)?
                    .into_iter()
                    .take(take)
                    .map(|value| if value { "TRUE".into() } else { "FALSE".into() })
                    .collect())
            } else {
                Ok(read_numeric_dataset(&dataset)?
                    .into_iter()
                    .take(take)
                    .map(format_number)
                    .collect())
            }
        }
        other => bail!("unsupported HDF5 object type {other:?} for obs column '{name}'"),
    }
}

fn read_groups(group: &Group, name: &str) -> Result<(Vec<String>, Vec<Option<usize>>)> {
    ensure!(
        group.link_exists(name),
        "obs column '{name}' does not exist"
    );
    if group.loc_type_by_name(name)? == LocationType::Group {
        let column = group.group(name)?;
        ensure!(
            attr_string(&column, "encoding-type").unwrap_or_default() == "categorical",
            "grouping column '{name}' must be categorical or string"
        );
        let names = read_string_dataset(&column.dataset("categories")?)?;
        let codes = read_integer_dataset(&column.dataset("codes")?)?;
        let groups = codes
            .into_iter()
            .map(|code| {
                if code >= 0 && (code as usize) < names.len() {
                    Some(code as usize)
                } else {
                    None
                }
            })
            .collect();
        return Ok((names, groups));
    }
    let dataset = group.dataset(name)?;
    ensure!(
        is_string(&dataset)?,
        "grouping column '{name}' must be categorical or string"
    );
    let values = read_string_dataset(&dataset)?;
    let mut names = Vec::new();
    let mut indices = HashMap::new();
    let groups = values
        .into_iter()
        .map(|value| {
            if let Some(index) = indices.get(&value) {
                *index
            } else {
                let index = names.len();
                names.push(value.clone());
                indices.insert(value, index);
                index
            }
        })
        .map(Some)
        .collect();
    Ok((names, groups))
}

fn aggregate_matrix(
    file: &File,
    path: &str,
    matrix_name: &str,
    cell_groups: &[Option<usize>],
    group_count: usize,
    feature_count: usize,
    scale: InputScale,
) -> Result<Vec<f64>> {
    let mut group_sizes = vec![0usize; group_count];
    for group in cell_groups.iter().flatten() {
        group_sizes[*group] += 1;
    }
    ensure!(
        group_sizes.iter().all(|size| *size > 0),
        "annotation contains an empty category"
    );
    let mut sums = vec![0.0; feature_count * group_count];
    match file.loc_type_by_name(path)? {
        LocationType::Dataset => aggregate_dense(
            &file.dataset(path)?,
            matrix_name,
            cell_groups,
            group_count,
            feature_count,
            scale,
            &mut sums,
        )?,
        LocationType::Group => {
            let matrix = file.group(path)?;
            let encoding = attr_string(&matrix, "encoding-type")?;
            let shape = read_shape(&matrix)?;
            ensure!(
                shape == [cell_groups.len(), feature_count],
                "matrix shape {:?} does not match {} cells x {} features",
                shape,
                cell_groups.len(),
                feature_count
            );
            match encoding.as_str() {
                "csr_matrix" => aggregate_csr(
                    &matrix,
                    matrix_name,
                    cell_groups,
                    group_count,
                    feature_count,
                    scale,
                    &mut sums,
                )?,
                "csc_matrix" => aggregate_csc(
                    &matrix,
                    matrix_name,
                    cell_groups,
                    group_count,
                    feature_count,
                    scale,
                    &mut sums,
                )?,
                _ => bail!("unsupported AnnData matrix encoding '{encoding}'"),
            }
        }
        other => bail!("unsupported HDF5 matrix object {other:?}"),
    }
    finish_group_means(&mut sums, &group_sizes);
    Ok(sums)
}

fn aggregate_dense(
    dataset: &Dataset,
    name: &str,
    groups: &[Option<usize>],
    group_count: usize,
    features: usize,
    scale: InputScale,
    sums: &mut [f64],
) -> Result<()> {
    let shape = dataset.shape();
    ensure!(
        shape == [groups.len(), features],
        "dense matrix shape {:?} does not match {} cells x {} features",
        shape,
        groups.len(),
        features
    );
    for start in (0..groups.len()).step_by(512) {
        let end = (start + 512).min(groups.len());
        let values = read_dense_slice(dataset, start, end)?;
        for (local_cell, row) in values.chunks_exact(features).enumerate() {
            if let Some(group) = groups[start + local_cell] {
                for (feature, value) in row.iter().enumerate() {
                    sums[feature * group_count + group] += to_linear(*value, scale, name);
                }
            }
        }
    }
    Ok(())
}

fn aggregate_csr(
    matrix: &Group,
    name: &str,
    groups: &[Option<usize>],
    group_count: usize,
    features: usize,
    scale: InputScale,
    sums: &mut [f64],
) -> Result<()> {
    let indptr = read_integer_dataset(&matrix.dataset("indptr")?)?;
    let indices_dataset = matrix.dataset("indices")?;
    let data_dataset = matrix.dataset("data")?;
    let nnz = usize::try_from(*indptr.last().unwrap_or(&0)).context("negative CSR pointer")?;
    ensure!(
        indptr.len() == groups.len() + 1
            && indices_dataset.size() == nnz
            && data_dataset.size() == nnz,
        "invalid CSR arrays"
    );
    let mut cell = 0usize;
    for start in (0..nnz).step_by(1_000_000) {
        let end = (start + 1_000_000).min(nnz);
        let indices = read_integer_slice(&indices_dataset, start, end)?;
        let values = read_numeric_slice(&data_dataset, start, end)?;
        for (offset, (feature, value)) in indices.into_iter().zip(values).enumerate() {
            let position = start + offset;
            while cell + 1 < indptr.len() && position >= indptr[cell + 1] as usize {
                cell += 1;
            }
            let feature = usize::try_from(feature).context("negative CSR feature index")?;
            ensure!(feature < features, "CSR feature index out of bounds");
            if let Some(group) = groups[cell] {
                sums[feature * group_count + group] += to_linear(value, scale, name);
            }
        }
    }
    Ok(())
}

fn aggregate_csc(
    matrix: &Group,
    name: &str,
    groups: &[Option<usize>],
    group_count: usize,
    features: usize,
    scale: InputScale,
    sums: &mut [f64],
) -> Result<()> {
    let indptr = read_integer_dataset(&matrix.dataset("indptr")?)?;
    let indices_dataset = matrix.dataset("indices")?;
    let data_dataset = matrix.dataset("data")?;
    let nnz = usize::try_from(*indptr.last().unwrap_or(&0)).context("negative CSC pointer")?;
    ensure!(
        indptr.len() == features + 1 && indices_dataset.size() == nnz && data_dataset.size() == nnz,
        "invalid CSC arrays"
    );
    let mut feature = 0usize;
    for start in (0..nnz).step_by(1_000_000) {
        let end = (start + 1_000_000).min(nnz);
        let indices = read_integer_slice(&indices_dataset, start, end)?;
        let values = read_numeric_slice(&data_dataset, start, end)?;
        for (offset, (cell, value)) in indices.into_iter().zip(values).enumerate() {
            let position = start + offset;
            while feature + 1 < indptr.len() && position >= indptr[feature + 1] as usize {
                feature += 1;
            }
            let cell = usize::try_from(cell).context("negative CSC cell index")?;
            ensure!(cell < groups.len(), "CSC cell index out of bounds");
            if let Some(group) = groups[cell] {
                sums[feature * group_count + group] += to_linear(value, scale, name);
            }
        }
    }
    Ok(())
}

fn read_shape(group: &Group) -> Result<[usize; 2]> {
    let attr = group
        .attr("shape")
        .context("sparse matrix has no shape attribute")?;
    let values = if let Ok(values) = attr.read_raw::<i64>() {
        values
    } else {
        attr.read_raw::<u64>()?
            .into_iter()
            .map(|x| x as i64)
            .collect()
    };
    ensure!(
        values.len() == 2 && values.iter().all(|x| *x >= 0),
        "invalid sparse matrix shape"
    );
    Ok([values[0] as usize, values[1] as usize])
}

fn read_integer_dataset(dataset: &Dataset) -> Result<Vec<i64>> {
    read_integer_slice(dataset, 0, dataset.size())
}

fn read_integer_slice(dataset: &Dataset, start: usize, end: usize) -> Result<Vec<i64>> {
    let dtype = dataset.dtype()?;
    macro_rules! read_int {
        ($ty:ty) => {
            if dtype.is::<$ty>() {
                return Ok(dataset
                    .read_slice_1d::<$ty, _>(start..end)?
                    .into_iter()
                    .map(i64::from)
                    .collect());
            }
        };
    }
    read_int!(i8);
    read_int!(u8);
    read_int!(i16);
    read_int!(u16);
    read_int!(i32);
    read_int!(u32);
    if dtype.is::<i64>() {
        return Ok(dataset.read_slice_1d::<i64, _>(start..end)?.to_vec());
    }
    if dtype.is::<u64>() {
        return dataset
            .read_slice_1d::<u64, _>(start..end)?
            .into_iter()
            .map(|x| i64::try_from(x).context("integer exceeds i64"))
            .collect();
    }
    bail!("dataset '{}' is not an integer array", dataset.name())
}

fn read_numeric_dataset(dataset: &Dataset) -> Result<Vec<f64>> {
    read_numeric_slice(dataset, 0, dataset.size())
}

fn read_numeric_slice(dataset: &Dataset, start: usize, end: usize) -> Result<Vec<f64>> {
    let dtype = dataset.dtype()?;
    if dtype.is::<f64>() {
        return Ok(dataset.read_slice_1d::<f64, _>(start..end)?.to_vec());
    }
    if dtype.is::<f32>() {
        return Ok(dataset
            .read_slice_1d::<f32, _>(start..end)?
            .into_iter()
            .map(f64::from)
            .collect());
    }
    Ok(read_integer_slice(dataset, start, end)?
        .into_iter()
        .map(|x| x as f64)
        .collect())
}

fn read_dense_slice(dataset: &Dataset, start: usize, end: usize) -> Result<Vec<f64>> {
    let dtype = dataset.dtype()?;
    if dtype.is::<f64>() {
        return Ok(dataset
            .read_slice_2d::<f64, _>((start..end, ..))?
            .iter()
            .copied()
            .collect());
    }
    if dtype.is::<f32>() {
        return Ok(dataset
            .read_slice_2d::<f32, _>((start..end, ..))?
            .iter()
            .copied()
            .map(f64::from)
            .collect());
    }
    macro_rules! dense_int {
        ($ty:ty) => {
            if dtype.is::<$ty>() {
                return Ok(dataset
                    .read_slice_2d::<$ty, _>((start..end, ..))?
                    .iter()
                    .copied()
                    .map(|x| x as f64)
                    .collect());
            }
        };
    }
    dense_int!(i8);
    dense_int!(u8);
    dense_int!(i16);
    dense_int!(u16);
    dense_int!(i32);
    dense_int!(u32);
    dense_int!(i64);
    dense_int!(u64);
    bail!("dense matrix has unsupported datatype")
}

fn read_bool_dataset(dataset: &Dataset) -> Result<Vec<bool>> {
    if dataset.dtype()?.is::<bool>() {
        return Ok(dataset.read_raw::<bool>()?);
    }
    Ok(read_integer_dataset(dataset)?
        .into_iter()
        .map(|value| value != 0)
        .collect())
}

fn is_string(dataset: &Dataset) -> Result<bool> {
    Ok(matches!(
        dataset.dtype()?.to_descriptor()?,
        TypeDescriptor::VarLenUnicode
            | TypeDescriptor::VarLenAscii
            | TypeDescriptor::FixedUnicode(_)
            | TypeDescriptor::FixedAscii(_)
    ))
}

fn format_number(value: f64) -> String {
    if value.is_nan() {
        "NA".into()
    } else {
        value.to_string()
    }
}

fn column_description(group: &Group, name: &str) -> Result<String> {
    if group.loc_type_by_name(name)? == LocationType::Group {
        let column = group.group(name)?;
        return Ok(attr_string(&column, "encoding-type").unwrap_or_else(|_| "group".into()));
    }
    let dataset = group.dataset(name)?;
    Ok(format!(
        "{:?}[{}]",
        dataset.dtype()?.to_descriptor()?,
        dataset.size()
    ))
}

fn matrix_description(file: &File, path: &str) -> Result<String> {
    if file.loc_type_by_name(path)? == LocationType::Dataset {
        let dataset = file.dataset(path)?;
        return Ok(format!("dense {:?}", dataset.shape()));
    }
    let matrix = file.group(path)?;
    Ok(format!(
        "{} {:?}",
        attr_string(&matrix, "encoding-type")?,
        read_shape(&matrix)?
    ))
}
