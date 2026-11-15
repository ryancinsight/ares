//! The README crate table names every workspace member.
//!
//! The workspace is `members = ["crates/*"]`, so a crate added under `crates/`
//! joins the build without touching any list a reader sees first. This test
//! reads the manifests and the README table and compares their (directory,
//! registry name) pairs, so a missing row, a stale row, a duplicated row and a
//! wrong registry name each fail with both sides printed. Every line of the
//! table must parse as a crate row, so a stale row in any other shape fails
//! too, naming the line.
//!
//! The comparison reads the workspace README and the sibling manifests above
//! this package. It is a repository-level check and its own test binary, kept
//! out of the published package by the `exclude` key in this crate's manifest:
//! an unpacked registry tarball or a vendored copy has no such layout, and the
//! analytical `oracles` harness is not where a documentation check belongs.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};

/// A workspace member: `(directory under crates/, registry name)`.
type Member = (String, String);

/// Why the README table cannot be read as a list of crate rows.
#[derive(Debug, PartialEq, Eq)]
enum TableError {
    /// There is no `## Crates` heading, or no table under it.
    MissingTable,
    /// The line after the table header is not a delimiter row; it is carried.
    MissingSeparator(String),
    /// A table line that is not a crate row; it is carried.
    MalformedRow(String),
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn invalid_data(path: &Path, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("{}: {message}", path.display()),
    )
}

fn with_path(path: &Path, source: &io::Error) -> io::Error {
    io::Error::new(source.kind(), format!("{}: {source}", path.display()))
}

fn read(path: &Path) -> io::Result<String> {
    fs::read_to_string(path).map_err(|source| with_path(path, &source))
}

/// The `name` key of the `[package]` table, which is the registry name.
///
/// `[lib]` and `[[bin]]` tables carry their own `name`, which is an import
/// path and not what a registry resolves.
fn package_name(manifest: &str) -> Option<String> {
    let mut in_package = false;
    for line in manifest.lines().map(str::trim) {
        if line.starts_with('[') {
            in_package = line == "[package]";
        } else if in_package
            && let Some((key, value)) = line.split_once('=')
            && key.trim() == "name"
        {
            return Some(value.trim().trim_matches('"').to_owned());
        }
    }
    None
}

/// Whether `line` is a table delimiter row such as `| --- | :-: |`.
fn is_separator(line: &str) -> bool {
    let line = line.trim();
    let cells: Vec<&str> = line
        .split('|')
        .map(str::trim)
        .filter(|cell| !cell.is_empty())
        .collect();
    line.starts_with('|')
        && !cells.is_empty()
        && cells.iter().all(|cell| {
            let dashes = cell.trim_matches(':');
            !dashes.is_empty() && dashes.chars().all(|c| c == '-')
        })
}

/// One crate row: ``| [`crates/<dir>`](crates/<dir>) | `<registry name>` | ... |``.
///
/// The link text and the link target must name the same directory, and the
/// registry name must be backticked, so a row in any other shape is an error.
fn crate_row(line: &str) -> Result<Member, TableError> {
    let malformed = || TableError::MalformedRow(line.to_owned());
    let cells: Vec<&str> = line.trim().split('|').map(str::trim).collect();
    let ["", crate_cell, registry_cell, ..] = cells.as_slice() else {
        return Err(malformed());
    };
    let (text, target) = crate_cell
        .strip_prefix("[`crates/")
        .and_then(|rest| rest.split_once("`](crates/"))
        .ok_or_else(malformed)?;
    let directory = target
        .strip_suffix(')')
        .filter(|directory| !directory.is_empty() && *directory == text)
        .ok_or_else(malformed)?;
    let name = registry_cell
        .strip_prefix('`')
        .and_then(|rest| rest.strip_suffix('`'))
        .filter(|name| !name.is_empty())
        .ok_or_else(malformed)?;
    Ok((directory.to_owned(), name.to_owned()))
}

/// The rows of the table under the `## Crates` heading.
///
/// The table runs from its header line to the first blank line; the header and
/// the delimiter row are not crate rows, and every line after them must be.
/// A table under any other heading is not read.
fn crate_table_rows(readme: &str) -> Result<Vec<Member>, TableError> {
    let mut table = readme
        .lines()
        .skip_while(|line| line.trim_end() != "## Crates")
        .skip(1)
        .take_while(|line| !line.starts_with("## "))
        .skip_while(|line| !line.starts_with('|'));
    table.next().ok_or(TableError::MissingTable)?;
    let separator = table.next().unwrap_or_default();
    if !is_separator(separator) {
        return Err(TableError::MissingSeparator(separator.to_owned()));
    }
    table
        .take_while(|line| !line.trim().is_empty())
        .map(crate_row)
        .collect()
}

fn workspace_members(root: &Path) -> io::Result<Vec<Member>> {
    let crates = root.join("crates");
    let mut members = Vec::new();
    for entry in fs::read_dir(&crates).map_err(|source| with_path(&crates, &source))? {
        let entry = entry.map_err(|source| with_path(&crates, &source))?;
        let path = entry.path();
        if !entry
            .file_type()
            .map_err(|source| with_path(&path, &source))?
            .is_dir()
        {
            continue;
        }
        let directory = entry
            .file_name()
            .into_string()
            .map_err(|_| invalid_data(&path, "the directory name is not UTF-8"))?;
        let manifest_path = path.join("Cargo.toml");
        let name = package_name(&read(&manifest_path)?)
            .ok_or_else(|| invalid_data(&manifest_path, "no `name` in the `[package]` table"))?;
        members.push((directory, name));
    }
    members.sort();
    Ok(members)
}

#[test]
fn readme_crate_table_matches_the_workspace_members() -> io::Result<()> {
    let root = workspace_root();
    let readme_path = root.join("README.md");
    let mut rows = crate_table_rows(&read(&readme_path)?)
        .map_err(|error| invalid_data(&readme_path, &format!("{error:?}")))?;
    rows.sort();
    assert_eq!(rows, workspace_members(&root)?);
    Ok(())
}

#[test]
fn package_name_reads_the_package_table_not_the_library_table() {
    let package_first = "[package]\nname = \"ares-solid\"\n\n[lib]\nname = \"ares\"\n";
    let library_first = "[lib]\nname = \"ares\"\n\n[package]\nname = \"ares-solid\"\n";
    assert_eq!(package_name(package_first).as_deref(), Some("ares-solid"));
    assert_eq!(package_name(library_first).as_deref(), Some("ares-solid"));
    assert_eq!(package_name("[lib]\nname = \"ares\"\n"), None);
}

const HEADER: &str = "| Crate | Registry name | What it is |\n| --- | --- | --- |\n";
const CORE_ROW: &str = "| [`crates/ares`](crates/ares) | `ares-solid` | The core. |\n";

fn readme_with_rows(rows: &str) -> String {
    format!("# Title\n\n## Crates\n\n{HEADER}{rows}\nProse after the table.\n")
}

#[test]
fn crate_table_rows_reads_the_rows_and_stops_at_the_blank_line() {
    let readme = format!(
        "{}\n## Decisions\n\n| [`crates/elsewhere`](crates/elsewhere) | `elsewhere` | Not a row. |\n",
        readme_with_rows(&format!(
            "{CORE_ROW}| [`crates/ares-operator`](crates/ares-operator) | `ares-operator` | The seam. |\n"
        ))
    );
    assert_eq!(
        crate_table_rows(&readme),
        Ok(vec![
            ("ares".to_owned(), "ares-solid".to_owned()),
            ("ares-operator".to_owned(), "ares-operator".to_owned()),
        ])
    );
}

#[test]
fn crate_table_rows_of_a_header_only_table_is_empty() {
    assert_eq!(crate_table_rows(&readme_with_rows("")), Ok(Vec::new()));
}

#[test]
fn crate_table_rows_does_not_read_a_table_under_another_heading() {
    let readme =
        format!("# Title\n\n## Crates\n\nNo table here.\n\n## Scope\n\n{HEADER}{CORE_ROW}");
    assert_eq!(crate_table_rows(&readme), Err(TableError::MissingTable));
    assert_eq!(
        crate_table_rows("# Title\n\n## Scope\n"),
        Err(TableError::MissingTable)
    );
}

#[test]
fn crate_table_rows_requires_the_delimiter_row() {
    let readme = format!("## Crates\n\n| Crate | Registry name |\n{CORE_ROW}");
    assert_eq!(
        crate_table_rows(&readme),
        Err(TableError::MissingSeparator(CORE_ROW.trim_end().to_owned()))
    );
}

#[test]
fn crate_table_rows_rejects_a_row_in_any_other_shape() {
    let stale_rows = [
        "| `ares-gone` | `ares-gone` | Stale. |",
        "| [`crates/ares-gone`](./crates/ares-gone) | `ares-gone` | Stale. |",
        "| [`crates/ares-gone`](crates/ares-other) | `ares-gone` | Text and target differ. |",
        "| [`crates/ares-gone`](crates/ares-gone) | ares-gone | Registry name not backticked. |",
        "| [`crates/ares-gone`](crates/ares-gone) |",
        "[`crates/ares-gone`](crates/ares-gone) | `ares-gone` | No leading pipe. |",
    ];
    for stale in stale_rows {
        let readme = readme_with_rows(&format!("{CORE_ROW}{stale}\n"));
        assert_eq!(
            crate_table_rows(&readme),
            Err(TableError::MalformedRow(stale.to_owned())),
            "{stale}"
        );
    }
}
