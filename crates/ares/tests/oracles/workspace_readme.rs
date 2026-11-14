//! The README crate table names every workspace member.
//!
//! The workspace is `members = ["crates/*"]`, so a crate added under `crates/`
//! joins the build without touching any list a reader sees first. This module
//! reads the manifests and the README table and compares their (directory,
//! registry name) pairs, so a missing row, a stale row, a duplicated row and a
//! wrong registry name each fail with both sides printed.
//!
//! The comparison reads files above this package, so it runs from the
//! workspace checkout and not from an unpacked registry tarball.

use std::fs;
use std::path::{Path, PathBuf};

/// A workspace member: `(directory under crates/, registry name)`.
type Member = (String, String);

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
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

/// The rows of the `## Crates` table, each `(directory, registry name)`.
///
/// A row is a line whose first cell links `crates/<directory>` and whose
/// second cell is the backticked registry name. The header, the separator and
/// any table under another heading are not rows.
fn crate_table_rows(readme: &str) -> Vec<Member> {
    readme
        .lines()
        .skip_while(|line| line.trim_end() != "## Crates")
        .skip(1)
        .take_while(|line| !line.starts_with("## "))
        .filter_map(|line| {
            let cells: Vec<&str> = line.split('|').map(str::trim).collect();
            let [_, link, registry, ..] = cells.as_slice() else {
                return None;
            };
            let directory = link.split_once("](crates/")?.1.strip_suffix(')')?;
            Some((directory.to_owned(), registry.trim_matches('`').to_owned()))
        })
        .collect()
}

fn workspace_members(root: &Path) -> Vec<Member> {
    let mut members: Vec<Member> = fs::read_dir(root.join("crates"))
        .expect("the workspace has a crates directory")
        .map(|entry| {
            let entry = entry.expect("a crates entry is readable");
            let directory = entry
                .file_name()
                .into_string()
                .expect("crate directory names are UTF-8");
            let manifest = fs::read_to_string(entry.path().join("Cargo.toml"))
                .expect("every crates entry holds a manifest");
            let name = package_name(&manifest).expect("every manifest names its package");
            (directory, name)
        })
        .collect();
    members.sort();
    members
}

#[test]
fn readme_crate_table_matches_the_workspace_members() {
    let root = workspace_root();
    let readme = fs::read_to_string(root.join("README.md")).expect("the workspace has a README");
    let mut rows = crate_table_rows(&readme);
    rows.sort();
    assert_eq!(rows, workspace_members(&root));
}

#[test]
fn package_name_reads_the_package_table_not_the_library_table() {
    let package_first = "[package]\nname = \"ares-solid\"\n\n[lib]\nname = \"ares\"\n";
    let library_first = "[lib]\nname = \"ares\"\n\n[package]\nname = \"ares-solid\"\n";
    assert_eq!(package_name(package_first).as_deref(), Some("ares-solid"));
    assert_eq!(package_name(library_first).as_deref(), Some("ares-solid"));
    assert_eq!(package_name("[lib]\nname = \"ares\"\n"), None);
}

#[test]
fn crate_table_rows_skips_the_header_the_separator_and_other_sections() {
    let readme = "\
# Title

## Crates

| Crate | Registry name | What it is |
| --- | --- | --- |
| [`crates/ares`](crates/ares) | `ares-solid` | The core. |
| [`crates/ares-operator`](crates/ares-operator) | `ares-operator` | The seam. |

## Decisions

| [`crates/elsewhere`](crates/elsewhere) | `elsewhere` | Not a member row. |
";
    assert_eq!(
        crate_table_rows(readme),
        vec![
            ("ares".to_owned(), "ares-solid".to_owned()),
            ("ares-operator".to_owned(), "ares-operator".to_owned()),
        ]
    );
}
