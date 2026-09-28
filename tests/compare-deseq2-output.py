#!/usr/bin/env python3
"""
Compare the first n lines of two DESeq2 output files.
Uses pytest.approx-like functionality to compare floating-point numbers.
"""

import argparse
import csv
import math
import os.path
import sys
from pathlib import Path


def compare_floats(a, b, rel=1e-4, abs_tol=1e-12):
    """Compare two floats, equivalent to `a == pytest.approx(b, nan_ok=True, rel=1e-4)`."""
    if a == b:
        return True
    if math.isnan(b):
        return math.isnan(a)
    if math.isinf(b):
        return False
    tolerance = max(rel * abs(b), abs_tol)
    return abs(b - a) <= tolerance

def read_top_n_lines(file_path, n=10):
    """Read the first n lines from a TSV file."""
    with open(file_path, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        lines = []
        for i, row in enumerate(reader):
            if i >= n:
                break
            lines.append(row)
    return lines

def parse_numeric(value):
    """Convert a string to a number (float or int) if possible, otherwise return the string."""
    if value == '' or value.lower() in ('na', 'nan', 'null', 'none'):
        return None
    try:
        # Handle very small numbers like 0.000...0005
        if 'e' in value.lower():
            return float(value)
        return float(value)
    except ValueError:
        return value

def compare_deseq2_files(file1, file2, rel_tol=1e-6, abs_tol=1e-12, top_n=9):
    """
    Compare the first n lines of two DESeq2 output files.

    Args:
        file1: Path to the first file
        file2: Path to the second file
        rel_tol: Relative tolerance for float comparison
        abs_tol: Absolute tolerance for float comparison
        top_n: Number of lines to compare

    Returns:
        tuple: (are_equal, report) where are_equal is a boolean and report is a string
    """
    # Read files
    try:
        lines1 = read_top_n_lines(file1, top_n)
        lines2 = read_top_n_lines(file2, top_n)
    except Exception as e:
        return False, f"Error reading files: {e}"

    if not lines1 or not lines2:
        return False, "One or both files are empty"

    # Find the ID column
    id_key1 = None
    id_key2 = None

    for key in lines1[0].keys():
        if key == 'Id':
            id_key1 = key
            break

    for key in lines2[0].keys():
        if key == 'Id':
            id_key2 = key
            break

    if not id_key1 or not id_key2:
        return False, "ID column not found in one or both files"

    # Create dictionaries indexed by ID
    dict1 = {row[id_key1]: row for row in lines1}
    dict2 = {row[id_key2]: row for row in lines2}

    # Find common IDs
    common_ids = set(dict1.keys()) & set(dict2.keys())
    all_ids = set(dict1.keys()) | set(dict2.keys())

    missing_in_1 = all_ids - set(dict1.keys())
    missing_in_2 = all_ids - set(dict2.keys())

    report = []
    are_equal = True

    # Report missing IDs
    if missing_in_1:
        report.append(f"IDs present in {file2} but not in {file1}: {sorted(missing_in_1)}")
        are_equal = False
    if missing_in_2:
        report.append(f"IDs present in {file1} but not in {file2}: {sorted(missing_in_2)}")
        are_equal = False

    # Compare common entries
    differences = []

    for id_ in sorted(common_ids):
        row1 = dict1[id_]
        row2 = dict2[id_]

        # Find common numeric columns (excluding the second "ID" column)
        all_keys = set(row1.keys()) | set(row2.keys())
        numeric_keys = set()

        # Determine which columns are numeric (excluding "ID")
        for key in all_keys:
            if key == 'ID' or key == 'dispersions.dds.':  # Skip the second ID and dispersions.dds. columns
                continue
            if key in row1 and key in row2:
                val1 = parse_numeric(row1[key])
                val2 = parse_numeric(row2[key])
                if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                    numeric_keys.add(key)

        # Compare each numeric column
        for key in sorted(numeric_keys):
            val1 = parse_numeric(row1.get(key, None))
            val2 = parse_numeric(row2.get(key, None))

            if val1 is None and val2 is None:
                continue

            # Comparison with tolerance
            if not compare_floats(val1, val2):
                differences.append(f"  {id_}.{key}: {val1} vs {val2}")
                are_equal = False

    if differences:
        report.append("\nNumeric value differences:")
        report.append("\n".join(differences))

    if are_equal:
        report.append(f"The {len(common_ids)} common entries are equal (within the specified tolerance).")

    return are_equal, "\n".join(report)

def main():

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description="Compare two DESeq2 output files.",
    )
    parser.add_argument(
        "-c", "--line-count", default=10, help="Number of lines to compare (default: 10)", type=int
    )
    parser.add_argument("file1", help="First DESeq2 output file")
    parser.add_argument("file2", help="Second DESeq2 output file")

    args = parser.parse_args()

    file1 = Path(args.file1)
    file2 = Path(args.file1)

    if not file1.exists() or not file2.exists():
        print(f"Error: One or both files do not exist ({file1} or {file2}).")
        sys.exit(1)

    are_equal, report = compare_deseq2_files(file1, file2, args.line_count)

    if not are_equal:
        print(f"\n❌ Differences detected between the two files ({file1} vs {file2}): ")
        print(report)

    sys.exit(0 if are_equal else 1)

if __name__ == "__main__":
    main()
