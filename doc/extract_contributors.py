#!/usr/bin/env python3
"""
extract_contributors.py

Scan source files for RooUnfold-style copyright headers and build a CONTRIBUTORS.md.

Usage:
  python extract_contributors.py [options] FILE [FILE ...]

Examples:
  python extract_contributors.py -o CONTRIBUTORS.md src/*.cxx python/*.py

Enhancements:
- Deduplicates and normalizes years/ranges across files per author.
- Computes involvement span per author (earliest..latest, inclusive).
- Sorts authors by involvement span (desc), then by #files (desc), then by name.
- Always lists explicit file paths each author appears in.

Notes:
- The script looks for a header block delimited by
    "BEGIN ROOUNFOLD COPYRIGHT" ... "END ROOUNFOLD COPYRIGHT"
  and then parses lines like:
    " - Full Name (2005-2007, 2009–2011, 2014, 2017)"
    " - Full Name <email@host> (2019–2020)"
- If the block is missing but a line starting with "Authors" is present,
  it will still try to parse following bullet lines.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

HEADER_BEGIN = "BEGIN ROOUNFOLD COPYRIGHT"
HEADER_END = "END ROOUNFOLD COPYRIGHT"

AUTHOR_LINE_RE = re.compile(
    r"""
    ^\s*                # leading space
    (?:[*#/]+\s*)?      # possible comment leader chars like *, #, /
    -\s*                # dash bullet
    (?P<name>[^<(]+?)   # author name (no '<' or '('), non-greedy
    \s*                 # optional whitespace
    (?:<(?P<email>[^>]+)>)?  # optional <email>
    \s*                 # optional whitespace
    \((?P<years>[^)]*)\)     # parenthesized years/ranges
    \s*$
    """,
    re.VERBOSE,
)

def extract_header_block(text: str) -> Optional[str]:
    """Return the ROOUNFOLD COPYRIGHT block if present, else None."""
    m = re.search(
        rf"{re.escape(HEADER_BEGIN)}(.*?){re.escape(HEADER_END)}",
        text,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if m:
        return m.group(1)
    # Fallback: try to find 'Authors' line and take following chunk as candidates
    m2 = re.search(r"Authors\s*\(.*?\):?(.*)$", text, flags=re.IGNORECASE | re.MULTILINE)
    if not m2:
        return None
    start = m2.end()
    pseudo = text[start : start + 4000]
    return pseudo

def parse_authors_from_text(text: str) -> List[Tuple[str, Optional[str], str]]:
    """
    Parse authors from the given text chunk.
    Returns a list of tuples: (name, email_or_None, years_string)
    """
    authors: List[Tuple[str, Optional[str], str]] = []
    for line in text.splitlines():
        m = AUTHOR_LINE_RE.match(line)
        if m:
            name = m.group("name").strip()
            email = m.group("email").strip() if m.group("email") else None
            years = m.group("years").strip()
            years = re.sub(r"\s+", " ", years)
            authors.append((name, email, years))
    return authors

def _tokenize_years(ystr: str) -> Set[int]:
    """
    Convert a mixed string of years/ranges into a set of int years.
    Supports separators ',', ';'. Supports hyphen '-' and en-dash '–' ranges.
    Ignores non-4-digit tokens.
    """
    years: Set[int] = set()
    # split on comma/semicolon
    for token in re.split(r"[;,]", ystr):
        token = token.strip()
        if not token:
            continue
        # normalize dash types to en-dash for detection, but handle '-' as well
        # Try range first
        m = re.match(r"^\s*(\d{4})\s*[\-–]\s*(\d{4})\s*$", token)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            if a <= b:
                years.update(range(a, b + 1))
            else:
                years.update(range(b, a + 1))
            continue
        # single year
        m = re.match(r"^\s*(\d{4})\s*$", token)
        if m:
            years.add(int(m.group(1)))
            continue
        # catch weird forms with extra words; find all 4-digit numbers
        for y in re.findall(r"\d{4}", token):
            years.add(int(y))
    return years

def normalize_years_union(years_strings: Set[str]) -> Tuple[str, Optional[int], Optional[int], int]:
    """
    Take a set of years strings from multiple files, union them, and recompact
    into non-overlapping sorted ranges with en-dashes.

    Returns (pretty, earliest, latest, span_years_inclusive).
    If no valid years found, pretty == "" and span == 0.
    """
    all_years: Set[int] = set()
    for ystr in years_strings:
        all_years.update(_tokenize_years(ystr))

    if not all_years:
        return ("", None, None, 0)

    ys = sorted(all_years)
    # compress into ranges
    ranges: List[Tuple[int, int]] = []
    start = prev = ys[0]
    for y in ys[1:]:
        if y == prev + 1:
            prev = y
        else:
            ranges.append((start, prev))
            start = prev = y
    ranges.append((start, prev))

    # pretty print
    parts = []
    for a, b in ranges:
        if a == b:
            parts.append(f"{a}")
        else:
            parts.append(f"{a}–{b}")
    pretty = ", ".join(parts)
    earliest, latest = ys[0], ys[-1]
    span = latest - earliest + 1  # inclusive span
    return (pretty, earliest, latest, span)

def collect_from_files(files: List[Path]):
    """
    Collect authors across files.
    Returns dictionaries and a derived per-author info map.
    """
    author_to_files: Dict[str, Set[str]] = defaultdict(set)
    author_to_emails: Dict[str, Set[str]] = defaultdict(set)
    author_to_years_raw: Dict[str, Set[str]] = defaultdict(set)

    for p in files:
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            print(f"[WARN] Could not read {p}: {e}")
            continue

        block = extract_header_block(text)
        if not block:
            continue

        for name, email, years in parse_authors_from_text(block):
            author_to_files[name].add(str(p))
            if email:
                author_to_emails[name].add(email)
            if years:
                author_to_years_raw[name].add(years)

    # Build derived info
    info = {}
    for name in author_to_files:
        emails = author_to_emails.get(name, set())
        years_raw = author_to_years_raw.get(name, set())
        pretty_years, earliest, latest, span = normalize_years_union(years_raw)
        info[name] = {
            "emails": sorted(emails),
            "files": sorted(author_to_files[name]),
            "years_pretty": pretty_years,
            "earliest": earliest,
            "latest": latest,
            "span": span,
        }
    return info

def write_contributors_md(out_path: Path, info: Dict[str, dict]) -> None:
    """Write a CONTRIBUTORS.md file summarizing authors across files, sorted by span."""
    # Sort by span desc, then by #files desc, then name
    ordering = sorted(
        info.items(),
        key=lambda kv: (- (kv[1]["span"] or 0), -len(kv[1]["files"]), kv[0].split()[-1].lower(), kv[0].lower()),
    )

    lines: List[str] = []
    lines.append("# Contributors\n")
    lines.append(f"Total: **{len(ordering)}**\n")
    lines.append("_Sorted by involvement span (earliest↔latest), then by number of files._\n")

    for name, meta in ordering:
        emails = meta["emails"]
        files = meta["files"]
        years_pretty = meta["years_pretty"]
        earliest = meta["earliest"]
        latest = meta["latest"]
        span = meta["span"]

        header = f"- **{name}**"
        if emails:
            header += f" — {', '.join(emails)}"
        lines.append(header)

        # Years + span
        if years_pretty:
            if earliest is not None and latest is not None:
                lines.append(f"  - Years: {years_pretty}  *(span {earliest}–{latest}, {span} year{'s' if span!=1 else ''})*")
            else:
                lines.append(f"  - Years: {years_pretty}")
        else:
            lines.append(f"  - Years: *(not available)*")

        # Files (explicit list)
        lines.append(f"  - Files ({len(files)}):")
        for f in files:
            lines.append(f"    - `{f}`")
        lines.append("")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out_path}")

def main():
    ap = argparse.ArgumentParser(description="Extract contributors from RooUnfold-style headers.")
    ap.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="List of source files to scan (Python, C/C++).",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("CONTRIBUTORS.md"),
        help="Output file path (default: CONTRIBUTORS.md).",
    )
    args = ap.parse_args()

    info = collect_from_files(args.files)
    if not info:
        print("No contributors found in the provided files.")
    write_contributors_md(args.output, info)

if __name__ == "__main__":
    main()
