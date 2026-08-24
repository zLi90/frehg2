#!/usr/bin/env python3
"""check_parameter_docs.py — keep docs/user-guide/parameters.md in lockstep
with the schema in src/core/ConfigSchema.cpp (plan §5.4: the key reference
is CI-checked against the schema; §11.1 rule 3: new keys land in schema,
docs, and tests in the same PR).

Checks two directions:
  1. every quoted lowercase token in buildRootSchema() (key names and
     allowed enum values) appears somewhere in the parameters table;
  2. every dotted-path segment documented in the table exists as a schema
     token (catches stale docs after a schema rename).
Exits nonzero with the offending names on any mismatch.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCHEMA = ROOT / "src/core/ConfigSchema.cpp"
DOCS = ROOT / "docs/user-guide/parameters.md"

# Tokens that appear in the docs' key paths but are structural markers or
# YAML literals rather than schema key names.
DOC_ONLY_TOKENS = {"auto", "true", "constant,file"}


def schema_tokens():
    text = SCHEMA.read_text(encoding="utf-8")
    begin = text.index("Spec buildRootSchema()")
    end = text.index("// ------", begin)
    body = text[begin:end]
    return set(re.findall(r'"([a-z][a-z0-9_]*)"', body))


def doc_tokens():
    """Two token sets from the parameter tables:
    - key-path segments from the first column (used for both directions);
    - any backticked lowercase token in a table row (used only to prove
      schema coverage: enum values are documented in the Type column, and
      description columns may mention legacy names)."""
    text = DOCS.read_text(encoding="utf-8")
    path_tokens = set()
    any_tokens = set()
    for line in text.splitlines():
        if not line.startswith("| `"):
            continue
        first_cell = line[1:].split(" | ")[0].strip()
        match = re.match(r"`([^`]+)`", first_cell)
        if match:
            for piece in re.split(r"[.\[\]{}|,]", match.group(1)):
                piece = piece.strip()
                if piece and re.fullmatch(r"[a-z][a-z0-9_]*", piece):
                    path_tokens.add(piece)
        for span in re.findall(r"`([^`]+)`", line):
            for piece in re.split(r"[.\[\]{}|,]", span):
                piece = piece.strip()
                if piece and re.fullmatch(r"[a-z][a-z0-9_]*", piece):
                    any_tokens.add(piece)
    return path_tokens, any_tokens


def main():
    in_schema = schema_tokens()
    path_tokens, any_tokens = doc_tokens()

    undocumented = sorted(in_schema - any_tokens)
    stale = sorted(path_tokens - in_schema - DOC_ONLY_TOKENS)

    status = 0
    if undocumented:
        print("schema tokens missing from docs/user-guide/parameters.md:")
        for name in undocumented:
            print(f"  - {name}")
        status = 1
    if stale:
        print("documented tokens absent from the schema (stale docs?):")
        for name in stale:
            print(f"  - {name}")
        status = 1
    if status == 0:
        print(f"check_parameter_docs: OK ({len(in_schema)} schema tokens documented)")
    return status


if __name__ == "__main__":
    sys.exit(main())
