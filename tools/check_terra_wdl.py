#!/usr/bin/env python3
"""Reject workflow-scope file writers in WDL documents and their imports."""
import argparse
import sys

import WDL


def workflow_writers(document):
    findings = []

    def visit(node):
        if isinstance(node, WDL.Expr.Apply) and node.function_name.startswith('write_'):
            findings.append((node.pos.line, node.function_name))
        for child in node.children:
            visit(child)

    if document.workflow:
        visit(document.workflow)
    return findings


def check_document(document, visited):
    if document.pos.abspath in visited:
        return []
    visited.add(document.pos.abspath)
    findings = [(document.pos.uri, line, function) for line, function in workflow_writers(document)]
    for imported in document.imports:
        findings.extend(check_document(imported.doc, visited))
    return findings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wdl', nargs='+')
    args = parser.parse_args()
    failures = []
    for path in args.wdl:
        document = WDL.load(path)
        if document.effective_wdl_version != '1.0':
            failures.append((path, 1, 'WDL version must be 1.0'))
        failures.extend(check_document(document, set()))
    for path, line, function in failures:
        print(f'{path}:{line}: workflow-scope file writer is not allowed: {function}', file=sys.stderr)
    if failures:
        return 1
    print('Terra static check passed. No workflow-scope file writers found.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
