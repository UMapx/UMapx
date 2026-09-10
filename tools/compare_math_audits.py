"""Compare every test ID in two TRX files; retain evidence for resolved defects."""
import argparse
import collections
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path


def read(path):
    root = ET.parse(path).getroot()
    ns = {'t': root.tag.split('}')[0][1:]}
    suites = {test.attrib['id']: test.find('t:TestMethod', ns).attrib['className'].split('.')[-1]
              for test in root.findall('./t:TestDefinitions/t:UnitTest', ns)}
    return {row.attrib['testId']: dict(id=row.attrib['testId'], suite=suites[row.attrib['testId']],
                                     case=row.attrib['testName'], outcome=row.attrib['outcome'])
            for row in root.findall('./t:Results/t:UnitTestResult', ns)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--current', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    before, after = read(args.baseline), read(args.current)
    common = before.keys() & after.keys()
    resolved = [after[k] for k in common if before[k]['outcome'] == 'Failed' and after[k]['outcome'] == 'Passed']
    regressed = [after[k] for k in common if before[k]['outcome'] == 'Passed' and after[k]['outcome'] != 'Passed']
    added = [after[k] for k in after.keys() - before.keys()]
    removed = [before[k] for k in before.keys() - after.keys()]
    result = dict(baseline_cases=len(before), current_cases=len(after), common_cases=len(common),
                  baseline_totals=dict(collections.Counter(row['outcome'] for row in before.values())),
                  current_totals=dict(collections.Counter(row['outcome'] for row in after.values())),
                  added_totals=dict(collections.Counter(row['outcome'] for row in added)),
                  resolved_by_suite=dict(sorted(collections.Counter(row['suite'] for row in resolved).items())),
                  resolved_count=len(resolved), regression_count=len(regressed), added_count=len(added), removed_count=len(removed),
                  resolved=sorted(resolved, key=lambda row: (row['suite'], row['case'])),
                  regressions=regressed, removed=removed,
                  evidence_sha256={name: hashlib.sha256(path.read_bytes()).hexdigest()
                                   for name, path in [('baseline_trx', args.baseline), ('current_trx', args.current)]})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ['resolved', 'regressions', 'removed', 'evidence_sha256']}, indent=2))


if __name__ == '__main__':
    main()
