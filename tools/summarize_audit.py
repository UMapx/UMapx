"""Export reproducible audit evidence from TRX and coverlet; Python standard library only.

Coverage means execution, not mathematical correctness. Files without sequence
points are listed explicitly rather than assigned a fictitious 100% score.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trx', type=Path, required=True)
    parser.add_argument('--coverage', type=Path, required=True, help='coverage.cobertura.xml')
    parser.add_argument('--output', type=Path, default=ROOT / 'docs' / 'audit')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    coverage = ET.parse(args.coverage).getroot()
    lines = collections.defaultdict(dict)
    for cls in coverage.findall('.//class'):
        name = cls.attrib['filename'].replace('\\', '/')
        for line in cls.findall('./lines/line'):
            number, hits = int(line.attrib['number']), int(line.attrib['hits'])
            lines[name][number] = max(lines[name].get(number, 0), hits)
    files = sorted(p for p in (ROOT / 'sources').rglob('*.cs') if not {'bin', 'obj'} & set(p.relative_to(ROOT).parts))
    inventory = []
    for path in files:
        name = path.relative_to(ROOT / 'sources').as_posix()
        points = lines.get(name, {})
        covered = sum(hits > 0 for hits in points.values())
        role = 'algorithm or executable support'
        if not points:
            role = 'no instrumented sequence points (interface, enum, delegate, or conditional code)'
        elif name.startswith('Video/'):
            role = 'in-memory video checks; live HTTP/device/capture integration excluded'
        elif name == 'Core/Debugger.cs':
            role = 'console diagnostic utility; no numerical tests'
        elif name == 'Core/Xml.cs':
            role = 'serialization contract; no mathematical algorithm'
        inventory.append(dict(file='sources/' + name, area=name.split('/')[0],
                              executable_lines=len(points), covered_lines=covered,
                              line_percent=round(100*covered/len(points), 2) if points else None,
                              uncovered_lines=sorted(n for n,hits in points.items() if hits == 0),
                              scope=role, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    dump(args.output / 'source-inventory.json', inventory)
    with (args.output / 'source-inventory.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(inventory[0]))
        writer.writeheader()
        writer.writerows({**row, 'uncovered_lines': ' '.join(map(str,row['uncovered_lines']))} for row in inventory)

    coverlet = json.loads(args.coverage.with_name('coverage.json').read_text(encoding='utf-8'))
    methods = []
    for module in coverlet.values():
        for filename, classes in module.items():
            relative = filename.replace('\\', '/').split('/sources/')[-1]
            for classname, members in classes.items():
                for signature, detail in members.items():
                    points = detail.get('Lines', {})
                    if points and not any(points.values()):
                        methods.append(dict(file='sources/'+relative, type=classname, method=signature, line=min(map(int,points))))
    dump(args.output / 'uncovered-methods.json', sorted(methods,key=lambda row:(row['file'],row['line'],row['method'])))

    trx = ET.parse(args.trx).getroot()
    ns = {'t': trx.tag.split('}')[0][1:]}
    definitions = {}
    for test in trx.findall('./t:TestDefinitions/t:UnitTest', ns):
        method = test.find('t:TestMethod', ns)
        definitions[test.attrib['id']] = (method.attrib['className'].split('.')[-1], method.attrib['name'])
    suites = collections.defaultdict(collections.Counter)
    failures = []
    for result in trx.findall('./t:Results/t:UnitTestResult', ns):
        suite, method = definitions[result.attrib['testId']]
        outcome = result.attrib['outcome']
        suites[suite][outcome] += 1
        if outcome != 'Passed':
            message = result.findtext('./t:Output/t:ErrorInfo/t:Message', '', ns).replace('\r','')
            stack = result.findtext('./t:Output/t:ErrorInfo/t:StackTrace', '', ns)
            locations = re.findall(r'(tests[\\/]UMapx.Tests[\\/][^:]+):line (\d+)', stack)
            location = next((item for item in locations if item[0].endswith(suite+'.cs')), locations[0] if locations else None)
            failures.append(dict(id=result.attrib['testId'], suite=suite, method=method,
                                 case=result.attrib['testName'], outcome=outcome, message=message[:3000],
                                 test_file=location[0].replace('\\','/') if location else '',
                                 test_line=int(location[1]) if location else None))
    failures.sort(key=lambda row:(row['suite'],row['method'],row['case']))
    dump(args.output / 'failures.json', failures)
    grouped = collections.defaultdict(list)
    for row in failures:
        # Reference parameters are retained in failures.json and the embedded data files.
        grouped[row['suite']+'.'+row['method']].append(row)
    md = ['# Failed audit checks', '',
          'Generated from the recorded TRX. These are failing test families, not independent root causes.',
          'The explanatory repair register is in [the expanded report](../math-audit-expanded-2026-09-10.md).', '',
          '| Test family | Failed cases | Example |', '| --- | ---: | --- |']
    for name, rows in grouped.items():
        example = rows[0]['message'].split('\n')[0].replace('|','\\|')
        file = rows[0]['test_file'] or 'tests/UMapx.Tests/'+rows[0]['suite']+'.cs'
        md.append(f'| [{name}](../../{file}) | {len(rows)} | {example} |')
    (args.output / 'failures.md').write_text('\n'.join(md)+'\n',encoding='utf-8')

    areas = collections.defaultdict(lambda:dict(files=0, instrumented_files=0, executed_files=0, covered_lines=0, executable_lines=0))
    for row in inventory:
        area=areas[row['area']];area['files']+=1;area['instrumented_files']+=row['executable_lines']>0
        area['executed_files']+=row['covered_lines']>0;area['covered_lines']+=row['covered_lines'];area['executable_lines']+=row['executable_lines']
    total = collections.Counter()
    for counts in suites.values(): total.update(counts)
    summary = dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                   source_digest=hashlib.sha256('\n'.join(row['file']+':'+row['sha256'] for row in inventory).encode()).hexdigest(),
                   source_files=len(files), totals=dict(total), suites=dict(sorted(suites.items())), areas=dict(sorted(areas.items())),
                   coverage={k:coverage.attrib[k] for k in ['lines-covered','lines-valid','branches-covered','branches-valid','line-rate','branch-rate']},
                   uncovered_methods=len(methods), failing_test_families=len(grouped),
                   evidence_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.trx,args.coverage,args.coverage.with_name('coverage.json')]})
    dump(args.output / 'summary.json', summary)
    md=['# Source coverage inventory', '',
        'Each repository source file is listed. Execution coverage is not proof of mathematical correctness.',
        'Missing methods and exact uncovered lines are available in the adjacent JSON and CSV files.', '',
        '| Source file | Covered / executable lines | Line coverage | Scope |', '| --- | ---: | ---: | --- |']
    for row in inventory:
        percent='n/a' if row['line_percent'] is None else str(row['line_percent'])+'%'
        md.append(f"| [{row['file']}](../../{row['file']}) | {row['covered_lines']} / {row['executable_lines']} | {percent} | {row['scope']} |")
    (args.output / 'source-inventory.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ['source_files','totals','coverage','uncovered_methods','failing_test_families']},indent=2))


if __name__ == '__main__':
    main()
