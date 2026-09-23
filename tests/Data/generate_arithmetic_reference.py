"""Generate independent 100-digit mpmath references at exact binary32 inputs."""
import json
from pathlib import Path
import struct
import mpmath as mp

mp.mp.dps = 100

def single(value):
    return struct.unpack('f', struct.pack('f', value))[0]

rows = []
real = [-1e30, -10000, -2, -1, -.5, -1e-20, 0, 1e-20, .5, 1, 2, 10000, 1e30]
imag = [-1e30, -20, -1.0000001, -1, -.99999994, -.5, -1e-20, -1e-40, 0, 1e-40, 1e-20, .5, .99999994, 1, 1.0000001, 20, 1e30]
for name, function in [('Asinh', mp.asinh), ('Acosh', mp.acosh), ('Actan', lambda z: mp.atan(1/z))]:
    for re in real:
        for im in imag:
            re, im = single(re), single(im)
            if name == 'Actan' and (re, im) in [(0, 0), (0, 1), (0, -1)]:
                continue
            value = function(mp.mpc(re, im))
            rows.append([name, re, im, 0, 0, float(value.real), float(value.imag)])
for base in [-1e30, -10000, -2, -1, -1e-20, 1e-20, 1, 2, 1e30]:
    for re, im in [(0, 0), (1, 0), (.5, 0), (-.5, 0), (.5, .25), (-.5, -.25), (0, 1)]:
        base, re, im = single(base), single(re), single(im)
        value = mp.power(mp.mpc(base, 0), mp.mpc(re, im))
        rows.append(['Pow', base, 0, re, im, float(value.real), float(value.imag)])
for name, function in [('Tanh', mp.tanh), ('Ctanh', mp.coth), ('Sech', mp.sech), ('Cosch', mp.csch)]:
    for re in [-1e30, -100, -89, -20, -1e-30, 0, 1e-30, 20, 89, 100, 1e30]:
        for im in [-20, -1.5707963267948966, -.5, 0, .5, 1.5707963267948966, 20]:
            re, im = single(re), single(im)
            if name in ['Ctanh', 'Cosch'] and re == 0 and im == 0:
                continue
            value = function(mp.mpc(re, im))
            rows.append([name, re, im, 0, 0, float(value.real), float(value.imag)])
output = Path(__file__).with_name('arithmetic-repair.json')
output.write_text(json.dumps({'generator': 'mpmath ' + mp.__version__, 'precision': mp.mp.dps, 'cases': rows}, indent=2) + '\n', encoding='utf-8')
print(f'{len(rows)} arithmetic reference cases -> {output}')
