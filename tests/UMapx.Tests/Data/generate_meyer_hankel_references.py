"""Generate independent Meyer and Hankel fixtures with mpmath 1.3.0 at 70 decimal digits.

Run from the repository root with mpmath on PYTHONPATH. Production code is not
loaded: Meyer is integrated in the frequency domain, and Hankel uses mpmath's
ordered positive Bessel zeros. Output arguments are rounded to float32 first.
"""
import json
import struct
from pathlib import Path
import mpmath as mp

mp.mp.dps = 70

def single(x):
    """Round a finite real argument to the library's IEEE float32 input format."""
    return struct.unpack('f', struct.pack('f', x))[0]

def adjacent(x, direction):
    """Return an adjacent positive float32 value in the chosen bit direction."""
    bits = struct.unpack('I', struct.pack('f', x))[0]
    return struct.unpack('f', struct.pack('I', bits + direction))[0]

meyer = []
for wavelet, centers in [(False, [0, .75]), (True, [.125, .5, .875, 1.25])]:
    points = {-3.0, -1.0, 0.0, .1, .37, 2.0, 4.0}
    for c in centers:
        points.add(c)
        if c:
            points.update(adjacent(c, d) for d in [-1, 1])
    if not wavelet:
        points.update(-v for v in list(points))
    for x in sorted({single(v) for v in points}):
        t = mp.mpf(x) - (mp.mpf('.5') if wavelet else 0)
        a = 2 * mp.pi / 3
        if wavelet:
            value = (mp.quad(lambda w: -mp.cos(3*w/4)*mp.cos(w*t), [a, 2*a])
                     + mp.quad(lambda w: mp.sin(3*w/8)*mp.cos(w*t), [2*a, 4*a])) / mp.pi
        else:
            value = (mp.quad(lambda w: mp.cos(w*t), [0, a])
                     + mp.quad(lambda w: mp.sin(3*w/4)*mp.cos(w*t), [a, 2*a])) / mp.pi
        meyer.append(dict(wavelet=wavelet, x=x, expected=float(value)))

hankel = []
for order in [0, 1, 5, 10, 20, 50]:
    roots = [mp.besseljzero(order, k) for k in range(1, 22)]
    for size in [1, 3, 9, 20]:
        denominator = roots[size]
        values = [[float(2*mp.besselj(order, roots[i]*roots[j]/denominator)
                         / (denominator*mp.besselj(order+1, roots[i])*mp.besselj(order+1, roots[j])))
                   for j in range(size)] for i in range(size)]
        hankel.append(dict(order=order, size=size, values=values))

Path('tests/UMapx.Tests/Data/meyer-hankel.json').write_text(
    json.dumps(dict(meyer=meyer, hankel=hankel), indent=2) + '\n', encoding='utf-8')
