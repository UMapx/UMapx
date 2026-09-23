"""Generate reference values with mpmath 1.3.0, 60 decimal digits.
All non-integer inputs are rounded to IEEE float32 before reference evaluation.
"""
import json, struct
from pathlib import Path
root=Path(__file__).resolve().parent

import mpmath as mp
mp.mp.dps=60
cases=[]
def f32(x): return struct.unpack('f',struct.pack('f',float(x)))[0]
def add(name,fn,args,kinds):
    values=[]; encoded=[]
    for x,kind in zip(args,kinds):
        if kind=='c':
            z=complex(x); re,im=f32(z.real),f32(z.imag); values.append(mp.mpc(re,im));encoded.append([re,im])
        elif kind=='i': values.append(int(x));encoded.append([int(x),0])
        else: val=f32(x);values.append(mp.mpf(val));encoded.append([val,0])
    try: expected=fn(*values)
    except (ValueError,ZeroDivisionError): return
    expected=mp.mpc(expected)
    if not mp.isfinite(expected) or abs(expected)>mp.mpf('1e35'): return
    if 'c' not in kinds and abs(expected.imag)>mp.mpf('1e-40'): return
    cases.append(dict(name=name,kinds=kinds,args=encoded,re=float(expected.real),im=float(expected.imag)))
functions={
 'Gamma':mp.gamma,'LogGamma':mp.loggamma,'DiGamma':lambda x:mp.digamma(x),
 'TriGamma':lambda x:mp.polygamma(1,x),'Erf':mp.erf,'Erfc':mp.erfc,'Erfi':mp.erfi,
 # UMapx uses integrands cos(t^2), sin(t^2), without the pi/2 scaling.
 'Fresnelc':lambda x:mp.sqrt(mp.pi/2)*mp.fresnelc(x*mp.sqrt(2/mp.pi)),
 'Fresnels':lambda x:mp.sqrt(mp.pi/2)*mp.fresnels(x*mp.sqrt(2/mp.pi)),'Zeta':mp.zeta,
 'Ci':mp.ci,'Si':mp.si,'Ei':mp.ei,'Li':mp.li,
}
for name,fn in functions.items():
    for x in [-5.25,-1.25,-.3,.1,.5,1,2,3,5,8,10,20,30]: add(name,fn,[x],['f'])
    for z in [.5+.25j,2+1j,-2+.5j,3-2j,1+5j,5+1j]: add(name,fn,[z],['c'])
    if name in ['Fresnelc','Fresnels']:
        for x in [5.9,6,6.1]:add(name,fn,[x],['f'])
for name,fn in [('J',mp.besselj),('Y',mp.bessely),('I',mp.besseli),('K',mp.besselk),('H',mp.struveh),('L',mp.struvel)]:
    for n in [0,1,2,5,10]:
        for x in [.1,1,3,10,19,20,21,30]: add(name,lambda x,k:fn(k,x),[x,n],['f','i'])
        for z in [1+.5j,3+2j,1+5j,20+5j,30j]: add(name,lambda x,k:fn(k,x),[z,n],['c','i'])
for x in [-.3,-.1,.1,1,3,10]:
    for k in [-1,0,1]:
        add('LambertW',mp.lambertw,[x,k],['f','i']);add('LambertW',mp.lambertw,[complex(x,.5),k],['c','i'])
for a,b in [(.5,.5),(2,3),(5,10),(20,30)]:
    add('Beta',mp.beta,[a,b],['f','f'])
    for x in [.01,.1,.5,.9,.99]:
        add('BetaIncomplete',lambda a,b,x:mp.betainc(a,b,0,x),[a,b,x],['f','f','f'])
        add('BetaIncompleteRegularized',lambda a,b,x:mp.betainc(a,b,0,x,regularized=True),[a,b,x],['f','f','f'])
for a,b in [(1,2),(2,3),(.5,1),(-2,1)]:
    for z in [-2,.5,2,5]:add('Hypergeom',mp.hyp1f1,[a,b,z],['f','f','f'])
for a,b,c in [(1,1,2),(.5,.5,1),(2,3,4),(-2,3,1)]:
    for z in [-2,-.5,.5,.9]: add('Hypergeom',mp.hyp2f1,[a,b,c,z],['f','f','f','f'])
(root/'special-functions.json').write_text('[' + chr(10) + (',' + chr(10)).join(json.dumps(case,separators=(',', ':')) for case in cases) + chr(10) + ']' + chr(10),encoding='utf-8')
print(f'{len(cases)} reference values generated with mpmath {mp.__version__}; dps={mp.mp.dps}')

