"""Independent mpmath 1.3.0 references, evaluated at exact float32 inputs."""
import json
import struct
from pathlib import Path
import mpmath as mp

mp.mp.dps = 60
root = Path(__file__).resolve().parent
cases = []

def f32(value):
    return struct.unpack('f', struct.pack('f', float(value)))[0]

def add(name, function, arguments, kinds):
    values, encoded = [], []
    for value, kind in zip(arguments, kinds):
        if kind == 'c':
            z = complex(value)
            re, im = f32(z.real), f32(z.imag)
            values.append(mp.mpc(re, im)); encoded.append([re, im])
        elif kind in ('i', 'b'):
            values.append(int(value)); encoded.append([int(value), 0])
        else:
            x = f32(value); values.append(mp.mpf(x)); encoded.append([x, 0])
    try:
        expected = mp.mpc(function(*values))
    except (ValueError, ZeroDivisionError):
        return
    if not mp.isfinite(expected) or abs(expected) > mp.mpf('1e35'):
        return
    # Faddeeva is complex even for a real argument.
    if 'c' not in kinds and name != 'Faddeeva' and abs(expected.imag) > mp.mpf('1e-40'):
        return
    cases.append(dict(name=name, kinds=list(kinds), args=encoded, re=float(expected.real), im=float(expected.imag)))

for kind, xs in [('f', [-2, -1, -.5, 0, .5, 1, 2]), ('c', [.5+.25j, -1+.5j, 2-1j])]:
    for name, fn in [('ChebyshevT',mp.chebyt), ('ChebyshevU',mp.chebyu), ('Legendre',mp.legendre), ('Hermite',mp.hermite)]:
        for n in [0,1,2,3,5,8,12]:
            for x in xs: add(name, lambda x,n:fn(n,x), [x,n], kind+'i')
    for name, fn in [('Abel',lambda x,a,n:1 if n==0 else x*(x-n*a)**(n-1)), ('Laguerre',lambda x,a,n:mp.laguerre(n,a,x)), ('Gegenbauer',lambda x,a,n:mp.gegenbauer(n,a,x))]:
        for n in [0,1,2,5,8]:
            for x in xs: add(name,fn,[x,.7,n],kind*2+'i')
    for name, fn in [('Gd',lambda x:2*mp.atan(mp.tanh(x/2))), ('Agd',lambda x:mp.atanh(mp.sin(x))), ('Cas',lambda x:mp.cos(x)+mp.sin(x)), ('Sinc',lambda x:mp.sinc(mp.pi*x)), ('Logistic',lambda x:1/(1+mp.exp(-x))), ('Faddeeva',lambda x:mp.exp(-x*x)*mp.erfc(-1j*x))]:
        for x in xs:add(name,fn,[x],kind)
    for x in xs:
        for positive in [False,True]:
            add('Dawson',lambda x,p:mp.sqrt(mp.pi)/2*(mp.exp(-x*x)*mp.erfi(x) if p else mp.exp(x*x)*mp.erf(x)),[x,positive],kind+'b')
        for n in [0,1,2,5,20]:
            if (complex(x).real>=0):add('Erlang',lambda y,n:(y**n/mp.factorial(n))/sum(y**k/mp.factorial(k) for k in range(n+1)),[x,n],kind+'i')
        for a in [.3,2]:add('Sinc',lambda x,a:mp.sinc(a*x),[x,a],kind*2)
    for a in [.5,2,5,10]:
        for x in ([.1,1,3,5,10,20] if kind=='f' else [.5+.25j,3+1j]):
            for name, fn in [('GammaP',lambda a,x:mp.gammainc(a,0,x)/mp.gamma(a)),('GammaQ',lambda a,x:mp.gammainc(a,x,mp.inf)/mp.gamma(a)),('GammaIncomplete',lambda a,x:mp.gammainc(a,0,x)),('GammaIncompleteComplemented',lambda a,x:mp.gammainc(a,x,mp.inf))]:
                add(name,fn,[a,x],kind*2)
    for x in ([.1,.5,1] if kind=='f' else [.1+.1j,.5+.25j]):
        for n in [0,1,2,3,4,5]:add('Gerf',lambda x,n:mp.factorial(n)/mp.sqrt(mp.pi)*mp.quad(lambda t:mp.exp(-t**n),[0,x]),[x,n],kind+'i')
    for h in [.1,1,3]:
        for a in [.1,1,3]:add('Owen',lambda h,a:mp.quad(lambda t:mp.exp(-h*h*(1+t*t)/2)/(1+t*t),[0,a])/(2*mp.pi),[h,a],kind*2)

for n in [0,1,2,5,10,30]:
    add('Factorial',mp.factorial,[n],'f');add('LogFactorial',lambda n:mp.loggamma(n+1),[n],'f')
    for k in range(0,7):
        add('FactorialDown',mp.ff,[n,k],'ff');add('FactorialUp',mp.rf,[n,k],'ff');add('Binomial',mp.binomial,[n,k],'ff')
        if k<=n:add('LogBinomial',lambda n,k:mp.log(mp.binomial(n,k)),[n,k],'ff')
for p in [.001,.1,.25,.5,.9,.999]:
    add('Erf',lambda p,b:mp.erfinv(p),[p,True],'fb')
    add('Q',lambda p,b:mp.sqrt(2)*mp.erfinv(1-2*p),[p,True],'fb')

(root/'special-functions-extended.json').write_text('[\n'+',\n'.join(json.dumps(c,separators=(',',':')) for c in cases)+'\n]\n',encoding='utf-8')

hankel=[]
for order in [0,1,2,5]:
    for size in [2,4,8,16]:
        zeros=[mp.besseljzero(order,k) for k in range(1,size+2)]
        z=zeros[-1]
        values=[[float(2*mp.besselj(order,zeros[i]*zeros[j]/z)/(z*mp.besselj(order+1,zeros[i])*mp.besselj(order+1,zeros[j]))) for j in range(size)] for i in range(size)]
        hankel.append(dict(order=order,size=size,values=values))
(root/'hankel.json').write_text(json.dumps(hankel,separators=(',',':'))+'\n',encoding='utf-8')
print(f'{len(cases)} special-function cases and {len(hankel)} Hankel matrices; mpmath {mp.__version__}; dps={mp.mp.dps}')
