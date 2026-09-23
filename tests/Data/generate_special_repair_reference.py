"""Additional independent reference cases for the special-function repairs.

Run with mpmath 1.3.0. Inputs are rounded to float32 before the 80-digit
reference is evaluated. Existing audit fixtures are intentionally unchanged.
"""
import json
import struct
from pathlib import Path
import mpmath as mp

mp.mp.dps = 80
cases = []

def add(name, function, args, kinds):
    encoded, values = [], []
    for x, kind in zip(args, kinds):
        if kind in 'ib':
            values.append(int(x)); encoded.append([int(x), 0])
        else:
            z = complex(x)
            re, im = [struct.unpack('f', struct.pack('f', v))[0] for v in (z.real, z.imag)]
            values.append(mp.mpc(re, im) if kind == 'c' else mp.mpf(re))
            encoded.append([re, im])
    try:
        expected = mp.mpc(function(*values))
    except (ValueError, ZeroDivisionError):
        return
    maximum = mp.mpf('1e300') if name in ['Factorial', 'FactorialDown', 'Binomial', 'Euler', 'Bernoulli'] and 'c' not in kinds else mp.mpf('3e38')
    if not mp.isfinite(expected) or abs(expected) > maximum:
        return
    if 'c' not in kinds and name != 'Faddeeva' and abs(expected.imag) > mp.mpf('1e-50'):
        return
    cases.append(dict(name=name, kinds=list(kinds), args=encoded,
                      re=float(expected.real), im=float(expected.imag)))

for name, fn in [('J',mp.besselj),('Y',mp.bessely),('I',mp.besseli),('K',mp.besselk),('H',mp.struveh),('L',mp.struvel)]:
    print(name, flush=True)
    for n in [-5,-2,-1,0,1,2,5,10,20,50]:
        for x in [.01,.5,8.99,9,9.01,11.99,12,12.01,15,35,60,100]:
            add(name,lambda x,n:fn(n,x),[x,n],'fi')
        for z in [.25+12j,-8+1j,-12-2j,10+10j,30j,60+.5j]:
            add(name,lambda z,n:fn(n,z),[z,n],'ci')

for name,fn in [('Gamma',mp.gamma),('LogGamma',mp.loggamma),('DiGamma',mp.digamma),('TriGamma',lambda z:mp.polygamma(1,z)),('Zeta',mp.zeta)]:
    for z in [-50.2,-20.001,-8.999,.00001,.99999,1.00001,40,100]:
        add(name,(lambda x:mp.log(abs(mp.gamma(x)))) if name=='LogGamma' else fn,[z],'f')
    for z in [.001+.001j,-.5+30j,-130+.5j,40+3j,-3.2-2j,.5+100j,1+9.06472028365j]:
        add(name,fn,[z],'c')

def lower_gamma(a,x):
    # mpmath 1.3.0 can recurse indefinitely while ordering complex endpoints
    # with negative real part. Use its independent 1F1 implementation there.
    return x**a/a*mp.hyp1f1(a,a+1,-x) if mp.re(x)<0 else mp.gammainc(a,0,x)

def upper_gamma(a,x):
    return mp.gamma(a)-lower_gamma(a,x) if mp.re(x)<0 else mp.gammainc(a,x,mp.inf)

for name,fn in [('GammaP',lambda a,x:lower_gamma(a,x)/mp.gamma(a)),('GammaQ',lambda a,x:upper_gamma(a,x)/mp.gamma(a)),('GammaIncomplete',lower_gamma),('GammaIncompleteComplemented',upper_gamma)]:
    print(name, flush=True)
    for a in [-2.3,-.5,.01,.3,5,30,100,500]:
        for x in [.001,.5,3,20,100,500,1000]:add(name,fn,[a,x],'ff')
    for a in [.3+.7j,-.5+.3j,3+12j]:
        for x in [.01+.02j,5+3j,-4+.1j,20-5j]:add(name,fn,[a,x],'cc')

for a,b in [(.01,.2),(.3,120),(30,120),(500,500)]:
    add('Beta',mp.beta,[a,b],'ff')
    add('LogBeta',lambda a,b:mp.log(mp.beta(a,b)),[a,b],'ff')
    add('BetaDerivative',lambda a,b:mp.diff(lambda t:mp.beta(t,b),a),[a,b],'ff')
    for x in [1e-10,.01,.3,.5,.9,.99999]:
        add('BetaIncomplete',lambda a,b,x:mp.betainc(a,b,0,x),[a,b,x],'fff')
        add('BetaIncompleteRegularized',lambda a,b,x:mp.betainc(a,b,0,x,regularized=True),[a,b,x],'fff')
for a,b in [(.3+.2j,2-.1j),(20+3j,30-2j)]:
    for name,fn in [('Beta',mp.beta),('LogBeta',lambda a,b:mp.loggamma(a)+mp.loggamma(b)-mp.loggamma(a+b))]:add(name,fn,[a,b],'cc')
    for x in [.1+.2j,.9+.1j,-2+.5j,2+.1j]:add('BetaIncomplete',lambda a,b,x:mp.betainc(a,b,0,x),[a,b,x],'ccc')

for a,b,c in [(1,1,2),(.5,.5,1),(2,3,4),(-5,3,2),(.3,2.5,4),(2.1,3.2,1.3)]:
    for z in [-50,-2,.999,.99999]:add('Hypergeom',mp.hyp2f1,[a,b,c,z],'ffff')
    for z in [.8+.9j,2+.1j,2-.1j,-3+2j]:add('Hypergeom',mp.hyp2f1,[a,b,c,z],'cccc')
for a,b in [(1,2),(2.3,.7),(-5,2),(.5,1)]:
    for z in [-100,-30,30,60]:add('Hypergeom',mp.hyp1f1,[a,b,z],'fff')
    for z in [30j,-30+20j,60j]:add('Hypergeom',mp.hyp1f1,[a,b,z],'ccc')

for name,fn in [('Erf',mp.erf),('Erfc',mp.erfc),('Erfi',mp.erfi),('Faddeeva',lambda z:mp.exp(-z*z)*mp.erfc(-1j*z)),('Fresnelc',lambda z:mp.sqrt(mp.pi/2)*mp.fresnelc(z*mp.sqrt(2/mp.pi))),('Fresnels',lambda z:mp.sqrt(mp.pi/2)*mp.fresnels(z*mp.sqrt(2/mp.pi))),('Ci',mp.ci),('Si',mp.si),('Ei',mp.ei),('Li',mp.li)]:
    for z in [-100,-30,-10,-6,-1e-10,1e-10,5.99,6,6.01,8,10,11.99,12,12.01,30,100]:add(name,fn,[z],'f')
    for z in [-3+2j,3-2j,6+6j,12+.01j,.01+12j,-15-2j,30j,30+20j]:add(name,fn,[z],'c')
for z in [-50,-10,-.001,.001,10,50,3+2j,-3-2j,12+.1j,1+12j]:
    for sign in [0,1]:add('Dawson',lambda z,p:mp.sqrt(mp.pi)/2*(mp.exp(-z*z)*mp.erfi(z) if p else mp.exp(z*z)*mp.erf(z)),[z,sign],'cb' if isinstance(z,complex) else 'fb')
for p in [1e-38,1e-20,1e-8,.25,.99999994]:add('Q',lambda p,_:mp.sqrt(2)*mp.erfinv(1-2*p),[p,1],'fb')
for p in [-.99999994,-.5,1e-10,.99999994]:add('Erf',lambda p,_:mp.erfinv(p),[p,1],'fb')

for z in [-.36787942,-.36,-.001,1e-30,1001,1e20,-.4+.01j,-.4-.01j,-3+.01j,-3-.01j,.01+1j]:
    for k in [-3,-1,0,1,3]:add('LambertW',mp.lambertw,[z,k],'ci' if isinstance(z,complex) else 'fi')
for z in [-1,.1,.5,1,2,.2+.3j,1+1j,-1+.2j]:
    for n in [0,1,2,3,4,5]:add('Gerf',lambda z,n:mp.factorial(n)/mp.sqrt(mp.pi)*mp.quad(lambda t:mp.exp(-t**n),[0,z]),[z,n],'ci' if isinstance(z,complex) else 'fi')

for n in [0,.3,35,100,170]:
    add('Factorial',mp.factorial,[n],'f')
    for k in [0,2,10,40]:
        for name,fn in [('FactorialUp',mp.rf),('FactorialDown',mp.ff),('Binomial',mp.binomial)]:add(name,fn,[n,k],'ff')
for y,n in [(100,100),(1000,1000),(80,200)]:add('Erlang',lambda y,n:y**n/mp.factorial(n)/mp.fsum(y**k/mp.factorial(k) for k in range(n+1)),[y,n],'fi')

# Recurrences at higher degrees also guard against exponentially recursive implementations.
for name,fn in [('ChebyshevT',mp.chebyt),('ChebyshevU',mp.chebyu),('Legendre',mp.legendre),('Hermite',mp.hermite)]:
    for n in [20,50]:
        for z in [-1,0,.3,1,1.1,.3+.2j]:add(name,lambda z,n:fn(n,z),[z,n],'ci' if isinstance(z,complex) else 'fi')
for name,fn in [('Laguerre',lambda z,a,n:mp.laguerre(n,a,z)),('Gegenbauer',lambda z,a,n:mp.gegenbauer(n,a,z))]:
    for z in [.3,.3+.2j]:add(name,fn,[z,.7,50],'cci' if isinstance(z,complex) else 'ffi')

path = Path(__file__).with_name('special-functions-repair.json')
for n in [0,1,2,5,10,20]:
    add('Euler',lambda n:mp.eulernum(n),[n],'i')
    add('Bernoulli',mp.bernoulli,[n],'i')
    for x in [0,.3,.5,1]:
        add('Euler',mp.eulerpoly,[n,x],'if')
        add('Bernoulli',mp.bernpoly,[n,x],'if')
for n in [-40,-10,-1,0,1,10,30,44]:
    add('Fibonacci',mp.fibonacci,[n],'i')
    add('Lucas',lambda n:mp.fibonacci(n-1)+mp.fibonacci(n+1),[n],'i')
for h,a in [(.1,100),(10,1),(1,-2),(1+.3j,.2+.1j),(-2+.1j,.5+.2j)]:
    add('Owen',lambda h,a:mp.quad(lambda t:mp.exp(-h*h*(1+t*t)/2)/(1+t*t),[0,a])/(2*mp.pi),[h,a],'cc' if isinstance(h,complex) else 'ff')
for x in [-100,100,-2,.3,.3+.2j]:
    kind='c' if isinstance(x,complex) else 'f'
    for name,fn in [('Gd',lambda z:2*mp.atan(mp.tanh(z/2))),('Logistic',lambda z:1/(1+mp.exp(-z)))]:add(name,fn,[x],kind)
    add('Heaviside',lambda x,k:(1+mp.tanh(x*k))/2,[x,1],kind*2)
    add('Gompertz',lambda t,a,b,c:a*mp.exp(-b*mp.exp(-c*t)),[x,2,.5,.3],kind*4)
    add('Dirac',lambda x,a:mp.exp(-(x/a)**2)/(abs(a)*mp.sqrt(mp.pi)),[x,.3],kind*2)
    add('Mahler',lambda x,t:mp.exp(x*(1+t-mp.exp(t))),[x,.1],kind*2)
for x in [1,.8,2,100]:add('Ssqrt',lambda x,k:mp.exp(mp.lambertw(mp.log(x),k)),[x,0],'fi')
for n in [0,1,10,1000]:
    add('Harm',mp.harmonic,[n],'i')
    add('Harm',lambda n,m:mp.fsum(mp.mpf(k)**(-m) for k in range(1,n+1)),[n,.7],'if')
path.write_text('[\n'+',\n'.join(json.dumps(c,separators=(',',':')) for c in cases)+'\n]\n',encoding='utf-8')
print(f'{len(cases)} additional cases; mpmath {mp.__version__}; dps={mp.mp.dps}')
