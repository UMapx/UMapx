"""Independent probability references; requires mpmath 1.3.0.

Parameters and evaluation points are rounded to binary32 before evaluation.
Moments and entropy are integrated from the reference density, not from UMapx
formulas. PowerNormal and PowerLognormal follow the NIST definitions cited by
their API documentation. Generated files are sufficient to run the C# tests.
"""
import json
import math
import re
import struct
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40
ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parents[2] / 'sources' / 'Distribution'
CASES = []
CATALOG = []
phi = lambda x: mp.exp(-x*x/2) / mp.sqrt(2*mp.pi)
Phi = lambda x: mp.erfc(-x/mp.sqrt(2))/2


def f32(x):
    return struct.unpack('f', struct.pack('f', float(x)))[0]


def emit(name, args, signature, member, expected, x=None, tolerance=3e-4):
    expected = mp.mpf(expected)
    if mp.isnan(expected):
        value = 'NaN'
    elif mp.isinf(expected):
        value = 'Infinity' if expected > 0 else '-Infinity'
    else:
        value = float(expected)
    CASES.append(dict(name=name, constructor=args, signature=signature,
                      member=member, x=x, expected=value, tolerance=tolerance))


def law(name, parameters, factory, signature=None, points=None, statistics=True):
    args = [f32(x) for x in parameters]
    signature = signature or ['f'] * len(args)
    values = [int(x) if kind == 'i' else mp.mpf(x) for x, kind in zip(args, signature)]
    cfg = factory(*values)
    lo, hi = cfg['support']
    discrete = cfg.get('discrete', False)
    raw_pdf = cfg['pdf']
    pdf = lambda x: mp.mpf(0) if not discrete and (x <= lo or x >= hi) else raw_pdf(x)
    cuts = sorted(set([lo, hi] + [mp.mpf(x) for x in cfg.get('cuts', [-8,-2,0,1,2,4,8,32]) if lo < x < hi]))

    def integrate(fn, upper=None):
        if discrete:
            # The selected parameters make the omitted tail smaller than 1e-25.
            low = int(lo) if mp.isfinite(lo) else -500
            high = int(hi) if mp.isfinite(hi) else 1000
            if upper is not None: high = min(high, math.floor(upper))
            return mp.fsum(fn(mp.mpf(k)) for k in range(low, high+1))
        lower, maximum = cfg.get('integration_limits', (lo, hi))
        bound = maximum if upper is None else min(maximum, upper)
        if bound <= lower: return mp.mpf(0)
        split = [lower] + [v for v in cuts if lower < v < bound] + [bound]
        return mp.quad(fn, split)

    cdf = cfg.get('cdf', lambda x: integrate(pdf, x))
    points = points or [-8,-2,-.5,0,.125,.5,1,2,3,5,10,30]
    code = (SOURCE / (name+'.cs')).read_text(encoding='utf-8-sig')
    unsupported = []
    for member in ['Mean','Variance','Median','Mode','Skewness','Excess','Entropy','Distribution']:
        match = re.search(r'public float(?:\[\])? '+member+r'\b[\s\S]*?(?=\n\s*///|\n\s*#endregion)', code)
        if match and ('throw new NotSupportedException' in match[0] or '=> float.NaN' in match[0] or 'return float.NaN;' in match[0] and '{' not in match[0].split('return float.NaN;')[0].split('get')[-1].strip('{} \n\r')):
            unsupported.append(member)
    CATALOG.append(dict(name=name, constructor=args, signature=signature,
                        support=[str(lo),str(hi)], discrete=discrete,
                        unsupported=unsupported))
    emit(name,args,signature,'Support.Min',lo)
    emit(name,args,signature,'Support.Max',hi)
    for point in sorted(set(points)):
        x = mp.mpf(f32(point))
        if ((x < lo or x > hi) and not cfg.get('periodic')) or (discrete and x != mp.floor(x)):
            density = mp.mpf(0)
        elif not discrete and (x == lo or x == hi):
            # Endpoint density values are arbitrary on a set of measure zero.
            density = None
        else:
            density = raw_pdf(x) if cfg.get('periodic') else pdf(x)
        if density is not None and mp.isfinite(density):
            emit(name,args,signature,'Function',density,float(x))
        if cfg.get('cdf_supported', True):
            cumulative = 0 if x < lo or (x == lo and not discrete) else 1 if x >= hi else cdf(x)
            emit(name,args,signature,'Distribution',cumulative,float(x))
    for member, value in cfg.get('properties', {}).items():
        if member not in unsupported:
            emit(name,args,signature,member,value)
    if statistics:
        requested = [m for m in ['Mean','Variance','Skewness','Excess','Entropy']
                     if m not in unsupported and m not in cfg.get('properties', {}) and m not in cfg.get('omit', [])]
        if any(m in requested for m in ['Mean','Variance','Skewness','Excess']):
            mass = integrate(pdf)
            assert abs(mass-1) < mp.mpf('1e-12'), (name, args, 'reference mass', mass)
            mean = integrate(lambda x: x*pdf(x))
            variance = integrate(lambda x: (x-mean)**2*pdf(x))
            stats = dict(Mean=mean,Variance=variance)
            if variance > 0:
                if 'Skewness' in requested: stats['Skewness'] = integrate(lambda x: (x-mean)**3*pdf(x))/variance**mp.mpf('1.5')
                if 'Excess' in requested: stats['Excess'] = integrate(lambda x: (x-mean)**4*pdf(x))/variance**2-3
            for member in requested:
                if member in stats: emit(name,args,signature,member,stats[member],tolerance=2e-3)
        if 'Entropy' in requested:
            entropy = integrate(lambda x: -pdf(x)*mp.log(pdf(x)) if pdf(x)>0 else mp.mpf(0))
            if cfg.get('entropy_bits'): entropy /= mp.log(2)
            emit(name,args,signature,'Entropy',entropy,tolerance=2e-3)
    print(name, args, len(CASES), flush=True)


def gamma_law(scale, shape):
    return dict(support=(0,mp.inf),pdf=lambda x: x**(shape-1)*mp.exp(-x/scale)/(mp.gamma(shape)*scale**shape),
                cdf=lambda x: mp.gammainc(shape,0,x/scale,regularized=True))


def beta_law(a,b):
    return dict(support=(0,1),pdf=lambda x:x**(a-1)*(1-x)**(b-1)/mp.beta(a,b),
                cdf=lambda x:mp.betainc(a,b,0,x,regularized=True),cuts=[mp.mpf('.5')])


def invgamma(a,b):
    return dict(support=(0,mp.inf),pdf=lambda x:b**a/mp.gamma(a)*x**(-a-1)*mp.exp(-b/x),
                cdf=lambda x:mp.gammainc(a,b/x,mp.inf,regularized=True))


def normal(s,m):
    return dict(support=(-mp.inf,mp.inf),pdf=lambda x:phi((x-m)/s)/s,cdf=lambda x:Phi((x-m)/s),cuts=[m-8*s,m,m+8*s],properties={'Median':m})


def poisson(l):
    return dict(support=(0,mp.inf),discrete=True,pdf=lambda k:mp.exp(-l)*l**k/mp.factorial(k),
                cdf=lambda x:mp.gammainc(mp.floor(x)+1,l,mp.inf,regularized=True))


def bs(m,b,g):
    z=lambda x:(mp.sqrt((x-m)/b)-mp.sqrt(b/(x-m)))/g
    return dict(support=(m,mp.inf),pdf=lambda x:phi(z(x))*(mp.sqrt((x-m)/b)+mp.sqrt(b/(x-m)))/(2*g*(x-m)),
                cdf=lambda x:Phi(z(x)),cuts=[m+b/4,m+b,m+4*b],properties={'Median':m+b})


def fisher(d1,d2):
    return dict(support=(0,mp.inf),pdf=lambda x:(d1/d2)**(d1/2)*x**(d1/2-1)/(mp.beta(d1/2,d2/2)*(1+d1*x/d2)**((d1+d2)/2)),
                cdf=lambda x:mp.betainc(d1/2,d2/2,0,d1*x/(d1*x+d2),regularized=True))


def lognormal(s,m):
    return dict(support=(0,mp.inf),pdf=lambda x:phi((mp.log(x)-m)/s)/(s*x),cdf=lambda x:Phi((mp.log(x)-m)/s),properties={'Median':mp.exp(m)})


def trapezoid(a,b,c,d):
    height=2/(d+c-b-a)
    pdf=lambda x:height*(x-a)/(b-a) if x<b else height if x<c else height*(d-x)/(d-c)
    return dict(support=(a,d),pdf=pdf,cuts=[b,c])


def shifted(m,s,k):
    if k==0:
        return dict(support=(-mp.inf,mp.inf),pdf=lambda x:mp.exp(-(x-m)/s)/(s*(1+mp.exp(-(x-m)/s))**2),
                    cdf=lambda x:1/(1+mp.exp(-(x-m)/s)),properties={'Median':m})
    cdf=lambda x:1/(1+(1+k*(x-m)/s)**(-1/k))
    return dict(support=(m-s/k,mp.inf) if k>0 else (-mp.inf,m-s/k),
                pdf=lambda x:(1+k*(x-m)/s)**(-1/k-1)/(s*(1+(1+k*(x-m)/s)**(-1/k))**2),cdf=cdf,
                properties={'Median':m},omit=['Skewness','Excess'])


def tukey(l):
    q=lambda p:mp.log(p/(1-p)) if l==0 else (p**l-(1-p)**l)/l
    def cdf(x):
        a,b=mp.mpf(0),mp.mpf(1)
        for _ in range(160):
            mid=(a+b)/2
            if q(mid)<x:a=mid
            else:b=mid
        return (a+b)/2
    def pdf(x):
        p=cdf(x)
        return p*(1-p) if l==0 else 1/(p**(l-1)+(1-p)**(l-1))
    variance=mp.pi**2/3 if l==0 else mp.inf if l<=mp.mpf('-.5') else 2/l**2*(1/(1+2*l)-mp.beta(l+1,l+1))
    return dict(support=(-1/l,1/l) if l>0 else (-mp.inf,mp.inf),pdf=pdf,cdf=cdf,
                properties={'Mean':0 if l>-1 else mp.nan,'Variance':variance,'Median':0})


law('Arcsine',[],lambda:dict(support=(0,1),pdf=lambda x:1/(mp.pi*mp.sqrt(x*(1-x))),cdf=lambda x:2/mp.pi*mp.asin(mp.sqrt(x)),properties={'Median':mp.mpf('.5')}))
for p in [.25,.7]:
    law('Bernoulli',[p],lambda p:dict(support=(0,1),discrete=True,pdf=lambda k:p if k==1 else 1-p))
for a,b in [(2,3),(.5,.75),(20,30)]:law('Beta',[a,b],beta_law)
law('BetaPrime',[2,8],lambda a,b:dict(support=(0,mp.inf),pdf=lambda x:x**(a-1)/(mp.beta(a,b)*(1+x)**(a+b)),cdf=lambda x:mp.betainc(a,b,0,x/(1+x),regularized=True)))
for n,p in [(5,.3),(2,.7),(10,.8)]:law('Binomial',[n,p],lambda n,p:dict(support=(0,n),discrete=True,pdf=lambda k:mp.binomial(n,k)*p**k*(1-p)**(n-k)),['i','f'])
for args in [(0,1,1),(.5,2,.7)]:law('BirnbaumSaunders',args,bs)
law('Burr',[2,5],lambda c,k:dict(support=(0,mp.inf),pdf=lambda x:c*k*x**(c-1)/(1+x**c)**(k+1),cdf=lambda x:1-(1+x**c)**(-k)))
law('Cauchy',[1.5,.25],lambda g,m:dict(support=(-mp.inf,mp.inf),pdf=lambda x:1/(mp.pi*g*(1+((x-m)/g)**2)),cdf=lambda x:mp.mpf('.5')+mp.atan((x-m)/g)/mp.pi,properties={'Median':m},omit=['Mean','Variance','Skewness','Excess']))
law('ChiSquare',[6],lambda n:gamma_law(mp.mpf(2),mp.mpf(n)/2),['i'])
law('Degenerate',[2],lambda n:dict(support=(n,n),discrete=True,pdf=lambda k:mp.mpf(1)),['i'])
law('Erlang',[3,2],lambda n,r:gamma_law(1/r,n),['i','f'])
law('Exponential',[1.25],lambda r:gamma_law(1/r,mp.mpf(1)))
law('FisherSnedecor',[4,12],lambda a,b:fisher(mp.mpf(a),mp.mpf(b)),['i','i'])
law('FisherZ',[4,12],lambda a,b:dict(support=(-mp.inf,mp.inf),pdf=lambda x:2*mp.exp(2*x)*fisher(a,b)['pdf'](mp.exp(2*x)),cdf=lambda x:fisher(a,b)['cdf'](mp.exp(2*x))))
law('FoldedNormal',[.75,1.25],lambda m,s:dict(support=(0,mp.inf),pdf=lambda x:(phi((x-m)/s)+phi((x+m)/s))/s,cdf=lambda x:Phi((x-m)/s)-Phi((-x-m)/s)))
law('Gamma',[1.5,4],gamma_law)
law('Gaussian',[1.25,.5],normal)
law('GeneralizedNormal',[.5,1.5,2.5],lambda m,a,b:dict(support=(-mp.inf,mp.inf),pdf=lambda x:b/(2*a*mp.gamma(1/b))*mp.exp(-(abs(x-m)/a)**b),cdf=lambda x:mp.mpf('.5')+mp.sign(x-m)*mp.gammainc(1/b,0,(abs(x-m)/a)**b,regularized=True)/2,cuts=[m],properties={'Median':m}))
law('Geometric',[.25],lambda p:dict(support=(0,mp.inf),discrete=True,pdf=lambda k:p*(1-p)**k,cdf=lambda x:1-(1-p)**(mp.floor(x)+1),entropy_bits=True))
# Truncate double-exponential tails before quadrature samples astronomically
# large exponents. Omitted mass is below 1e-60 for these selected parameters.
law('Gompertz',[.75,1.25],lambda e,b:dict(support=(0,mp.inf),pdf=lambda x:b*e*mp.exp(b*x-e*mp.expm1(b*x)),cdf=lambda x:-mp.expm1(-e*mp.expm1(b*x)),integration_limits=(0,20)))
law('Gumbel',[.5,1.5],lambda m,b:dict(support=(-mp.inf,mp.inf),pdf=lambda x:mp.exp(-(x-m)/b-mp.exp(-(x-m)/b))/b,cdf=lambda x:mp.exp(-mp.exp(-(x-m)/b)),cuts=[m],integration_limits=(m-10*b,m+150*b)))
law('HyperbolicSecant',[],lambda:dict(support=(-mp.inf,mp.inf),pdf=lambda x:1/(2*mp.cosh(mp.pi*x/2)),cdf=lambda x:2/mp.pi*mp.atan(mp.exp(mp.pi*x/2)),properties={'Median':0}))
law('Hypergeometric',[30,12,8],lambda n,k,d:dict(support=(max(0,k+d-n),min(k,d)),discrete=True,pdf=lambda x:mp.binomial(k,x)*mp.binomial(n-k,d-x)/mp.binomial(n,d)))
law('InverseChiSquare',[12],lambda n:invgamma(mp.mpf(n)/2,mp.mpf('.5')),['i'])
law('InverseGamma',[6,2],invgamma)
law('InverseGaussian',[1.25,2],lambda m,l:dict(support=(0,mp.inf),pdf=lambda x:mp.sqrt(l/(2*mp.pi*x**3))*mp.exp(-l*(x-m)**2/(2*m*m*x)),cdf=lambda x:Phi(mp.sqrt(l/x)*(x/m-1))+mp.exp(2*l/m)*Phi(-mp.sqrt(l/x)*(x/m+1))))
law('Kumaraswamy',[2,3],lambda a,b:dict(support=(0,1),pdf=lambda x:a*b*x**(a-1)*(1-x**a)**(b-1),cdf=lambda x:1-(1-x**a)**b))
law('Laplace',[1.25,.5],lambda a,m:dict(support=(-mp.inf,mp.inf),pdf=lambda x:a/2*mp.exp(-a*abs(x-m)),cdf=lambda x:mp.exp(a*(x-m))/2 if x<=m else 1-mp.exp(-a*(x-m))/2,cuts=[m],properties={'Median':m}))
law('Levy',[.5,1.5],lambda m,c:dict(support=(m,mp.inf),pdf=lambda x:mp.sqrt(c/(2*mp.pi))*mp.exp(-c/(2*(x-m)))/(x-m)**mp.mpf('1.5'),cdf=lambda x:mp.erfc(mp.sqrt(c/(2*(x-m)))),omit=['Mean','Variance','Skewness','Excess']))
law('Logarithmic',[.6],lambda p:dict(support=(1,mp.inf),discrete=True,pdf=lambda k:-p**k/(k*mp.log(1-p))))
law('LogGaussian',[.6,.2],lognormal)
law('Logistic',[.5,1.5],lambda m,s:shifted(m,s,mp.mpf(0)))
law('LogLogistic',[1.5,6],lambda a,b:dict(support=(0,mp.inf),pdf=lambda x:b/a*(x/a)**(b-1)/(1+(x/a)**b)**2,cdf=lambda x:1/(1+(x/a)**(-b)),properties={'Median':a}))
law('Nakagami',[2,1.5],lambda m,w:dict(support=(0,mp.inf),pdf=lambda x:2*m**m/(mp.gamma(m)*w**m)*x**(2*m-1)*mp.exp(-m*x*x/w),cdf=lambda x:mp.gammainc(m,0,m*x*x/w,regularized=True)))
law('NegativeBinomial',[4,.6],lambda r,p:dict(support=(0,mp.inf),discrete=True,pdf=lambda k:mp.binomial(k+r-1,k)*p**r*(1-p)**k,cdf=lambda x:mp.betainc(r,mp.floor(x)+1,0,p,regularized=True)),['i','f'])
law('Pareto',[1.25,6],lambda m,k:dict(support=(m,mp.inf),pdf=lambda x:k*m**k/x**(k+1),cdf=lambda x:1-(m/x)**k))
for l in [3.5,100]:law('Poisson',[l],poisson)
law('PowerNormal',[2],lambda p:dict(support=(-mp.inf,mp.inf),pdf=lambda x:p*phi(x)*Phi(-x)**(p-1),cdf=lambda x:1-Phi(-x)**p),statistics=False)
law('PowerLognormal',[2,1],lambda p,s:dict(support=(0,mp.inf),pdf=lambda x:p*phi(mp.log(x)/s)*Phi(-mp.log(x)/s)**(p-1)/(s*x),cdf=lambda x:1-Phi(-mp.log(x)/s)**p),statistics=False)
law('Rademacher',[],lambda:dict(support=(-1,1),discrete=True,pdf=lambda k:mp.mpf('.5') if abs(k)==1 else mp.mpf(0)))
law('Rayleigh',[1.5],lambda s:dict(support=(0,mp.inf),pdf=lambda x:x/s**2*mp.exp(-x*x/(2*s*s)),cdf=lambda x:1-mp.exp(-x*x/(2*s*s))))
for k in [0,.2,-.2]:law('ShiftedLogLogistic',[.5,1.5,k],shifted)
law('Student',[8],lambda n:dict(support=(-mp.inf,mp.inf),pdf=lambda x:mp.gamma((n+1)/2)/(mp.sqrt(n*mp.pi)*mp.gamma(n/2))*(1+x*x/n)**(-(n+1)/2),cdf=lambda x:mp.mpf('.5')+mp.sign(x)*(1-mp.betainc(n/2,mp.mpf('.5'),0,n/(n+x*x),regularized=True))/2,properties={'Median':0}))
law('SymmetricGeometric',[.4],lambda p:dict(support=(-mp.inf,mp.inf),discrete=True,pdf=lambda k:p*((1-p)/(1+p))**abs(k)))
for args in [(0,1,2,4),(1,2,4,8)]:law('Trapezoidal',args,trapezoid)
law('Triangular',[0,4,1],lambda a,b,c:dict(support=(a,b),pdf=lambda x:2*(x-a)/((b-a)*(c-a)) if x<c else 2*(b-x)/((b-a)*(b-c)),cdf=lambda x:(x-a)**2/((b-a)*(c-a)) if x<c else 1-(b-x)**2/((b-a)*(b-c)),cuts=[c]))
for l in [0,.25,1,2,-.25,-.75,-2]:law('TukeyLambda',[l],tukey,statistics=False)
law('Uniform',[-1,3],lambda a,b:dict(support=(a,b),pdf=lambda x:1/(b-a),cdf=lambda x:(x-a)/(b-a)))
law('UniformDiscrete',[-2,3],lambda a,b:dict(support=(a,b),discrete=True,pdf=lambda k:mp.mpf(1)/(b-a+1)),['i','i'])
law('UQuadratic',[1,5],lambda a,b:dict(support=(a,b),pdf=lambda x:12*(x-(a+b)/2)**2/(b-a)**3,cdf=lambda x:mp.mpf('.5')+4*(x-(a+b)/2)**3/(b-a)**3,cuts=[(a+b)/2]))
law('Weibull',[1.5,2],lambda l,k:dict(support=(0,mp.inf),pdf=lambda x:k/l*(x/l)**(k-1)*mp.exp(-(x/l)**k),cdf=lambda x:1-mp.exp(-(x/l)**k)))
law('Wigner',[2],lambda r:dict(support=(-r,r),pdf=lambda x:2*mp.sqrt(r*r-x*x)/(mp.pi*r*r),cdf=lambda x:mp.mpf('.5')+(x*mp.sqrt(r*r-x*x)+r*r*mp.asin(x/r))/(mp.pi*r*r)))
law('WrappedCauchy',[.5,1.2],lambda m,g:dict(support=(-mp.pi,mp.pi),pdf=lambda x:mp.sinh(g)/(2*mp.pi*(mp.cosh(g)-mp.cos(x-m))),cdf_supported=False,periodic=True,properties={'Mean':m,'Variance':1-mp.exp(-g)},omit=['Skewness','Excess']))

(ROOT/'distributions.json').write_text('[\n'+',\n'.join(json.dumps(c,separators=(',',':')) for c in CASES)+'\n]\n',encoding='utf-8')
(ROOT/'distribution-catalog.json').write_text(json.dumps(CATALOG,indent=2)+'\n',encoding='utf-8')
print('Generated',len(CASES),'cases for',len(set(c['name'] for c in CASES)),'probability distributions.')
