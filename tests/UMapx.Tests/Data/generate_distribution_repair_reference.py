"""Independent probability references at exact binary32 inputs."""
import json
import math
from pathlib import Path
import struct
import mpmath as mp

mp.mp.dps = 70
rows = []
def single(x): return struct.unpack('f', struct.pack('f', float(x)))[0]
def emit(name, args, member, expected, x=None, signature=None):
    value = 'NaN' if mp.isnan(expected) else 'Infinity' if expected == mp.inf else '-Infinity' if expected == -mp.inf else float(expected)
    rows.append(dict(name=name, constructor=[single(a) for a in args], signature=signature or ['f']*len(args), member=member,
                     x=single(x) if x is not None else None, expected=value))
def args32(*args): return tuple(mp.mpf(single(x)) for x in args)
phi=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
tail=lambda x:mp.erfc(x/mp.sqrt(2))/2

for power in [.0001,.05,.5,1,2,10,100,1e-30]:
    p,=args32(power)
    for xx in [-12,-8,-2,0,.5,2,8,12,26,50,100,1e15]:
        x=mp.mpf(single(xx)); logTail=mp.log(tail(x))
        emit('PowerNormal',[p],'Function',p*phi(x)*mp.exp((p-1)*logTail),x)
        emit('PowerNormal',[p],'Distribution',-mp.expm1(p*logTail),x)
for power in [.05,.5,1,2,10]:
    for sigma in [.1,.7,3]:
        p,s=args32(power,sigma)
        for xx in [1e-30,1e-5,.1,.5,1,2,10,1e5,1e30]:
            x=mp.mpf(single(xx)); z=mp.log(x)/s; logTail=mp.log(tail(z))
            emit('PowerLognormal',[p,s],'Function',p*phi(z)*mp.exp((p-1)*logTail)/(s*x),x)
            emit('PowerLognormal',[p,s],'Distribution',-mp.expm1(p*logTail),x)
print('Power laws:',len(rows),flush=True)

for shape in [.0001,.01,.1,.7,1,2,10,100]:
    g,=args32(shape)
    # Integrate -f log f directly in log(x/beta), including the Jacobian.
    limit=2*mp.asinh(14*g/2)
    def entropy_integrand(t):
        u=mp.exp(t); z=(mp.sqrt(u)-1/mp.sqrt(u))/g
        density=phi(z)*(mp.sqrt(u)+1/mp.sqrt(u))/(2*g*u)
        return -density*mp.log(density)*u
    baseEntropy=mp.quad(entropy_integrand,[-limit,0,limit])
    for scale in [.1,1,4]:
        m,b,g=args32(-3.25,scale,shape)
        emit('BirnbaumSaunders',[m,b,g],'Median',m+b)
        emit('BirnbaumSaunders',[m,b,g],'Entropy',baseEntropy+mp.log(b))
print('Birnbaum-Saunders:',len(rows),flush=True)

for eta in [1e-20,1e-6,.001,.1,.75,1,10,100,10000]:
    e,=args32(eta); scale=max(1,e); end=scale*mp.log1p(200/e)
    survival=lambda t:mp.exp(-e*mp.expm1(t/scale))
    mean=mp.quad(survival,[0,end/2,end])/scale
    second=2*mp.quad(lambda t:t*survival(t),[0,end/2,end])/(scale*scale)
    for rate in [.5,2]:
        e,b=args32(eta,rate)
        emit('Gompertz',[e,b],'Mean',mean/b)
        emit('Gompertz',[e,b],'Variance',(second-mean*mean)/(b*b))
print('Gompertz:',len(rows),flush=True)

for value in [1e-6,.1,.5,.69314718,.9,1,1.0000001,4.25,20,100,1000,10000]:
    lam,=args32(value)
    end=int(mp.ceil(lam+14*mp.sqrt(lam)+50)); probabilities=[mp.exp(-lam)]
    for k in range(1,end+1): probabilities.append(probabilities[-1]*lam/k)
    cumulative=mp.mpf(0); median=None; entropy=mp.mpf(0)
    for k,prob in enumerate(probabilities):
        cumulative+=prob
        if median is None and cumulative>=mp.mpf('.5'):median=k
        if prob: entropy-=prob*mp.log(prob)
    emit('Poisson',[lam],'Median',median)
    emit('Poisson',[lam],'Entropy',entropy)
    for xx in sorted(set([-.5,0,.5,1,2,30,max(0,int(lam)-2),int(lam),int(lam)+2])):
        x=mp.mpf(single(xx)); k=int(mp.floor(x))
        emit('Poisson',[lam],'Function',mp.exp(-lam)*lam**k/mp.factorial(k) if k>=0 and k==x else 0,x)
        emit('Poisson',[lam],'Distribution',mp.fsum(probabilities[:k+1]) if k>=0 else 0,x)
print('Poisson:',len(rows),flush=True)

for n in [0,1,2,5,30,100,1000]:
    for pp in [0,1e-6,.01,.3,.5,.7,.99,1]:
        p=mp.mpf(single(pp)); probabilities=[mp.binomial(n,k)*p**k*(1-p)**(n-k) for k in range(n+1)]
        cumulative=mp.mpf(0)
        for k,prob in enumerate(probabilities):
            cumulative+=prob
            if cumulative>=mp.mpf('.5'):median=k;break
        emit('Binomial',[n,p],'Median',median,signature=['i','f'])
        for xx in sorted(set([-1,0,.5,1,n//2,n,n+1])):
            x=mp.mpf(single(xx));k=int(mp.floor(x))
            emit('Binomial',[n,p],'Function',probabilities[k] if k==x and 0<=k<=n else 0,x,['i','f'])
            emit('Binomial',[n,p],'Distribution',mp.fsum(probabilities[:k+1]) if k>=0 else 0,x,['i','f'])
print('Binomial:',len(rows),flush=True)

for a0,b0 in [(.25,4.1),(1,5),(2,8),(10,20),(100,1000)]:
    a,b=args32(a0,b0)
    raw=[mp.rf(a,k)/mp.rf(b-k,k) for k in range(1,5)]
    variance=raw[1]-raw[0]**2
    excess=(raw[3]-4*raw[0]*raw[2]+6*raw[0]**2*raw[1]-3*raw[0]**4)/variance**2-3
    emit('BetaPrime',[a,b],'Excess',excess)
for d1,d2 in [(1,9),(4,12),(10,30),(100,1000)]:
    a=mp.mpf(d1)/2;b=mp.mpf(d2)/2
    raw=[mp.rf(a,k)/mp.rf(b-k,k) for k in range(1,5)];variance=raw[1]-raw[0]**2
    excess=(raw[3]-4*raw[0]*raw[2]+6*raw[0]**2*raw[1]-3*raw[0]**4)/variance**2-3
    emit('FisherSnedecor',[d1,d2],'Excess',excess,signature=['i','i'])
for d1,d2 in [(1,2),(4,12),(100,200)]:
    a,b=args32(d1/2,d2/2)
    for xx in [-100,-30,-3,0,.25,3,30,100]:
        x=mp.mpf(single(xx)); y=mp.exp(2*x)*a/b
        emit('FisherZ',[d1,d2],'Function',2*y**a/(mp.beta(a,b)*(1+y)**(a+b)),x)
        emit('FisherZ',[d1,d2],'Distribution',mp.betainc(a,b,0,y/(1+y),regularized=True),x)

for name,parameters in [('InverseChiSquare',[(1,),(2,),(5,),(12,),(100,)]),('Levy',[(0,.01),(3,1),(-2,10)]),
                        ('Kumaraswamy',[(.5,.7),(1,1),(2,3),(10,4)]),('Wigner',[(.01,),(1,),(10,)])]:
    for parameters0 in parameters:
        params=args32(*parameters0)
        if name=='InverseChiSquare':
            a=params[0]/2;b=mp.mpf('.5')
            density=lambda y:b**a/mp.gamma(a)*mp.exp(-a*y-b*mp.exp(-y))
            entropy=mp.quad(lambda y:-density(y)*(mp.log(density(y))-y),[-100,-10,0,10,200])
        elif name=='Levy':
            a=mp.mpf('.5');b=params[1]/2
            density=lambda y:b**a/mp.gamma(a)*mp.exp(-a*y-b*mp.exp(-y))
            entropy=mp.quad(lambda y:-density(y)*(mp.log(density(y))-y),[-100,-10,0,10,200])
        elif name=='Kumaraswamy':
            a,b=params;pdf=lambda x:a*b*x**(a-1)*(1-x**a)**(b-1)
            entropy=mp.quad(lambda x:-pdf(x)*mp.log(pdf(x)),[0,.5,1])
        else:
            r=params[0];pdf=lambda x:2*mp.sqrt(r*r-x*x)/(mp.pi*r*r)
            entropy=mp.quad(lambda x:-pdf(x)*mp.log(pdf(x)),[-r,0,r])
        emit(name,params,'Entropy',entropy,signature=['i'] if name=='InverseChiSquare' else None)
print('Entropy and moments:',len(rows),flush=True)

for values in [(0,1,2,4,2,2,1),(1,2,5,7,2,2,1),(0,0,1,1,2,2,1),(0,1,1,3,2,2,1),
               (-3,-1,2,4,.5,3,2),(10000,10001,10003,10004,4,.5,.3)]:
    a,b,c,d,n1,n3,alpha=args32(*values)
    norm=1/(alpha*(b-a)/n1+(alpha+1)*(c-b)/2+(d-c)/n3)
    def pdf(x):
        if x<=a or x>=d:return mp.mpf(0)
        if x<b:return norm*alpha*((x-a)/(b-a))**(n1-1)
        if x<c:return norm*(1+(alpha-1)*(c-x)/(c-b))
        return norm*((d-x)/(d-c))**(n3-1)
    cuts=sorted(set([a,b,c,d])); mean=mp.quad(lambda x:x*pdf(x),cuts)
    variance=mp.quad(lambda x:(x-mean)**2*pdf(x),cuts)
    emit('Trapezoidal',[a,b,c,d,n1,n3,alpha],'Mean',mean)
    emit('Trapezoidal',[a,b,c,d,n1,n3,alpha],'Variance',variance)
for value in [-2,-1,-.75,-.5,-.499,-.1,-.001,-1e-6,0,1e-6,.001,.1,1,2,10]:
    l,=args32(value)
    variance=mp.inf if l<=-.5 else mp.pi**2/3 if l==0 else 2/l**2*(1/(1+2*l)-mp.beta(l+1,l+1))
    emit('TukeyLambda',[l],'Mean',0 if l>-1 else mp.nan)
    emit('TukeyLambda',[l],'Variance',variance)
for k in [1,2,3,6,10,100,1000]:
    lo=mp.mpf(0);hi=mp.mpf(k)
    for i in range(100):
        mid=(lo+hi)/2
        if mp.gammainc(mp.mpf(k)/2,0,mid/2,regularized=True)<.5:lo=mid
        else:hi=mid
    emit('ChiSquare',[k],'Median',(lo+hi)/2,signature=['i'])

output=Path(__file__).with_name('distribution-repair.json')
output.write_text(json.dumps({'generator':'mpmath '+mp.__version__,'precision':mp.mp.dps,'cases':rows},indent=2)+'\n',encoding='utf-8')
print(len(rows),'reference cases ->',output,flush=True)
