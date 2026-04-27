import numpy as np
import math
from scipy.stats import poisson, nbinom, norm, binom, betabinom, bernoulli
from scipy.stats import beta as betadis
from scipy.special import kv
import warnings

def sichel(max, mu, sigma, v):
    warnings.filterwarnings("ignore")

    pmf = np.empty(max)
    alpha = math.sqrt(1 / math.pow(sigma, 2) + 2 * mu / sigma) 
    w = math.sqrt(math.pow(mu, 2) + math.pow(alpha, 2)) - mu
    for x in range(max):
        if x == 0:
            pmf[0] = math.pow(w / alpha, v) * kv(v, alpha) / kv(v, w)
        elif x == 1:
            pmf[1] = pmf[0] * (mu * w / alpha) * kv(v+1, alpha) / kv(v, alpha)
        else:
            pmf[x] = (2 * mu * w / math.pow(alpha, 2)) * ((x+v -1) / x) * pmf[x-1] + math.pow(mu * w / alpha, 2) / (x * (x-1)) * pmf[x-2]
    return pmf

def pmf(max, HMM, state, marker):
    if (HMM['emission'][state][marker]['distribution'] == 'NBI'):
        p = float(HMM['emission'][state][marker]['parameters']['p'])
        r = float(HMM['emission'][state][marker]['parameters']['r'])
        return nbinom.pmf(np.arange(0, max), p=p, n=r)

    elif (HMM['emission'][state][marker]['distribution'] == 'PO'):
        lam = float(HMM['emission'][state][marker]['parameters']['lambda'])
        return np.array([poisson.pmf(v, lam) for v in np.arange(0, max)])

    elif (HMM['emission'][state][marker]['distribution'] == 'BI'):
        p = float(HMM['emission'][state][marker]['parameters']['p'])
        n = int(HMM['emission'][state][marker]['parameters']['n'])
        return binom.pmf(np.arange(0, max), p=p, n=n)

    elif (HMM['emission'][state][marker]['distribution'] == 'ZAP'):
        lam = float(HMM['emission'][state][marker]['parameters']['lambda'])
        pi = float(HMM['emission'][state][marker]['parameters']['pi'])
        return np.array([pi if v == 0 else (1-pi) *  poisson.pmf(v, lam) / (1-poisson.pmf(0, lam)) for v in np.arange(0, max)])

    elif (HMM['emission'][state][marker]['distribution'] == 'BB'):
        alpha = float(HMM['emission'][state][marker]['parameters']['alpha'])
        beta = float(HMM['emission'][state][marker]['parameters']['beta'])
        if 'n' in HMM['emission'][state][marker]['parameters']:
            n = int(HMM['emission'][state][marker]['parameters']['n'])
            return betabinom.pmf(np.arange(0, max), a=alpha, b=beta, n=n)
        x = np.linspace(0, 1, 200)
        return np.array([betadis.cdf(x[i], a=alpha, b=beta) - betadis.cdf(x[i-1], a=alpha, b=beta) for i in range(1, 200)])

    elif (HMM['emission'][state][marker]['distribution'] == 'BNB'):
        alpha = float(HMM['emission'][state][marker]['parameters']['alpha'])
        beta = float(HMM['emission'][state][marker]['parameters']['beta'])
        r = float(HMM['emission'][state][marker]['parameters']['r'])
        return [math.exp((math.lgamma(v+r) + math.lgamma(alpha+r) + math.lgamma(beta+v) + math.lgamma(alpha+beta)) - (math.lgamma(v+1) + math.lgamma(r) + math.lgamma(alpha + r + beta + v) + math.lgamma(alpha) + math.lgamma(beta))) for v in np.arange(0, max)]

    elif (HMM['emission'][state][marker]['distribution'] == 'ZANBI'):
        p = float(HMM['emission'][state][marker]['parameters']['p'])
        r = float(HMM['emission'][state][marker]['parameters']['r'])
        pi = float(HMM['emission'][state][marker]['parameters']['pi'])
        return np.array([pi if v == 0 else (1-pi) *  nbinom.pmf(v, p=p, n=r) / (1-nbinom.pmf(0, p=p, n=r)) for v in np.arange(0, max)])

    elif (HMM['emission'][state][marker]['distribution'] == 'ZABNB'):
        alpha = float(HMM['emission'][state][marker]['parameters']['alpha'])
        beta = float(HMM['emission'][state][marker]['parameters']['beta'])
        r = float(HMM['emission'][state][marker]['parameters']['r'])
        pi = float(HMM['emission'][state][marker]['parameters']['pi'])
        v = 0
        zeroP = math.exp((math.lgamma(v+r) + math.lgamma(alpha+r) + math.lgamma(beta+v) + math.lgamma(alpha+beta)) - (math.lgamma(v+1) + math.lgamma(r) + math.lgamma(alpha + r + beta + v) + math.lgamma(alpha) + math.lgamma(beta)))
        return [pi if v == 0 else (1-pi) * math.exp((math.lgamma(v+r) + math.lgamma(alpha+r) + math.lgamma(beta+v) + math.lgamma(alpha+beta)) - (math.lgamma(v+1) + math.lgamma(r) + math.lgamma(alpha + r + beta + v) + math.lgamma(alpha) + math.lgamma(beta))) / (1-zeroP) for v in np.arange(0, max)]
    
    elif (HMM['emission'][state][marker]['distribution'] == 'GA'):
        mean = float(HMM['emission'][state][marker]['parameters']['mean'])
        std = float(HMM['emission'][state][marker]['parameters']['std'])
        return norm.pdf(np.arange(0, max), loc=mean, scale=std)

    elif (HMM['emission'][state][marker]['distribution'] == 'SI'):
        mu = float(HMM['emission'][state][marker]['parameters']['mu'])
        sigma = float(HMM['emission'][state][marker]['parameters']['sigma'])
        v = float(HMM['emission'][state][marker]['parameters']['v'])
        return sichel(max, mu, sigma, v)

    elif (HMM['emission'][state][marker]['distribution'] == 'ZASI'):
        mu = float(HMM['emission'][state][marker]['parameters']['mu'])
        sigma = float(HMM['emission'][state][marker]['parameters']['sigma'])
        v = float(HMM['emission'][state][marker]['parameters']['v'])
        pi = float(HMM['emission'][state][marker]['parameters']['pi'])
        pmf = sichel(max, mu, sigma, v)
        zeroP = pmf[0]
        pmf[0] = pi
        for i in range(1, len(pmf)):
            pmf[i] = (1.0-pi) * pmf[i] / (1-zeroP)
        return pmf
    
    elif (HMM['emission'][state][marker]['distribution'] == 'AB'):
        alpha = float(HMM['emission'][state][marker]['parameters']['alpha'])
        beta = float(HMM['emission'][state][marker]['parameters']['beta'])
        pi = float(HMM['emission'][state][marker]['parameters']['pi'])
        x = np.linspace(0, 1, 200)
        return np.array([pi] + [(1 - pi) * (betadis.cdf(x[i], a=alpha, b=beta) - betadis.cdf(x[i-1], a=alpha, b=beta)) for i in range(1, 200)])
    
    elif (HMM['emission'][state][marker]['distribution'] == 'B'):
        p = float(HMM['emission'][state][marker]['parameters']['p'])
        return np.array([bernoulli.pmf(v, p) for v in np.arange(0, max)])

    else:
        print("Unknown distribution")
        exit(1)
        
def param(distribution, mean, std, zeroFreq, n, max):
    mean = np.max([mean, 0.00001])
    std = np.max([std, 0.00001])

    if (distribution == 'NBI'):
        var = std**2
        if mean >= var:
            var = mean * 1.25
        p = mean / var
        return {"p": p, "r": mean**2 / (var-mean)}

    elif (distribution == 'PO'):
        return {"lambda": mean}
    
    elif (distribution == 'BI'):
        return {"p": mean / n, "n": int(n)}

    elif (distribution == 'BB'):
        pi = mean / max
        theta = np.max([0.001, (std**2 - max * pi * (1-pi)) / (max**2 * pi * (1-pi) - std**2)])
        return {"alpha": pi / theta, "beta": (1-pi) / theta, "n": int(max)}

    elif (distribution == 'BNB'):
        var = std**2
        if mean >= var:
            var = mean * 1.25
        r = np.max([1.0, mean**2 / (var-mean)])
        p = np.min([mean / var , 1.0])
        alpha = np.max([1.01, p * (p * (1-p) / (0.05**2) - 1)])
        beta = np.max([0.01, (1-p) * (p * (1-p) / (0.05**2) - 1)])
        return {"alpha": alpha, "beta": beta, "r": r}
    
    elif (distribution == 'ZAP'):
        return {"lambda": mean, "pi": zeroFreq}
    
    elif (distribution == 'ZANBI'):
        var = std**2
        if mean >= var:
            var = mean * 1.25
        p = mean / var
        r = mean**2 / (var-mean)
        return {"p": p, "r": r, "pi": zeroFreq}

    elif (distribution == 'ZABNB'):
        r = np.max([1.0, mean**2 / (std**2-mean)])
        p = mean / std**2 
        alpha = np.max([1.01, p * (p * (1-p) / (0.05**2) - 1)])
        beta = np.max([0.01, (1-p) * (p * (1-p) / (0.05**2) - 1)])
        return {"alpha": alpha, "beta": beta, "r": r, "pi": zeroFreq}
    
    elif (distribution == 'GA'):
        return {"mean": mean, "std": std}

    elif (distribution == 'SI'):
        w = mean / (std**2 / mean -1)
        alpha2 = np.max([0.0001, (w + mean)**2 - mean**2])
        s = np.max([0.1, (math.sqrt(alpha2 + mean**2) + mean) / alpha2])
        return {"mu": mean, "sigma": s, 'v': -0.5}

    elif (distribution == 'ZASI'):
        w = mean / (std**2 / mean -1)
        alpha2 = np.max([0.0001, (w + mean)**2 - mean**2])
        s = np.max([0.1, (math.sqrt(alpha2 + mean**2) + mean) / alpha2])
        return {"mu": mean, "sigma": s, 'v': -0.5, "pi": zeroFreq}
    
    elif (distribution == 'B'):
        return {"p": mean}
    
    else:
        print("Unknown distribution")
        exit(1)

def meth_param(distribution, mean, std, zeros):
    mean = np.max([mean, 0.00001])
    std = np.max([std, 0.00001])
    if (distribution == 'BI'):
        return {"p": mean}

    elif (distribution == 'BB'):
        #only valid if alpha, beta > 0 or in in terms of the mean and std: std**2 < mean * (1.0 - mean)
        x = (mean * (1.0 - mean) / std**2) - 1.0
        alpha = np.max([0.0001, mean * x])
        beta = np.max([0.0001, (1.0 - mean) * x])
        return {"alpha": alpha, "beta": beta}

    elif (distribution == 'AB'):
        #only valid if alpha, beta > 0 or in in terms of the mean and std: std**2 < mean * (1.0 - mean)
        x = (mean * (1.0 - mean) / std**2) - 1.0
        alpha = mean * x
        beta = (1.0 - mean) * x
        return {"alpha": alpha, "beta": beta, "pi": max(zeros, 0.0001)}

