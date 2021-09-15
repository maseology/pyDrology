
from scipy.signal import lfilter
import numpy as np
from scipy.signal.signaltools import lfilter
from scipy.stats import linregress
from datetime import timedelta
# import matplotlib.pyplot as plt


########################################################
# automated extraction of the baseflow recession coefficient k as in Linsley, Kohler, Paulhus (1975) pg.230
# ref: Linsley, R.K., M.A. Kohler, J.L.H. Paulhus, 1975. Hydrology for Engineers 2nd ed. McGraw-Hill. 482pp.
########################################################
def recessionCoef(df):
    # collect recession dates
    d = df.to_dict('index')
    x, y = [], []
    for k,v in d.items():
        k1 = k + timedelta(days=1)
        if k1 in d: 
            if v['Val'] > d[k1]['Val']: 
                x.append(d[k1]['Val'])
                y.append(v['Val'])   
    # xt = x.copy()
    # yt = y.copy()

    while True:
        lnreg = linregress(x,y)
        # print(lnreg)
        if lnreg.rvalue > 0.995: break
        rem = []
        for i in range(len(x)): 
            if y[i] > lnreg.slope*x[i]: rem.append(i)
        for i in sorted(rem, reverse = True):
            del x[i]
            del y[i]

    x2 = np.vstack([x, np.zeros(len(x))]).T # x needs to be a column vector instead of a 1D vector for this https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
    a, _,_,_ = np.linalg.lstsq(x2, y, rcond=None) 
    # print(1/a[0]) # =k
    # plt.scatter(xt, yt, alpha=0.5)    
    # plt.plot(x2,a*x2,"r", alpha=0.75)
    # plt.show()
    return 1/a[0]


########################################################
# recession days
# time (in days) after peak discharge from which quick flow ceases and total flow is entirely slow flow
# ref: Linsley, R.K., M.A. Kohler, J.L.H. Paulhus, 1975. Hydrology for Engineers 2nd ed. McGraw-Hill. 482pp.
########################################################
def Ndays(careaKM2): 0.827 * careaKM2 ** 0.2


########################################################
# digital filter methods of automatic baseflow separaion
########################################################
def digitalFilter(v, a1, b0, b1, nPasses=1):
    if nPasses <= 1: return np.minimum(v, lfilter([b0,b1], [1.0,-a1], v, axis=0))
    f = v
    for i in range(nPasses):
        if (i+1) % 2 == 0: f = np.flip(f)
        f = lfilter([b0,b1], [1.0,-a1], f, axis=0)
        # f = np.minimum(v, f)
    return np.minimum(v,f)


########################################################
# the "Clarifica" technique (a.k.a. Graham method); named in (Clarifica, 2002) as a "forward and backward-step averaging approach."
########################################################
def clarifica(v):
    # ref: Clarifica Inc., 2002. Water Budget in Urbanizing Watersheds: Duffins Creek Watershed. Report prepared for the Toronto and Region Conservation Authority.
    # Clarifica method baseflow, 5-day avg running, 6-day min running   
    s = v.rolling('6D').min() # 6-day running minimum discharge   
    s = s.rolling('5D').mean().shift(-1, "D") # 5-day running average (3 days previous, 1 day ahead)
    return np.minimum(v,s)


########################################################
# returns grand estimate of baseflow
########################################################
def estimateBaseflow(df, k):
    nPass = 1
    # # Lyne, V. and M. Hollick, 1979. Stochastic time-variable rainfall-runoff modelling. Hydrology and Water Resources Symposium, Institution of Engineers Australia, Perth: 89-92.
    # #  k <- 0.925 # Ranges from 0.9 to 0.95 (Nathan and McMahon, 1990).  
    # #  nPasses = 3 # Commonly used (Chapman, 1999)
    # a = k
    # b = (1-k)/2
    # c = (1-k)/2
    # nPass = 3

    # # Chapman, T.G., 1991. Comment on the evaluation of automated techniques for base flow and recession analyses, by R.J. Nathan and T.A. McMahon. Water Resource Research 27(7): 1783-1784
    # a = (3*k-1)/(3-k)
    # b = (1-k)/(3-k)
    # c = (1-k)/(3-k)

    # # Chapman, T.G. and A.I. Maxwell, 1996. Baseflow separation - comparison of numerical methods with tracer experiments. Institute Engineers Australia National Conference. Publ. 96/05, 539-545.
    # a = k/(2-k)
    # b = (1-k)/(2-k)
    # c = 0

    # df['baseflow'] = digitalFilter(df[['Val']].to_numpy(),a,b,c,nPass)

    # Clarifica Inc., 2002. Water Budget in Urbanizing Watersheds: Duffins Creek Watershed. Report prepared for the Toronto and Region Conservation Authority.
    df['baseflow'] = clarifica(df)
