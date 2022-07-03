import numpy as np
import pandas as pd
from scipy import interpolate
import math
e = math.e

# %%
class MBS:
    def __init__(self):
        self.init_balance = 1
        self.cpnRate = 0.045
        self.CPR = [0.1, 0.2, 0.3, 0.4]
        self.num_cpn_remained = 150

    def Refi(self, two_year_rate, ten_year_rate):
        mortgage_rate = 0.024 + 0.2 * two_year_rate + 0.6 * ten_year_rate
        refi = 0.2406 - 0.1389 * np.arctan(5.952 * (1.089 - self.cpnRate / mortgage_rate))
        return refi

    def Cashflow(self):
        cf_dt = pd.DataFrame(['Month, Balance', 'Interest', 'Scheduled Principal', 'Prepayment'])
        for t in month:
            print('hi')
        return 1


# %%

rate = [0.033, 0.034, 0.035, 0.040, 0.042, 0.044, 0.048, 0.0475]
tenor = [0.25, 0.5, 1, 5, 7, 10, 20, 30]


class Hull_White_2_Factor:
    def __init__(self):
        self.zero_curve = None
        self.type = 'IRS'
        self.instrument = {'Notion': 100}
        self.parameters = [1, 2, 3, 4, 5] # a, b, sigma, eta, rho

    def zc_input(self, tenor, rate): # input zero curve data
        self.zero_curve = {'tenor': tenor, 'rate': rate}

    def zc_interp(self): # zero curve interpolation
        return interpolate.interp1d(self.zero_curve['tenor'], self.zero_curve['rate'])

    def instrument_input(self, maturity, strike, bsprice): # input the data of the instrument for calibration
        self.instrument.update({'Type': 'IRS', 'Strike': strike, 'Price': bsprice})

    def Sigma(self, T, t):
        [a, b, s1, s2, r] = self.parameters
        res = s1**2 / (2*a**3) * (1-e**(-a*t))**2 * (1-e**(-2*a*T))\
              + s2**2 / (2*b**3) * (1-e**(-b*t))**2 * (1-e**(-2*b*T))\
              + 2*r*s1*s2 / (a*b*(a+b)) * (1-e**(-a*t)) * (1-e**(-b*t)) * (1-e**(-(a-b)*T))

        # Discount factor with IR interpolation
        def P(T):
            R = interpolate.interp1d(self.zero_curve['tenor'], self.zero_curve['rate'])(T)
            res = e**(-R*T)
            return res

        def Cap(t, T, N, X):


        return res*P(T)



m = Hull_White_2_Factor()
m.zc_input(tenor, rate)
m.instrument_input(maturity=[1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30],
                   strike=[0.0353, 0.0366, 0.0378, 0.0390, 0.0402, 0.0421, 0.0439, 0.0456, 0.0471, 0.0471, 0.0471],
                   bsprice=[0.1532, 0.6416, 1.3366, 2.0290, 2.7366, 4.2960, 6.5992, 9.6787, 12.2580, 14.0969, 15.7873])

