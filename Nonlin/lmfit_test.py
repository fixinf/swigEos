__author__ = 'const'
import lmfit

def f(params=lmfit.Parameters()):
    x = params['x'].value
    return [x**2 - 2]

fit_params = lmfit.Parameters()
fit_params.add('x', 2)
lmfit.minimize(f, fit_params)
print(fit_params['x'].value)
lmfit.report_fit(fit_params)