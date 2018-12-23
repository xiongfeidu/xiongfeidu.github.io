### Xiong-Fei Du
### Andrew ID: xiongfed
### 15-112 Section I, Spring 2016

### MATHLAB Term Project 15-112

### updated version 1.05
### release date 29 May 2018

import math

from math import *

from tkinter import *

####################################
# edit math module to adjust for domain errors
####################################

arcsin = asin
arccos = acos
arctan = atan
arcsinh = asinh
arccosh = acosh
arctanh = atanh

####################################
# other useful functions
####################################

ln = lambda x: math.log(x)
log = lambda x, base = 10: math.log(x, base)
sec = lambda x: 1/math.cos(x)
csc = lambda x: 1/math.sin(x)
cot = lambda x: 1/math.tan(x)

def roundif(x, epsilon = 1e-10):
    if round(x) == x or abs((round(x) - x) / x) < epsilon: x = round(x)
    return x

def roundsf(x, sf):
    if x == 0: return x
    elif abs(x) < 1: return round(x, sf - int(math.log(abs(x), 10)))
    else: return round(x, sf - int(math.log(abs(x), 10)) - 1)

def root(x, a):
    epsilon = 0.0001
    if x >= 0: return x**(1/a)
    elif abs(round(a) - a) > epsilon: return None
    elif a % 2 == 0: return None
    else: return -(-x)**(1/a)

####################################
# discontinuous functions
####################################

def heaviside(x):
    if x < 0: return 0
    elif x == 0: return 0.5
    elif x > 0: return 1

def sgn(x):
    if x < 0: return -1
    elif x == 0: return 0
    elif x > 0: return 1

def delta(x, h = 100):
    if x < 0: return 0
    elif x < 1/h: return h
    elif x > 1/h: return 0

def boxcar(x, a, b):
    if x == a or x == b: return 0.5
    elif x > a and x < b: return 1
    elif x < a or x > b: return 0

def rectangular(x):
    if equal(abs(x), 0.5): return 0.5
    elif abs(x) < 0.5: return 1
    elif abs(x) > 0.5: return 0

def ramp(x):
    if x < 0: return 0
    elif x >= 0: return x

def square(x, period = 1):
    if period == 0: return 0
    return sgn(sin(2*math.pi*x/period))

def triangle(x, period = 1):
    if period == 0: return 0
    return 2/math.pi*math.asin(sin(math.pi*x*2/period))

def sawtooth(x, period = 1):
    if period == 0: return 0
    if equal(x % 1, 0): return 1
    return -2/math.pi*math.atan(cot(math.pi*x/period))

####################################
# calculus module
####################################

def derivative(f, x):
    # use limit definition of a derivative
    sf = 8
    h = 0.0000001
    try:
        return roundsf((f(x + h/2) - f(x - h/2))/h, sf)
    except: return None

def secondDerivative(f, x):
    # use limit definition of second derivative
    sf = 5
    h = 0.00001
    try:
        return roundsf((f(x+h) - 2*f(x) + f(x-h))/h**2, sf)
    except: return None

def integral(f, a, b):
    #use Simpson's rule (approximates functions with piecewise parabolas)
    intervals = 4000
    sf = 10
    divide = 3
    runningSum = 0
    dx = (b - a)/intervals
    try:
        for i in range(intervals + 1):
            if i == 0: factor = 1
            elif i == intervals: factor = 1
            elif i % 2 == 0: factor = 2
            elif i % 2 == 1: factor = 4
            runningSum += factor * f(i*dx + a) * dx
        factor = 3
        return roundsf((runningSum/factor), sf)
    except: return None

def limit(f, a):
    digits = 6
    # use delta-epsilon definition of limit
    delta = 10**(-10**2)
    epsilon = 0.01
    try: # first try to plug it in
        return roundsf(f(a), digits)
    except: # approach from left and right
        try:
            left = f(a + delta)
            right = f(a - delta)
            assert(abs((left - right) / left) < epsilon)
            return roundsf((left + right)/2, digits)
        except: return "limit does not exist"

####################################
# monte carlo integration
####################################

def doubleIntegral(f, constraints, xMin, xMax, yMin, yMax,
                   totalPoints = 20000, string = True):
    # use monte carlo method to estimate a double integral on complex domain
    points = []
    maxCount = totalPoints * 100
    count = 0
    sf = 6
    # generate random points within bounding box and check if in constraints
    while len(points) < totalPoints:
        count += 1
        if count > maxCount: return None
        x, y = random.uniform(xMin, xMax), random.uniform(yMin, yMax)
        good = True
        for constraint in constraints:
            if constraint(x, y) == False: good = False; break
        if good: points.append((x, y))
    # area = area of bounding box
    area = (xMax - xMin)*(yMax - yMin)*totalPoints/count
    total = 0
    for x, y in points: total += f(x, y)
    # apply mean value theorem of integrals
    if string: return str(roundsf(total/totalPoints*area, sf)) + " +/- 1%"
    else: return total/totalPoints*area

def tripleIntegral(f, constraints, xMin, xMax, yMin, yMax, zMin, zMax,
                   totalPoints = 20000, string = True):
    # use monte carlo method to estimate a triple integral on complex domain
    points = []; count = 0; maxCount = totalPoints * 100
    sf = 6
    # generate random points within bounding box and check if in constraints
    while len(points) < totalPoints:
        count += 1
        if count > maxCount: return None
        x, y, z = (random.uniform(xMin, xMax), random.uniform(yMin, yMax),
                   random.uniform(zMin, zMax))
        good = True
        for constraint in constraints:
            if constraint(x, y, z) == False: good = False; break
        if good: points.append((x, y, z))
    # volume = volume of bounding box
    volume = (xMax - xMin)*(yMax - yMin)*(zMax - zMin)*totalPoints/count
    total = 0
    for x, y, z in points: total += f(x, y, z)
    # apply mean value theorem of integrals
    if string: return str(roundsf(total/totalPoints*volume, sf)) + " +/- 1%"
    else: return total/totalPoints*volume

####################################
# Newton's method
####################################

def zero(f, guess):
    # finds zero using Newton's method
    epsilon = 10**-10
    maxCount = 10000
    sf = 7
    try:
        x = guess
        count = 0
        y = f(x)
        prev = x
        while abs(y) > epsilon or abs(prev - x) > 10**-(sf + 1):
            count += 1
            assert(count < maxCount)
            der = derivative(f, x)
            y = f(x)
            prev = x
            x -= y/der
        return roundsf(x, sf)
    except:
        return roundsf(x, sf//2) if abs(y) < epsilon else None

####################################
# partial derivatives
####################################

def partialDerivativeX(f, x, y, z = None):
    # apply limit definition of a derivative
    sf = 7
    h = 0.0000001
    try:
        fx_minus = f(x - h/2, y, z)
        fx_plus = f(x + h/2, y, z)
        return roundsf((fx_plus - fx_minus)/h, sf)
    except: return None

def partialDerivativeY(f, x, y, z = None):
    # apply limit definition of a derivative
    sf = 7
    h = 0.0000001
    try:
        fx_minus = f(x, y - h/2, z)
        fx_plus = f(x, y + h/2, z)
        return roundsf((fx_plus - fx_minus)/h, sf)
    except: return None

def partialDerivativeZ(f, x, y, z):
    # apply limit definition of a derivative
    sf = 7
    h = 0.0000001
    try:
        fx_minus = f(x, y, z - h/2)
        fx_plus = f(x, y, z + h/2)
        return roundsf((fx_plus - fx_minus)/h, sf)
    except: return None

def secondPartialDerivativeX(f, x, y, z = None):
    # apply limit definition of second derivative
    sf = 4
    h = 0.00001
    try:
        fx = f(x, y, z)
        fx_plus = f(x + h, y, z)
        fx_minus = f(x - h, y, z)
        return roundsf((fx_plus - 2*fx + fx_minus)/h**2, sf)
    except: return None

def secondPartialDerivativeY(f, x, y, z = None):
    # apply limit definition of second derivative
    sf = 4
    h = 0.00001
    try:
        fx = f(x, y, z)
        fx_plus = f(x, y + h, z)
        fx_minus = f(x, y - h, z)
        return roundsf((fx_plus - 2*fx + fx_minus)/h**2, sf)
    except: return None

def secondPartialDerivativeZ(f, x, y, z):
    # apply limit definition of second derivative
    sf = 4
    h = 0.00001
    try:
        fx = f(x, y, z)
        fx_plus = f(x, y, z + h)
        fx_minus = f(x, y, z - h)
        return roundsf((fx_plus - 2*fx + fx_minus)/h**2, sf)
    except: return None

####################################
# vector calculus
####################################

def gradient(f, x, y):
    # gradient of f(x, y)
    answer = partialDerivativeX(f, x, y), partialDerivativeY(f, x, y)
    return answer if None not in answer else None

def gradient3(f, x, y, z):
    # gradient of f(x, y, z)
    answer = (partialDerivativeX(f, x, y, z), partialDerivativeY(f, x, y, z),
              partialDerivativeZ(f, x, y, z))
    return answer if None not in answer else None

def minimize3D(f, point):
    # gradient descent algorithm
    x, y = point
    gamma = 0.01 # step
    epsilon = 0.0001 # small number near zero
    sf = 4
    big = 10**10 # serves as an upper bound so that weird things don't happen
    try:
        maxCount = 1000000
        count = 0
        while (abs(partialDerivativeX(f, x, y)) > epsilon and
               abs(partialDerivativeY(f, x, y)) > epsilon):
            derX = partialDerivativeX(f, x, y)
            derY = partialDerivativeY(f, x, y)
            x -= gamma*derX
            y -= gamma*derY
            count += 1
            assert(count < maxCount and abs(derX) < big and abs(derY) < big)
        val = f(x, y, None)
        return roundsf(x, sf), roundsf(y, sf), roundsf(val, sf)
    except: return None

def maximize3D(f, point):
    # gradient ascent algorithm
    x, y = point
    gamma = 0.01
    epsilon = 0.0001
    sf = 4
    big = 10**10
    try:
        maxCount = 1000000
        count = 0
        while (abs(partialDerivativeX(f, x, y)) > epsilon and
               abs(partialDerivativeY(f, x, y)) > epsilon):
            derX = partialDerivativeX(f, x, y)
            derY = partialDerivativeY(f, x, y)
            x += gamma*derX
            y += gamma*derY
            count += 1
            assert(count < maxCount and abs(derX) < big and abs(derY) < big)
        val = f(x, y, None)
        return roundsf(x, sf), roundsf(y, sf), roundsf(val, sf)
    except: return None

def minimize(f, x):
    # gradient descent algorithm on one variable
    gamma = 0.01
    epsilon = 0.00001
    sf = 5
    big = 10**10
    try:
        maxCount = 100000
        count = 0
        der = 1
        while abs(der) > epsilon:
            der = derivative(f, x)
            x -= gamma*der
            count += 1
            assert(count < maxCount and abs(der) < big)
        return roundsf(x, sf), roundsf(f(x), sf)
    except: return None

def maximize(f, x):
    # gradient ascent algorithm on one variable
    gamma = 0.01
    epsilon = 0.00001
    sf = 5
    big = 10**10
    try:
        maxCount = 100000
        count = 0
        der = 1
        while abs(der) > epsilon:
            der = derivative(f, x)
            x += gamma*der
            count += 1
            assert(count < maxCount and abs(der) < big)
        return roundsf(x, sf), roundsf(f(x), sf)
    except: return None

def curl(P, Q, R, x, y, z):
    # P = x component, Q = y component, R = z component
    # finds curl of 3D vector field at a point
    try:
        i = partialDerivativeY(R, x, y, z) - partialDerivativeZ(Q, x, y, z)
        j = partialDerivativeZ(P, x, y, z) - partialDerivativeX(R, x, y, z)
        k = partialDerivativeX(Q, x, y, z) - partialDerivativeY(P, x, y, z)
        return i, j, k
    except: return None

def divergence(P, Q, R, x, y, z):
    # finds divergence of 3D vector field at a point
    try:
        return (partialDerivativeX(P, x, y, z) + partialDerivativeY(Q, x, y, z)
                + partialDerivativeZ(R, x, y, z))
    except: return None

def laplacian(f, x, y, z):
    # finds Laplacian of f(x, y, z) at a point
    try:
        return (secondPartialDerivativeX(f, x, y, z) +
                secondPartialDerivativeY(f, x, y, z) +
                secondPartialDerivativeZ(f, x, y, z))
    except: return None

####################################
# continuous probability
####################################

def normalDistribution(mean, stdev, a, b):
    # normal probability density function formula
    return (0.5*(1+math.erf((b - mean)/(stdev*math.sqrt(2)))) -
            0.5*(1+math.erf((a - mean)/(stdev*math.sqrt(2)))))

def normalPDF(mean, stdev, x):
    # normal probability density function formula
    return math.exp(-0.5*((x - mean)/stdev)**2)/stdev/sqrt(2*math.pi)

def inverseNormal(probability, mean, stdev):
    # use binary search divide and conquer algorithm
    if probability >= 1 or probability <= 0 or stdev <= 0: return None
    factor = 100
    lowerGuess = mean - factor*stdev
    upperGuess = mean + factor*stdev
    epsilon = 0.00000001 # degree of accuracy
    sf = 7
    if probability < epsilon: return -math.inf
    elif (1 - probability) < epsilon: return inf
    guess = (upperGuess + lowerGuess)/2
    # binary search while error is greater than epsilon
    while (abs(normalDistribution(mean, stdev, mean - factor*stdev, guess)
               - probability) > epsilon):
        if normalDistribution(mean,stdev,mean-factor*stdev,guess)>probability:
            upperGuess = guess
        else: lowerGuess = guess
        guess = (upperGuess + lowerGuess)/2
    return roundsf(guess, sf)

def tPDF(value, degreesOfFreedom):
    # student's t probability density function formula
    v = degreesOfFreedom
    try: return (math.gamma(v/2+0.5)/math.sqrt(v*math.pi)/math.gamma(v/2)*
            (1 + value**2/v)**(-v/2-0.5))
    except:
        try: return (math.exp(math.lgamma(v/2+0.5) - math.lgamma(v/2))/
            math.sqrt(v*math.pi)*(1 + value**2/v)**(-v/2-0.5))
        except: return None
        
def tDistribution(lower, upper, degreesOfFreedom):
    # student's t distribution is just an integral. Use integration
    maxBound = 30
    if lower < -maxBound: lower = -maxBound
    if upper > maxBound: upper = maxBound
    try:
        return integral(lambda x: tPDF(x, degreesOfFreedom), lower, upper)
    except: return None

def inverseT(probability, degreesOfFreedom):
    # use binary search divide and conquer algorithm
    if probability >= 1 or probability <= 0: return None
    factor = 200
    lowerGuess = -factor
    upperGuess = factor
    epsilon = 0.00001 # degree of accuracy
    sf = 5
    if probability < epsilon: return -math.inf
    elif (1 - probability) < epsilon: return inf
    guess = (upperGuess + lowerGuess)/2
    # binary search while error is greater than epsilon
    while (abs(tDistribution(-factor, guess, degreesOfFreedom)
               - probability) > epsilon):
        if tDistribution(-factor, guess, degreesOfFreedom)>probability:
            upperGuess = guess
        else: lowerGuess = guess
        guess = (upperGuess + lowerGuess)/2
    return roundsf(guess, sf)

def exponentialPDF(l, x):
    return l * math.exp(-l * x)

def exponentialDistribution(l, lower, upper):
    return math.exp(- l * lower) - math.exp(- l * upper)

def gammaPDF(alpha, beta, x):
    return (x**(alpha-1) * math.exp(-x/beta))/(math.gamma(alpha) * beta**alpha)

def gammaDistribution(alpha, beta, lower, upper):
    return integral(lambda x: gammaPDF(alpha, beta, x), lower, upper)

def betaPDF(alpha, beta, x):
    B = math.exp(math.lgamma(alpha+beta)-math.lgamma(alpha)-math.lgamma(beta))
    return B * x**(alpha-1) * (1-x)**(beta-1)

def betaDistribution(alpha, beta, lower, upper):
    return integral(lambda x: betaPDF(alpha, beta, x), lower, upper)

####################################
# discrete mathematics
####################################

### edit math module to prevent crashing

def nPr(n, r):
    # permutation formula
    result = 1
    for _ in range(r):
        result *= n
        n -= 1
    return result

def nCr(n, r):
    # combination formula
    if n - r < r: r = n - r
    result = 1
    L = [i + 1 for i in range(r)]
    for _ in range(r):
        result *= n
        for i in L:
            if result % i == 0:
                result //= i
                L.remove(i)
                break
        n -= 1
    for i in L:
        result //= i
    return result

def series(expression, start, end):
    # explicit series f(i) from i = start to i = end
    try:
        currentSum = 0
        for i in range(start, end + 1):
            currentSum += expression(i)
        return currentSum
    except: return None

def sequenceRecursive(expression, initial, iterations):
    # recursive sequence formula f(i) (i = previous term)
    try:
        i = initial
        result = [i]
        for _ in range(iterations):
            i = expression(i)
            result.append(i)
        return result
    except: return None

def sequenceExplicit(expression, start, end):
    # explicit sequence formula f(i) from i = start to i = end
    try:
        result = []
        for i in range(start, end + 1):
            result.append(expression(i))
        return result
    except: return None

####################################
# discrete probability
####################################

def binomial(trials, probability, successes):
    # binomial distribution formula
    try:
        assert(probability >= 0 and probability <= 1)
        combinations = nCr(trials, successes)
        probabilitySuccesses = probability**successes
        probabilityFailures = (1-probability)**(trials - successes)
        return combinations*probabilitySuccesses*probabilityFailures
    except: return None

def negativeBinomial(successes, probability, trials):
    # negative binomial distribution formula
    try:
        assert(probability >= 0 and probability <= 1)
        combinations = nCr(trials - 1, successes - 1)
        probabilitySuccesses = probability**successes
        probabilityFailures = (1-probability)**(trials - successes)
        return combinations*probabilitySuccesses*probabilityFailures
    except: return None

def geometric(probability, trials):
    # geometric distribution formula
    try:
        assert(probability >= 0 and probability <= 1)
        assert(type(trials) == int or equal(round(trials), trials))
        if trials == 0: return 0
        return probability*(1-probability)**(trials - 1)
    except: return None

def hypergeometric(N, K, n, k):
    # hypergeometric distribution formula
    # N = population, K = possible successes, n = draws, k = successes
    try: return nCr(K, k)*(nCr(N - K, n - k))/(nCr(N, n))
    except:
        if (type(N) == int and N >= 0 and type(K) == int and K >= 0 and
            type(n) == int and n >= 0 and type(k) == int and k >= 0): return 0
        else: return None

def poisson(expected, value):
    # poisson distribution formula
    try: return expected**value*math.exp(-expected)/math.factorial(value)
    except: return None

####################################
# regression
####################################

def linearRegression(L):
    avgX = avg([x for (x, y) in L]) # get averages
    avgY = avg([y for (x, y) in L])
    SSxx, SSxy, SSyy = 0, 0, 0
    for (x, y) in L: # find aggregate statistics
        SSxx += (x - avgX)**2
        SSxy += (x - avgX)*(y - avgY)
        SSyy += (y - avgY)**2
    slope = SSxy/SSxx
    yInt = avgY - slope*avgX
    # find correlation coefficient squared
    R2 = SSxy**2/SSxx/SSyy
    return slope, yInt, R2

def exponentialRegression(L):
    newL = [(x, math.log(y)) for (x, y) in L]
    b, a, R = linearRegression(newL)
    return math.exp(a), b, R

def logarithmicRegression(L):
    newL = [(math.log(x), y) for (x, y) in L]
    return linearRegression(newL)

def powerRegression(L):
    newL = [(math.log(x), math.log(y)) for (x, y) in L]
    b, a, R = linearRegression(newL)
    return math.exp(a), b, R

def isFlat(L):
    s = set()
    for (x, y) in L:
        s.add(y)
    return len(s) == 1
   
def polynomial(coeffs, x):
    result = 0
    for i in range(len(coeffs)):
        result += coeffs[i]*x**(len(coeffs) - i - 1)
    return result

def polynomialRegression(points, m):
    xs = [x for (x, y) in points]
    ys = [y for (x, y) in points]
    X = Matrix([[x**j for j in range(m+1)] for x in xs])
    XT = X.transpose()
    coeffsMatrix = (XT * X).inverse() * XT * Vector(ys)
    coeffs = coeffsMatrix.L
    yAvg = avg(ys)
    poly = lambda x: sum([coeffs[i]*x**i for i in range(len(coeffs))])
    SSreg = sum([(poly(x) - yAvg)**2 for x in xs])
    SStot = sum([(y - yAvg)**2 for y in ys])
    R2 = SSreg/SStot
    # check that R^2 is within reason
    assert(R2 <= 2), "Polynomial degree numerical error: R^2 = %f" % R2
    return coeffs[::-1], R2

####################################
# statistics
####################################

def avg(L):
    # average = sum/len
    if len(L) == 0: return None
    for element in L:
        if type(element) != int and type(element) != float: return None
    return sum(L)/len(L)

def standardDeviation(L):
    # sample standard deviation formula
    average = avg(L)
    if average == None or len(L) < 2: return None
    runningSum = 0
    for element in L: runningSum += (element - average)**2
    return (runningSum / (len(L) - 1))**0.5

def popStandardDeviation(L):
    # population standard deviation formula
    average = avg(L)
    if average == None or len(L) < 2: return None
    runningSum = 0
    for element in L: runningSum += (element - average)**2
    return (runningSum /len(L))**0.5

def median(L):
    # compute median
    L = sorted(L)
    if len(L) == 0: return None
    elif len(L) == 1: return L[0]
    elif len(L) % 2 == 1:
        return L[len(L) // 2]
    else: return (L[len(L) // 2 - 1] + L[len(L) // 2]) / 2

def firstQuartile(L):
    # compute first quartile (median of lower half)
    L = sorted(L)
    return median(L[:len(L)//2])

def thirdQuartile(L):
    # compute third quartile (median of upper half)
    L = sorted(L)
    return median(L[len(L)//2:]) if len(L)%2==0 else median(L[len(L)//2 + 1:])

def sumSquared(L):
    # compute sum of squares
    total = 0
    for element in L:
        if type(element) != int and type(element) != float: return None
        else: total += element**2
    return total

####################################
# statistical confidence and hypothesis testing
####################################

def zIntervalStats(mean, stdev, n, confidence):
    # finds the z interval given aggregate stats
    sf = 5
    z = inverseNormal((confidence + 1)/2, 0, 1)
    return (roundsf(mean - z*stdev/math.sqrt(n), sf),
            roundsf(mean + z*stdev/math.sqrt(n), sf))

def zIntervalData(L, stdev, confidence):
    # finds the z interval given a list
    mean = avg(L)
    n = len(L)
    return zIntervalStats(mean, stdev, n, confidence)

def zIntervalProportion(successes, trials, confidence):
    p = successes / trials
    return zIntervalStats(p, math.sqrt(p * (1 - p)), trials, confidence)

def zIntervalStatsTwoSample(mean1, stdev1, n1, mean2, stdev2, n2, confidence):
    sf = 5
    z = inverseNormal((confidence + 1)/2, 0, 1)
    mean = mean1 - mean2
    stdev = math.sqrt(stdev1**2/n1 + stdev2**2/n2)
    return (roundsf(mean - z*stdev, sf), roundsf(mean + z*stdev, sf))

def zIntervalDataTwoSample(L1, stdev1, L2, stdev2, confidence):
    m1 = avg(L1)
    n1 = len(L1)
    m2 = avg(L2)
    n2 = len(L2)
    return zIntervalStatsTwoSample(m1, stdev1, n1, m2, stdev2, n2, confidence)

def zIntervalProportionTwo(succ1, trials1, succ2, trials2, confidence):
    sf = 5
    p1 = succ1 / trials1
    p2 = succ2 / trials2
    diff = p1 - p2
    stdev = math.sqrt(p1 * (1 - p1) / trials1 + p2 * (1 - p2) / trials2)
    z = inverseNormal((confidence + 1)/2, 0, 1)
    return (roundsf(diff - z*stdev, sf), roundsf(diff + z*stdev, sf))

def tIntervalStats(mean, stdev, n, confidence):
    # finds the t interval given aggregate stats
    sf = 4
    t = inverseT((confidence + 1)/2, n - 1)
    return (roundsf(mean - t*stdev/math.sqrt(n), sf),
            roundsf(mean + t*stdev/math.sqrt(n), sf))

def tIntervalData(L, confidence):
    # finds the t interval given a list
    mean = avg(L)
    stdev = standardDeviation(L)
    n = len(L)
    return tIntervalStats(mean, stdev, n, confidence)

def tIntervalStatsTwoSample(mean1, stdev1, n1, mean2, stdev2, n2, confidence):
    sf = 4
    s12n1 = stdev1**2/n1
    s22n2 = stdev2**2/n2
    df = (s12n1 + s22n2)**2 / ((s12n1)**2/(n1 - 1) + (s22n2)**2/(n2 - 1))
    t = inverseT((confidence + 1)/2, df)
    mean = mean1 - mean2
    stdev = math.sqrt(s12n1 + s22n2)
    return (roundsf(mean - t*stdev, sf), roundsf(mean + t*stdev, sf))

def tIntervalDataTwoSample(L1, L2, confidence):
    mean1 = avg(L1)
    s1 = standardDeviation(L1)
    n1 = len(L1)
    mean2 = avg(L2)
    s2 = standardDeviation(L2)
    n2 = len(L2)
    return tIntervalStatsTwoSample(mean1, s1, n1, mean2, s2, n2, confidence)

####################################
# see this https://docs.python.org/3/reference/datamodel.html
####################################

####################################
# vectors
####################################

import copy

class Vector(object):

    def __init__(self, L):
        self.L = copy.copy(L)
        self.n = len(L)

    def mag(self):
        return math.sqrt(sum([self.L[i]**2 for i in range(self.n)]))

    def __abs__(self):
        return self.mag()

    def __len__(self):
        return self.n

    def __iter__(self):
        return iter(self.L)

    def __reversed__(self):
        return Vector(list(reversed(self.L)))

    def __contains__(self, item):
        return item in self.L

    def __getitem__(self, key):
        return self.L[key]

    def __setitem__(self, key, value):
        self.L[key] = value

    def __add__(self, other):
        assert(type(other) == Vector), "Second argument not a vector"
        assert(self.n == other.n), "Dimensional inconsistency"
        return Vector([self.L[i] + other.L[i] for i in range(self.n)])

    def __sub__(self, other):
        assert(type(other) == Vector), "Second argument not a vector"
        assert(self.n == other.n), "Dimensional inconsistency"
        return Vector([self.L[i] - other.L[i] for i in range(self.n)])

    def __mul__(self, other):
        assert(type(other) in (int, float)), "Second argument not scalar"
        return Vector([self.L[i] * other for i in range(self.n)])

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        assert(type(other) in (int, float)), "Second argument not scalar"
        return Vector([self.L[i] / other for i in range(self.n)])

    def __floordiv__(self, other):
        assert(type(other) in (int, float)), "Second argument not scalar"
        return Vector([self.L[i] // other for i in range(self.n)])

    def dot(self, other):
        assert(type(other) == Vector), "Second argument not a vector"
        assert(self.n == other.n), "Dimensional inconsistency"
        return sum([self.L[i] * other.L[i] for i in range(self.n)])

    def cross(self, other):
        assert(type(other) == Vector), "Second argument not a vector"
        assert(self.n == other.n), "Dimensional inconsistency"
        assert(self.n == 3), "Cross product only defined in three dimensions"
        a1, a2, a3 = tuple(self.L)
        b1, b2, b3 = tuple(other.L)
        return Vector([a2*b3-a3*b2, a3*b1-a1*b3, a1*b2 - a2*b1])

    def __round__(self, ndigits = None):
        if ndigits == None:
            return Vector(list(map(roundif, self.L)))
        else:
            return Vector(list(map(lambda x: round(x, ndigits), self.L)))

    def __eq__(self, other):
        if (self.n != other.n):
            return False
        for i in range(self.n):
            if self.L[i] != other.L[i]:
                return False
        return True

    def __neg__(self):
        return Vector([-x for x in self.L])

    def __pos__(self):
        return Vector([+x for x in self.L])

    def comp(self, other):
        return self.dot(other) / other.mag()

    def proj(self, other):
        return other * other.dot(self) / other.dot(other)

    def rej(self, other):
        return self - self.proj(other)

    def unit(self):
        mag = self.mag()
        return Vector([x/mag for x in self.L])

    def latex(self):
        s = str(self.L)
        s = s.replace("(", "")
        s = s.replace(")", "")
        s = s.replace("[", "")
        s = s.replace("]", "")
        return "$\\vec{v} = \langle " + s + " \\rangle$\n"

    def __repr__(self):
        return str(self.L)

def gs(L):
    # Gram-Schmidt process
    # Disadvantage: unstable, Advantage: can be parallel
    result = []
    for i in range(len(L)):
        v = L[i]
        for j in range(i):
            v -= L[i].proj(result[j])
        result.append(v)
    return result

def mgs(L):
    # modified Gram-Schmidt process
    # Advantage: stable, Disadvantage: inherently sequential
    result = []
    for i in range(len(L)):
        v = L[i]
        for j in range(i):
            v -= v.proj(result[j])
        result.append(v)
    return result

####################################
# matrices
####################################

class Matrix(object):

    def __init__(self, L):
        self.L = copy.deepcopy(L)
        self.m = len(L)
        self.n = len(L[0])
        for row in L:
            assert(len(row) == self.n), "Dimensional inconsistency"

    def row(self, i):
        return Vector(self.L[i])

    def col(self, j):
        return Vector([self.L[i][j] for i in range(self.m)])

    def dim(self):
        return self.m, self.n

    def __getitem__(self, key):
        return self.L[key]

    def isSquare(self):
        return self.m == self.n

    def tr(self):
        assert(self.isSquare()), "Matrix is not square"
        return sum([self.L[i][i] for i in range(self.n)])

    def minor(self, i, j):
        L = self.L
        m, n = self.m, self.n
        M = [[L[a][b] for b in range(n) if b != j] for a in range(m) if a != i]
        return Matrix(M).det()

    def cofactor(self, i, j):
        sgn = 1 if (i + j) % 2 == 0 else -1
        return sgn * self.minor(i, j)

    def det(self):
        assert(self.isSquare()), "Matrix is not square"
        L = self.L
        n = self.n
        if n == 1: return L[0][0]
        else: return sum([L[i][0] * self.cofactor(i, 0) for i in range(n)])

    def transpose(self):
        L = self.L
        m, n = self.n, self.m
        return Matrix([[L[j][i] for j in range(n)] for i in range(m)])

    def __add__(self, other):
        assert(type(other) == Matrix), "Second argument is not a matrix"
        assert(self.m == other.m), "Dimensional inconsistency"
        assert(self.n == other.n), "Dimensional inconsistency"
        m = self.m
        n = self.n
        A = self.L
        B = other.L
        return Matrix([[A[i][j] + B[i][j] for j in range(n)] for i in range(m)])

    def __sub__(self, other):
        assert(type(other) == Matrix), "Second argument is not a matrix"
        assert(self.m == other.m), "Dimensional inconsistency"
        assert(self.n == other.n), "Dimensional inconsistency"
        m = self.m
        n = self.n
        A = self.L
        B = other.L
        return Matrix([[A[i][j] - B[i][j] for j in range(n)] for i in range(m)])

    def __mul__(self, other):
        if type(other) in (int, float):
            c = other
            A = self.L
            m = self.m
            n = self.n
            return Matrix([[c*A[i][j] for j in range(n)] for i in range(m)])
        elif type(other) == Vector:
            assert(self.n == other.n), "Dimensional inconsistency"
            A = self
            B = other
            m = self.m
            n = other.n
            AB = [(A.row(i)).dot(B) for i in range(m)]
            return Vector(AB)
        assert(type(other) == Matrix), "Second argument is not a matrix"
        assert(self.n == other.m), "Dimensional inconsistency"
        A = self
        B = other
        m = self.m
        n = other.n
        AB = [[(A.row(i)).dot(B.col(j)) for j in range(n)] for i in range(m)]
        return Matrix(AB)

    def __rmul__(self, other):
        if type(other) in (int, float):
            return self * other
        elif type(other) == Vector:
            m = self.m
            A = other
            B = self
            AB = [A.dot(B.col(i)) for i in range(m)]
            return Matrix([AB])
        else:
            assert(False), "Invalid argument"

    def __truediv__(self, other):
        assert(type(other) in (int, float)), "Second argument not an int or float"
        c = other
        A = self.L
        m = self.m
        n = self.n
        return Matrix([[A[i][j]/c for j in range(n)] for i in range(m)])
    
    def __floordiv__(self, other):
        assert(type(other) in (int, float)), "Second argument not an int or float"
        c = other
        A = self.L
        m = self.m
        n = self.n
        return Matrix([[A[i][j]//c for j in range(n)] for i in range(m)])

    def __pow__(self, k):
        assert(type(k) == int), "Exponent is not an int"
        assert(self.isSquare()), "Matrix is not square"
        if k == 0: return identity(self.n)
        elif k > 0: return self ** (k - 1) * self
        else: return self.inverse() ** (-k)

    def __neg__(self):
        A = self.L
        m = self.m
        n = self.n
        return Matrix([[-A[i][j] for j in range(n)] for i in range(m)])

    def __pos__(self):
        A = self.L
        m = self.m
        n = self.n
        return Matrix([[+A[i][j] for j in range(n)] for i in range(m)])

    def __abs__(self):
        return self.det()

    def ref(self):
        m, n = self.m, self.n
        L = [[self.L[i][j] for j in range(n)] for i in range(m)]
        i, j = 0, 0
        while i < m:
            # eliminate previous rows
            for k in range(i):
                l = k
                while l < self.n and L[k][l] == 0: l += 1
                if l == n: break
                c = L[i][l] / L[k][l]
                L[i] = [L[i][l] - c * L[k][l] for l in range(n)]
            # make necessary swaps
            swap = False
            while j < self.n and L[i][j] == 0:
                for k in range(i + 1, m):
                    if L[k][j] != 0:
                        L[i], L[k] = L[k], L[i]
                        swap = True
                        break
                if not swap: j += 1
            if j == n: break
            if not swap: # set pivot to 1
                assert(L[i][j] != 0)
                L[i] = [L[i][k] / L[i][j] if k >= j else 0 for k in range(n)]
                i += 1
                j += 1
        return Matrix(L)

    def rref(self):
        L = self.ref().L
        for i in range(self.m):
            for k in range(i + 1, self.m):
                j = i
                while j < self.n and L[k][j] == 0: j += 1
                if j == self.n: continue
                c = L[i][j] / L[k][j]
                L[i] = [L[i][l] - c * L[k][l] for l in range(self.n)]
        return Matrix(L)

    def __round__(self, ndigits = None):
        L = self.L
        m, n = self.m, self.n
        if ndigits == None:
            L = [[roundif(L[i][j]) for j in range(n)] for i in range(m)]
        else:
            L = [[round(L[i][j], ndigits) for j in range(n)] for i in range(m)]
        return Matrix(L)

    def __eq__(self, other):
        assert(type(other) == Matrix), "Second argument is not a matrix"
        return self.L == other.L

    def inverse(self, check = True):
        if check: assert(self.isSquare()), "Matrix is not square"
        n = self.n
        f = lambda i, j, n: self.L[i][j] if j < n else 1 if i + n == j else 0
        L = [[f(i,j,n) for j in range(2*n)] for i in range(n)]
        B = Matrix(L).rref().L
        result = Matrix([[B[i][j + n] for j in range(n)] for i in range(n)])
        if check:
            inv = result.inverse(False)
            eps = [(self.L[i][j] - inv.L[i][j])/self.L[i][j]
                       for j in range(n) for i in range(n)]
            for x in eps:
                assert(x < 1e-4), "Matrix is not invertible with stability"
        return result

    def qr(self):
        detATA = (self.transpose() * self).det()
        assert(detATA != 0), "Matrix columns not linearly independent"
        QT = Matrix([
                u.unit().L
                    for u in mgs([Vector(v) for v in self.transpose().L])])
        R = QT * self
        return QT.transpose(), R

    def eigen(self):
        return None # Not yet implemented
        # http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
        m, n = self.m, self.n
        A = self
        T = identity(n)
        for _ in range(10):
            Q, R = A.qr()
            A = R * Q
            T = T * Q
        D = Matrix([[A.L[i][j] if i == j else 0
                     for j in range(n)] for i in range(m)])
        return A, T

    def latex(self):
        s = str(self)
        s = s.replace("],", "\\\\")
        s = s.replace("]]", "\\\\")
        s = s.replace("),", "\\\\")
        s = s.replace("))", "\\\\")
        s = s.replace(",", " &")
        s = s.replace(" [", "")
        s = s.replace("[[", "")
        s = s.replace(" (", "")
        s = s.replace("((", "")
        return "$\n\\begin{pmatrix}\n" + s + "\n\\end{pmatrix}\n$\n"

    def __repr__(self):
        if type(self.L) == list: return str(self.L).replace("],", "],\n")
        if type(self.L) == tuple: return str(self.L).replace("),", "),\n")

def identity(n):
    return Matrix([[1 if i == j else 0 for j in range(n)] for i in range(n)])

####################################
# other math functions
####################################

def inRange(x, a, b):
    # is a <= x <= b
    return x>=a and x<=b

def equal(a, b):
    # float1 == float2
    epsilon = 0.00000001
    return abs(a-b) < epsilon

def minus(a, b):
    epsilon = 0.0000000001
    if a == 0: return a
    elif abs(b - a) / abs(a) < epsilon: return 0
    else: return a - b

import random

####################################
# color stuff
####################################

def colorRGB(z, maxZ):
    # generates an RGB color based on the height of a point on a graph
    RGB = 255
    red = max(min(round(math.sqrt(2)*RGB*math.sin(math.pi*z/maxZ)), RGB), 0)
    green = max(min(round(math.sqrt(2)*RGB*math.cos(math.pi*z/maxZ)), RGB), 0)
    blue = max(min(round(math.sqrt(2)*RGB*math.sin(-math.pi*z/maxZ)), RGB), 0)
    return red, green, blue

def isLightColor(red, green, blue):
    # is RGB a light color - i.e. a combination is greater than 256//2
    upper = 128
    return (green > upper or (red > upper and green > upper) or
            (red > upper and blue > upper) or (green > upper and blue > upper))

def randomColor():
    # generate random color
    RGB = 255
    red = random.randint(0, RGB)
    green = random.randint(0, RGB)
    blue = random.randint(0, RGB)
    if isLightColor(red, green, blue): return randomColor()
    return red, green, blue

####################################
# function object
####################################

from types import FunctionType

class Function(object):
    def __init__(self, fxn):
        if type(fxn) == str:
            self.fxn = fxn
            self.generateLambda()
        elif isinstance(fxn, FunctionType):
            self.fxn = ""
            self.f = fxn
        else: assert(False), "Function input must be string or function"
        self.vectors = []
        self.color = "#%02x%02x%02x" % (randomColor())

    def generateLambda(self):
        try:
            assert(self.fxn.replace(" ", "") != "")
            self.f = eval("lambda %s: (%s)" % (self.strFunc(), self.fxn))
        except:
            self.f = eval("lambda %s: None" % self.strFunc())

    def evalWrapper(self, answer):
        # evaluate a function (string) at a given point
        try:
            assert(type(answer) in (int, float, tuple))
            if type(answer) == tuple:
                for a in answer:
                    assert(type(a) in (int, float))
            self.isValid = True
            return answer
        except:
            self.isValid = False
            return None

    def __repr__(self):
        return self.fxn

    def generateVectors(self, data):
        pass

    def plot(self, canvas, data):
        self.generateVectors(data)
        self.draw(canvas, data)

####################################
# 3D
####################################

class ThreeD(Function):

    def draw(self, canvas, data):
        try:
            # loop through one direction creating parallel line segments
            for row in self.vectors:
                prevVector = row[0]
                for i in range(len(row) - 1):
                    vector = row[i+1]
                    if prevVector == None or vector == None:
                        prevVector = vector; continue
                    try:
                        z, maxZ = vector[2], data.maxZ
                        color = "#%02x%02x%02x" % colorRGB(z, maxZ)
                        drawLine3D(prevVector, vector, canvas, data, color)
                    except: pass
                    prevVector = vector
            # loop through another direction creating parallel line segments
            for col in range(len(self.vectors[0])):
                if len(self.vectors) == 1: break
                prevVector = self.vectors[0][col]
                for row in range(len(self.vectors) - 1):
                    vector = self.vectors[row+1][col]
                    if prevVector == None or vector == None:
                        prevVector = vector; continue
                    try:
                        z, maxZ = vector[2], data.maxZ
                        color = "#%02x%02x%02x" % colorRGB(z, maxZ)
                        drawLine3D(prevVector, vector, canvas, data, color)
                    except: pass
                    prevVector = vector
        except: pass

def drawLine3D(prevVector, vector, canvas, data, color):
    # draws a line between two points given 2 vectors in 3D
    if (prevVector != None and vector != None and
        checkBounds(prevVector, vector, data)):
        prevPoint = vectorToPoint(prevVector, data)
        point = vectorToPoint(vector, data)
        canvas.create_line(prevPoint, point, fill = color)

class Cartesian3D(ThreeD):

    def strFunc(self):
        return "x, y"

    def evaluate(self, x, y):
        try: return self.evalWrapper(self.f(x = x, y = y))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a 2D list of vectors given a function
        self.vectors = []
        stepX = (data.maxX - data.minX)/data.increment
        stepY = (data.maxY - data.minY)/data.increment
        for i in range(data.increment + 1):
            row = []
            x = data.minX + i*stepX
            for j in range(data.increment + 1):
                y = data.minY + j*stepY
                z = self.evaluate(x = x, y = y)
                if z != None: vector = x, y, z
                else: vector = None
                row.append(vector)
            self.vectors.append(row)

class Parametric3D(ThreeD):

    def __init__(self, fxn, t_min = None, t_max = None):
        super().__init__(fxn)
        self.minT = t_min
        self.maxT = t_max
        self.isFuncInput = isinstance(fxn, FunctionType)
        self.isValid = self.isFuncInput or self.checkValidity()

    def checkValidity(self):
        # check if input is valid. If so, get min t and max t
        errorMessage = """Please type to enter a valid function
in the form 'x(t), y(t), z(t), t_min, t_max'"""
        functionLength = 5
        split = self.fxn.split(",")
        try:
            assert(len(split) >= functionLength)
            self.minT = eval(split[-2])
            self.maxT = eval(split[-1])
            for num in (self.minT, self.maxT):
                assert(type(num) == int or type(num) == float)
            return True
        except:
            self.errorMessage = errorMessage
            return False

    def strFunc(self):
        return "t"

    def generateVectors(self, data):
        # generate a 2D list of vectors given a parametric function
        self.isValid = self.isFuncInput or self.checkValidity()
        if not self.isValid: return
        self.vectors = []
        stepT = (self.maxT - self.minT)/data.incrementCurve
        row = []
        for i in range(data.incrementCurve+1):
            t = self.minT + i*stepT
            vector = self.evaluate(t)
            row.append(vector)
        self.vectors.append(row)

    def evaluate(self, t):
        try: tempAns = self.f(t = t)
        except: self.isValid = False; return
        if len(tempAns) == 3:
            return self.evalWrapper(tempAns)
        elif len(tempAns) == 5:
            x, y, z, _, _ = tempAns
            return self.evalWrapper((x, y, z))
        else:
            return self.evalWrapper(None)

class CylindricalRDependent(ThreeD):

    def strFunc(self):
        return "z, theta"

    def evaluate(self, z, theta):
        try: return self.evalWrapper(self.f(z = z, theta = theta))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a 2D list of vectors given a function
        self.vectors = []
        stepZ = (data.maxZ - data.minZ)/data.increment
        stepTheta = 2*math.pi/data.increment
        for i in range(data.increment + 1):
            row = []
            z = i*stepZ + data.minZ
            for j in range(data.increment+1):
                theta = j*stepTheta
                r = self.evaluate(z, theta)
                if r != None:
                    try: vector = cylindricalToCartesian(r, theta, z)
                    except: vector = None
                else: vector = None
                row.append(vector)
            self.vectors.append(row)

class CylindricalZDependent(ThreeD):

    def strFunc(self):
        return "r, theta"

    def evaluate(self, r, theta):
        try: return self.evalWrapper(self.f(r = r, theta = theta))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a 2D list of vectors given a function
        self.vectors = []
        maxR = (data.maxX**2 + data.maxY**2)**0.5
        stepR = maxR/data.increment
        stepTheta = 2*math.pi/data.increment
        for i in range(data.increment+1):
            row = []
            r = i*stepR
            for j in range(data.increment+1):
                theta = j*stepTheta
                z = self.evaluate(r = r, theta = theta)
                if z != None:
                    try: vector = cylindricalToCartesian(r, theta, z)
                    except: vector = None
                else: vector = None
                row.append(vector)
            self.vectors.append(row)

class Spherical(ThreeD):

    def strFunc(self):
        return "theta, phi"

    def evaluate(self, theta, phi):
        try: return self.evalWrapper(self.f(theta = theta, phi = phi))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a 2D list of vectors given a function
        self.vectors = []
        stepTheta = 2*math.pi/data.increment
        stepPhi = math.pi/data.increment
        for i in range(data.increment+1):
            row = []
            theta = i*stepTheta
            for j in range(data.increment+1):
                phi = j*stepPhi
                rho = self.evaluate(theta = theta, phi = phi)
                if rho != None:
                    try: vector = sphericalToCartesian(rho, theta, phi)
                    except: vector = None
                else: vector = None
                row.append(vector)
            self.vectors.append(row)

class VectorField3D(ThreeD):

    def strFunc(self):
        return "x, y, z"

    def evaluate(self, x, y, z):
        try: return self.evalWrapper(self.f(x = x, y = y, z = z))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a list of vectors given a vector field in 3D
        self.vectors = []
        stepX = (data.maxX - data.minX)/data.incrementField3D
        stepY = (data.maxY - data.minY)/data.incrementField3D
        stepZ = (data.maxZ - data.minZ)/data.incrementField3D
        for i in range(data.incrementField3D+1):
            x = stepX*i + data.minX
            for j in range(data.incrementField3D+1):
                y = stepY*j + data.minY
                for k in range(data.incrementField3D+1):
                    z = stepZ*k + data.minZ
                    tail = x, y, z
                    vector = self.evaluate(x = x, y = y, z = z)
                    try:
                        head = (tail[0]+vector[0], tail[1]+vector[1],
                            tail[2]+vector[2])
                        self.vectors.append((head, tail))
                    except: pass

    def draw(self, canvas, data):
        # draws vector field (from head to tail)
        for head, tail in self.vectors:
            tailPoint = vectorToPoint(tail, data)
            headPoint = vectorToPoint(head, data)
            canvas.create_line(tailPoint, headPoint, fill = self.color)

def distance(head, tail):
    # finds distance between two points (head and tail)
    return ((head[0] - tail[0])**2 + (head[1] - tail[1])**2 +
            (head[2] - tail[2])**2)**0.5

####################################
# 2D
####################################

class TwoD(Function):

    def draw(self, canvas, data):
        # draws a 2D graph given a list of position vectors
        if self.vectors == []: return
        prevVector = self.vectors[0]
        for i in range(len(self.vectors) - 1):
            vector = self.vectors[i + 1]
            drawLine2D(prevVector, vector, canvas, data, self.color)
            prevVector = vector

def drawLine2D(prevVector, vector, canvas, data, color):
    # check if vectors are valid and are in bounds
    if (prevVector != None and vector != None and
        checkBounds2D(prevVector, vector, data)):
        # draw line from previous vector to current vector
        prevPoint = (data.originX+prevVector[0]*data.scaleX,
                     data.originY-prevVector[1]*data.scaleY)
        point = (data.originX+vector[0]*data.scaleX,
                 data.originY-vector[1]*data.scaleY)
        canvas.create_line(prevPoint, point, fill = color)

class Cartesian2DyDep(TwoD):

    def strFunc(self):
        return "x"

    def evaluate(self, x):
        try: return self.evalWrapper(self.f(x = x))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a list of vectors given a function
        self.vectors = []
        stepX = (data.maxX - data.minX)/data.increment
        for i in range(data.increment + 1):
            x = i*stepX + data.minX
            y = self.evaluate(x = x)
            vector = (x, y) if y != None else None
            self.vectors.append(vector)

class Cartesian2DxDep(TwoD):

    def strFunc(self):
        return "y"

    def evaluate(self, y):
        try: return self.evalWrapper(self.f(y = y))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a list of vectors given a function
        self.vectors = []
        stepY = (data.maxY - data.minY)/data.increment
        for i in range(data.increment + 1):
            y = i*stepY + data.minY
            x = self.evaluate(y = y)
            vector = (x, y) if x != None else None
            self.vectors.append(vector)

class Polar(TwoD):

    def strFunc(self):
        return "theta"

    def evaluate(self, theta):
        try:
            r = self.evalWrapper(self.f(theta = theta))
            assert(type(r) in (int, float))
            return r
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a list of vectors given a function
        self.vectors = []
        stepTheta = 2*math.pi/data.increment
        for i in range(data.increment + 1):
            theta = i*stepTheta
            r = self.evaluate(theta = theta)
            vector = polarToCartesian(r, theta) if r != None else None
            self.vectors.append(vector)

class Parametric2D(TwoD):

    def __init__(self, fxn, t_min = None, t_max = None):
        super().__init__(fxn)
        self.minT = t_min
        self.maxT = t_max
        self.isFuncInput = isinstance(fxn, FunctionType)
        self.isValid = self.isFuncInput or self.checkValidity()

    def strFunc(self):
        return "t"

    def evaluate(self, t):
        try: tempAns = self.f(t = t)
        except: self.isValid = False; return
        if len(tempAns) == 2:
            return self.evalWrapper(tempAns)
        elif len(tempAns) == 4:
            x, y, _, _ = tempAns
            return self.evalWrapper((x, y))
        else:
            return self.evalWrapper(None)

    def checkValidity(self):
        # checks that an input is valid. If so, find min t and max t
        errorMessage = """Please type to enter a valid function
in the form 'x(t), y(t), t_min, t_max'"""
        functionLength = 4
        split = self.fxn.split(",")
        try:
            assert(len(split) >= functionLength)
            self.minT = eval(split[-2])
            self.maxT = eval(split[-1])
            assert(type(self.minT) in (int, float) and
                   type(self.maxT) in (int, float))
            return True
        except:
            self.errorMessage = errorMessage
            return False

    def generateVectors(self, data):
        # generate a list of vectors given a function
        self.isValid = self.isFuncInput or self.checkValidity()
        if not self.isValid: return
        self.vectors = []
        stepT = (self.maxT - self.minT)/data.incrementCurve
        for i in range(data.incrementCurve + 1):
            t = self.minT + i*stepT
            vector = self.evaluate(t)
            self.vectors.append(vector)

    def __repr__(self):
        return str(self.fxn)

class VectorField2D(TwoD):

    def strFunc(self):
        return "x, y"

    def evaluate(self, x, y):
        try: return self.evalWrapper(self.f(x = x, y = y))
        except: self.isValid = False

    def generateVectors(self, data):
        # generate a list of vectors given a vector field in 2D
        self.vectors = []
        stepX = (data.maxX - data.minX)/data.incrementField2D
        stepY = (data.maxY - data.minY)/data.incrementField2D
        for i in range(data.incrementField2D+1):
            x = stepX*i + data.minX
            for j in range(data.incrementField2D+1):
                y = stepY*j + data.minY
                tail = x, y
                try:
                    dx, dy = self.evaluate(x = x, y = y)
                    dx /= data.incrementField2D
                    dy /= data.incrementField2D
                    head = x + dx, y + dy
                    self.vectors.append((head, tail))
                except: pass

    def draw(self, canvas, data):
        # draws vectors in vector field
        scaleX, scaleY = data.scaleX, data.scaleY
        for (x0, y0), (x1, y1) in self.vectors:
            tailPt = (data.originX + scaleX * x1, data.originY - scaleY * y1)
            headPt = (data.originX + scaleX * x0, data.originY - scaleY * y0)
            canvas.create_line(tailPt, headPt, fill = self.color)

class Point(TwoD):

    def __init__(self, x, y):
        # creates a point that can be drawn on a canvas
        self.x = x
        self.y = y
        self.radius = 2

    def generateVectors(self, data):
        pass

    def draw(self, canvas, data, label = True, color = "black"):
        # draws point on canvas
        scaleX, scaleY = data.scaleX, data.scaleY
        pointX, pointY = data.originX+scaleX*self.x,data.originY-scaleY*self.y
        x, y = self.x, self.y
        # check if point is on the plane
        if inRange(x, data.minX, data.maxX) and inRange(y,data.minY,data.maxY):
            canvas.create_oval(pointX - self.radius, pointY - self.radius,
                               pointX + self.radius, pointY + self.radius,
                               fill = color)
            if type(label) == str:
                canvas.create_text(pointX, pointY, anchor = SW, text = label)
            elif label:
                canvas.create_text(pointX, pointY, anchor = SW,
                                   text = "(%s, %s)" % (self.x, self.y))


####################################
# differential equations
####################################

class DiffEq(VectorField2D):

    def __init__(self, fxn, x0 = None, y0 = None):
        super().__init__(fxn)
        self.x0 = x0
        self.y0 = y0
        self.solution = []

    def strFunc(self):
        return "x, y"

    def evaluate(self, x, y):
        try: tempSlope = self.f(x = x, y = y)
        except: self.isValid = False; return
        if tempSlope == None:
            return self.evalWrapper(None)
        elif type(tempSlope) == tuple:
            slope = tempSlope[0]
            if type(slope) in (float, int):
                return self.evalWrapper(slope)
            else:
                return self.evalWrapper(None)
        else:
            return self.evalWrapper(tempSlope)

    def generateVectors(self, data):
        # creates a set of slope vectors given a diff eq
        self.vectors = []
        stepX = (data.maxX - data.minX)/data.incrementField2D
        stepY = (data.maxY - data.minY)/data.incrementField2D
        for i in range(data.incrementField2D+1):
            x = stepX*i + data.minX
            for j in range(data.incrementField2D+1):
                y = stepY*j + data.minY
                tail = x, y
                slope = self.evaluate(x = x, y = y)
                if slope == None: continue
                angle = math.atan(slope*stepX/stepY)
                head = (x + math.cos(angle) * stepX * data.scaleFactor,
                        y + math.sin(angle) * stepY * data.scaleFactor)
                self.vectors.append((head, tail))
        self.getSolution(data)
        return self.vectors

    def getSolution(self, data):
        # gets initial condition of solution
        self.solution = []
        if self.fxn != "":
            split = self.fxn.split(",")
            try:
                assert(len(split) > 2)
                self.x0 = eval(split[-2])
                self.y0 = eval(split[-1])
                self.solveDiffEq(data)
            except:
                self.x0 = None
                self.y0 = None
        elif type(self.x0) in (int, float) and type(self.y0) in (int, float):
            self.solveDiffEq(data)
        else: pass

    def solveDiffEq(self, data):
        # solves diff eq using Runge-Kutta 4 method
        factor = 4
        step = (data.maxX - data.minX)/data.stepFactor/factor
        x, y = self.x0, self.y0
        xNext, yNext = x, y
        prevSlope = 0
        self.solution = [(x, y)]
        self.goForward(data, x, y, xNext, yNext, step, prevSlope, factor)
        self.goBackward(data, x, y, xNext, yNext, step, prevSlope, factor)
        return self.solution

    def goForward(self, data, x, y, xNext, yNext, step, prevSlope, factor):
        # go forward from initial condition
        while checkBounds2D((x,y),(xNext,yNext),data):
            for i in range(factor):
                slope1 = self.evaluate(x = x, y = y)
                slope2 = self.evaluate(x = x + step/2, y = y + slope1*step/2)
                slope3 = self.evaluate(x = x + step/2, y = y + slope2*step/2)
                slope4 = self.evaluate(x = x + step, y = y + slope3*step)
                num = 6
                try: slope = (slope1 + 2*slope2 + 2*slope3 + slope4)/num
                except: return
                if checkDifferentiability(slope, prevSlope) == False:
                    self.solution.append((x, y))
                    return
                x += step
                y += slope*step
                prevSlope = slope
            self.solution.append((x, y))

    def goBackward(self, data, x, y, xNext, yNext, step, prevSlope, factor):
        # go backwards from initial condition
        while checkBounds2D((x,y),(xNext,yNext), data):
            for i in range(factor):
                slope1 = self.evaluate(x = x, y = y)
                slope2 = self.evaluate(x = x - step/2, y = y - slope1*step/2)
                slope3 = self.evaluate(x = x - step/2, y = y - slope2*step/2)
                slope4 = self.evaluate(x = x - step, y = y - slope3*step)
                num = 6
                try: slope = (slope1 + 2*slope2 + 2*slope3 + slope4)/num
                except: return
                if checkDifferentiability(slope, prevSlope) == False:
                    self.solution.insert(0, (x, y))
                    return
                x -= step
                y -= slope*step
                prevSlope = slope
            self.solution.insert(0, (x, y))

    def draw(self, canvas, data):
        # draw slope field onto plane as a vector field and a solution
        super().draw(canvas, data)
        if self.solution != [] and self.vectors != []:
            prevVector = self.solution[0]
            for i in range(len(self.solution) - 1):
                vector = self.solution[i + 1]
                drawLine2D(prevVector, vector, canvas, data, self.color)
                prevVector = vector            

class SecondOrderDE(TwoD):

    def __init__(self, fxn, x0 = None, y0 = None, yPrime0 = None):
        super().__init__(fxn)
        if type(fxn) == str:
            self.isStr = True
            self.getInitialValues()
        else:
            self.x0 = x0
            self.y0 = y0
            self.yPrime0 = yPrime0
            self.isStr = False
            self.solve = True
        self.vectors = []

    def strFunc(self):
        return "x, y, Dy"

    def evaluate(self, x, y, Dy):
        try: tempSlope = self.f(x = x, y = y, Dy = Dy)
        except: self.isValid = False; return
        if tempSlope == None:
            return self.evalWrapper(None)
        elif type(tempSlope) == tuple:
            slope = tempSlope[0]
            if type(slope) in (float, int):
                return self.evalWrapper(slope)
            else:
                return self.evalWrapper(None)
        else:
            return self.evalWrapper(tempSlope)

    def generateVectors(self, data):
        # solves 2nd order initial value problem with Runge Kutta 4 algorithm
        if self.isStr: self.getInitialValues()
        if not self.solve:
            self.vectors = []
            return
        x, y, yPrime = self.x0, self.y0, self.yPrime0
        factor = 4
        step = (data.maxX - data.minX)/data.stepFactor/factor
        self.vectors = [(x, y)]
        self.goForward(x, y, yPrime, data, factor, step)
        self.goBackwards(x, y, yPrime, data, factor, step)
        return self.vectors

    def goForward(self, x, y, yPrime, data, factor, step):
        # go forward from the initial value
        while checkBounds2D((x, y), (x, y), data):
            try:
                for i in range(factor):
                    Dy, DDy = self.goForwardRK4(x, y, yPrime, step)
                    x += step
                    y += Dy
                    yPrime += DDy
            except: break
            self.vectors.append((x, y))
        if self.vectors != []: self.vectors.pop()

    def goBackwards(self, x, y, yPrime, data, factor, step):
        # go backwards from the initial value
        while checkBounds2D((x, y), (x, y), data):
            try:
                for i in range(factor):
                    Dy, DDy = self.goBackwardRK4(x, y, yPrime, step)
                    x -= step
                    y -= Dy
                    yPrime -= DDy
            except: break
            self.vectors.insert(0, (x, y))
        if self.vectors != []: self.vectors.pop(0)

    def goForwardRK4(self, x, y, yPrime, step):
        # find delta y and delta squared y using Runge Kutta 4 algorithm
        Dy1 = yPrime*step
        DDy1 = self.evaluate(x = x, y = y, Dy = yPrime)*step
        x += step/2
        y2 = y + Dy1/2; yPrime2 = yPrime + DDy1/2
        Dy2 = yPrime2*step
        DDy2 = self.evaluate(x = x, y = y2, Dy = yPrime2)*step
        y3 = y + Dy2/2; yPrime3 = yPrime + DDy2/2
        Dy3 = yPrime3*step
        DDy3 = self.evaluate(x = x, y = y3, Dy = yPrime3)*step
        x += step/2
        y4 = y + Dy3; yPrime4 = yPrime + DDy3
        Dy4 = yPrime4*step
        DDy4 = self.evaluate(x = x, y = y4, Dy = yPrime4)*step
        num = 6
        Dy = (Dy1 + 2*Dy2 + 2*Dy3 + Dy4)/num
        DDy = (DDy1 + 2*DDy2 + 2*DDy3 + DDy4)/num
        return Dy, DDy

    def goBackwardRK4(self, x, y, yPrime, step):
        # find delta y and delta squared y using Runge Kutta 4 algorithm
        Dy1 = yPrime*step
        DDy1 = self.evaluate(x = x, y = y, Dy = yPrime)*step
        x -= step/2
        y2 = y - Dy1/2; yPrime2 = yPrime - DDy1/2
        Dy2 = yPrime2*step
        DDy2 = self.evaluate(x = x, y = y2, Dy = yPrime2)*step
        y3 = y - Dy2/2; yPrime3 = yPrime - DDy2/2
        Dy3 = yPrime3*step
        DDy3 = self.evaluate(x = x, y = y3, Dy = yPrime3)*step
        x -= step/2
        y4 = y - Dy3; yPrime4 = yPrime - DDy3
        Dy4 = yPrime4*step
        DDy4 = self.evaluate(x = x, y = y4, Dy = yPrime4)*step
        num = 6
        Dy = (Dy1 + 2*Dy2 + 2*Dy3 + Dy4)/num
        DDy = (DDy1 + 2*DDy2 + 2*DDy3 + DDy4)/num
        return Dy, DDy

    def getInitialValues(self):
        # gets initial condition of solution
        split = self.fxn.split(",")
        conditions = 3
        try:
            assert(len(split) == conditions + 1)
            self.x0 = eval(split[1])
            self.y0 = eval(split[2])
            self.yPrime0 = eval(split[3])
            self.solve = True
        except: self.solve = False

Order1ODE = DiffEq
Order2ODE = SecondOrderDE

class OrdDiffEq(object):

    def __init__(self, f, x0, y0, x_min, x_max, step):
        self.f = f
        self.x0 = x0
        self.y0 = y0
        self.minX = x_min
        self.maxX = x_max
        self.step = step

    def solveEuler(self):
        f = self.f
        x0, y0 = self.x0, self.y0
        x_min, x_max = self.minX, self.maxX
        step = self.step
        x, y = x0, y0
        solution = [(x, y)]
        while x <= x_max:
            slope = f(x = x, y = y)
            x += step
            y += slope*step
            solution.append((x, y))
        x, y = x0, y0
        while x >= x_min:
            slope = f(x = x, y = y)
            x -= step
            y -= slope*step
            solution.insert(0, (x, y))
        return solution

    def solveHeun(self):
        f = self.f
        x0, y0 = self.x0, self.y0
        x_min, x_max = self.minX, self.maxX
        step = self.step
        x, y = x0, y0
        solution = [(x, y)]
        while x <= x_max:
            slope1 = f(x = x, y = y)
            slope2 = f(x = x + step, y = y + slope1*step)
            slope = (slope1 + slope2)/2
            x += step
            y += slope*step
            solution.append((x, y))
        x, y = x0, y0
        while x >= x_min:
            slope1 = f(x = x, y = y)
            slope2 = f(x = x - step, y = y - slope1*step)
            slope = (slope1 + slope2)/2
            x -= step
            y -= slope*step
            solution.insert(0, (x, y))
        return solution

    def solveRK4(self):
        f = self.f
        x0, y0 = self.x0, self.y0
        x_min, x_max = self.minX, self.maxX
        step = self.step
        x, y = x0, y0
        solution = [(x, y)]
        while x <= x_max:
            slope1 = f(x = x, y = y)
            slope2 = f(x = x + step/2, y = y + slope1*step/2)
            slope3 = f(x = x + step/2, y = y + slope2*step/2)
            slope4 = f(x = x + step, y = y + slope3*step)
            slope = (slope1 + 2*slope2 + 2*slope3 + slope4)/6
            x += step
            y += slope*step
            solution.append((x, y))
        x, y = x0, y0
        while x >= x_min:
            slope1 = f(x = x, y = y)
            slope2 = f(x = x - step/2, y = y - slope1*step/2)
            slope3 = f(x = x - step/2, y = y - slope2*step/2)
            slope4 = f(x = x - step, y = y - slope3*step)
            slope = (slope1 + 2*slope2 + 2*slope3 + slope4)/6
            x -= step
            y -= slope*step
            solution.insert(0, (x, y))
        return solution

class HeatEq(ThreeD):

    def __init__(self, fxn, alpha, ic,
                 t0 = 0, t_min = -5, t_max = 5, x_min = -5, x_max = 5):
        super().__init__(fxn)
        self.incrementFactor = 64
        self.alpha = alpha
        self.ic = ic
        self.t0 = t0
        self.t_min = t_min
        self.t_max = t_max
        self.x_min = x_min
        self.x_max = x_max
        self.vectors = []

    def generateVectors(self, data, skip = True):
        # note: look into crank-nicolson method as a more accurate solution
        increment = data.increment * 2
        t_min, t_max = self.t_min, self.t_max
        x_min, x_max = self.x_min, self.x_max
        alpha = self.alpha
        dt = (t_max - t_min) / increment / self.incrementFactor
        epsilon = 1e-4
        dx_min = ((2 + epsilon) * alpha * dt) ** 0.5
        dx_main = (x_max - x_min) / increment / self.incrementFactor
        dx = max(dx_min, dx_main)
        truncate = math.ceil(((x_max - x_min) / increment) / dx)
        r = alpha * dt / (dx**2)
        self.vectors = []
        self.ftcs(self.ic, r,
                  self.t0, t_min, t_max, dt,
                  x_min, x_max, dx, skip)
        self.btcs(self.ic, r,
                  self.t0, t_min, t_max, dt,
                  x_min, x_max, dx, skip)
        vs = self.vectors
        if skip:
            for i in range(len(vs)):
                vs[i] = [vs[i][j]
                             for j in range(len(vs[i]))
                                 if j % truncate == 0]
        return vs

    def solve(self, dt, dx):
        # note: look into crank-nicolson method as a more accurate solution
        old = self.vectors
        t_min, t_max = self.t_min, self.t_max
        x_min, x_max = self.x_min, self.x_max
        alpha = self.alpha
        t0 = self.t0
        epsilon = 1e-4
        dx_min = ((2 + epsilon) * alpha * dt) ** 0.5
        dx = max(dx_min, dx)
        r = alpha * dt / (dx**2)
        self.vectors = []
        self.ftcs(self.ic, r, t0, t_min, t_max, dt, x_min, x_max, dx, False)
        self.btcs(self.ic, r, t0, t_min, t_max, dt, x_min, x_max, dx, False)
        solution = self.vectors
        self.vectors = old
        return solution

    def ftcs(self, f, r, t0, t_min, t_max, dt, x_min, x_max, dx, skip):
        r2 = 1 - 2*r
        t = t0
        count = 0
        inc = self.incrementFactor
        while t <= t_max:
            i = 0
            x = x_min
            step = []
            while x <= x_max:
                if t == t0:
                    step.append((x, t, f(x)))
                elif x == x_min:
                    step.append((x, t, self.vectors[-1][0][-1]))
                elif x + dx > x_max:
                    step.append((x, t, self.vectors[-1][-1][-1]))
                else:
                    current = r2 * self.vectors[-1][i][-1]
                    current += r * self.vectors[-1][i - 1][-1]
                    current += r * self.vectors[-1][i + 1][-1]
                    step.append((x, t, current))
                i += 1
                x += dx
            if skip and count % inc != 1 and len(self.vectors) > 0:
                self.vectors.pop()
            self.vectors.append(step)
            t += dt
            count += 1

    def btcs(self, f, r, t0, t_min, t_max, dt, x_min, x_max, dx, skip):
        r2 = 1 + 2*r
        t = t0 - dt
        count = 1
        inc = self.incrementFactor
        while t > t_min:
            i = 0
            x = x_min
            step = []
            while x < x_max:
                if t == t0:
                    step.append((x, t, f(x)))
                elif x == x_min:
                    step.append((x, t, self.vectors[0][0][-1]))
                elif x + dx >= x_max:
                    step.append((x, t, self.vectors[0][-1][-1]))
                else:
                    current = r2 * self.vectors[0][i][-1]
                    current -= r * self.vectors[0][i - 1][-1]
                    current -= r * self.vectors[0][i + 1][-1]
                    step.append((x, t, current))
                i += 1
                x += dx
            if skip and count % inc != 1 and len(self.vectors) > 0:
                self.vectors.pop(0)
            self.vectors.insert(0, step)
            count += 1
            t -= dt

    def getInitialValues(self):
        split = self.fxn.split(",")
        # TODO: INCOMPLETE

class WaveEq(ThreeD):

    def __init__(self, fxn, c2, ic, icPrime,
                 t0 = 0, t_min = -5, t_max = 5, x_min = -5, x_max = 5):
        super().__init__(fxn)
        self.incrementFactor = 64
        self.c2 = c2
        self.ic = ic
        self.icPrime = icPrime
        self.t0 = t0
        self.t_min = t_min
        self.t_max = t_max
        self.x_min = x_min
        self.x_max = x_max
        self.vectors = []
        self.derVectors = []

    def generateVectors(self, data, skip = True):
        # note: look into crank-nicolson method as a more accurate solution
        increment = data.increment * 2
        t_min, t_max = self.t_min, self.t_max
        x_min, x_max = self.x_min, self.x_max
        c2 = self.c2
        dt = (t_max - t_min) / increment / self.incrementFactor
        epsilon = 1e-4
        dx_min = ((2 + epsilon) * c2 * dt) ** 0.5
        dx_main = (x_max - x_min) / increment / self.incrementFactor
        dx = max(dx_min, dx_main)
        truncate = math.ceil(((x_max - x_min) / increment) / dx)
        self.vectors = []
        self.derVectors = []
        self.ftcs(self.ic, self.icPrime, self.c2,
                  self.t0, t_min, t_max, dt,
                  x_min, x_max, dx, skip)
        self.btcs(self.ic, self.icPrime, self.c2,
                  self.t0, t_min, t_max, dt,
                  x_min, x_max, dx, skip)
        vs = self.vectors
        if skip:
            for i in range(len(vs)):
                vs[i] = [vs[i][j]
                             for j in range(len(vs[i]))
                                 if j % truncate == 0]
        return vs

    def solve(self, dt, dx):
        # note: look into crank-nicolson method as a more accurate solution
        old = self.vectors
        oldDer = self.derVectors
        t_min, t_max = self.t_min, self.t_max
        x_min, x_max = self.x_min, self.x_max
        alpha = self.alpha
        t0 = self.t0
        epsilon = 1e-4
        dx_min = ((2 + epsilon) * alpha * dt) ** 0.5
        dx = max(dx_min, dx)
        r = alpha * dt / (dx**2)
        self.vectors = []
        self.derVectors = []
        self.ftcs(self.ic, self.icPrime, self.c2,
                  self.t0, t_min, t_max, dt, x_min, x_max, dx, False)
        self.btcs(self.ic, self.icPrime, self.c2,
                  self.t0, t_min, t_max, dt, x_min, x_max, dx, False)
        solution = self.vectors
        self.vectors = old
        self.derVectors = oldDer
        return solution

    def ftcs(self, f, deriv, c2, t0, t_min, t_max, dt, x_min, x_max, dx, skip):
        t = t0
        count = 0
        inc = self.incrementFactor
        while t <= t_max:
            i = 0
            x = x_min
            stepDer = []
            step = []
            while x <= x_max:
                if t == t0:
                    der = deriv(x)
                    pos = f(x)
                elif x == x_min:
                    der = self.derVectors[-1][0][-1]
                    pos = self.vectors[-1][0][-1] + dt * der
                elif x + dx > x_max:
                    der = self.derVectors[-1][-1][-1]
                    pos = self.vectors[-1][-1][-1] + dt * der
                else:
                    prevDer = self.derVectors[-1][i][-1]
                    deltaDer = (-2 * self.vectors[-1][i][-1] +
                                self.vectors[-1][i - 1][-1] +
                                self.vectors[-1][i + 1][-1]) / dx**2
                    der = self.derVectors[-1][i][-1] + deltaDer * c2 * dt
                    pos = self.vectors[-1][i][-1] + dt * (der + prevDer) / 2
                stepDer.append((x, t, der))
                step.append((x, t, pos))
                i += 1
                x += dx
            if skip and count % inc != 1 and len(self.vectors) > 0:
                self.derVectors.pop()
                self.vectors.pop()
            self.derVectors.append(stepDer)
            self.vectors.append(step)
            t += dt
            count += 1

    def btcs(self, f, deriv, c2, t0, t_min, t_max, dt, x_min, x_max, dx, skip):
        t = t0
        count = 0
        inc = self.incrementFactor
        while t >= t_min:
            i = 0
            x = x_min
            stepDer = []
            step = []
            while x <= x_max:
                if t == t0:
                    der = deriv(x)
                    pos = f(x)
                elif x == x_min:
                    der = self.derVectors[0][0][-1]
                    pos = self.vectors[0][0][-1] + dt * der
                elif x + dx > x_max:
                    der = self.derVectors[0][-1][-1]
                    pos = self.vectors[0][-1][-1] + dt * der
                else:
                    prevDer = self.derVectors[0][i][-1]
                    deltaDer = (-2 * self.vectors[0][i][-1] +
                                self.vectors[0][i - 1][-1] +
                                self.vectors[0][i + 1][-1]) / dx**2
                    der = self.derVectors[0][i][-1] - deltaDer * c2 * dt
                    pos = self.vectors[0][i][-1] - dt * (der + prevDer) / 2
                stepDer.append((x, t, der))
                step.append((x, t, pos))
                i += 1
                x += dx
            if skip and count % inc != 1 and len(self.vectors) > 0:
                self.derVectors.pop(0)
                self.vectors.pop(0)
            self.derVectors.insert(0, stepDer)
            self.vectors.insert(0, step)
            t -= dt
            count += 1

####################################
# physics module
####################################

class ForceField(VectorField3D):

    def __init__(self, fxn):
        super().__init__(fxn)
        try: self.i, self.j, self.k = tuple(self.fxn.split(","))
        except: pass

    def strFunc(self):
        return "x, y, z, t"

    def evaluate(self, x, y, z, t):
        try: return self.evalWrapper(self.f(x = x, y = y, z = z, t = t))
        except: self.isValid = False

class Particle(ThreeD):

    def __init__(self,
                 ef = ForceField("0,0,0"),
                 mf = ForceField("0,0,0"),
                 gf = ForceField("0,0,0"),
                 ff = ForceField("0,0,0"),
                 charge = 1e-9, mass = 1,
                 x0 = 0, y0 = 0, z0 = 0,
                 xPrime0 = 0, yPrime0 = 0, zPrime0 = 0,
                 t_min = 0, t_max = 10):
        super().__init__("")
        self.electricField = ef # units N/C
        self.magneticField = mf # units T
        self.gravitationalField = gf # units N/kg
        self.forceField = ff # units N
        self.charge = charge # 1 nanocoulomb
        self.mass = mass # 1 kg
        self.x0, self.y0, self.z0 = x0, y0, z0 # start at origin, start at rest
        self.xPrime0, self.yPrime0, self.zPrime0 = xPrime0, yPrime0, zPrime0
        self.minT = t_min
        self.maxT = t_max # timescale: 0 to 10 seconds
        self.solvable = True

    def strFunc(self):
        return ""

    def generateVectors(self, data):
        # simulate the motion of the particle
        if not self.solvable:
            self.vectors = [[]]
            return
        step = (self.maxT - self.minT)/data.incrementCurve
        self.vectors = self.solve(step)

    def solve(self, dt):
        vectors = []
        step = dt
        t = self.minT
        x, y, z = self.x0, self.y0, self.z0
        xPrime, yPrime, zPrime = self.xPrime0, self.yPrime0, self.zPrime0
        q, m = self.charge, self.mass
        row = [(x, y, z)]
        while t <= self.maxT: # increase time in small steps
            Dx, Dy, Dz, DDx, DDy, DDz = self.RungeKutta4(
                q, m, t, x, y, z, xPrime, yPrime, zPrime, step)
            t += step # change time, position, velocity
            x += Dx; y += Dy; z += Dz
            xPrime += DDx; yPrime += DDy; zPrime += DDz
            vector = x, y, z
            row.append(vector)
        row.pop()
        vectors.append(row)
        return vectors

    def RungeKutta4(self, q, m, t, x, y, z, xPrime, yPrime, zPrime, step):
        # calculate the change in position, velocity using Runge Kutta 4
        Dx1, Dy1, Dz1, DDx1, DDy1, DDz1 = self.RKStep(q, m, t,
                x, y, z, xPrime, yPrime, zPrime, step)
        Dx2, Dy2, Dz2, DDx2, DDy2, DDz2 = self.RKStep(q, m, t + step/2,
                x + Dx1/2, y + Dy1/2, z + Dz1/2,
                xPrime + DDx1/2, yPrime + DDy1/2, zPrime + DDz1/2, step)
        Dx3, Dy3, Dz3, DDx3, DDy3, DDz3 = self.RKStep(q, m, t + step/2,
                x + Dx2/2, y + Dy2/2, z + Dz2/2,
                xPrime + DDx2/2, yPrime + DDy2/2, zPrime + DDz2/2, step)
        Dx4, Dy4, Dz4, DDx4, DDy4, DDz4 = self.RKStep(q, m, t + step,
                x + Dx3, y + Dy3, z + Dz3,
                xPrime + DDx3, yPrime + DDy3, zPrime + DDz3, step)
        Dx = (Dx1 + 2*Dx2 + 2*Dx3 + Dx4)/6
        Dy = (Dy1 + 2*Dy2 + 2*Dy3 + Dy4)/6
        Dz = (Dz1 + 2*Dz2 + 2*Dz3 + Dz4)/6
        DDx = (DDx1 + 2*DDx2 + 2*DDx3 + DDx4)/6
        DDy = (DDy1 + 2*DDy2 + 2*DDy3 + DDy4)/6
        DDz = (DDz1 + 2*DDz2 + 2*DDz3 + DDz4)/6
        return Dx, Dy, Dz, DDx, DDy, DDz

    def RKStep(self, q, m, t, x, y, z, xPrime, yPrime, zPrime, step):
        # compute calculations involved in each Runge Kutta step
        Dx = xPrime*step
        Dy = yPrime*step
        Dz = zPrime*step
        F = self.forceField.evaluate(x = x, y = y, z = z, t = t)
        g = self.gravitationalField.evaluate(x = x, y = y, z = z, t = t)
        E = self.electricField.evaluate(x = x, y = y, z = z, t = t)
        B = self.magneticField.evaluate(x = x, y = y, z = z, t = t)
        DDx = (q/m*(E[0] + yPrime*B[2] - zPrime*B[1]) + F[0]/m + g[0])*step
        DDy = (q/m*(E[1] + zPrime*B[0] - xPrime*B[2]) + F[1]/m + g[1])*step
        DDz = (q/m*(E[2] + xPrime*B[1] - yPrime*B[0]) + F[2]/m + g[2])*step
        return Dx, Dy, Dz, DDx, DDy, DDz

####################################
# dialog boxes
####################################

class Button(object): # buttons can be clicked on by user

    def __init__(self, text, cx, cy, width, height):
        self.text = text
        self.cx = cx
        self.cy = cy
        self.width = width
        self.height = height

    def clickedOn(self, x, y):
        # did user click on button?
        return (x < self.cx + self.width//2 and x > self.cx - self.width//2 and
                y < self.cy + self.height//2 and y > self.cy - self.width//2)

    def draw(self, canvas, data):
        # draw button on canvas
        canvas.create_rectangle(self.cx - self.width//2,
                                self.cy - self.height//2,
                                self.cx + self.width//2,
                                self.cy + self.height//2, fill ="deep sky blue")
        canvas.create_text(self.cx, self.cy, text = self.text,
                           justify = CENTER)

class DialogBar(object): # dialog bars that users can type into to input data

    def __init__(self, text, cx, cy, width, height, color = "black"):
        self.text = text
        self.cx = cx
        self.cy = cy
        self.width = width
        self.height = height
        self.color = color
        self.allowEdit = False # is this bar selected

    def clickedOn(self, x, y):
        # did user click on bar? if yes - allow user to edit
        ans = (x < self.cx + self.width//2 and x > self.cx - self.width//2 and
               y < self.cy + self.height//2 and y > self.cy - self.height//2)
        if ans: self.allowEdit = True
        else: self.allowEdit = False
        return ans

    def edit(self, key):
        # edits the text given the keystroke
        if not self.allowEdit: return None
        if key == "BackSpace" and len(self.text) > 0:
            self.text = self.text[:-1]
        elif key == "asterisk": self.text += "*"
        elif key == "slash": self.text += "/"
        elif key == "plus": self.text += "+"
        elif key == "minus": self.text += "-"
        elif key == "asciicircum": self.text += "**"
        elif key == "parenleft": self.text +="("
        elif key == "period": self.text += "."
        elif key == "parenright": self.text +=")"
        elif key == "comma": self.text += ","
        elif key == "space": self.text += " "
        elif key == "quotedbl": self.text += '"'
        elif key == "quoteright": self.text += "'"
        elif key == "less": self.text += "<"
        elif key == "greater": self.text += ">"
        elif len(key) == 1: self.text += key

    def __repr__(self): return self.text

    def draw(self, canvas, data):
        # draws dialog bar on canvas
        canvas.create_rectangle(self.cx - self.width//2,
                                self.cy - self.height//2,
                                self.cx + self.width//2,
                                self.cy + self.height//2,
                                fill = "PaleTurquoise1")
        text = (self.text + "|") if self.allowEdit else self.text
        canvas.create_text(self.cx - self.width//2 + data.miniMargin,
                           self.cy, text = text, anchor = W, fill = self.color,
                           font = "MS 10")

class DialogBox(object): # dialog box allows user to input pieces of info

    def __init__(self, title, inputs, data):
        self.title = title
        self.inputs = inputs
        self.outputs = []
        self.edit = 0
        factor = 0.66
        self.cx = data.width//2
        self.cy = data.height//2
        self.width = data.width*factor
        height = 25
        self.nudge = 50
        self.height = (len(self.inputs) + 2)*(height + data.margin)
        self.initBarsAndButtons(data, factor, height)

    def initBarsAndButtons(self, data, factor, height):
        # adds dialog bars into box given inputs
        for i in range(len(self.inputs)):
            self.outputs.append(DialogBar("", self.cx + self.width//2//2
                - self.nudge//2,
                self.cy-self.height//2+(height+data.margin)*(i + 1)+data.margin,
                self.width//2 - data.margin*2 + self.nudge, height))
        self.outputs[0].allowEdit = True
        buttonWidth = 60
        # adds buttons "OK" and "Cancel"
        self.buttons = [Button("OK", self.cx - buttonWidth//2 - data.margin,
            self.cy + self.height//2 - data.margin - height//2, buttonWidth,
            height), Button("Cancel", self.cx + buttonWidth//2 + data.margin,
            self.cy + self.height//2 - data.margin - height//2, buttonWidth,
            height)]
            
    def clickedOn(self, x, y):
        # did user click on anything in the box?
        for bar in self.outputs:
            if bar.clickedOn(x, y):
                self.edit = self.outputs.index(bar)
                bar.allowEdit = True
            else: bar.allowEdit = False

    def draw(self, canvas, data):
        # draws dialog box onto canvas
        factor = 4
        height = 25
        margin = data.margin
        canvas.create_rectangle(self.cx - self.width//2,
            self.cy - self.height//2, self.cx + self.width//2,
            self.cy + self.height//2, fill = "light blue")
        canvas.create_rectangle(self.cx - self.width//2,
            self.cy - self.height//2, self.cx + self.width//2,
            self.cy + - self.height//2 + height, fill = "sky blue")
        canvas.create_text(self.cx, self.cy - self.height//2,
                           text = self.title, justify = CENTER,
                           anchor = N, font = "Arial 12")
        for i in range(len(self.inputs)):
            canvas.create_text(self.cx - self.nudge,
                self.cy - self.height//2 + (height + margin)*(i + 1) + margin,
                text = self.inputs[i], font = "Arial 10",
                justify = RIGHT, anchor = E)
            self.outputs[i].draw(canvas, data)
        for button in self.buttons: button.draw(canvas, data)

    def editOutputs(self, key):
        # if you press tab, you go to the next line. Or else, edit the bar
        self.outputs[self.edit].edit(key)
        if key == "Tab" and self.edit < len(self.inputs) - 1:
            self.outputs[self.edit].allowEdit = False
            self.edit += 1
            self.outputs[self.edit].allowEdit = True

class DialogFunctions(object): # displays a function similar to dialog box

    def __init__(self, data):
        self.edit = 0
        self.left = data.margin
        self.bottom = data.height - data.margin
        factor = 0.67
        self.width = data.width*factor
        height = 30
        self.height = (len(data.functions[data.mode])*(height + data.margin)
                       + data.margin)
        self.cx = self.left + self.width//2
        self.cy = self.bottom - self.height//2
        self.reinit(data)

    def reinit(self, data):
        # recalculate width, height
        height = 30
        self.height = (len(data.functions[data.mode])*(height + data.margin)
                       + data.margin)
        self.cy = self.bottom - self.height//2
        self.outputs = []
        for i in range(len(data.functions[data.mode])):
            # make a dialog bar for every function
            graph = data.functions[data.mode][i]
            modeText = getModeText(graph)
            self.outputs.append(DialogBar(modeText + " = " + graph.fxn,
                self.cx, self.cy - self.height//2 + i*(height + data.margin)
                + height//2 + data.margin, self.width - data.margin*2, height,
                color = graph.color))
        # change which bar is selected for editing
        try: self.outputs[self.edit].allowEdit = True
        except:
            self.edit = len(data.functions[data.mode]) - 1
            self.outputs[self.edit].allowEdit = True
        if data.mode == "Differential Equations":
            reinitDifferentialEquations(data)

    def regenerateText(self, data):
        try:
            i = self.edit
            graph = data.functions[data.mode][i]
            self.outputs[i].text = getModeText(graph) + " = " + graph.fxn
        except: self.reinit(data)

    def clickedOn(self, x, y):
        # what did the user click on in the box?
        for bar in self.outputs:
            if bar.clickedOn(x, y):
                self.edit = self.outputs.index(bar)
                bar.allowEdit = True
            else: bar.allowEdit = False

    def draw(self, canvas, data):
        # draw functions box
        factor = 4
        margin = data.margin
        canvas.create_rectangle(self.left, self.bottom - self.height,
                                self.left + self.width, self.bottom,
                                fill = "light blue")
        for output in self.outputs:
            output.draw(canvas, data)

def reinitDifferentialEquations(data):
    try: graph = data.functions[data.mode][data.dialogFunctions.edit]
    except:
        data.dialogFunctions.edit -= 1
        return reinitDifferentialEquations(data)
    if type(graph) == DiffEq: data.order = "First Order"
    elif type(graph) == SecondOrderDE: data.order = "Second Order"
    getMenuDiffEq(data)

def getModeText(graph):
    # mode text displays what the user should input
    if type(graph) == Cartesian3D: modeText = "z = f(x, y)"
    elif type(graph) == CylindricalZDependent: modeText = "z = f(r, theta)"
    elif type(graph) == CylindricalRDependent: modeText = "r = f(z, theta)"
    elif type(graph) == Spherical: modeText = "rho = f(theta, phi)"
    elif type(graph) == Parametric3D:
        modeText = "r(t) = x(t), y(t), z(t), t_min, t_max"
    elif type(graph) == VectorField3D:
        modeText = "F(x, y, z) = P(x,y,z), Q(x,y,z), R(x,y,z)"
    elif type(graph) == Cartesian2DyDep: modeText = "y = f(x)"
    elif type(graph) == Cartesian2DxDep: modeText = "x = f(y)"
    elif type(graph) == Polar: modeText = "r = f(theta)"
    elif type(graph) == Parametric2D:
        modeText = "r(t) = x(t), y(t), t_min, t_max"
    elif type(graph) == VectorField2D: modeText = "F(x, y) = P(x,y), Q(x,y)"
    elif type(graph) == DiffEq: modeText = "y' = f(x, y) [, x_0, y_0]"
    elif type(graph) == SecondOrderDE:
        modeText = "y\" = f(x, y, Dy), x_0, y_0, y'_0"
    else: modeText = ""
    return modeText

####################################
# graphics projection
####################################

def vectorToPoint(vector, data):
    # these few lines of code took me ages to map a 3D vector onto a plane!
    try:
        y = (data.originY - vector[2]*data.scaleZ*math.sin(data.phi) +
             vector[1]*data.scaleY*math.cos(data.phi)*math.sin(data.theta) +
             vector[0]*data.scaleX*math.cos(data.theta)*math.cos(data.phi))
        x = (data.originX + vector[1]*data.scaleY*math.cos(data.theta) -
             vector[0]*data.scaleX*math.sin(data.theta))
        return x, y
    except: return None, None

####################################
# execute buttons
####################################

def executeMenu(data):
    # execute whatever is selected in the menu by the user
    if data.menuSelected[0] == -1 or data.menuSelected[1] == -1: return
    option = data.dropdown[data.menuSelected[0]][data.menuSelected[1]]
    executeMenu1(data, option)

def executeMenu1(data, option):
    # executes options in the menu
    if option == "Clear all":
        if data.mode == "MATHLAB Statistics": initStats(data); return
        data.functions[data.mode] = []
        if data.mode == "MATHLAB 3D": addNewFunction3D(data, data.mode)
        elif data.mode == "MATHLAB 2D": addNewFunction2D(data, data.mode)
        data.exeMessage = "All functions cleared!"
    elif option == "Remove current" or option == "Remove plot":
        data.functions[data.mode].pop(data.dialogFunctions.edit)
        if data.functions[data.mode] == []:
            if data.mode == "MATHLAB 3D": addNewFunction3D(data, data.mode)
            elif data.mode == "MATHLAB 2D": addNewFunction2D(data, data.mode)
            elif data.mode == "Differential Equations":
                data.functions[data.mode] = [SecondOrderDE("")]
        data.exeMessage = "Function removed!"
        data.dialogFunctions.reinit(data)
    elif option == "New field": executeGraph(data, option)
    else: executeMenu2(data, option)

def executeMenu2(data, option):
    # executes options in the menu
    if option == "Clear field":
        for graph in data.functions[data.mode]:
            if type(graph) == DiffEq:
                graph.fxn = ""; graph.vectors = []
                data.dialogFunctions.regenerateText(data)
        data.exeMessage = "Slope field cleared!"
    elif option == "Add plot":
        graph = data.functions[data.mode][data.dialogFunctions.edit]
        if type(graph) != SecondOrderDE or graph.fxn != "":
            data.functions[data.mode].append(SecondOrderDE(""))
            data.dialogFunctions.edit = len(data.functions[data.mode]) - 1
        executeGraph(data, option)
    elif option.endswith("Order"):
        if option.startswith("First"): executeFirstOrder(data)
        elif option.startswith("Second"): executeSecondOrder(data)
        data.dialogFunctions.reinit(data)
        data.menu[-1] = data.order = option
    else: executeMenu3(data, option)

def executeMenu3(data, option):
    # executes options in the menu
    if option in data.coordinates:
        data.coordinate = option; data.dialogBox = None
        executeGraph(data, option)
        if data.mode == "MATHLAB 3D":
            addNewFunction3D(data, option)
        elif data.mode == "MATHLAB 2D":
            addNewFunction2D(data, option)
    elif option in data.modes:
        data.mode = option
        if data.mode == "MATHLAB 3D": init3D(data)
        elif data.mode == "MATHLAB 2D": init2D(data)
        elif data.mode == "MATHLAB Calculator": initCalc(data)
        elif data.mode == "MATHLAB Statistics": reinitStats(data)
        elif data.mode == "Differential Equations": initDiffEq(data)
        elif data.mode == "PHYSLAB": initPhys(data)
    else: executeMenu4(data, option)

def executeMenu4(data, option):
    # executes options in the menu
    if option == "Edit axes": changeAxes(data, option)
    elif data.mode == "MATHLAB Calculator": executeCalculator(data, option)
    elif data.mode == "MATHLAB Statistics": executeStatistics(data, option)
    elif data.mode == "PHYSLAB": executeMenuPhys(data, option)

def executeFirstOrder(data):
    data.dropdown[-1] = ["New field", "Clear field", "Edit axes"]
    addField = True
    for graph in data.functions[data.mode]:
        if type(graph) == DiffEq:
            addField = False
            data.dialogFunctions.edit = data.functions[data.mode].index(graph)
            break
    if addField:
        data.functions[data.mode] += [DiffEq("")]
        data.dialogFunctions.edit = len(data.functions[data.mode]) - 1

def executeSecondOrder(data):
    data.dropdown[-1] = ["Add plot", "Remove plot", "Edit axes"]
    if data.functions[data.mode][data.dialogFunctions.edit].fxn == "":
        data.functions[data.mode][
            data.dialogFunctions.edit] = SecondOrderDE("")
    elif type(data.functions[data.mode][-1]) != SecondOrderDE:
        data.functions[data.mode] += [SecondOrderDE("")]
    data.dialogFunctions.edit = len(data.functions[data.mode]) - 1

def executeMenuPhys(data, option):
    data.exeMessage = "Type to enter a value. Hit return when done."
    particle = data.functions[data.mode][0]
    if option.endswith("Force") or option.endswith("Field"):
        executeMenuForceField(data, option, particle)
    elif option == "Edit Particle Properties":
        executeMenuParticleProperties(data, option, particle)
    elif option == "Edit Initial Conditions":
        executeMenuInitialConditions(data, option, particle)
    elif option == "Edit Timescale":
        executeMenuTimescale(data, option, particle)

def executeMenuForceField(data, option, particle):
    units, field = getUnitsField(particle, option)
    inputs = ["(%s) Enter i = P(x, y, z, t):" % units,
              "(%s) Enter j = Q(x, y, z, t):" % units,
              "(%s) Enter k = R(x, y, z, t):" % units]
    data.dialogBox = DialogBox(option, inputs, data)
    i, j, k = range(len(inputs))
    try:
        data.dialogBox.outputs[i].text = field.i
        data.dialogBox.outputs[j].text = field.j
        data.dialogBox.outputs[k].text = field.k
    except: pass

def getUnitsField(particle, option):
    if option == "General Force":
        units = "F"
        field = particle.forceField
    elif option == "Gravitational Field":
        units = "N/kg"
        field = particle.gravitationalField
    elif option == "Electric Field":
        units = "N/C"
        field = particle.electricField
    elif option == "Magnetic Field":
        units = "T"
        field = particle.magneticField
    return units, field

def executeMenuParticleProperties(data, option, particle):
    inputs = ["Particle mass (kg):", "Particle charge (C):"]
    data.dialogBox = DialogBox(option, inputs, data)
    mass, charge = range(len(inputs))
    try:
        data.dialogBox.outputs[mass].text = str(particle.mass)
        data.dialogBox.outputs[charge].text = str(particle.charge)
    except: pass

def executeMenuInitialConditions(data, option, particle):
    inputs = ["(m) Enter initial x:", "(m) Enter initial y:",
              "(m) Enter initial z:", "(m/s) Enter initial x':",
              "(m/s) Enter initial y':", "(m/s) Enter initial z':"]
    data.dialogBox = DialogBox(option, inputs, data)
    x0, y0, z0, xPrime0, yPrime0, zPrime0 = range(len(inputs))
    try:
        data.dialogBox.outputs[x0].text = str(particle.x0)
        data.dialogBox.outputs[y0].text = str(particle.y0)
        data.dialogBox.outputs[z0].text = str(particle.z0)
        data.dialogBox.outputs[xPrime0].text = str(particle.xPrime0)
        data.dialogBox.outputs[yPrime0].text = str(particle.yPrime0)
        data.dialogBox.outputs[zPrime0].text = str(particle.zPrime0)
    except: pass

def executeMenuTimescale(data, option, particle):
    inputs = ["Start time (s):", "End time (s):"]
    data.dialogBox = DialogBox(option, inputs, data)
    minT, maxT = range(len(inputs))
    try:
        data.dialogBox.outputs[minT].text = str(particle.minT)
        data.dialogBox.outputs[maxT].text = str(particle.maxT)
    except: pass

def executeGraph(data, option):
    # open a dialog box if specific coordinate systems were picked
    data.exeMessage = "Type to enter a function."
    if data.mode == "PHYSLAB": executeGraphPhys(data, option); return
    elif option == "Parametric Curves":
        inputs = ["Enter x(t):", "Enter y(t):", "Enter z(t):",
                  "Enter minimum t:", "Enter maximum t:"]
    elif option == "Vector Field":
        inputs = ["Enter i = P(x, y, z):", "Enter j = Q(x, y, z):",
                  "Enter k = R(x, y, z):"]
    elif option == "Parametric":
        inputs = ["Enter x(t):", "Enter y(t):",
                  "Enter minimum t:", "Enter maximum t:"]
    elif option == "Vector Field 2D":
        inputs = ["Enter i = P(x, y):", "Enter j = Q(x, y):"]
    elif option == "New field":
        inputs = ["Enter y' = f(x, y):", "(Optional) initial x:",
                  "(Optional) initial y:"]
    elif option == "Add plot": inputs = ["(y' := Dy) y\" = f(x, y, Dy):",
        "Enter initial x:", "Enter initial y:", "Enter initial y\':"]
    else: return
    data.dialogBox = DialogBox(option, inputs, data)

def changeAxes(data, option):
    # open a dialog box to allow user to change axes
    data.exeMessage = "Type to enter a function."
    inputs = ["Enter minimum x:", "Enter maximum x:",
              "Enter minimum y:", "Enter maximum y:"]
    if data.mode == "MATHLAB 3D" or data.mode == "PHYSLAB":
        inputs.append("Enter minimum z:")
        inputs.append("Enter maximum z:")
    data.dialogBox = DialogBox(option, inputs, data)

def executeCalculator(data, option):
    # open dialog box for whatever is selected by the user in the menu
    data.exeMessage = "Type to enter a value. Hit return when done."
    inputs = None
    if option == "Derivative at point":
        inputs = ["Enter a function f(x):", "Enter a value for x:"]
    elif option == "Second derivative":
        inputs = ["Enter a function f(x):", "Enter a value for x:"]
    elif option == "Definite Integral":
        inputs = ["Enter a function f(x):", "Enter lower bound:",
                  "Enter upper bound:"]
        data.exeMessage += " Note that this function assumes continuity."
    elif option == "Limit":
        inputs = ["Enter a function f(x):", "Enter a value for x:"]
    else: inputs = executeCalculator1(data, option, inputs)
    data.dialogBox = DialogBox(option, inputs, data)

def executeCalculator1(data, option, inputs):
    # open dialog box for whatever is selected by the user in the menu
    if option == "Definite Double Integral":
        inputs = ["Enter a function f(x, y):", "(Optional) Constraint 1:",
                  "(Optional) Constraint 2:", "(Optional) Constraint 3:",
                  "(Optional) Constraint 4:", "(Optional) Constraint 5:",
                  "Enter minimum x:", "Enter maximum x:", "Enter minimum y:",
                  "Enter maximum y:"]
    elif option == "Definite Triple Integral":
        inputs = ["Enter a function f(x, y, z):", "(Optional) Constraint 1:",
                  "(Optional) Constraint 2:", "(Optional) Constraint 3:",
                  "(Optional) Constraint 4:", "(Optional) Constraint 5:",
                  "Enter minimum x:", "Enter maximum x:", "Enter minimum y:",
                  "Enter maximum y:", "Enter minimum z:", "Enter maximum z:"]
    elif option == "Solve for 0":
        inputs = ["Enter a function f(x):", "Enter an initial guess:"]
    elif option == "Local Minimum":
        inputs = ["Enter a function f(x):", "Enter an x near minimum:"]
    else: inputs = executeCalculator2(data, option, inputs)
    return inputs

def executeCalculator2(data, option, inputs):
    # open dialog box for whatever is selected by the user in the menu
    if option == "Local Maximum":
        inputs = ["Enter a function f(x):", "Enter an x near maximum:"]
    elif option == "Local Minimum":
        inputs = ["Enter a function f(x):", "Enter an x near minimum:"]
    elif option == "Local Maximum 3D":
        inputs = ["Enter a function f(x, y):",
                  "Enter an x near maximum:", "Enter a y near maximum:"]
    elif option == "Local Minimum 3D":
        inputs = ["Enter a function f(x, y):",
                  "Enter an x near minimum:", "Enter a y near minimum:"]
    elif option == "Gradient at a point":
        inputs = ["Enter a function f(x, y, z):",
                  "Enter x:", "Enter y:", "Enter z:"]
    elif option == "Curl at a point":
        inputs = ["Enter i = P(x, y, z):", "Enter j = Q(x, y, z):",
                  "Enter k = R(x, y, z):", "Enter x:", "Enter y:", "Enter z:"]
    else: inputs = executeCalculator3(data, option, inputs)
    return inputs

def executeCalculator3(data, option, inputs):
    # open dialog box for whatever is selected by the user in the menu
    if option == "Divergence at a point":
        inputs = ["Enter i = P(x, y, z):", "Enter j = Q(x, y, z):",
                  "Enter k = R(x, y, z):", "Enter x:", "Enter y:", "Enter z:"]
    elif option == "Laplacian at a point":
        inputs = ["Enter a function f(x, y, z):",
                  "Enter x:", "Enter y:", "Enter z:"]
    elif option == "Normal Distribution":
        inputs = ["Mean:", "Standard Deviation:",
                  "Enter lower bound:", "Enter upper bound:"]
    elif option == "Inverse Normal":
        inputs = ["Probability less than value:","Mean:","Standard Deviation:"]
    elif option == "Student's t-distribution":
        inputs = ["Enter lower bound:","Enter upper bound:",
                  "Enter sample size:"]
        data.exeMessage += (
        " Note that small sample sizes and large intervals cause large errors.")
    else: inputs = executeCalculator4(data, option, inputs)
    return inputs

def executeCalculator4(data, option, inputs):
    # open dialog box for whatever is selected by the user in the menu
    if option == "Inverse Student's t":
        inputs = ["Probability less than value:","Enter sample size:"]
    elif option == "Exponential Distribution":
        inputs = ["Inverse of mean (lambda):","Enter lower bound:",
                  "Enter upper bound:"]
    elif option == "Gamma Distribution":
        inputs = ["Enter alpha:","Enter beta:","Enter lower bound:",
                  "Enter upper bound:"]
    elif option == "Beta Distribution":
        inputs = ["Enter alpha:","Enter beta:","Enter lower bound:",
                  "Enter upper bound:"]
    else: inputs = executeCalculator5(data, option, inputs)
    return inputs

def executeCalculator5(data, option, inputs):
    # open dialog box for whatever is selected by the user in the menu
    if option == "Binomial Distribution":
        inputs = ["Number of trials:","Probability of success:",
                  "Number of successes:"]
    elif option == "Negative Binomial":
        inputs = ["Number of successes:","Probability of success:",
                  "Number of trials:"]
    elif option == "Geometric Distribution":
        inputs = ["Probability of success:", "Number of trials:"]
    elif option == "Hypergeometric":
        inputs = ["Total population:", "Total possible successes:",
                  "Number of draws:", "Actual successes:"]
    elif option == "Poisson Distribution":
        inputs = ["Expected value (Average):", "Actual value:"]
    else: inputs = executeCalculator6(data, option, inputs)
    return inputs

def executeCalculator6(data, option, inputs):
    # open dialog box for whatever is selected by the user in the menu
    if option == "Factorial (n!)":
        inputs = ["Enter integer:"]
    elif option == "Permutation (nPr)":
        inputs = ["Enter n:", "Enter r:"]
    elif option == "Combination (nCr)":
        inputs = ["Enter n:", "Enter r:"]
    elif option == "Series":
        inputs = ["Enter expression f(i):","Enter start:","Enter end:"]
    elif option == "Sequence: Recursive":
        inputs = ["Enter recursive formula f(i):","Enter initial term:",
                  "Enter number of iterations:"]
        data.exeMessage += " 'i' refers to the previous term."
    elif option == "Sequence: Explicit":
        inputs = ["Enter explicit formula f(i):","Enter start:","Enter end:"]
        data.exeMessage += " 'i' refers to the index of the term."
    else: return
    return inputs

def executeStatistics(data, option):
    # open dialog box for whatever is selected by the user in the menu
    if option.endswith("Regression"):
        try:
            assert(len(data.points) > 1)
            points = [(point.x, point.y) for point in data.points]
            executeStatisticsRegression(data, option, points)
        except: pass
    elif option == "Statistics: Two-Variable":
        executeStatisticsTwoVar(data)
    elif option.startswith("Statistics"):
        executeStatisticsData(data, option)
    else: executeStatisticsDialog(data, option)

def executeStatisticsRegression(data, option, points):
    # get regression function
    if option.startswith("Linear"):
        m, b, r = linearRegression(points)
        if b >= 0: data.bestFit = Cartesian2DyDep("%.4f*x + %.4f" % (m, b))
        else: data.bestFit = Cartesian2DyDep("%.4f*x %.4f" % (m, b))
    elif option.startswith("Exponential"):
        a, b, r = exponentialRegression(points)
        data.bestFit = Cartesian2DyDep("%.4f*exp(%.4f*x)" % (a, b))
    elif option.startswith("Logarithmic"):
        a, b, r = logarithmicRegression(points)
        if b >= 0: data.bestFit = Cartesian2DyDep("%.4f*ln(x) + %.4f" % (a, b))
        else: data.bestFit = Cartesian2DyDep("%.4f*ln(x) %.4f" % (a, b))
    elif option.startswith("Power"):
        a, b, r = powerRegression(points)
        data.bestFit = Cartesian2DyDep("%.4f*x**%.4f" % (a, b))
    elif option.startswith("Polynomial"):
        data.tempPoints = points
        executeStatisticsDialog(data, option)
    data.bestFit.generateVectors(data)
    if not option.startswith("Polynomial"): data.correlationCoefficient = r

def executeStatisticsData(data, option):
    # get statistics for whatever column selected
    col = 0 if option.endswith("1") else 1
    text = []
    decimal = 5
    length = 11
    try:
        # build aggregate statistics for column
        column = buildColumn(data, col); assert(column != [])
        result = [len(column), sum(column), sumSquared(column),
            avg(column),standardDeviation(column),popStandardDeviation(column),
            min(column), firstQuartile(column), median(column),
            thirdQuartile(column), max(column)]
        data.inputMessage = """One Variable Statistics\n\nCount:\nSum:
Sum of Squares:\nAverage:\nSample Standard Deviation:
Population Standard Deviation:\nMinimum:\nQuartile 1:\nMedian:\nQuartile 3:
Maximum:"""
        for element in result:
            if type(element) == float: element = round(element, decimal)
            text.append(str(element))
        data.outputMessage = "Column %s\n\n" % option[-1] + "\n".join(text)
    except: pass

def buildColumn(data, col):
    # gets all numbers in a column
    column = []
    for row in range(data.rows):
        if type(data.evalTable[row][col]) in (int, float):
            column.append(data.evalTable[row][col])
    return column

def executeStatisticsTwoVar(data):
    decimal = 5
    decimalSmall = 8
    try:
        assert(len(data.points) > 1)
        points = [(point.x, point.y) for point in data.points]
        X = [x for (x, _) in points]
        Y = [y for (_, y) in points]
        XY = [x * y for (x, y) in points]
        avgX = avg(X)
        avgY = avg(Y)
        SSxx, SSxy, SSyy = 0, 0, 0
        sumXY = sum(XY)
        covariance = avg(XY) - avgX * avgY
        for (x, y) in points:
            SSxx += (x - avgX)**2
            SSxy += (x - avgX)*(y - avgY)
            SSyy += (y - avgY)**2
        m = SSxy/SSxx
        b = avgY - m*avgX
        r = math.sqrt(SSxy**2/SSxx/SSyy)
        result = [len(points), sumXY, covariance, SSxx, SSyy, SSxy, m, b, r]
        data.inputMessage = """Two Variable Statistics\n\nCount:
Sum of Products:\nCovariance:\nSS_xx:\nSS_yy:\nSS_xy:\nSlope:\ny-Intercept:
|Correlation Coefficient|:"""
        text = []
        for element in result:
            if type(element) == float:
                if abs(element) < 1: element = round(element, decimalSmall)
                else: element = round(element, decimal)
            text.append(str(element))
        data.outputMessage = "\n\n" + "\n".join(text)
    except: pass

def executeStatisticsDialog(data, option):
    # open dialog box for whatever is selected by the user in the menu
    if option == "z-interval: 1-Var Statistics":
        inputs = ["Enter observed mean:","Known standard deviation:",
                  "Enter sample size:","Enter confidence level:"]
    elif option.startswith("z-interval: Column"):
        inputs = ["Known standard deviation:","Enter confidence level:"]
    elif option == "z-interval: 2-Var Statistics":
        inputs = ["Enter observed mean 1:","Known standard deviation 1:",
                  "Enter sample size 1:","Enter observed mean 2:",
                  "Known standard deviation 2:",
                  "Enter sample size 2:","Enter confidence level:"]
    elif option == "z-interval: 2-Var Data":
        inputs = ["Known standard deviation 1:",
                  "Known standard deviation 2:", "Enter confidence level:"]
    elif option == "z-interval: 1-Var Prop":
        inputs = ["Enter observed successes:","Enter total trials:",
                  "Enter confidence level:"]
    elif option == "z-interval: 2-Var Prop":
        inputs = ["Enter observed successes 1:","Enter total trials 1:",
                  "Enter observed successes 2:","Enter total trials 2:",
                  "Enter confidence level:"]
    elif option == "t-interval: 1-Var Statistics":
        inputs = ["Enter observed mean:","Sample standard deviation:",
                  "Enter sample size:","Enter confidence level:"]
    elif option.startswith("t-interval: Column"):
        inputs = ["Enter confidence level:"]
    elif option == "t-interval: 2-Var Statistics":
        inputs = ["Enter observed mean 1:","Sample standard deviation 1:",
                  "Enter sample size 1:","Enter observed mean 2:",
                  "Sample standard deviation 2:",
                  "Enter sample size 2:","Enter confidence level:"]
    elif option == "t-interval: 2-Var Data":
        inputs = ["Enter confidence level:"]
    elif option == "Polynomial Regression":
        inputs = ["Enter order:"]
    else: return
    data.dialogBox = DialogBox(option, inputs, data)

def executeInput(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Derivative at point":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'derivative(lambda x: %s, %s)\n\n\n' % (f,x))
    elif data.inputTitle == "Second derivative":
        data.functions[data.mode][0] += (
            'secondDerivative(lambda x: %s, %s)\n\n\n' % tuple(data.inputs))
    elif data.inputTitle == "Definite Integral":
        f, a, b = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'integral(lambda x: %s, %s, %s)\n\n\n' % (f,a,b))
    elif data.inputTitle == "Limit":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'limit(lambda x: %s, %s)\n\n\n' % (f,x)
    else: executeInput1(data)
    evaluatePage(data)
    data.inputs = None
    data.exeMessage = ""

def executeInput1(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Solve for 0":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'zero(lambda x: %s, %s)\n\n\n' % (f,x)
    elif data.inputTitle == "Local Minimum":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'minimize(lambda x: %s, %s)\n\n\n' % (f,x))
    elif data.inputTitle == "Local Maximum":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'maximize(lambda x: %s, %s)\n\n\n' % (f,x))
    elif data.inputTitle == "Local Minimum 3D":
        f, x, y = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'minimize3D(lambda x, y, z: %s, (%s, %s))\n\n\n' % (f,x,y))
    elif data.inputTitle == "Local Maximum 3D":
        f, x, y = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'maximize3D(lambda x, y, z: %s, (%s, %s))\n\n\n' % (f,x,y))
    else: executeInput2(data)

def executeInput2(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Definite Double Integral":
        num = 4
        xMin, xMax, yMin, yMax = tuple(data.inputs[-num:])
        for i in range(num): data.inputs.pop()
        f = data.inputs.pop(0)
        constraints = ", ".join(["lambda x, y: %s" % cons
                                 for cons in data.inputs if cons.strip() != ""])
        data.functions[data.mode][0] += (
            'doubleIntegral(lambda x, y: %s, [%s], %s, %s, %s, %s)\n\n\n' % (
                f,constraints,xMin,xMax,yMin,yMax))
    else: executeInput3(data)

def executeInput3(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Definite Triple Integral":
        num = 6
        xMin, xMax, yMin, yMax, zMin, zMax = tuple(data.inputs[-num:])
        for i in range(num): data.inputs.pop()
        f = data.inputs.pop(0)
        constraints = ", ".join(["lambda x, y, z: %s" % cons
                                 for cons in data.inputs if cons.strip() != ""])
        data.functions[data.mode][0] += (
            'tripleIntegral(%s: %s, [%s], %s, %s, %s, %s, %s, %s)\n\n\n' % (
                "lambda x, y, z", f, constraints,
                xMin, xMax, yMin, yMax, zMin, zMax))
    elif data.inputTitle == "Gradient at a point":
        f, x, y, z = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'gradient3(lambda x, y, z: %s, %s, %s, %s)\n\n\n' % (f, x, y, z))
    elif data.inputTitle == "Curl at a point":
        P, Q, R, x, y, z = tuple(data.inputs)
        f = "lambda x, y, z: "
        data.functions[data.mode][0] += (
            'curl(%s%s, %s%s, %s%s, %s, %s, %s)\n\n\n' % (
                f, P, f, Q, f, R, x, y, z))
    else: executeInput4(data)

def executeInput4(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Divergence at a point":
        P, Q, R, x, y, z = tuple(data.inputs)
        f = "lambda x, y, z: "
        data.functions[data.mode][0] += (
            'divergence(%s%s, %s%s, %s%s, %s, %s, %s)\n\n\n' % (
                f, P, f, Q, f, R, x, y, z))
    elif data.inputTitle == "Laplacian at a point":
        f, x, y, z = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'laplacian(lambda x, y, z: %s, %s, %s, %s)\n\n\n' % (f, x, y, z))
    elif data.inputTitle == "Normal Distribution":
        mean, stdev, lower, upper = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'normalDistribution(%s, %s, %s, %s)\n\n\n'%(mean,stdev,lower,upper))
    elif data.inputTitle == "Inverse Normal":
        probability, mean, stdev = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'inverseNormal(%s, %s, %s)\n\n\n' % (probability, mean, stdev))
    else: executeInput5(data)

def executeInput5(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Student's t-distribution":
        lower, upper, sample = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'tDistribution(%s, %s, %s - 1)\n\n\n'%(lower,upper,sample))
    elif data.inputTitle == "Inverse Student's t":
        probability, sample = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'inverseT(%s, %s - 1)\n\n\n'%(probability,sample))
    elif data.inputTitle == "Exponential Distribution":
        lambd, lower, upper = tuple(data.inputs)
        data.functions[data.mode][0]+=(
            'exponentialDistribution(%s, %s, %s)\n\n\n' % (lambd, lower, upper))
    else: executeInput6(data)

def executeInput6(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Gamma Distribution":
        a, b, lower, upper = tuple(data.inputs)
        data.functions[data.mode][0]+=(
            'gammaDistribution(%s, %s, %s, %s)\n\n\n' % (a, b, lower, upper))
    elif data.inputTitle == "Beta Distribution":
        a, b, lower, upper = tuple(data.inputs)
        data.functions[data.mode][0]+=(
            'betaDistribution(%s, %s, %s, %s)\n\n\n' % (a, b, lower, upper))
    elif data.inputTitle == "Binomial Distribution":
        trials, probability, successes = tuple(data.inputs)
        data.functions[data.mode][0]+=('binomial(%s, %s, %s)\n\n\n'%
                                         (trials, probability, successes))
    elif data.inputTitle == "Negative Binomial":
        successes, probability, trials = tuple(data.inputs)
        data.functions[data.mode][0]+=('negativeBinomial(%s, %s, %s)\n\n\n'%
                                         (successes, probability, trials))
    else: executeInput7(data)
    
def executeInput7(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Hypergeometric":
        N, K, n, k = tuple(data.inputs)
        data.functions[data.mode][0]+=('hypergeometric(%s, %s, %s, %s)\n\n\n'%
                                       (N, K, n, k))
    elif data.inputTitle == "Geometric Distribution":
        probability, trials = tuple(data.inputs)
        data.functions[data.mode][0]+=('geometric(%s, %s)\n\n\n'%
                                         (probability, trials))
    elif data.inputTitle == "Poisson Distribution":
        expected, value = tuple(data.inputs)
        data.functions[data.mode][0]+=('poisson(%s, %s)\n\n\n'%
                                         (expected, value))
    elif data.inputTitle == "Factorial (n!)":
        n = data.inputs[0]
        data.functions[data.mode][0] += 'factorial(%s)\n\n\n' % n
    elif data.inputTitle == "Permutation (nPr)":
        n, r = tuple(data.inputs)
        data.functions[data.mode][0] += 'nPr(%s, %s)\n\n\n' % (n, r)
    else: executeInput8(data)

def executeInput8(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Combination (nCr)":
        n, r = tuple(data.inputs)
        data.functions[data.mode][0] += 'nCr(%s, %s)\n\n\n' % (n, r)
    elif data.inputTitle == "Series":
        expression, start, end = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'series(lambda i: %s, %s, %s)\n\n\n' % (expression, start, end))
    elif data.inputTitle == "Sequence: Explicit":
        expression, start, end = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'sequenceExplicit(lambda i: %s, %s, %s)\n\n\n' % (
                expression, start, end))
    elif data.inputTitle == "Sequence: Recursive":
        expression, initial, iterations = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'sequenceRecursive(lambda i: %s, %s, %s)\n\n\n' %
                (expression, initial, iterations))

def executeInputStats(data):
    # executes whatever the user inputs into the dialog box
    dec = 6
    if data.inputTitle == "z-interval: 1-Var Statistics":
        executeInputStatsZStats(data, dec)
    elif data.inputTitle.startswith("z-interval: Column"):
        executeInputStatsZData(data, dec)
    elif data.inputTitle == "z-interval: 2-Var Statistics":
        executeInputStatsZStatsTwo(data, dec)
    elif data.inputTitle == "z-interval: 2-Var Data":
        executeInputStatsZDataTwo(data, dec)
    elif data.inputTitle == "z-interval: 1-Var Prop":
        executeInputStatsZProp(data, dec)
    elif data.inputTitle == "z-interval: 2-Var Prop":
        executeInputStatsZPropTwo(data, dec)
    elif data.inputTitle == "t-interval: 1-Var Statistics":
        executeInputStatsTStats(data, dec)
    elif data.inputTitle.startswith("t-interval: Column"):
        executeInputStatsTData(data, dec)
    elif data.inputTitle == "t-interval: 2-Var Statistics":
        executeInputStatsTStatsTwo(data, dec)
    elif data.inputTitle == "t-interval: 2-Var Data":
        executeInputStatsTDataTwo(data, dec)
    elif data.inputTitle == "Polynomial Regression":
        degree = int(data.inputs[0])
        points = data.tempPoints
        coeffs, r = polynomialRegression(points, degree)
        fxn = generateFunctionPolynomialRegression(coeffs, degree)
        data.bestFit = Cartesian2DyDep(fxn)
        data.bestFit.generateVectors(data)
        data.correlationCoefficient = r

def generateFunctionPolynomialRegression(coeffs, degree):
    # generate a function line (string) from a list of coefficients
    fxn = ""
    for i in range(len(coeffs)):
        coeff = coeffs[i]
        if i == 0:
            if degree == 0: fxn += "%.4f" % coeffs[i]
            elif degree == 1: fxn += "%.4f*x" % coeffs[i]
            else: fxn += "%.4f*x**%d" % (coeffs[i], degree)
        elif degree - i == 0:
            if coeffs[i] >= 0: fxn += " + %.4f" % coeffs[i]
            else: fxn += " %.4f" % coeffs[i]
        elif degree - i == 1:
            if coeffs[i] >= 0: fxn += " + %.4f*x" % coeffs[i]
            else: fxn += " %.4f*x" % coeffs[i]
        else:
            if coeffs[i] >= 0: fxn += " + %.4f*x**%d" % (coeffs[i], degree - i)
            else: fxn += " %.4f*x**%d" % (coeffs[i], degree - i)
    return fxn

def executeInputStatsZStats(data, dec):
    # executes input for z-interval stats
    mean, stdev, n, confidence = tuple([eval(elem) for elem in data.inputs])
    z = inverseNormal((confidence + 1)/2, 0, 1)
    output = zIntervalStats(mean, stdev, n, confidence)
    data.inputMessage = """Confidence z-Interval\n\nObserved Mean:
Population Standard Deviation:\nSample Size:\nConfidence level:\n
z:\nConfidence Interval:"""
    data.outputMessage = "Statistics\n\n%s\n\n%f\n%s" % (
        "\n".join(data.inputs), z, str(output))

def executeInputStatsZStatsTwo(data, dec):
    # executes input for z-interval stats
    mean1, stdev1, n1, mean2, stdev2, n2, conf = tuple(
        [eval(elem) for elem in data.inputs])
    z = inverseNormal((conf + 1)/2, 0, 1)
    output = zIntervalStatsTwoSample(mean1, stdev1, n1, mean2, stdev2, n2, conf)
    data.inputMessage = """Confidence z-Interval\n\nObserved Mean 1:
Population Standard Deviation 1:\nSample Size 1:\nObserved Mean 2:
Population Standard Deviation 2:\nSample Size 2:\nConfidence level:\n
z:\nConfidence Interval:"""
    data.outputMessage = "Two-Variable Statistics\n\n%s\n\n%f\n%s" % (
        "\n".join(data.inputs), z, str(output))

def executeInputStatsZData(data, dec):
    # executes input for z-interval data
    col = int(data.inputTitle[-1]) - 1
    column = buildColumn(data, col)
    stdev, confidence = tuple([eval(elem) for elem in data.inputs])
    mean = avg(column)
    z = inverseNormal((confidence + 1)/2, 0, 1)
    output = zIntervalData(column, stdev, confidence)
    data.inputMessage = """Confidence z-Interval\n\nObserved Mean:
Population Standard Deviation:\nSample Size:\nConfidence level:\n
z:\nConfidence Interval:"""
    data.outputMessage = "Column %s\n\n%s\n%s\n%d\n%s\n\n%f\n%s" % (
        data.inputTitle[-1], str(round(mean, dec)), str(stdev),
        len(column), str(confidence), z, str(output))

def executeInputStatsZDataTwo(data, dec):
    col1 = buildColumn(data, 0)
    col2 = buildColumn(data, 1)
    mean1 = avg(col1)
    mean2 = avg(col2)
    n1 = len(col1)
    n2 = len(col2)
    stdev1, stdev2, confidence = tuple([eval(elem) for elem in data.inputs])
    z = inverseNormal((confidence + 1)/2, 0, 1)
    output = zIntervalDataTwoSample(col1, stdev1, col2, stdev2, confidence)
    data.inputMessage = """Confidence z-Interval\n\nObserved Mean 1:
Population Standard Deviation 1:\nSample Size 1:\nObserved Mean 2:
Population Standard Deviation 2:\nSample Size 2:\nConfidence level:\n
z:\nConfidence Interval:"""
    data.outputMessage = "%s\n\n%s\n%s\n%d\n%s\n%s\n%d\n%s\n\n%f\n%s" % (
        "Two-Sample Data",
        str(round(mean1, dec)), str(stdev1), len(col1),
        str(round(mean2, dec)), str(stdev2), len(col2),
        str(confidence), z, str(output))

def executeInputStatsZProp(data, dec):
    # executes input for z-interval proportions
    successes, trials, confidence = tuple([eval(elem) for elem in data.inputs])
    z = inverseNormal((confidence + 1)/2, 0, 1)
    output = zIntervalProportion(successes, trials, confidence)
    data.inputMessage = """Confidence z-Interval\n\nObserved Successes:
Total Trials:\nConfidence level:\n\nz:\nConfidence Interval:"""
    data.outputMessage = "1-Variable Proportion\n\n%s\n\n%f\n%s" % (
        "\n".join(data.inputs), z, str(output))

def executeInputStatsZPropTwo(data, dec):
    # executes input for z-interval proportions
    succ1, trials1, succ2, trials2, confidence = tuple(
        [eval(elem) for elem in data.inputs])
    z = inverseNormal((confidence + 1)/2, 0, 1)
    output = zIntervalProportionTwo(succ1, trials1, succ2, trials2, confidence)
    data.inputMessage = """Confidence z-Interval\n\nObserved Successes 1:
Total Trials 1:\nObserved Successes 2:\nTotal Trials 2:
Confidence level:\n\nz:\nConfidence Interval:"""
    data.outputMessage = "2-Variable Proportion\n\n%s\n\n%f\n%s" % (
        "\n".join(data.inputs), z, str(output))

def executeInputStatsTStats(data, dec):
    # executes input for t-interval stats
    mean, stdev, n, confidence = tuple([eval(elem) for elem in data.inputs])
    t = inverseT((confidence + 1)/2, n - 1)
    output = tIntervalStats(mean, stdev, n, confidence)
    data.inputMessage = """Confidence t-Interval\n\nObserved Mean:
Sample Standard Deviation:\nSample Size:\nConfidence level:\n
t:\nConfidence Interval:"""
    data.outputMessage = "Statistics\n\n%s\n\n%.3f\n%s" % (
        "\n".join(data.inputs), t, str(output))

def executeInputStatsTStatsTwo(data, dec):
    # executes input for t-interval stats - 2 vars
    mean1, stdev1, n1, mean2, stdev2, n2, conf = tuple(
        [eval(elem) for elem in data.inputs])
    s12n1 = stdev1**2/n1
    s22n2 = stdev2**2/n2
    df = (s12n1 + s22n2)**2 / ((s12n1)**2/(n1 - 1) + (s22n2)**2/(n2 - 1))
    t = inverseT((conf + 1)/2, df)
    output = tIntervalStatsTwoSample(mean1, stdev1, n1, mean2, stdev2, n2, conf)
    data.inputMessage = """Confidence t-Interval\n\nObserved Mean:
Sample Standard Deviation:\nSample Size:\nConfidence level:\n
Degrees of Freedom:\nt:\nConfidence Interval:"""
    data.outputMessage = "Statistics\n\n%s\n\n%f\n%.3f\n%s" % (
        "\n".join(data.inputs), df, t, str(output))

def executeInputStatsTData(data, dec):
    # executes input for t-interval data
    col = int(data.inputTitle[-1]) - 1
    column = buildColumn(data, col)
    confidence = eval(data.inputs[0])
    n = len(column)
    t = inverseT((confidence + 1)/2, n - 1)
    output = tIntervalData(column, confidence)
    data.inputMessage = """Confidence t-Interval\n\nObserved Mean:
Sample Standard Deviation:\nSample Size:\nConfidence level:\n
t:\nConfidence Interval:"""
    data.outputMessage = "Column %s\n\n%s\n%s\n%d\n%s\n\n%.3f\n%s" % (
        data.inputTitle[-1], str(avg(column)),
        str(round(standardDeviation(column), dec)), n, str(confidence), t,
        str(output))

def executeInputStatsTDataTwo(data, dec):
    # executes input for t-interval data - 2 vars
    col1 = buildColumn(data, 0)
    col2 = buildColumn(data, 1)
    mean1 = avg(col1)
    mean2 = avg(col2)
    n1 = len(col1)
    n2 = len(col2)
    stdev1 = standardDeviation(col1)
    stdev2 = standardDeviation(col2)
    confidence = eval(data.inputs[0])
    s12n1 = stdev1**2/n1
    s22n2 = stdev2**2/n2
    df = (s12n1 + s22n2)**2 / ((s12n1)**2/(n1 - 1) + (s22n2)**2/(n2 - 1))
    t = inverseT((confidence + 1)/2, df)
    output = tIntervalDataTwoSample(col1, col2, confidence)
    data.inputMessage = """Confidence t-Interval\n\nObserved Mean 1:
Sample Standard Deviation 1:\nSample Size 1:\nObserved Mean 2:
Sample Standard Deviation 2:\nSample Size 2:\nConfidence level:\n
Degrees of Freedom:\nt:\nConfidence Interval:"""
    data.outputMessage = "%s\n\n%s\n%f\n%d\n%s\n%f\n%d\n%s\n\n%f\n%.3f\n%s" % (
        "Two-Sample Data", str(mean1), stdev1, n1, str(mean2), stdev2, n2,
        str(confidence), df, t, str(output))

def executeInputGraph(data):
    # executes input for certain types of graphs
    mode = data.mode
    if data.inputTitle == "Parametric Curves":
        data.functions[mode][data.dialogFunctions.edit].fxn = (
            "%s, %s, %s, %s, %s" % tuple(data.inputs))
    elif data.inputTitle == "Vector Field":
        data.functions[mode][data.dialogFunctions.edit].fxn = (
            "%s, %s, %s" % tuple(data.inputs))
    elif data.inputTitle == "Parametric":
        data.functions[mode][data.dialogFunctions.edit].fxn = (
            "%s, %s, %s, %s" % tuple(data.inputs))
    else: executeInputGraph1(data)
    data.functions[mode][data.dialogFunctions.edit].generateLambda()
    data.functions[mode][data.dialogFunctions.edit].generateVectors(data)
    data.dialogFunctions.reinit(data)
    data.exeMessage = ""

def executeInputGraph1(data):
    # executes input for certain types of graphs
    mode = data.mode
    if data.inputTitle == "Vector Field 2D":
        P, Q = tuple(data.inputs)
        data.functions[mode][data.dialogFunctions.edit].fxn = ("%s, %s" % (P,Q))
    elif data.inputTitle == "New field":
        dy_dx, x_0, y_0 = tuple(data.inputs)
        if x_0 == "" or y_0 == "":
            data.functions[mode][data.dialogFunctions.edit].fxn = dy_dx
        else: data.functions[mode][data.dialogFunctions.edit].fxn = (
            "%s, %s, %s" % (dy_dx, x_0, y_0))
    elif data.inputTitle == "Add plot":
        fxn, x0, y0, yPrime0 = tuple(data.inputs)
        data.functions[mode][data.dialogFunctions.edit].fxn="%s, %s, %s, %s" % (
            (fxn, x0, y0, yPrime0))
    elif data.inputTitle == "Edit axes": editAxes(data)
    else: return
    data.exeMessage = ""

def executeInputPhys(data):
    # executes whatever the user inputs into the dialog box in physics module
    particle = data.functions[data.mode][0]
    if data.inputTitle == "General Force":
        particle.forceField = ForceField("%s,%s,%s" % tuple(data.inputs))
    elif data.inputTitle == "Gravitational Field":
        particle.gravitationalField=ForceField("%s,%s,%s" % tuple(data.inputs))
    elif data.inputTitle == "Electric Field":
        particle.electricField = ForceField("%s,%s,%s" % tuple(data.inputs))
    elif data.inputTitle == "Magnetic Field":
        particle.magneticField = ForceField("%s,%s,%s" % tuple(data.inputs))
    else: executeInputPhys1(data)
    data.inputs = None
    data.exeMessage = ""
    try:
        particle.solvable = True
        particle.generateVectors(data)
    except: particle.solvable = False

def executeInputPhys1(data):
    # executes whatever the user inputs into the dialog box in physics module
    particle = data.functions[data.mode][0]
    if data.inputTitle == "Edit Particle Properties":
        mass, charge = tuple(data.inputs)
        particle.mass = eval(mass)
        particle.charge = eval(charge)
    elif data.inputTitle == "Edit Initial Conditions":
        x0, y0, z0, xPrime0, yPrime0, zPrime0 = tuple(data.inputs)
        particle.x0 = eval(x0)
        particle.y0 = eval(y0)
        particle.z0 = eval(z0)
        particle.xPrime0 = eval(xPrime0)
        particle.yPrime0 = eval(yPrime0)
        particle.zPrime0 = eval(zPrime0)
    elif data.inputTitle == "Edit Timescale":
        minT, maxT = tuple(data.inputs)
        particle.minT = eval(minT)
        particle.maxT = eval(maxT)
    elif data.inputTitle == "Edit axes": editAxes(data)

def editAxes(data):
    # scales the axes given user input
    scale = data.scaleFactor
    try:
        for element in data.inputs:
            if type(eval(element)) not in (int, float): return
    except: return
    if data.mode == "MATHLAB 3D" or data.mode == "PHYSLAB":
        data.minX,data.maxX,data.minY,data.maxY,data.minZ,data.maxZ=tuple(
            [eval(element) for element in data.inputs])
        if data.minZ > data.maxZ: data.minZ,data.maxZ=data.maxZ,data.minZ
        data.scaleZ = data.width/max(abs(data.maxZ),abs(data.minZ))*scale
    else: data.minX,data.maxX,data.minY,data.maxY = tuple(
            [eval(element) for element in data.inputs])
    adjustAxes(data, scale)
    if data.mode in ("MATHLAB 2D", "Differential Equations"):
        try: recalculateAxes(data)
        except: pass
    for graph in data.functions[data.mode]: graph.generateVectors(data)

def adjustAxes(data, scale):
    # adjust max and min, recalculate the scale
    if data.minX > data.maxX: data.minX,data.maxX=data.maxX,data.minX
    if data.minY > data.maxY: data.minY,data.maxY=data.maxY,data.minY
    try:
        data.scaleX = data.width/max(abs(data.maxX),abs(data.minX))*scale
        data.scaleY = data.width/max(abs(data.maxY),abs(data.minY))*scale
    except: pass

def recalculateAxes(data):
    # reposition the center of the graph and recalculate the scale
    if data.minX > 0: data.minX = 0
    if data.maxX < 0: data.maxX = 0
    if data.minY > 0: data.minY = 0
    if data.maxY < 0: data.maxY = 0
    factor = 3/16
    startX, startY = data.width*factor, data.height*factor
    factor = 5/8
    maxWidth = data.width*factor
    maxHeight = data.height*factor
    data.originX = startX - data.minX/(data.maxX - data.minX)*maxWidth
    data.originY = startY + data.maxY/(data.maxY - data.minY)*maxHeight
    data.scaleX = maxWidth/(data.maxX - data.minX)
    data.scaleY = maxHeight/(data.maxY - data.minY)

def executeButtons(data, button):
    # executes any button clicked in the top right
    try: option = data.buttons[button]
    except: return
    if option == "H":
        if data.mode == "MATHLAB 3D": initGraph3D(data)
        elif data.mode == "MATHLAB 2D": initGraph2D(data)
        elif data.mode == "Differential Equations": initGraphDiffEq(data)
        elif data.mode == "PHYSLAB": initGraph3D(data)
    elif option == "+": zoomIn(data)
    elif option == "-": zoomOut(data)
    elif option == "?": data.helpPressed = not data.helpPressed
    if data.mode in ("MATHLAB Calculator", "MATHLAB Statistics"): return
    for graph in data.functions[data.mode]: graph.generateVectors(data)

def zoomIn(data):
    # zooms in on the graph
    data.maxX /= data.zoomFactor
    data.minX /= data.zoomFactor
    data.maxY /= data.zoomFactor
    data.minY /= data.zoomFactor
    data.scaleX *= data.zoomFactor
    data.scaleY *= data.zoomFactor
    if data.mode in ("MATHLAB 3D", "PHYSLAB"):
        data.maxZ /= data.zoomFactor
        data.minZ /= data.zoomFactor
        data.scaleZ *= data.zoomFactor

def zoomOut(data):
    # zooms out on the graph
    data.maxX *= data.zoomFactor
    data.minX *= data.zoomFactor
    data.maxY *= data.zoomFactor
    data.minY *= data.zoomFactor
    data.scaleX /= data.zoomFactor
    data.scaleY /= data.zoomFactor
    if data.mode in ("MATHLAB 3D", "PHYSLAB"):
        data.maxZ *= data.zoomFactor
        data.minZ *= data.zoomFactor
        data.scaleZ /= data.zoomFactor

def drawHelpMessage(data, canvas):
    # draws the help message onto canvas
    margin = 50
    canvas.create_rectangle(margin, margin, data.width - margin,
                            data.height - margin, fill = "sky blue")
    canvas.create_text(data.width//2, data.height//2, text = data.helpMessage,
                       justify = CENTER)

####################################
# coordinate system conversions
####################################

def sphericalToCartesian(rho, theta, phi):
    # convert spherical to cartesian coordinates
    return (rho*math.cos(theta)*math.sin(phi),
            rho*math.sin(theta)*math.sin(phi),
            rho*math.cos(phi))

def cylindricalToCartesian(r, theta, z):
    # convert cylindrical to cartesian coordinates
    return r*math.cos(theta), r*math.sin(theta), z

def polarToCartesian(r, theta):
    # convert polar to cartesian coordinates (2D)
    return r*math.cos(theta), r*math.sin(theta)

def checkBounds2D(prevPoint, point, data):
    # make sure both points are in bounds
    for i in prevPoint:
        if type(i) not in (int, float): return False
    for i in point:
        if type(i) not in (int, float): return False
    return (inRange(prevPoint[0], data.minX, data.maxX) and
            inRange(prevPoint[1], data.minY, data.maxY) and
            inRange(point[0], data.minX, data.maxX) and
            inRange(point[1], data.minY, data.maxY))

def checkBounds(prevVector, vector, data):
    # make sure both points are in bounds
    return (inRange(vector[0], data.minX, data.maxX) and
            inRange(prevVector[0], data.minX, data.maxX) and
            inRange(vector[1], data.minY, data.maxY) and
            inRange(prevVector[1], data.minY, data.maxY) and
            inRange(vector[2], data.minZ, data.maxZ) and
            inRange(prevVector[2], data.minZ, data.maxZ))

def checkDifferentiability(slope, prevSlope):
    # checks if something is differentiable (i.e. no big change in slope)
    smallNumber = 1
    if abs(slope) < smallNumber or abs(prevSlope) < smallNumber: return True
    return abs(slope - prevSlope)/abs(prevSlope) < smallNumber

####################################
# drawAxes
####################################

def drawAxes2D(canvas, data):
    # draws 2D axes
    canvas.create_line(data.originX+data.scaleX*data.minX, data.originY,
                       data.originX, data.originY)
    canvas.create_line(data.originX+data.scaleX*data.maxX, data.originY,
                       data.originX, data.originY)
    canvas.create_line(data.originX, data.originY-data.scaleY*data.maxY,
                       data.originX, data.originY)
    canvas.create_line(data.originX, data.originY-data.scaleY*data.minY,
                       data.originX, data.originY)
    canvas.create_text(
        data.originX+data.scaleX*data.maxX,data.originY+data.room,
                       text = "x = %.1f" % data.maxX, font = "Arial 7 bold")
    canvas.create_text(
        data.originX,data.originY-data.scaleY*data.maxY-data.room,
                       text = "y = %.1f" % data.maxY, font = "Arial 7 bold")

def drawAxes3D(canvas, data):
    # draws 3D axes
    axes = [(data.maxX,0,0),(0,data.maxY,0),(0,0,data.maxZ)]
    negAxes = [(data.minX,0,0),(0,data.minY,0),(0,0,data.minZ)]
    axesNames = ("x","y","z")
    maxValues = data.maxX, data.maxY, data.maxZ
    for i in range(len(axes)):
        vector = axes[i]
        point = vectorToPoint(vector,data)
        canvas.create_line(data.originX,data.originY,point)
        if data.mode == "PHYSLAB":
            text = "%s = %.3f m" % (axesNames[i], maxValues[i])
        else: text = "%s = %.1f" % (axesNames[i], maxValues[i])
        canvas.create_text(point[0]-data.room,point[1]-data.room,
                           text = text, font = "Arial 7 bold")
        vector = negAxes[i]
        point = vectorToPoint(vector,data)
        canvas.create_line(data.originX,data.originY,point)

####################################
# evaluate calculator page
####################################

import string

def evaluatePage(data):
    # evaluates all inputs on calculator page
    try: ans = eval(data.results[-1])
    except: pass
    skiplines = 3; lines = data.functions[data.mode][0].splitlines()
    for lineIndex in reversed(range(len(lines))):
        line = lines[lineIndex]
        try:
            if line == "": assert(lineIndex % skiplines != 0); continue
            answer = eval(line)
            assert(answer != None)
            assert(type(answer) in (int, float, str, tuple, list, bool))
            if type(answer) == float: answer = float("%.12g" % answer) # 12 sf
            try: generatePicture(data, line)
            except: pass
            data.results.append(str(answer))
            break
        except: data.results.append("Error!"); break

def generatePicture(data, line):
    if line.startswith("integral"):
        generateIntegralPicture(data, line)
    elif line.startswith("derivative"):
        generateDerivativePicture(data, line)
    elif line.startswith("zero"):
        generateZeroPicture(data, line)
    elif line.startswith("maximize(") or line.startswith("minimize("):
        generateMaxMinPicture(data, line)

class Struct(object): pass

def generateIntegralPicture(data, line):
    # generates a temporary struct to lay the groundwork to draw the integral
    line = line.strip().replace("integral(", "(")
    (fxn, start, end) = eval(line)
    tempData = Struct()
    startX, startY = data.width//2, data.height//2
    lengthFactor = 0.75
    maxWidth, maxHeight = startX*lengthFactor, startY*lengthFactor
    tempData.type = "integral"
    tempData.maxX = max(end, 0); tempData.minX = min(start, 0)
    tempData.start, tempData.end = start, end
    tempData.originX=startX-tempData.minX/(tempData.maxX-tempData.minX)*maxWidth
    tempData.increment = 100
    tempData.room = data.room
    graph = Cartesian2DyDep(fxn)
    graph.generateVectors(tempData)
    tempData.minY = min([y for x, y in graph.vectors] + [0])
    tempData.maxY = max([y for x, y in graph.vectors] + [0])
    assert(tempData.maxY != tempData.minY and tempData.minX != tempData.maxX)
    tempData.originY = (
        startY + tempData.minY/(tempData.maxY - tempData.minY)*maxHeight)
    tempData.scaleX = maxWidth/(tempData.maxX - tempData.minX)
    tempData.scaleY = maxHeight/(tempData.maxY - tempData.minY)
    tempData.graph = graph
    data.tempData = tempData

def drawIntegralPicture(canvas, data):
    # draws integral picture
    assert(data.type == "integral")
    drawAxes2D(canvas, data)
    data.graph.draw(canvas, data)
    vectors = data.graph.vectors
    prevPoints = None
    epsilon = 0.01
    for x, y in vectors:
        points = ((data.originX + x*data.scaleX, data.originY - y*data.scaleY),
                  (data.originX + x*data.scaleX, data.originY))
        if (prevPoints != None and
            x > data.start - epsilon and x < data.end + epsilon):
            canvas.create_polygon(prevPoints, points, fill = data.graph.color)
        prevPoints = tuple(reversed(list(points)))

def generateDerivativePicture(data, line):
    # generates a temporary struct to lay the groundwork to draw tangent line
    line = line.strip().replace("derivative(", "(")
    (fxn, x) = eval(line)
    tempData = Struct()
    tempData.fxn, tempData.x = fxn, x
    startX, startY = data.width//2, data.height//2
    lengthFactor = 0.75
    maxWidth, maxHeight = startX*lengthFactor, startY*lengthFactor
    tempData.type = "derivative"; tempData.factor = 1.2
    tempData.maxX = max(x*tempData.factor, 1/2)
    tempData.minX = min(x*tempData.factor, -1/2)
    tempData.start, tempData.end = tempData.minX, tempData.maxX
    tempData.originX=startX-tempData.minX/(tempData.maxX-tempData.minX)*maxWidth
    tempData.increment = 100
    tempData.room = data.room
    graph = Cartesian2DyDep(fxn); graph.generateVectors(tempData)
    tempData.minY = min([y for x, y in graph.vectors] + [0])
    tempData.maxY = max([y for x, y in graph.vectors] + [0])
    assert(tempData.maxY != tempData.minY and tempData.minX != tempData.maxX)
    tempData.originY = (
        startY + tempData.minY/(tempData.maxY - tempData.minY)*maxHeight)
    tempData.scaleX = maxWidth/(tempData.maxX - tempData.minX)
    tempData.scaleY = maxHeight/(tempData.maxY - tempData.minY)
    tempData.graph = graph
    data.tempData = tempData

def drawDerivativePicture(canvas, data):
    # draws derivative picture
    assert(data.type == "derivative")
    drawAxes2D(canvas, data)
    data.graph.draw(canvas, data)
    y0 = data.fxn(data.x)
    x = data.x/data.factor if abs(data.x) > 1/2 else -1/2
    der = derivative(data.fxn, data.x)
    y = -der*(data.x - x) + y0
    start=(data.originX+x*data.scaleX,data.originY-y*data.scaleY)
    x = data.x*data.factor if abs(data.x) > 1/2 else 1/2
    y = -der*(data.x - x) + y0
    end =(data.originX +x*data.scaleX,data.originY-y*data.scaleY)
    canvas.create_line(start, end)

def generateZeroPicture(data, line):
    # generates a temporary struct to lay the groundwork to draw zero point
    x = eval(line); y = 0
    line = line.replace("zero(", "(")
    (fxn, guess) = eval(line)
    generatePointPicture(data, fxn, x, y)

def generateMaxMinPicture(data, line):
    # generates a temporary struct to lay the groundwork to draw max min point
    x, y = eval(line)
    line = line.replace("minimize(", "(")
    line = line.replace("maximize(", "(")
    (fxn, guess) = eval(line)
    generatePointPicture(data, fxn, x, y)

def generatePointPicture(data, fxn, x, y):
    # generates a temporary struct to lay the groundwork to draw point
    tempData = Struct()
    startX, startY = data.width//2, data.height//2
    lengthFactor = 0.75
    maxWidth, maxHeight = startX*lengthFactor, startY*lengthFactor
    tempData.type = "point"; tempData.factor = 2
    tempData.maxX = max(x*tempData.factor, 1/2)
    tempData.minX = min(x*tempData.factor, -1/2)
    tempData.start, tempData.end = tempData.minX, tempData.maxX
    tempData.originX=startX-tempData.minX/(tempData.maxX-tempData.minX)*maxWidth
    tempData.increment = 100
    tempData.room = data.room
    graph = Cartesian2DyDep(fxn); graph.generateVectors(tempData)
    tempData.minY = min([y for x, y in graph.vectors] + [0])
    tempData.maxY = max([y for x, y in graph.vectors] + [0])
    assert(tempData.maxY != tempData.minY and tempData.minX != tempData.maxX)
    tempData.originY = (
        startY + tempData.minY/(tempData.maxY - tempData.minY)*maxHeight)
    tempData.scaleX = maxWidth/(tempData.maxX - tempData.minX)
    tempData.scaleY = maxHeight/(tempData.maxY - tempData.minY)
    tempData.graph = graph
    tempData.point = Point(x, y)
    data.tempData = tempData

def drawPointPicture(canvas, data):
    # draws derivative picture
    assert(data.type == "point")
    drawAxes2D(canvas, data)
    data.graph.draw(canvas, data)
    data.point.draw(canvas, data)

def evaluateStatsTable(data):
    # evaluates whatever is in the table in statistics
    data.evalTable = [[None]*data.cols for i in range(data.rows)]
    data.points = []
    decimal = 6
    for row in range(data.rows):
        isPoint = True
        for col in range(data.cols):
            try:
                ans = eval(data.table[row][col])
                if type(ans) == float: ans = round(ans, decimal)
                elif type(ans) != int: isPoint = False
                data.evalTable[row][col] = ans
            except:
                data.evalTable[row][col] = data.table[row][col]
                isPoint = False
        if isPoint:
            x, y = tuple(data.evalTable[row])
            data.points.append(Point(x, y))
    generateStatsGraph(data)

def generateStatsGraph(data):
    # creates the scale of the graph in stats mode
    data.maxX = max([data.defaultMax] + [point.x for point in data.points])
    data.minX = min([data.defaultMin] + [point.x for point in data.points])
    data.maxY = max([data.defaultMax] + [point.y for point in data.points])
    data.minY = min([data.defaultMin] + [point.y for point in data.points])
    data.scaleX = data.width/data.maxX*data.scaleFactor
    data.scaleY = data.height/data.maxY*data.scaleFactor
    factor = 0.45
    startX, startY = data.width*factor, data.height*(1-factor)
    maxWidth, maxHeight = data.width*factor, data.height*factor
    data.originX = startX - data.minX/(data.maxX - data.minX)*maxWidth
    data.originY = startY + data.minY/(data.maxY - data.minY)*maxHeight
    data.scaleX = maxWidth/(data.maxX - data.minX)
    data.scaleY = maxHeight/(data.maxY - data.minY)

####################################
# add messages
####################################

def addMessage(data, canvas):
    # adds messages to canvas
    point = data.width - data.margin, data.height - data.margin
    text = "%s\n%s" % (data.exeMessage, data.message)
    canvas.create_text(point, text = text, anchor = SE, justify = RIGHT)

####################################
# add menus
####################################

def addMenu(data, canvas):
    # draws menu on canvas
    menuWidth = data.menuWidth
    menuHeight = data.menuHeight
    for i in range(len(data.menu)):
        text = data.menu[i]
        start = i*menuWidth
        end = (i+1)*menuWidth
        color = "sky blue" if data.menuActivated[i] else "light blue"
        canvas.create_rectangle(start, 0, end, menuHeight, fill = color)
        canvas.create_text((end+start)/2, menuHeight/2, text = text)
        if data.menuActivated[i]: addDropdownMenu(data, canvas, i)

def addDropdownMenu(data, canvas, i):
    # draws dropdown menu on canvas
    menu = data.dropdown[i]
    for j in range(len(menu)):
        x1, x2 = i*data.menuWidth, (i+1)*data.menuWidth
        y1, y2 = (j+1)*data.menuHeight, (j+2)*data.menuHeight
        text = menu[j]
        canvas.create_rectangle(x1, y1, x2, y2, fill = "sky blue")
        canvas.create_text(x1 + data.miniMargin, (y1+y2)/2, anchor = W,
                           text = text)

####################################
# draw buttons
####################################

def drawButtons(data, canvas):
    # draws buttons on canvas
    for i in range(len(data.buttons)):
        x1,y1=data.width-data.room-data.buttonSize,i*data.buttonSize+data.room
        x2,y2=data.width - data.room, (i+1)*data.buttonSize+data.room
        canvas.create_rectangle(x1,y1,x2,y2, fill = "deep sky blue")
        canvas.create_text((x1+x2)/2,(y1+y2)/2, text = data.buttons[i])
    
####################################
# customize these functions
####################################

def init(data):
    data.mode = "MATHLAB 3D"
    data.modes = ["MATHLAB 3D", "MATHLAB 2D", "MATHLAB Calculator",
                  "MATHLAB Statistics", "Differential Equations", "PHYSLAB"]
    data.keyTable = {"asterisk": "*",
                     "slash": "/",
                     "plus": "+",
                     "minus": "-",
                     "asciicircum": "**",
                     "parenleft": "(",
                     "parenright": ")",
                     "period": ".",
                     "comma": ",",
                     "space": " ",
                     "quotedbl": '"',
                     "quoteright": "'",
                     "less": "<",
                     "greater": ">",
                     "percent": "%",
                     "colon": ":",
                     "semicolon": ";",
                     "bracketleft": "[",
                     "bracketright" : "]",
                     "exclam": "!",
                     "at": "@",
                     "numbersign": "#",
                     "dollar": "$",
                     "ambersand": "&",
                     "underscore": "_",
                     "equal": "=",
                     "backslash": "\\",
                     "bar": "|",
                     "braceleft": "{",
                     "braceright": "}",
                     "question": "?",
                     "asciitilde": "~",
                     "quoteleft": "`"
                     }
    # init functions as a dictionary
    data.functions = {}
    data.functions["MATHLAB 3D"] = [Cartesian3D("")]
    data.functions["MATHLAB 2D"] = [Cartesian2DyDep("")]
    data.functions["Differential Equations"] = [DiffEq("")]
    data.functions["MATHLAB Calculator"] = [""]
    data.functions["PHYSLAB"] = [Particle()]
    data.margin = 10
    data.results = []
    # initial mode is 3D
    initStats(data)
    init3D(data)

def initCalc(data):
    # initialize calculator mode
    data.message = "Type to enter an expression. Hit return to evaluate."
    data.exeMessage = ""
    generateHelpMessageCalc(data)
    data.helpPressed = False
    data.dialogBox = None
    data.inputData = []
    data.bigMargin = 25
    data.biggerMargin = 50
    initMenuCalc(data)
    data.buttons = ["?"]
    data.buttonSize = 25

def generateHelpMessageCalc(data):
    data.helpMessage = """Help: MATHLAB Calculator\n
Simply type an expression and hit return to evaluate.
For example: If you type < '1' '+' '1' 'Return' >, it will display '2'.\n
You can try this for various built-in functions.
For example: An input of 'sinh(1)' will return 1.175...
Or, an input of sqrt(2) will return 1.414...\n
You can also go into the menu to click on a built-in function.
When the dialog box appears, type in the dialog box to input data.
Hit 'tab' to proceed to the next input and hit 'return' to submit your input.\n
Note that pi and e can be easily referenced as 'pi' and 'e' respectively.
The last answer can always be referenced by calling 'ans'\n
Have fun!\n\nClick on '?' to return.\n\nMATHLAB (Mathematics Laboratory) was \
developed by Xiong-Fei Du in Spring 2016
as a term project for 15-112 at Carnegie Mellon University."""

def initMenuCalc(data):
    # initialize calculator menu
    data.calculus = ["Limit", "Derivative at point", "Second derivative",
                     "Definite Integral", "Solve for 0", "Local Maximum",
                     "Local Minimum", "Definite Double Integral",
                     "Definite Triple Integral", "Gradient at a point",
                     "Curl at a point", "Divergence at a point",
                     "Laplacian at a point", "Local Maximum 3D",
                     "Local Minimum 3D"]
    data.probability = ["Normal Distribution", "Inverse Normal",
                        "Student's t-distribution",
                        "Inverse Student's t", "Exponential Distribution",
                        "Gamma Distribution", "Beta Distribution",
                        "Binomial Distribution", "Negative Binomial",
                        "Geometric Distribution", "Hypergeometric",
                        "Poisson Distribution"]
    data.discrete = ["Factorial (n!)","Permutation (nPr)","Combination (nCr)",
                     "Series", "Sequence: Explicit", "Sequence: Recursive"]
    data.menu = ["MATHLAB", "Calculus", "Probability", "Discrete"]
    data.dropdown = [data.modes,data.calculus,data.probability,data.discrete]
    data.menuWidth = 180
    data.menuHeight = 25
    reinitMenu(data)

def initStats(data):
    # restart stats mode
    data.rows = 28
    data.cols = 2
    data.table = [[""]*data.cols for i in range(data.rows)]
    data.bestFit = None
    data.bestFitNums = None
    data.inputMessage = ""
    data.outputMessage = ""
    initGraphStats(data)
    evaluateStatsTable(data)

def reinitStats(data):
    # initialize stats mode
    data.message = "Type to enter a value. Click on cell to navigate."
    data.exeMessage = ""
    data.cellHeight = 25
    data.cellWidth = 120
    data.cellSelected = (0, 0)
    initGraphStats(data)
    evaluateStatsTable(data)
    initMenuStats(data)
    data.dialogBox = None
    data.buttons = ["?"]
    data.buttonSize = 25
    generateHelpMessageStats(data)
    data.helpPressed = False

def generateHelpMessageStats(data):
    data.helpMessage = """Help: MATHLAB Statistics\n
Click on a cell to select a cell, then type to give it a value.\n
You can type any expression, number, or combination of letters inside a cell.
For example, if you type cos(0) into a cell, that cell will evaluate to 1.\n
If you would like to plot points on the graph,
simply treat the left column as your x values and the right column as your y \
values.\nThese inputs will be evaluated and plotted automatically.
You can also find the line of best fit of those points using the regression \
options in the menu.\nYou can access other built-in analysis tools as well, \
such as confidence intervals.\n
Note that pi and e can be easily referenced as 'pi' and 'e' respectively.\n
Have fun!\n\nClick on '?' to return.\n\nMATHLAB (Mathematics Laboratory) was \
developed by Xiong-Fei Du in Spring 2016
as a term project for 15-112 at Carnegie Mellon University."""

def initMenuStats(data):
    # initialize stats menu
    data.stats = ["Statistics: Column 1", "Statistics: Column 2",
                  "Statistics: Two-Variable", "Linear Regression",
                  "Exponential Regression", "Logarithmic Regression",
                  "Power Regression", "Polynomial Regression"]
    data.tests = ["z-interval: 1-Var Statistics", "z-interval: Column 1",
                  "z-interval: Column 2", "z-interval: 2-Var Statistics",
                  "z-interval: 2-Var Data", "z-interval: 1-Var Prop",
                  "z-interval: 2-Var Prop",
                  "t-interval: 1-Var Statistics",
                  "t-interval: Column 1", "t-interval: Column 2",
                  "t-interval: 2-Var Statistics", "t-interval: 2-Var Data"]
    data.tableOptions = ["Clear all"]
    data.menu = ["MATHLAB", "Statistics", "Confidence", "Data"]
    data.dropdown = [data.modes,data.stats,data.tests, data.tableOptions]
    data.menuWidth = 180
    data.menuHeight = 25
    data.border = 50
    reinitMenu(data)

def initGraphStats(data):
    # initialize the graph in stats mode
    factor = 2
    data.originX = data.width*(1-1/factor)
    data.originY = data.height/factor
    data.defaultMin, data.defaultMax = -1, 5
    data.minX,data.maxX = data.defaultMin, data.defaultMax
    data.minY,data.maxY = data.defaultMin, data.defaultMax
    data.scaleFactor = 0.4
    data.scaleX = data.width/data.maxX*data.scaleFactor
    data.scaleY = data.height/data.maxY*data.scaleFactor

def init2D(data):
    # initialize 2D mode
    data.dialogFunctions = DialogFunctions(data)
    data.dialogBox = None
    data.increment = 200
    data.incrementCurve = 500
    data.incrementField2D = 10
    initGraph2D(data)
    data.message = ""
    data.exeMessage = ""
    data.buttons = ["?","H","+","-"]
    data.buttonSize = 25
    data.room = 10
    initMenu2D(data)
    generateHelpMessage2D(data)
    data.helpPressed = False

def generateHelpMessage2D(data):
    data.helpMessage = """Help: MATHLAB 2D\n
Type to enter a function.
Once you have given a valid function, the graph will automatically display.\n
Be sure to give a valid function. For example, if the command line begins with \
'y = f(x)',\n be sure to input a function in terms of x.
For example: 'sqrt(1 - x**2)' is a valid input. 'a - 7*b' is not valid.
Or, if the command line reads 'r = f(theta)', '2*sin(theta) - 3' is valid.
Or, if the command line reads 'r(t) = x(t), y(t), t_min, t_max',
an input of 'cos(t), sin(t), 0, 7' is a valid input.\n
Also, be sure to include '*' when implying multiplication and '**' when \
implying powers.\nGood: '7*x**2 + 2'. Bad: '7x2 + 2' (This is not valid!)\n
If you ever want to change the coordinate system, simply click on what you'd \
like in the menu.\n\nHit 'Return' to add a new line.\n
Note that pi and e can be easily referenced as 'pi' and 'e' respectively.\n
Press '+' or '-' to zoom in and out. Press 'H' to return to default view.\n
You can also clear functions in the menu to get rid of ones you don't need.\n
Have fun!\n\nClick on '?' to return.\n\nMATHLAB (Mathematics Laboratory) was \
developed by Xiong-Fei Du in Spring 2016
as a term project for 15-112 at Carnegie Mellon University."""

def init3D(data):
    # load data.xyz as appropriate
    # initialize 3D mode
    data.dialogFunctions = DialogFunctions(data)
    data.dialogBox = None
    data.room = 10
    data.dAngleFactor = 12
    data.dAngle = math.pi/data.dAngleFactor
    data.increment = 36
    data.incrementCurve = 400
    data.incrementField3D = 4
    data.miniMargin = 5
    initGraph3D(data)
    data.exeMessage = "Press arrows to rotate."
    data.message = ""
    data.buttons = ["?","H","+","-"]
    data.buttonSize = 25
    initMenu3D(data)
    generateHelpMessage3D(data)
    data.helpPressed = False

def generateHelpMessage3D(data):
    data.helpMessage = """Help: MATHLAB 3D\n\nType to enter a function.
Once you have given a valid function, the graph will automatically display.\n
Be sure to give a valid function. For example, if the command line begins with \
'z = f(x, y)',\n be sure to input a function in terms of x and y.
For example: 'sinh(x**2) + sqrt(x*y)' is a valid input. 'a - 7*b' is not valid.
Or, if the command line reads 'z = f(r, theta)', '2*r*sin(theta)' is valid.
Or, if the command line reads 'r(t) = x(t), y(t), z(t), t_min, t_max',
an input of 'cos(t), sin(t), t, -10, 10' is a valid input.\n
Also, be sure to include '*' when implying multiplication and '**' when \
implying powers.\nGood: '7*x**2 + 2'. Bad: '7x2 + 2' (This is not valid!)\n
If you ever want to change the coordinate system, simply click on what you'd \
like in the menu.\n\nHit 'Return' to add a new line.\n
Note that pi and e can be easily referenced as 'pi' and 'e' respectively.\n
Press '+' or '-' to zoom in and out. Press 'H' to return to default view.
Press the arrow keys to rotate the angle of view.\n
You can also clear functions in the menu to get rid of ones you don't need.\n
Have fun!\n\nClick on '?' to return.\n\nMATHLAB (Mathematics Laboratory) was \
developed by Xiong-Fei Du in Spring 2016
as a term project for 15-112 at Carnegie Mellon University."""

def initGraph2D(data):
    # initialize 2D graph
    data.originX = data.width//2
    data.originY = data.height//2
    data.minX,data.maxX = -5,5
    data.minY,data.maxY = -5,5
    data.scaleFactor = 5/16
    data.scaleX = data.width/data.maxX*data.scaleFactor
    data.scaleY = data.height/data.maxY*data.scaleFactor
    data.zoomFactor = 1.2
    for graph in data.functions[data.mode]: graph.generateVectors(data)

def initMenu2D(data):
    # initialize menu in 2D mode
    data.coordinate = "Cartesian: y-Dependent"
    data.coordinates = ["Cartesian: y-Dependent","Cartesian: x-Dependent",
                        "Polar","Parametric","Vector Field 2D"]
    data.functionOptions = ["Clear all", "Remove current", "Edit axes"]
    data.menu = ["MATHLAB", "Coordinates", "Functions"]
    data.dropdown = [data.modes,data.coordinates,data.functionOptions]
    data.menuWidth = 180
    data.menuHeight = 25
    reinitMenu(data)

def initDiffEq(data):
    # initialize diff eq mode
    data.orders = ["First Order", "Second Order"]
    data.dialogFunctions = DialogFunctions(data)
    data.dialogBox = None
    data.incrementField2D = 20
    data.stepFactor = 500
    data.message = ""
    data.exeMessage = ""
    data.buttons = ["?","H","+","-"]
    data.buttonSize = 25
    data.room = 10
    initGraphDiffEq(data)
    initMenuDiffEq(data)
    generateHelpMessageDiffEq(data)
    data.helpPressed = False
    reinitDifferentialEquations(data)

def generateHelpMessageDiffEq(data):
    data.helpMessage = """Help: MATHLAB First-Order and Second-Order Ordinary \
Differential Equations\n\nType to enter a first-order or second-order \
ordinary differential equation.\nOnce you have given a valid function, \
the slope field and/or solution will automatically display.\n
Be sure to give a valid function. For first-order equations, the command line
begins with "y' = f(x, y)". Therefore, be sure to input a function in terms of x and y.
For example: 'cos(x*y)' is a valid input. 'sin(t - a)' is not valid.
For second-order, the command line begins with 'y" = f(x, y, Dy)', where Dy
represents y'. For example, 'y**2 - 2*Dy' is a valid input.\n
Optional for first-order equations: You can also specify an initial condition
by entering an initial point separated by commas.
For example, 'cos(x*y), -2, 3' is a valid input.\n\nHowever, note that an \
initial condition must be provided for second-order equations.\n
Also, be sure to include '*' when implying multiplication and '**' when \
implying powers.\nGood: '7*x**2 + 2'. Bad: '7x2 + 2' (This is not valid!)\n
Note that pi and e can be easily referenced as 'pi' and 'e' respectively.\n
Press '+' or '-' to zoom in and out. Press 'H' to return to default view.
You can also clear the slope field in the menu.\n
Have fun!\n\nClick on '?' to return.\n\nMATHLAB (Mathematics Laboratory) was \
developed by Xiong-Fei Du in Spring 2016
as a term project for 15-112 at Carnegie Mellon University."""

def initGraphDiffEq(data):
    initGraph2D(data)

def initMenuDiffEq(data):
    # initialize diff eq menu
    try: data.order
    except: data.order = "First Order"
    getMenuDiffEq(data)
    data.menuWidth = 180
    data.menuHeight = 25
    reinitMenu(data)

def getMenuDiffEq(data):
    if data.order.startswith("First"):
        data.functionOptions = ["New field", "Clear field", "Edit axes"]
    elif data.order.startswith("Second"):
        data.functionOptions = ["Add plot", "Remove plot", "Edit axes"]
    data.menu = ["MATHLAB", "Differential Equations", data.order]
    data.dropdown = [data.modes, data.orders, data.functionOptions]

def initGraph3D(data):
    # initialize 3D graph
    data.scaleFactor = 0.33333333333
    data.theta = math.pi*data.scaleFactor/2
    data.phi = math.pi*data.scaleFactor
    data.originX = data.width//2
    data.originY = data.height//2
    data.minX,data.maxX = -5,5
    data.minY,data.maxY = -5,5
    data.minZ,data.maxZ = -5,5
    data.scaleX = data.width/(data.maxX - data.minX)*data.scaleFactor*2
    data.scaleY = data.width/(data.maxY - data.minY)*data.scaleFactor*2
    data.scaleZ = data.width/(data.maxZ - data.minZ)*data.scaleFactor*2
    data.zoomFactor = 1.2
    for graph in data.functions[data.mode]: graph.generateVectors(data)

def initPhys(data):
    # initialize physics
    init3D(data)
    data.incrementCurve = 1000
    initMenuPhys(data)
    generateHelpMessagePhys(data)

def generateHelpMessagePhys(data):
    data.helpMessage = """Help: PHYSLAB\n\nMATHLAB Physics Simulator\n
Go into the menu to make changes to the default parameters and forces.\n
In the menu, 'General Force' can describe either contact or field forces.\n
Be sure to follow the units given.\n
The initial conditions describe the particle at the starting time.
The default particle starts from rest at the origin.\n
This simulator is for point masses and point charges.\n
Uncharacteristically sharp turns in the plot may be signs of either the time \
span being set\ntoo large or chaotic motion, which are subject \
to computational imprecision.\n
Also, be sure to include '*' when implying multiplication and '**' when \
implying powers.\nGood: '7*x**2 + 2'. Bad: '7x2 + 2' (This is not valid!)\n
Press '+' or '-' to zoom in and out. Press 'H' to return to default view.
Press the arrow keys to rotate the angle of view.\n
Have fun!\n\nClick on '?' to return.\n\nMATHLAB (Mathematics Laboratory) was \
developed by Xiong-Fei Du in Spring 2016
as a term project for 15-112 at Carnegie Mellon University."""

def initMenuPhys(data):
    # initialize physics graph menu
    data.fields = ["General Force", "Gravitational Field",
                   "Electric Field", "Magnetic Field"]
    data.particle = ["Edit Particle Properties", "Edit Initial Conditions",
                     "Edit Timescale", "Edit axes"]
    data.menu = ["MATHLAB", "Forces & Fields", "Simulation"]
    data.dropdown = [data.modes, data.fields, data.particle]
    data.menuWidth = 180
    data.menuHeight = 25
    reinitMenu(data)

def initMenu3D(data):
    # initialize 3D graph menu
    data.coordinate = "Cartesian"
    data.coordinates = ["Cartesian","Cylindrical: Z-Dependent",
                        "Cylindrical: R-Dependent","Spherical",
                        "Parametric Curves","Vector Field"]
    data.functionOptions = ["Clear all", "Remove current", "Edit axes"]
    data.menu = ["MATHLAB", "Coordinates", "Functions"]
    data.dropdown = [data.modes,data.coordinates,data.functionOptions]
    data.menuWidth = 180
    data.menuHeight = 25
    reinitMenu(data)

def reinitMenu(data):
    # restart the menu
    data.dropdownHeight = 0
    data.menuSelected  = [-1, -1]
    data.menuActivated = [False]*len(data.menu)

def mouseMotion(event, data):
    # mouse moves over menu
    if event.y < data.menuHeight and event.x < len(data.menu)*data.menuWidth:
        data.menuActivated = [False]*len(data.menu)
        selected = event.x // data.menuWidth
        data.menuSelected[0] = selected         
        data.menuActivated[selected] = True
        data.dropdownHeight = data.menuHeight*(len(data.dropdown[selected]) + 1)
    elif (event.x > data.menuSelected[0]*data.menuWidth and
          event.x < (data.menuSelected[0] + 1)*data.menuWidth and
          event.y < data.dropdownHeight): pass
    else: reinitMenu(data)        

def mousePressed(event, data):
    # use event.x and event.y
    data.exeMessage = ""
    if data.mode == "MATHLAB Statistics": mousePressedStats(event, data)
    if data.dialogBox != None: dialogBoxMousePressed(event, data)
    elif data.mode in ("MATHLAB 3D", "MATHLAB 2D", "Differential Equations"):
        data.dialogFunctions.clickedOn(event.x, event.y)
        if data.mode == "Differential Equations":
            reinitDifferentialEquations(data)
    # check if menu is pressed
    if (event.x > data.menuSelected[0]*data.menuWidth and
        event.x < (data.menuSelected[0]+1)*data.menuWidth and
        event.y < data.dropdownHeight):
        data.menuSelected[1] = event.y // data.menuHeight -1
        executeMenu(data)
    # check if mouse is pressed
    if (event.x + data.room + data.buttonSize > data.width and
        event.x + data.room < data.width and
        event.y - data.room > 0 and
        event.y - data.room - data.buttonSize*len(data.buttons)):
        button = (event.y - data.room)//data.buttonSize
        executeButtons(data, button)

def dialogBoxMousePressed(event, data):
    # mouse pressed for dialog box
    data.dialogBox.clickedOn(event.x, event.y)
    for button in data.dialogBox.buttons:
        if button.clickedOn(event.x, event.y):
            if button.text == "OK": OKPressed(data)
            elif button.text == "Cancel": data.dialogBox = None

def keyPressed(event, data):
    # use event.char and event.keysym
    mode = data.mode
    if mode == "MATHLAB 3D" or mode == "PHYSLAB":
        if event.keysym == "Left": data.theta -= data.dAngle
        elif event.keysym == "Right": data.theta += data.dAngle
        elif event.keysym == "Up": data.phi -= data.dAngle
        elif event.keysym == "Down": data.phi += data.dAngle
        if data.phi < 0: data.phi = 0
        elif data.phi > math.pi: data.phi = math.pi
    if data.dialogBox != None:
        dialogBoxKeyPressed(event, data); return
    if event.keysym == "BackSpace": backspacePressed(event, data)
    elif event.keysym == "Return": returnPressed(event, data)
    modifyFunction(data, event.keysym)
    data.exeMessage = ""
    if data.mode in ("MATHLAB 3D", "MATHLAB 2D", "Differential Equations"):
        data.dialogFunctions.regenerateText(data)

def dialogBoxKeyPressed(event, data):
    # key pressed when dialog box is activated
    mode = data.mode
    data.dialogBox.editOutputs(event.keysym)
    if event.keysym == "Return": OKPressed(data)

def OKPressed(data):
    # what to do once you press OK on a dialog box
    data.inputTitle = data.dialogBox.title
    data.inputs = [bar.text for bar in data.dialogBox.outputs]
    if data.mode == "MATHLAB Calculator": executeInput(data)
    elif data.mode == "MATHLAB Statistics":
        try: executeInputStats(data)
        except: pass
    elif data.mode == "PHYSLAB":
        try: executeInputPhys(data)
        except: pass
    else: executeInputGraph(data)
    data.dialogBox = None

def returnPressed(event, data):
    # what to do when you press return
    mode = data.mode
    if mode == "Differential Equations":
        if data.order.startswith("First"): return
        elif data.order.startswith("Second"):
            data.functions[data.mode].append(SecondOrderDE(""))
            data.dialogFunctions.edit = len(data.functions[mode]) - 1
    elif mode == "MATHLAB Calculator":
        if data.functions[mode][0][-1] != "\n":
            evaluatePage(data)
            data.functions[mode][0] += "\n\n\n"
        return
    elif mode == "MATHLAB Statistics":
        row, col = data.cellSelected
        if row < data.rows - 1: row += 1
        data.cellSelected = row, col
        return
    elif mode == "MATHLAB 3D": addNewFunction3D(data, event.keysym)
    elif mode == "MATHLAB 2D": addNewFunction2D(data, event.keysym)

def backspacePressed(event, data):
    # what to do when you press backspace
    mode = data.mode
    if data.mode == "MATHLAB Calculator":
        if (data.functions[mode][0] == "" or
            data.functions[mode][0][-1] == "\n"): return
        data.functions[mode][0]=data.functions[mode][0][:-1]
    elif data.mode == "MATHLAB Statistics":
        row, col = data.cellSelected
        text = data.table[row][col]
        text = text[:-1] if len(text) > 0 else text
        data.table[row][col] = text
    else:
        num = data.dialogFunctions.edit
        data.functions[mode][num].fxn=data.functions[mode][num].fxn[:-1]

def addNewFunction3D(data, selected):
    # add new function in 3D
    mode = data.mode; fxnNum = data.dialogFunctions.edit
    if data.functions[mode]!=[] and data.functions[mode][fxnNum].fxn == "":
        data.functions[mode].pop(fxnNum)
    if data.coordinate == "Cartesian":
        data.functions[mode].append(Cartesian3D(""))
    elif data.coordinate == "Cylindrical: Z-Dependent":
        data.functions[mode].append(CylindricalZDependent(""))
    elif data.coordinate == "Cylindrical: R-Dependent":
        data.functions[mode].append(CylindricalRDependent(""))
    elif data.coordinate == "Cylindrical: Z-Dependent":
        data.functions[mode].append(CylindricalZDependent(""))
    elif data.coordinate == "Spherical":
        data.functions[mode].append(Spherical(""))
    elif data.coordinate == "Parametric Curves":
        data.functions[mode].append(Parametric3D(""))
    elif data.coordinate == "Vector Field":
        data.functions[mode].append(VectorField3D(""))
    data.dialogFunctions.edit = len(data.functions[mode]) - 1
    data.dialogFunctions.reinit(data)

def addNewFunction2D(data, selected):
    # add new function in 2D
    mode = data.mode
    fxnNum = data.dialogFunctions.edit
    if data.functions[mode]!=[] and data.functions[mode][fxnNum].fxn == "":
        data.functions[mode].pop(fxnNum)
    if data.coordinate == "Cartesian: x-Dependent":
        data.functions[mode].append(Cartesian2DxDep(""))
    elif data.coordinate == "Cartesian: y-Dependent":
        data.functions[mode].append(Cartesian2DyDep(""))
    elif data.coordinate == "Polar":
        data.functions[mode].append(Polar(""))
    elif data.coordinate == "Parametric":
        data.functions[mode].append(Parametric2D(""))
    elif data.coordinate == "Vector Field 2D":
        data.functions[mode].append(VectorField2D(""))
    data.dialogFunctions.edit = len(data.functions[mode]) - 1
    data.dialogFunctions.reinit(data)

def modifyFunction(data, key):
    # modify the calculator input page
    mode = data.mode
    if mode == "MATHLAB Statistics": modifyStatsTable(data, key)
    elif mode != "MATHLAB Calculator": modifyFunctionGraph(data, key)
    elif key in data.keyTable: data.functions[mode][0] += data.keyTable[key]
    elif len(key) == 1: data.functions[mode][0]+=key

def modifyFunctionGraph(data, key):
    # modify the function text
    fun = data.functions[data.mode][data.dialogFunctions.edit]
    if key in data.keyTable: fun.fxn += data.keyTable[key]
    elif len(key) == 1: fun.fxn += key
    elif key != "BackSpace": return
    fun.generateLambda()
    fun.generateVectors(data)

def modifyStatsTable(data, key):
    # modify the stats table
    row, col = data.cellSelected
    if key == "asterisk": data.table[row][col] += "*"
    elif key == "slash": data.table[row][col] += "/"
    elif key == "plus": data.table[row][col] += "+"
    elif key == "minus": data.table[row][col] += "-"
    elif key == "asciicircum": data.table[row][col] += "**"
    elif key == "parenleft": data.table[row][col] += "("
    elif key == "period": data.table[row][col] += "."
    elif key == "parenright": data.table[row][col] += ")"
    elif key == "comma": data.table[row][col] += ","
    elif key == "space": data.table[row][col] += " "
    elif len(key) == 1: data.table[row][col] += key
    evaluateStatsTable(data)

def timerFired(data):
    pass

def mousePressedStats(event, data):
    # if you click on a cell in the stats table
    x = event.x
    y = event.y
    col = (x - data.border)//data.cellWidth
    row = (y - data.border)//data.cellHeight
    if row in range(data.rows) and col in range(data.cols):
        data.cellSelected = row, col

def redrawAll(canvas, data):
    # draw on canvas
    if data.mode == "MATHLAB 3D": redrawAll3D(canvas, data)
    elif data.mode == "MATHLAB 2D": redrawAll2D(canvas, data)
    elif data.mode == "MATHLAB Calculator": redrawAllCalc(canvas, data)
    elif data.mode == "MATHLAB Statistics": redrawAllStats(canvas, data)
    elif data.mode == "Differential Equations": redrawAllDiffEq(canvas, data)
    elif data.mode == "PHYSLAB": redrawAllPhys(canvas, data)
    if data.dialogBox != None: data.dialogBox.draw(canvas, data)
    if data.helpPressed: drawHelpMessage(data, canvas)
    drawButtons(data, canvas)
    addMessage(data, canvas)
    addMenu(data, canvas)

def redrawAll2D(canvas, data):
    # draw the 2D graph
    drawAxes2D(canvas, data)
    for graph in data.functions[data.mode]:
        graph.draw(canvas, data)
    data.dialogFunctions.draw(canvas, data)

def redrawAll3D(canvas, data):
    # draw the 3D graph
    drawAxes3D(canvas, data)
    for graph in data.functions[data.mode]:
        graph.draw(canvas, data)
    data.dialogFunctions.draw(canvas, data)

def redrawAllDiffEq(canvas, data):
    # draw the diff eq graph
    drawAxes2D(canvas, data)
    if data.order.startswith("First"):
        for graph in data.functions[data.mode]:
            if type(graph) == DiffEq: graph.draw(canvas, data); break
    elif data.order.startswith("Second"):
        for graph in data.functions[data.mode]:
            if type(graph) == SecondOrderDE: graph.draw(canvas, data)
    data.dialogFunctions.draw(canvas, data)

def redrawAllCalc(canvas, data):
    font = "Arial 10"
    # draw picture
    try: drawIntegralPicture(canvas, data.tempData)
    except: pass
    try: drawDerivativePicture(canvas, data.tempData)
    except: pass
    try: drawPointPicture(canvas, data.tempData)
    except: pass
    # draw the calculator page input and output
    canvas.create_text(data.bigMargin, data.height - data.biggerMargin,
                       font = font, text = "%s|"% data.functions[data.mode][0],
                       anchor = SW)
    canvas.create_text(data.biggerMargin, data.height - data.biggerMargin,
                       font = font, text = "%s\n\n"%"\n\n\n".join(data.results),
                       anchor = SW)

def redrawAllStats(canvas, data):
    # draw table, graph, data in stats mode
    drawTable(canvas, data)
    drawAxes2D(canvas, data)
    for point in data.points: point.draw(canvas, data)
    if data.bestFit != None:
        data.bestFit.draw(canvas, data)
        r = data.correlationCoefficient
        text = "y = %s; R^2 = %.4f" % (data.bestFit.fxn, r)
        canvas.create_text(data.originX + data.room, data.originY, anchor = NW,
            text = text, fill = data.bestFit.color)
    drawGraphAndBox(canvas, data)

def drawTable(canvas, data):
    # draw table in stats mode
    for row in range(len(data.table)):
        for col in range(len(data.table[0])):
            xLeft = data.border + col*data.cellWidth
            xRight= data.border + (col+1)*data.cellWidth
            yTop = data.border + row*data.cellHeight
            yBottom = data.border + (row+1)*data.cellHeight
            color = "red" if data.cellSelected == (row, col) else "white"
            canvas.create_rectangle(xLeft, yTop, xRight, yBottom, fill = color)
            if (row, col) == data.cellSelected: text = data.table[row][col]
            else: text = data.evalTable[row][col]
            canvas.create_text(xRight, yBottom, anchor = SE, justify = RIGHT,
                               text = text)

def drawGraphAndBox(canvas, data):
    # draws both the graph and box in stats mode
    centerFactor = 0.42
    edgeFactor = 0.95
    inputFactor = 0.72
    font = "Arial 10"
    canvas.create_rectangle(data.width*centerFactor,
                            data.height*(1 - centerFactor),
                            data.width *edgeFactor,
                            data.height*edgeFactor, fill = "sky blue")
    if data.inputMessage != "":
        canvas.create_text(data.width*inputFactor+data.margin,
                           data.height*(1 - centerFactor) + data.margin,
                           font = font, text = data.outputMessage,
                           justify = LEFT, anchor = NW)
        canvas.create_text(data.width*inputFactor,
                           data.height*(1 - centerFactor) + data.margin,
                           font = font, text = data.inputMessage,
                           justify = RIGHT, anchor = NE)

def redrawAllPhys(canvas, data):
    # draw physics
    drawAxes3D(canvas, data)
    for graph in data.functions[data.mode]:
        graph.draw(canvas, data)

####################################
# use the run function as-is
####################################

def run(width=300, height=300):
    def redrawAllWrapper(canvas, data):
        canvas.delete(ALL)
        redrawAll(canvas, data)
        canvas.update()    

    def mousePressedWrapper(event, canvas, data):
        mousePressed(event, data)
        redrawAllWrapper(canvas, data)

    def mouseMotionWrapper(event, canvas, data):
        mouseMotion(event, data)
        redrawAllWrapper(canvas, data)

    def keyPressedWrapper(event, canvas, data):
        keyPressed(event, data)
        redrawAllWrapper(canvas, data)

    def timerFiredWrapper(canvas, data):
        timerFired(data)
        redrawAllWrapper(canvas, data)
        # pause, then call timerFired again
        canvas.after(data.timerDelay, timerFiredWrapper, canvas, data)
    # Set up data and call init
    class Struct(object): pass
    data = Struct()
    data.width = width
    data.height = height
    data.timerDelay = 10**10 #1000 # milliseconds
    init(data)
    # create the root and the canvas
    root = Tk()
    root.title("MATHLAB 1.05 - Xiong-Fei Du")
    canvas = Canvas(root, width=data.width, height=data.height, bg = "azure")
    canvas.pack()
    # set up events
    root.bind("<Button-1>", lambda event:
                            mousePressedWrapper(event, canvas, data))
    root.bind("<Key>", lambda event:
                            keyPressedWrapper(event, canvas, data))
    root.bind("<Motion>", lambda event:
                            mouseMotionWrapper(event, canvas, data))
    timerFiredWrapper(canvas, data)
    # and launch the app
    root.mainloop()  # blocks until window is closed
    print("bye!")

####################################
# test functions
####################################

def testDerivative():
    assert(equal(derivative(lambda x: x**2, 6), 12))
    assert(equal(derivative(lambda x: x**3, 6), 108))
    assert(equal(derivative(sin, 0), 1))

def testIntegral():
    assert(equal(integral(lambda x: x**3, 0, 0), 0))
    assert(equal(integral(lambda x: x**2, 0, 1), 1/3))
    assert(equal(integral(sin, 0, 2*math.pi), 0))

def testZero():
    assert(equal(zero(lambda x: x**3 - 1, 2), 1))
    assert(equal(zero(lambda x: x**2 - 16, 91), 4))
    assert(equal(zero(sin, 1), 0))

def testAll():
    testDerivative()
    testIntegral()
    testZero()
    print("Passed!")

#testAll()

#run(width = 800, height = 800)

####################################
# non GUI handlers
####################################

def start(width = 800, height = 800):
    global data
    data = Struct()
    data.width = width
    data.height = height
    init(data)
    global canvas
    canvas = generateCanvas(data)
    print("Welcome to MATHLAB 1.05 non-GUI!")
    print("Check pop-up window for graphs.")
    return data, canvas

def setMode(select = None):
    if select == None:
        choices = {"A": "MATHLAB 3D", "B": "MATHLAB 2D",
                   "C": "MATHLAB Calculator", "D": "MATHLAB Statistics",
                   "E": "Differential Equations", "F": "PHYSLAB"}
        print("Select mode:")
        for choice in sorted(choices):
            print("\t%s) %s" % (choice, choices[choice]))
        select = input("Enter mode: ")
        if select in choices:
            executeMenu1(data, choices[select])
            print("Mode set to %s" % choices[select])
        else: print("Failure. Please try again.")
    elif select in data.modes:
        executeMenu1(data, select)
        print("Mode set to %s" % select)
    else: print("Failure. Please try again.")

def generateCanvas(data):
    root = Tk()
    root.title("MATHLAB 1.05 non-GUI Graph Window - Xiong-Fei Du")
    canvas = Canvas(
        root, width = data.width, height = data.height, bg = "azure")
    canvas.pack()
    data.plots = []
    
    def pressKey(event):
        dirs = {"Up", "Down", "Left", "Right"}
        if data.mode in {"MATHLAB 3D", "PHYSLAB"} and event.keysym in dirs:
            keyPressed(event, data)
            canvas.create_rectangle(0, 0, data.width, data.height,
                                    fill = "azure", width = 0)
            drawAxes()
            for f in data.plots:
                f.draw(canvas, data)

    root.bind("<Key>", pressKey)
    data.root = root
    drawLogo(canvas, data)
    drawWelcomeMessage(canvas, data)
    return canvas

def drawWelcomeMessage(canvas, data):
    canvas.create_text(data.width//2, data.height//2,
        text = "MATHLAB 1.05 non-GUI\n\nAll graphs will appear here!",
        justify = CENTER, font = "Cambria 12 bold")

def drawLogo(canvas, data):
    init3D(data)
    logo = CylindricalZDependent("cos(r*3)/3 - 1")
    logo.generateVectors(data)
    logo.draw(canvas, data)

def clear(color = "azure"):
    canvas.create_rectangle(0, 0, data.width, data.height, fill = color,
        width = 0)
    data.plots = []

def wait():
    data.root.mainloop()

def getData():
    return data

def getCanvas():
    return canvas

def axes(*args):
    if len(args) == 0:
        inputs = changeAxes(data, "")
        data.dialogBox = None
        data.inputs = []
        for i in inputs:
            data.inputs.append(input(i + " "))
        editAxes(data)
        print("Done!")
    elif len(args) == 4 or len(args) == 6:
        data.inputs = list(map(str, args))
        editAxes(data)
        print("Done!")
    else: print("Failed. Try again.")

def rotate(select = None):
    event = Struct()
    if select == None:
        choices = {"A": "Up", "B": "Down", "C": "Left", "D": "Right"}
        print("Select direction:")
        for choice in sorted(choices):
            print("\t%s) %s" % (choice, choices[choice]))
        select = input("Enter direction: ")
        if select in choices:
            event.keysym = choices[select]
            keyPressed(event, data)
            print("Success!")
        else: print("Failure. Please try again."); return
    elif select in {"Up", "Down", "Left", "Right"}:
        event.keysym = select
        keyPressed(event, data)
        print("Success!")
    else: print("Failure. Please try again."); return

def drawAxes():
    if data.mode in ("MATHLAB 3D", "PHYSLAB"): drawAxes3D(canvas, data)
    elif data.mode in ("MATHLAB 2D", "Differential Equations"):
        drawAxes2D(canvas, data)

def plot(*args):
    for f in args:
        f.plot(canvas, data)
        data.plots.append(f)

def scatterPlot(points, color = "black"):
    for x, y in points:
        p = Point(x, y)
        p.draw(canvas, data, label = False, color = color)
        data.plots.append(p)

def linePlot(points, color = None):
    f = Cartesian2DyDep(lambda x: None)
    if color != None: f.color = color
    f.vectors = points
    f.draw(canvas, data)
    data.plots.append(f)

def meshPlot(points):
    f = Cartesian3D(lambda x, y: None)
    f.vectors = points
    f.draw(canvas, data)
    data.plots.append(f)

def about():
    print("""\
MATHLAB (Mathematics Laboratory) was developed by Xiong-Fei Du in
Spring 2016 as a term project for 15-112 at Carnegie Mellon University.
This stand-alone non-GUI version was an addition to MATHLAB v. 1.05,
released 29 May 2018.""")
