### Xiong-Fei Du
### Andrew ID: xiongfed
### 15-112 Section I

### MATHLAB Term Project 15-112

### updated version 1.04
### release date 25 Aug 2016

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

ln = lambda x: math.log(x)
log = lambda x, base = 10: math.log(x, base)
sec = lambda x: 1/math.cos(x)
csc = lambda x: 1/math.sin(x)
cot = lambda x: 1/math.tan(x)

####################################
# discontinuous functions
####################################

def root(x, a):
    epsilon = 0.0001
    if x >= 0: return x**(1/a)
    elif abs(round(a) - a) > epsilon: return None
    elif a % 2 == 0: return None
    else: return -(-x)**(1/a)

def heaviside(x):
    if x < 0: return 0
    elif x == 0: return 0.5
    elif x > 0: return 1

def sgn(x):
    if x < 0: return -1
    elif x == 0: return 0
    elif x > 0: return 1

def delta(x):
    h = 100
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
    decimal = 12
    h = 0.0000001
    try:
        x -= h/2
        fx_minus = eval(f)
        x += h
        fx_plus = eval(f)
        return round(float("%.7g" % ((fx_plus - fx_minus)/h)), decimal) # 7 sf
    except: return None

def secondDerivative(f, x):
    # use limit definition of second derivative
    decimal = 12
    h = 0.00001
    try:
        fx = eval(f)
        x += h
        fx_plus = eval(f)
        x -= h*2
        fx_minus = eval(f)
        return round(float("%.4g"%((fx_plus - 2*fx + fx_minus)/h**2)), decimal)
        # 4 sf
    except: return None

def integral(f, a, b):
    #use Simpson's rule (approximates functions with piecewise parabolas)
    intervals = 4000
    decimal = 12
    divide = 3
    runningSum = 0
    dx = (b - a)/intervals
    try:
        for i in range(intervals+1):
            x = i*dx + a
            fx = eval(f)
            if i == 0: factor = 1
            elif i == intervals: factor = 1
            elif i % 2 == 0: factor = 2
            elif i % 2 == 1: factor = 4
            runningSum += factor*fx*dx
        factor = 3
        return round(float("%.10g" % (runningSum/factor)), decimal) # 10 sf
    except: return None

def limit(f, a):
    digits = 6
    # use delta-epsilon definition of limit
    delta = 10**(-10**2)
    epsilon = 0.01
    try: # first try to plug it in
        x = a
        answer = eval(f)
        return round(answer, digits)
    except: # approach from left and right
        try:
            x = a + delta
            left = eval(f)
            x = a - delta
            right = eval(f)
            assert(abs(left-right) < epsilon)
            return round((left+right)/2,digits)
        except: return "limit does not exist"

####################################
# monte carlo integration
####################################

def doubleIntegral(f, constraints, xMin, xMax, yMin, yMax):
    # use monte carlo method to estimate a double integral on complex domain
    points = []
    totalPoints = 20000; maxCount = 1000000
    count = 0
    decimal = 6
    # generate random points within bounding box and check if in constraints
    while len(points) < totalPoints:
        count += 1
        if count > maxCount: return None
        x, y = random.uniform(xMin, xMax), random.uniform(yMin, yMax)
        good = True
        for constraint in constraints:
            if eval(constraint) == False: good = False; break
        if good: points.append((x, y))
    # area = area of bounding box
    area = (xMax - xMin)*(yMax - yMin)*totalPoints/count
    if "x" not in f and "y" not in f and type(eval(f)) in (int, float):
        total = eval(f)*totalPoints
    else:
        total = 0
        for x, y in points: total += eval(f)
    # apply mean value theorem of integrals
    return str(round(total/totalPoints*area, decimal)) + " +/- 1%"

def tripleIntegral(f, constraints, xMin, xMax, yMin, yMax, zMin, zMax):
    # use monte carlo method to estimate a triple integral on complex domain
    points = []; count = 0; totalPoints = 20000; maxCount = 1000000
    decimal = 6
    # generate random points within bounding box and check if in constraints
    while len(points) < totalPoints:
        count += 1
        if count > maxCount: return None
        x, y, z = (random.uniform(xMin, xMax), random.uniform(yMin, yMax),
                   random.uniform(zMin, zMax))
        good = True
        for constraint in constraints:
            if eval(constraint) == False: good = False; break
        if good: points.append((x, y))
    # volume = volume of bounding box
    volume = (xMax - xMin)*(yMax - yMin)*(zMax - zMin)*totalPoints/count
    if ("x" not in f and "y" not in f and "z" not in f and
        type(eval(f))in(int, float)):
        total = eval(f)*totalPoints
    else:
        total = 0
        for x, y in points: total += eval(f)
    # apply mean value theorem of integrals
    return str(round(total/totalPoints*volume, decimal)) + " +/- 1%"

####################################
# Newton's method
####################################

def zero(f, guess):
    # finds zero using Newton's method
    epsilon = 10**-10
    maxCount = 10000
    decimal = 6
    try:
        x = guess
        count = 0
        y = eval(f)
        prev = x
        while abs(y) > epsilon or abs(prev - x) > 10**-(decimal + 1):
            count += 1
            assert(count < maxCount)
            der = derivative(f, x)
            y = eval(f)
            prev = x
            x -= y/der
        return round(x, decimal)
    except:
        return round(x, decimal//2) if abs(y) < epsilon else None
        #except: return None

####################################
# partial derivatives
####################################

def partialDerivativeX(f, x, y, z = None):
    # apply limit definition of a derivative
    decimal = 12
    h = 0.0000001
    try:
        x -= h/2
        fx_minus = eval(f)
        x += h
        fx_plus = eval(f)
        return round(float("%.7g" % ((fx_plus - fx_minus)/h)), decimal) # 7 sf
    except: return None

def partialDerivativeY(f, x, y, z = None):
    # apply limit definition of a derivative
    decimal = 12
    h = 0.0000001
    try:
        y -= h/2
        fx_minus = eval(f)
        y += h
        fx_plus = eval(f)
        return round(float("%.7g" % ((fx_plus - fx_minus)/h)), decimal) # 7 sf
    except: return None

def partialDerivativeZ(f, x, y, z):
    # apply limit definition of a derivative
    decimal = 12
    h = 0.0000001
    try:
        z -= h/2
        fx_minus = eval(f)
        z += h
        fx_plus = eval(f)
        return round(float("%.7g" % ((fx_plus - fx_minus)/h)), decimal) # 7 sf
    except: return None

def secondPartialDerivativeX(f, x, y, z = None):
    # apply limit definition of second derivative
    decimal = 12
    h = 0.00001
    try:
        fx = eval(f)
        x += h
        fx_plus = eval(f)
        x -= h*2
        fx_minus = eval(f)
        return round(float("%.4g"%((fx_plus - 2*fx + fx_minus)/h**2)), decimal)
        # 4 sf
    except: return None

def secondPartialDerivativeY(f, x, y, z = None):
    # apply limit definition of second derivative
    decimal = 12
    h = 0.00001
    try:
        fx = eval(f)
        y += h
        fx_plus = eval(f)
        y -= h*2
        fx_minus = eval(f)
        return round(float("%.4g"%((fx_plus - 2*fx + fx_minus)/h**2)), decimal)
        # 4 sf
    except: return None

def secondPartialDerivativeZ(f, x, y, z):
    # apply limit definition of second derivative
    decimal = 12
    h = 0.00001
    try:
        fx = eval(f)
        z += h
        fx_plus = eval(f)
        z -= h*2
        fx_minus = eval(f)
        return round(float("%.4g"%((fx_plus - 2*fx + fx_minus)/h**2)), decimal)
        # 4 sf
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
    decimal = 3
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
        return round(x, decimal), round(y, decimal), round(eval(f), decimal)
    except: return None

def maximize3D(f, point):
    # gradient ascent algorithm
    x, y = point
    gamma = 0.01
    epsilon = 0.0001
    decimal = 3
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
        return round(x, decimal), round(y, decimal), round(eval(f), decimal)
    except: return None

def minimize(f, x):
    # gradient descent algorithm on one variable
    gamma = 0.01
    epsilon = 0.00001
    decimal = 4
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
        return round(x, decimal), round(eval(f), decimal)
    except: return None

def maximize(f, x):
    # gradient ascent algorithm on one variable
    gamma = 0.01
    epsilon = 0.00001
    decimal = 4
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
        return round(x, decimal), round(eval(f), decimal)
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
    return math.e**(-0.5*((x - mean)/stdev)**2)/stdev/sqrt(2*math.pi)

def inverseNormal(probability, mean, stdev):
    # use binary search divide and conquer algorithm
    if probability >= 1 or probability <= 0 or stdev <= 0: return None
    factor = 100
    lowerGuess = mean - factor*stdev
    upperGuess = mean + factor*stdev
    epsilon = 0.00000001 # degree of accuracy
    decimal = 6
    if probability < epsilon: return "negative infinity"
    elif (1 - probability) < epsilon: return "infinity"
    guess = (upperGuess + lowerGuess)/2
    # binary search while error is greater than epsilon
    while (abs(normalDistribution(mean, stdev, mean - factor*stdev, guess)
               - probability) > epsilon):
        if normalDistribution(mean,stdev,mean-factor*stdev,guess)>probability:
            upperGuess = guess
        else: lowerGuess = guess
        guess = (upperGuess + lowerGuess)/2
    return round(guess, decimal)

def tPDF(value, degreesOfFreedom):
    # student's t probability density function formula
    v = degreesOfFreedom
    try: return (math.gamma(v/2+0.5)/math.sqrt(v*math.pi)/math.gamma(v/2)*
            (1 + value**2/v)**(-v/2-0.5))
    except:
        try: return (math.e**(math.lgamma(v/2+0.5) - math.lgamma(v/2))/
            math.sqrt(v*math.pi)*(1 + value**2/v)**(-v/2-0.5))
        except: return None
        
def tDistribution(lower, upper, degreesOfFreedom):
    # student's t distribution is just an integral. Use integration
    maxBound = 30
    if lower < -maxBound: lower = -maxBound
    if upper > maxBound: upper = maxBound
    try:
        return integral("tPDF(x, %d)" % degreesOfFreedom, lower, upper)
    except: return None

def inverseT(probability, degreesOfFreedom):
    # use binary search divide and conquer algorithm
    if probability >= 1 or probability <= 0: return None
    factor = 200
    lowerGuess = -factor
    upperGuess = factor
    epsilon = 0.00001 # degree of accuracy
    decimal = 4
    if probability < epsilon: return "negative infinity"
    elif (1 - probability) < epsilon: return "infinity"
    guess = (upperGuess + lowerGuess)/2
    # binary search while error is greater than epsilon
    while (abs(tDistribution(-factor, guess, degreesOfFreedom)
               - probability) > epsilon):
        if tDistribution(-factor, guess, degreesOfFreedom)>probability:
            upperGuess = guess
        else: lowerGuess = guess
        guess = (upperGuess + lowerGuess)/2
    return round(guess, decimal)

####################################
# discrete mathematics
####################################

### edit math module to prevent crashing

def factorial(n):
    try: return math.factorial(n)
    except: return None

def nPr(n, r):
    # permutation formula
    try: return round(math.factorial(n)/math.factorial(n-r))
    except: return None

def nCr(n, r):
    # combination formula
    try: return round(math.factorial(n)/math.factorial(r)/math.factorial(n-r))
    except: return None

def series(expression, start, end):
    # explicit series f(i) from i = start to i = end
    try:
        currentSum = 0
        for i in range(start, end + 1):
            currentSum += eval(expression)
        return currentSum
    except: return None

def sequenceRecursive(expression, initial, iterations):
    # recursive sequence formula f(i) (i = previous term)
    try:
        i = initial
        result = [i]
        for value in range(iterations):
            i = eval(expression)
            result.append(i)
        return result
    except: return None

def sequenceExplicit(expression, start, end):
    # explicit sequence formula f(i) from i = start to i = end
    try:
        result = []
        for i in range(start, end + 1):
            result.append(eval(expression))
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
    try: return expected**value*math.e**(-expected)/math.factorial(value)
    except: return None

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
    try: newL = [(x, math.log(y)) for (x, y) in L]
    except: return None
    b, a, R = linearRegression(newL)
    return math.e**a, b, R

def logarithmicRegression(L):
    try:
        newL = [(math.log(x), y) for (x, y) in L]
        return linearRegression(newL)
    except: return None

def powerRegression(L):
    try:
        newL = [(math.log(x), math.log(y)) for (x, y) in L]
        b, a, R = linearRegression(newL)
        return math.e**a, b, R
    except: return None

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
    decimal = 4
    z = inverseNormal((confidence + 1)/2, 0, 1)
    return (round(mean - z*stdev/math.sqrt(n), decimal),
            round(mean + z*stdev/math.sqrt(n), decimal))

def zIntervalData(L, stdev, confidence):
    # finds the z interval given a list
    mean = avg(L)
    n = len(L)
    return zIntervalStats(mean, stdev, n, confidence)

def tIntervalStats(mean, stdev, n, confidence):
    # finds the t interval given aggregate stats
    decimal = 3
    t = inverseT((confidence + 1)/2, n - 1)
    return (round(mean - t*stdev/math.sqrt(n), decimal),
            round(mean + t*stdev/math.sqrt(n), decimal))

def tIntervalData(L, confidence):
    # finds the t interval given a list
    mean = avg(L)
    stdev = standardDeviation(L)
    n = len(L)
    return tIntervalStats(mean, stdev, n, confidence)

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

class Function(object):
    def __init__(self, fxn):
        self.fxn = fxn
        self.vectors = []
        self.color = "#%02x%02x%02x" % (randomColor())

    def evaluate(self, x = None, y = None, z = None, r = None, theta = None,
                 rho = None, phi = None, t = None, Dy = None):
        # evaluate a function (string) at a given point
        try:
            answer = eval(self.fxn)
            assert(type(answer) == int or type(answer) == float
                   or type(answer) == tuple)
            assert(answer != None)
            self.isValid = True
            return answer
        except:
            self.isValid = False
            return None

    def __repr__(self):
        return self.fxn

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
                    z, maxZ = vector[2], data.maxZ
                    color = "#%02x%02x%02x" % colorRGB(z, maxZ)
                    drawLine3D(prevVector, vector, canvas, data, color)
                    prevVector = vector
            # loop through another direction creating parallel line segments
            for col in range(len(self.vectors[0])):
                if len(self.vectors) == 1: break
                prevVector = self.vectors[0][col]
                for row in range(len(self.vectors) - 1):
                    vector = self.vectors[row+1][col]
                    if prevVector == None or vector == None:
                        prevVector = vector; continue
                    z, maxZ = vector[2], data.maxZ
                    color = "#%02x%02x%02x" % colorRGB(z, maxZ)
                    drawLine3D(prevVector, vector, canvas, data, color)
                    prevVector = vector
            self.dialogBar.draw(canvas, data)
        except: pass

def drawLine3D(prevVector, vector, canvas, data, color):
    # draws a line between two points given 2 vectors in 3D
    if (prevVector != None and vector != None and
        checkBounds(prevVector, vector, data)):
        prevPoint = vectorToPoint(prevVector, data)
        point = vectorToPoint(vector, data)
        canvas.create_line(prevPoint, point, fill = color)

class Cartesian3D(ThreeD):
     
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

    def __init__(self, fxn):
        super().__init__(fxn)
        self.isValid = self.checkValidity()

    def checkValidity(self):
        # check if input is valid. If so, get min t and max t
        errorMessage = """Please type to enter a valid function
in the form 'x(t), y(t), z(t), t_min, t_max'"""
        functionLength = 5
        split = self.fxn.split(",")
        try:
            assert(len(split) == functionLength)
            self.minT = eval(split[-2])
            self.maxT = eval(split[-1])
            for num in (self.minT, self.maxT):
                assert(type(num) == int or type(num) == float)
            return True
        except:
            self.errorMessage = errorMessage
            return False

    def generateVectors(self, data):
        # generate a 2D list of vectors given a parametric function
        self.isValid = self.checkValidity()
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
        # evaluate a parametric function (string) at a given t
        try:
            split = self.fxn.split(",")
            x = eval(split[0])
            y = eval(split[1])
            z = eval(split[2])
            for answer in (x, y, z):
                assert(type(answer) == int or type(answer) == float)
                assert(answer != None)
            self.isValid = True
            return x, y, z
        except:
            self.isValid = False
            return None

    def __repr__(self):
        return str(self.fxn)

class CylindricalRDependent(ThreeD):
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
                r = self.evaluate(z = z, theta = theta)
                if r != None:
                    try: vector = cylindricalToCartesian(r, theta, z)
                    except: vector = None
                else: vector = None
                row.append(vector)
            self.vectors.append(row)

class CylindricalZDependent(ThreeD):
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
                    vector = self.evaluate(x = x, y = y, z = z, t = t)
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

    def __init__(self, fxn):
        super().__init__(fxn)
        self.isValid = self.checkValidity()

    def checkValidity(self):
        # checks that an input is valid. If so, find min t and max t
        errorMessage = """Please type to enter a valid function
in the form 'x(t), y(t), t_min, t_max'"""
        functionLength = 4
        split = self.fxn.split(",")
        try:
            assert(len(split) == functionLength)
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
        self.isValid = self.checkValidity()
        if not self.isValid: return
        self.vectors = []
        stepT = (self.maxT - self.minT)/data.incrementCurve
        for i in range(data.increment + 1):
            t = self.minT + i*stepT
            vector = self.evaluate(t)
            self.vectors.append(vector)

    def evaluate(self, t):
        # evaluate the function at a point t
        try:
            split = self.fxn.split(",")
            x = eval(split[0])
            y = eval(split[1])
            for answer in (x, y):
                assert(type(answer) == int or type(answer) == float)
                assert(answer != None)
            self.isValid = True
            return x, y
        except:
            self.isValid = False
            return None

    def __repr__(self):
        return str(self.fxn)

class VectorField2D(TwoD):

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
                vector = self.evaluate(x = x, y = y)
                try:
                    vector = (vector[0]/data.incrementField2D,
                              vector[1]/data.incrementField2D)
                    head = tail[0]+vector[0], tail[1]+vector[1]
                    self.vectors.append((head, tail))
                except: pass

    def draw(self, canvas, data):
        # draws vectors in vector field
        scaleX, scaleY = data.scaleX, data.scaleY
        for head, tail in self.vectors:
            tailPoint=(data.originX+scaleX*tail[0],data.originY-scaleY*tail[1])
            headPoint=(data.originX+scaleX*head[0],data.originY-scaleY*head[1])
            canvas.create_line(tailPoint, headPoint, fill = self.color)

class Point(TwoD):

    def __init__(self, x, y):
        # creates a point that can be drawn on a canvas
        self.x = x
        self.y = y
        self.radius = 2

    def draw(self, canvas, data):
        # draws point on canvas
        scaleX, scaleY = data.scaleX, data.scaleY
        pointX, pointY = data.originX+scaleX*self.x,data.originY-scaleY*self.y
        x, y = self.x, self.y
        # check if point is on the plane
        if inRange(x, data.minX, data.maxX) and inRange(y,data.minY,data.maxY):
            canvas.create_oval(pointX - self.radius, pointY - self.radius,
                               pointX + self.radius, pointY + self.radius,
                               fill = "black")
            canvas.create_text(pointX, pointY, anchor = NW,
                               text = "(%s, %s)" % (self.x, self.y))

####################################
# differential equations
####################################

class DiffEq(VectorField2D):

    def __init__(self, fxn):
        super().__init__(fxn)
        self.solve = False
            
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
                elif type(slope) == tuple:
                    slope = slope[0]
                    if type(slope) != float and type(slope) != int: return
                angle = math.atan(slope*stepX/stepY)
                head = (tail[0]+math.cos(angle)*stepX*data.scaleFactor,
                        tail[1]+math.sin(angle)*stepY*data.scaleFactor)
                self.vectors.append((head, tail))
        self.getSolution(data)

    def getSolution(self, data):
        # gets initial condition of solution
        split = self.fxn.split(",")
        if len(split) - 1 == 2:
            try:
                self.x0 = eval(split[1])
                self.y0 = eval(split[2])
                self.solve = True
                self.solveDiffEq(data)
            except: self.solve = False
        else: self.solve = False

    def solveDiffEq(self, data):
        # solves diff eq using Runge-Kutta 4 method
        if not self.solve: return
        self.solution = []
        factor = 4
        step = (data.maxX - data.minX)/data.stepFactor/factor
        x, y = self.x0, self.y0
        xNext, yNext = x, y
        prevSlope = 0
        self.solution = [(x, y)]
        self.goForward(data, x, y, xNext, yNext, step, prevSlope, factor)
        self.goBackward(data, x, y, xNext, yNext, step, prevSlope, factor)

    def goForward(self, data, x, y, xNext, yNext, step, prevSlope, factor):
        # go forward from initial condition
        while checkBounds2D((x,y),(xNext,yNext),data):
            for i in range(factor):
                slope1 = self.evaluate(x = x, y = y)[0]
                slope2 = self.evaluate(x = x + step/2, y = y + slope1*step/2)[0]
                slope3 = self.evaluate(x = x + step/2, y = y + slope2*step/2)[0]
                slope4 = self.evaluate(x = x + step, y = y + slope3*step)[0]
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
                slope1 = self.evaluate(x = x, y = y)[0]
                slope2 = self.evaluate(x = x - step/2, y = y - slope1*step/2)[0]
                slope3 = self.evaluate(x = x - step/2, y = y - slope2*step/2)[0]
                slope4 = self.evaluate(x = x - step, y = y - slope3*step)[0]
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
        if self.solve and self.vectors != []:
            prevVector = self.solution[0]
            for i in range(len(self.solution) - 1):
                vector = self.solution[i + 1]
                drawLine2D(prevVector, vector, canvas, data, self.color)
                prevVector = vector            

class SecondOrderDE(TwoD):

    def __init__(self, fxn):
        super().__init__(fxn)
        self.getInitialValues()
        self.vectors = []

    def generateVectors(self, data):
        # solves 2nd order initial value problem with Runge Kutta 4 algorithm
        self.getInitialValues()
        if not self.solve:
            self.vectors = []
            return
        x, y, yPrime = self.x0, self.y0, self.yPrime0
        factor = 4
        step = (data.maxX - data.minX)/data.stepFactor/factor
        self.vectors = [(x, y)]
        self.goForward(x, y, yPrime, data, factor, step)
        self.goBackwards(x, y, yPrime, data, factor, step)

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
        DDy1 = self.evaluate(x = x, y = y, Dy = yPrime)[0]*step
        x += step/2
        y2 = y + Dy1/2; yPrime2 = yPrime + DDy1/2
        Dy2 = yPrime2*step
        DDy2 = self.evaluate(x = x, y = y2, Dy = yPrime2)[0]*step
        y3 = y + Dy2/2; yPrime3 = yPrime + DDy2/2
        Dy3 = yPrime3*step
        DDy3 = self.evaluate(x = x, y = y3, Dy = yPrime3)[0]*step
        x += step/2
        y4 = y + Dy3; yPrime4 = yPrime + DDy3
        Dy4 = yPrime4*step
        DDy4 = self.evaluate(x = x, y = y4, Dy = yPrime4)[0]*step
        num = 6
        Dy = (Dy1 + 2*Dy2 + 2*Dy3 + Dy4)/num
        DDy = (DDy1 + 2*DDy2 + 2*DDy3 + DDy4)/num
        return Dy, DDy

    def goBackwardRK4(self, x, y, yPrime, step):
        # find delta y and delta squared y using Runge Kutta 4 algorithm
        Dy1 = yPrime*step
        DDy1 = self.evaluate(x = x, y = y, Dy = yPrime)[0]*step
        x -= step/2
        y2 = y - Dy1/2; yPrime2 = yPrime - DDy1/2
        Dy2 = yPrime2*step
        DDy2 = self.evaluate(x = x, y = y2, Dy = yPrime2)[0]*step
        y3 = y - Dy2/2; yPrime3 = yPrime - DDy2/2
        Dy3 = yPrime3*step
        DDy3 = self.evaluate(x = x, y = y3, Dy = yPrime3)[0]*step
        x -= step/2
        y4 = y - Dy3; yPrime4 = yPrime - DDy3
        Dy4 = yPrime4*step
        DDy4 = self.evaluate(x = x, y = y4, Dy = yPrime4)[0]*step
        num = 6
        Dy = (Dy1 + 2*Dy2 + 2*Dy3 + Dy4)/num
        DDy = (DDy1 + 2*DDy2 + 2*DDy3 + DDy4)/num
        return Dy, DDy

    def getInitialValues(self):
        # gets initial condition of solution
        split = self.fxn.split(",")
        conditions = 3
        try:
            self.x0, self.y0, self.yPrime0 = tuple(
                [eval(element) for element in split[-conditions:]])
            self.solve = True
        except: self.solve = False

####################################
# physics module
####################################

class ForceField(VectorField3D):

    def __init__(self, fxn):
        super().__init__(fxn)
        try: self.i, self.j, self.k = tuple(self.fxn.split(","))
        except: pass

class Particle(ThreeD):

    def __init__(self):
        super().__init__("")
        self.electricField = ForceField("0,0,0") # units N/C
        self.magneticField = ForceField("0,0,0") # units T
        self.gravitationalField = ForceField("0,0,0") # units N/kg
        self.forceField = ForceField("0,0,0") # units N
        self.charge = 1e-9 # 1 nanocoulomb
        self.mass = 1 # 1 kg
        self.x0, self.y0, self.z0 = 0, 0, 0 # start at origin
        self.xPrime0, self.yPrime0, self.zPrime0 = 0, 0, 0 # start at rest
        self.minT = 0
        self.maxT = 10 # timescale: 0 to 10 seconds
        self.solve = True

    def generateVectors(self, data):
        # simulate the motion of the particle
        if not self.solve:
            self.vectors = [[]]
            return
        self.vectors = []
        step = (self.maxT - self.minT)/data.incrementCurve
        row = []
        t = self.minT
        x, y, z = self.x0, self.y0, self.z0
        xPrime, yPrime, zPrime = self.xPrime0, self.yPrime0, self.zPrime0
        q, m = self.charge, self.mass
        for i in range(data.incrementCurve + 1): # increase time in small steps
            Dx, Dy, Dz, DDx, DDy, DDz = self.RungeKutta4(
                q, m, t, x, y, z, xPrime, yPrime, zPrime, step)
            t += step # change time, position, velocity
            x += Dx; y += Dy; z += Dz
            xPrime += DDx; yPrime += DDy; zPrime += DDz
            vector = x, y, z
            row.append(vector)
        self.vectors.append(row)

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
        num = 6
        Dx = (Dx1 + 2*Dx2 + 2*Dx3 + Dx4)/num
        Dy = (Dy1 + 2*Dy2 + 2*Dy3 + Dy4)/num
        Dz = (Dz1 + 2*Dz2 + 2*Dz3 + Dz4)/num
        DDx = (DDx1 + 2*DDx2 + 2*DDx3 + DDx4)/num
        DDy = (DDy1 + 2*DDy2 + 2*DDy3 + DDy4)/num
        DDz = (DDz1 + 2*DDz2 + 2*DDz3 + DDz4)/num
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
    #print(data.dialogFunctions.edit)
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
        data.fxnNum = len(data.functions[data.mode]) - 1
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
            data.fxnNum = len(data.functions[data.mode]) - 1
        elif data.mode == "MATHLAB 2D":
            addNewFunction2D(data, option)
            data.fxnNum = len(data.functions[data.mode]) - 1
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
    elif option == "Binomial Distribution":
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
    else: inputs = executeCalculator5(data, option, inputs)
    return inputs

def executeCalculator5(data, option, inputs):
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
        data.bestFit = Cartesian2DyDep("%.4f*e**(%.4f*x)" % (a, b))
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

def executeStatisticsDialog(data, option):
    # open dialog box for whatever is selected by the user in the menu
    if option == "z-interval: Statistics":
        inputs = ["Enter observed mean:","Known standard deviation:",
                  "Enter sample size:","Enter confidence level:"]
    elif option.startswith("z-interval: Column"):
        inputs = ["Known standard deviation:","Enter confidence level:"]
    elif option == "t-interval: Statistics":
        inputs = ["Enter observed mean:","Sample standard deviation:",
                  "Enter sample size:","Enter confidence level:"]
    elif option.startswith("t-interval: Column"):
        inputs = ["Enter confidence level:"]
    elif option == "Polynomial Regression":
        inputs = ["Enter order:"]
    else: return
    data.dialogBox = DialogBox(option, inputs, data)

def executeInput(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Derivative at point":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'derivative("%s", %s)\n\n\n' % (f,x)
    elif data.inputTitle == "Second derivative":
        data.functions[data.mode][0] += (
            'secondDerivative("%s", %s)\n\n\n' % tuple(data.inputs))
    elif data.inputTitle == "Definite Integral":
        f, a, b = tuple(data.inputs)
        data.functions[data.mode][0] += 'integral("%s", %s, %s)\n\n\n' % (f,a,b)
    elif data.inputTitle == "Limit":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'limit("%s", %s)\n\n\n' % (f,x)
    else: executeInput1(data)
    evaluatePage(data)
    data.inputs = None
    data.exeMessage = ""

def executeInput1(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Solve for 0":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'zero("%s", %s)\n\n\n' % (f,x)
    elif data.inputTitle == "Local Minimum":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'minimize("%s", %s)\n\n\n' % (f,x)
    elif data.inputTitle == "Local Maximum":
        f, x = tuple(data.inputs)
        data.functions[data.mode][0] += 'maximize("%s", %s)\n\n\n' % (f,x)
    elif data.inputTitle == "Local Minimum 3D":
        f, x, y = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'minimize3D("%s", (%s, %s))\n\n\n' % (f,x,y))
    elif data.inputTitle == "Local Maximum 3D":
        f, x, y = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'maximize3D("%s", (%s, %s))\n\n\n' % (f,x,y))
    else: executeInput2(data)

def executeInput2(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Definite Double Integral":
        num = 4
        xMin, xMax, yMin, yMax = tuple(data.inputs[-num:])
        for i in range(num): data.inputs.pop()
        f = data.inputs.pop(0)
        constraints = data.inputs
        while "" in constraints: constraints.remove("")
        data.functions[data.mode][0] += (
            'doubleIntegral("%s", %s, %s, %s, %s, %s)\n\n\n' % (
                f,str(constraints),xMin,xMax,yMin,yMax))
    else: executeInput3(data)

def executeInput3(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Definite Triple Integral":
        num = 6
        xMin, xMax, yMin, yMax, zMin, zMax = tuple(data.inputs[-num:])
        for i in range(num): data.inputs.pop()
        f = data.inputs.pop(0)
        constraints = data.inputs
        while "" in constraints: constraints.remove("")
        data.functions[data.mode][0] += (
            'tripleIntegral("%s", %s, %s, %s, %s, %s, %s, %s)\n\n\n' % (
                f,str(constraints),xMin,xMax,yMin,yMax,zMin,zMax))
    elif data.inputTitle == "Gradient at a point":
        f, x, y, z = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'gradient3("%s", %s, %s, %s)\n\n\n' % (f, x, y, z))
    elif data.inputTitle == "Curl at a point":
        P, Q, R, x, y, z = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'curl("%s", "%s", "%s", %s, %s, %s)\n\n\n' % (P, Q, R, x, y, z))
    else: executeInput4(data)

def executeInput4(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Divergence at a point":
        P, Q, R, x, y, z = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'divergence("%s", "%s", "%s", %s, %s, %s)\n\n\n' % (P, Q, R, x, y, z))
    elif data.inputTitle == "Laplacian at a point":
        f, x, y, z = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'laplacian("%s", %s, %s, %s)\n\n\n' % (f, x, y, z))
    elif data.inputTitle == "Normal Distribution":
        mean, stdev, lower, upper = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'normalDistribution(%s, %s, %s, %s)\n\n\n'%(mean,stdev,lower,upper))
    elif data.inputTitle == "Inverse Normal":
        probability, mean, stdev = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'inverseNormal(%s, %s, %s)\n\n\n'%(probability, mean, stdev))
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
    elif data.inputTitle == "Binomial Distribution":
        trials, probability, successes = tuple(data.inputs)
        data.functions[data.mode][0]+=('binomial(%s, %s, %s)\n\n\n'%
                                         (trials, probability, successes))
    elif data.inputTitle == "Negative Binomial":
        successes, probability, trials = tuple(data.inputs)
        data.functions[data.mode][0]+=('negativeBinomial(%s, %s, %s)\n\n\n'%
                                         (successes, probability, trials))
    else: executeInput6(data)

def executeInput6(data):
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
    else: executeInput7(data)

def executeInput7(data):
    # executes whatever the user inputs into the dialog box
    if data.inputTitle == "Combination (nCr)":
        n, r = tuple(data.inputs)
        data.functions[data.mode][0] += 'nCr(%s, %s)\n\n\n' % (n, r)
    elif data.inputTitle == "Series":
        expression, start, end = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'series("%s", %s, %s)\n\n\n' % (expression, start, end))
    elif data.inputTitle == "Sequence: Explicit":
        expression, start, end = tuple(data.inputs)
        data.functions[data.mode][0] += (
            'sequenceExplicit("%s", %s, %s)\n\n\n' % (expression, start, end))
    elif data.inputTitle == "Sequence: Recursive":
        expression, initial, iterations = tuple(data.inputs)
        data.functions[data.mode][0]+=('sequenceRecursive("%s", %s, %s)\n\n\n'%
                                         (expression, initial, iterations))

def executeInputStats(data):
    # executes whatever the user inputs into the dialog box
    dec = 6
    if data.inputTitle == "z-interval: Statistics":
        executeInputStatsZStats(data, dec)
    elif data.inputTitle.startswith("z-interval: Column"):
        executeInputStatsZData(data, dec)
    elif data.inputTitle == "t-interval: Statistics":
        executeInputStatsTStats(data, dec)
    elif data.inputTitle.startswith("t-interval: Column"):
        executeInputStatsTData(data, dec)
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

def executeInputStatsZStats(data, dec):
    # executes input for z-interval stats
    mean, stdev, n, confidence = tuple(data.inputs)
    output = zIntervalStats(eval(mean),eval(stdev),eval(n),eval(confidence))
    data.inputMessage = """Confidence z-Interval\n\nObserved Mean:
Population Standard Deviation:\nSample Size:\nConfidence level:\n
Confidence Interval:"""
    data.outputMessage = "Statistics\n\n%s\n\n%s" % (
        "\n".join([str(element) for element in [mean,stdev,n,confidence]]),
        str(output))

def executeInputStatsZData(data, dec):
    # executes input for z-interval data
    col = int(data.inputTitle[-1]) - 1
    column = buildColumn(data, col)
    stdev, confidence = tuple(data.inputs)
    output = zIntervalData(column, eval(stdev), eval(confidence))
    data.inputMessage = """Confidence z-Interval\n\nObserved Mean:
Population Standard Deviation:\nSample Size:\nConfidence level:\n
Confidence Interval:"""
    data.outputMessage = "Column %s\n\n%s\n\n%s" % (data.inputTitle[-1],
        "\n".join([str(element) for element in
        [round(avg(column),dec),round(eval(stdev),dec),len(column),confidence]]),
        str(output))

def executeInputStatsTStats(data, dec):
    # executes input for t-interval stats
    mean, stdev, n, confidence = tuple(data.inputs)
    output = tIntervalStats(eval(mean),eval(stdev),eval(n),eval(confidence))
    data.inputMessage = """Confidence t-Interval\n\nObserved Mean:
Sample Standard Deviation:\nSample Size:\nConfidence level:\n
Confidence Interval:"""
    data.outputMessage = "Statistics\n\n%s\n\n%s" % (
        "\n".join([str(element) for element in [mean,stdev,n,confidence]]),
        str(output))

def executeInputStatsTData(data, dec):
    # executes input for t-interval data
    col = int(data.inputTitle[-1]) - 1
    column = buildColumn(data, col)
    confidence = data.inputs[0]
    output = tIntervalData(column, eval(confidence))
    data.inputMessage = """Confidence t-Interval\n\nObserved Mean:
Sample Standard Deviation:\nSample Size:\nConfidence level:\n
Confidence Interval:"""
    data.outputMessage = "Column %s\n\n%s\n\n%s" % (data.inputTitle[-1],
        "\n".join([str(element) for element in
        [round(avg(column), dec), round(standardDeviation(column), dec),
         len(column), confidence]]), str(output))

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
        particle.solve = True
        particle.generateVectors(data)
    except: particle.solve = False

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
    skiplines = 3; lines = data.functions[data.mode][data.fxnNum].splitlines()
    for lineIndex in reversed(range(len(lines))):
        line = lines[lineIndex]
        try:
            if line == "": assert(lineIndex % skiplines != 0); continue
            answer = eval(line)
            assert(answer != None)
            assert(type(answer) in (int, float, str, tuple, list, bool))
            decimal = 12
            if type(answer) == float: answer = round(answer, decimal)
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
    line = line.replace("integral", "")
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
    line = line.replace("derivative", "")
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
    x = data.x
    y0 = eval(data.fxn)
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
    line = line.replace("zero", "")
    (fxn, guess) = eval(line)
    generatePointPicture(data, fxn, x, y)

def generateMaxMinPicture(data, line):
    # generates a temporary struct to lay the groundwork to draw max min point
    x, y = eval(line)
    line = line.replace("minimize", "")
    line = line.replace("maximize", "")
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
    data.fxnNum = 0
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
                        "Inverse Student's t", "Binomial Distribution",
                        "Negative Binomial", "Geometric Distribution",
                        "Hypergeometric", "Poisson Distribution"]
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
    data.stats = ["Statistics: Column 1",
                  "Statistics: Column 2", "Linear Regression",
                  "Exponential Regression", "Logarithmic Regression",
                  "Power Regression"]#, "Polynomial Regression"]
    data.tests = ["z-interval: Statistics", "z-interval: Column 1",
                  "z-interval: Column 2", "t-interval: Statistics",
                  "t-interval: Column 1", "t-interval: Column 2"]
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
    data.fxnNum = len(data.functions[data.mode]) - 1
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
    data.fxnNum = len(data.functions[data.mode]) - 1
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
    data.fxnNum = 0
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
            data.functions[mode][data.fxnNum] += "\n\n\n"
        return
    elif mode == "MATHLAB Statistics":
        row, col = data.cellSelected
        if row < data.rows - 1: row += 1
        data.cellSelected = row, col
        return
    elif mode == "MATHLAB 3D": addNewFunction3D(data, event.keysym)
    elif mode == "MATHLAB 2D": addNewFunction2D(data, event.keysym)
    data.fxnNum += 1

def backspacePressed(event, data):
    # what to do when you press backspace
    mode = data.mode
    num = data.fxnNum
    if data.mode == "MATHLAB Calculator":
        if (data.functions[mode][0] == "" or
            data.functions[mode][0][-1] == "\n"): return
        data.functions[mode][num]=data.functions[mode][num][:-1]
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
    elif key == "asterisk": data.functions[mode][data.fxnNum]+="*"
    elif key == "slash": data.functions[mode][data.fxnNum] += "/"
    elif key == "plus": data.functions[mode][data.fxnNum] += "+"
    elif key == "minus": data.functions[mode][data.fxnNum] += "-"
    elif key == "asciicircum": data.functions[mode][data.fxnNum]+="**"
    elif key == "parenleft": data.functions[mode][data.fxnNum]+="("
    elif key == "period": data.functions[mode][data.fxnNum] += "."
    elif key == "parenright":data.functions[mode][data.fxnNum]+=")"
    elif key == "comma": data.functions[mode][data.fxnNum] += ","
    elif key == "space": data.functions[mode][data.fxnNum] += " "
    elif key == "quotedbl": data.functions[mode][data.fxnNum] += '"'
    elif key == "quoteright": data.functions[mode][data.fxnNum] += "'"
    elif key == "less": data.functions[mode][data.fxnNum] += "<"
    elif key == "greater": data.functions[mode][data.fxnNum] += ">"
    elif key == "percent": data.functions[mode][data.fxnNum] += "%"
    elif len(key) == 1: data.functions[mode][data.fxnNum]+=key

def modifyFunctionGraph(data, key):
    # modify the function text
    mode = data.mode
    fxnNum = data.dialogFunctions.edit
    if key == "asterisk": data.functions[mode][fxnNum].fxn += "*"
    elif key == "slash": data.functions[mode][fxnNum].fxn += "/"
    elif key == "plus": data.functions[mode][fxnNum].fxn += "+"
    elif key == "minus": data.functions[mode][fxnNum].fxn += "-"
    elif key == "asciicircum": data.functions[mode][fxnNum].fxn += "**"
    elif key == "parenleft": data.functions[mode][fxnNum].fxn+="("
    elif key == "period": data.functions[mode][fxnNum].fxn += "."
    elif key == "parenright":data.functions[mode][fxnNum].fxn+=")"
    elif key == "comma": data.functions[mode][fxnNum].fxn += ","
    elif key == "space": data.functions[mode][fxnNum].fxn += " "
    elif key == "quotedbl": data.functions[mode][fxnNum].fxn += '"'
    elif key == "quoteright": data.functions[mode][fxnNum].fxn += "'"
    elif len(key) == 1: data.functions[mode][fxnNum].fxn += key
    elif key != "BackSpace": return
    data.functions[mode][fxnNum].generateVectors(data)

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
    root.title("MATHLAB 1.04 - Xiong-Fei Du")
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
    assert(equal(derivative("x**2", 6), 12))
    assert(equal(derivative("x**3", 6), 108))
    assert(equal(derivative("sin(x)", 0), 1))

def testIntegral():
    assert(equal(integral("x**3", 0, 0), 0))
    assert(equal(integral("x**2", 0, 1), 0.33333))
    assert(equal(integral("sin(x)", 0, 2*math.pi), 0))

def testZero():
    assert(equal(zero("x**3 - 1", 2), 1))
    assert(equal(zero("x**2 - 16", 91), 4))
    assert(equal(zero("sin(x)", 1), 0))

def testAll():
    testDerivative()
    testIntegral()
    testZero()
    print("Passed!")

#testAll()

run(width = 800, height = 800)
