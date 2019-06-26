import operator, math
import numpy as np
from numpy.polynomial.polynomial import Polynomial as NumpyPolynomial

class Polynomial(NumpyPolynomial):
    '''
    Coefficients [1,2,3] will correspond to 
    p = 1 + 2x + 3x^2
    '''

    def roots(self,method = "companion_laguerre", eps = 1e-16, r0 = None, unique = False):
        if unique:
            # Remove multiplicity
            g = gcd(self, self.deriv())
            if g.degree() >= 1:
                p = self // g
            else:
                # There are no multiple roots so use the original polynomial.
                p = self

        if method == "companion":
            return self.roots_companion()
        if method == "laguerre":
            return self.roots_laguerre(eps,r0)
        if method == "companion_laguerre":
            return self.roots_companion_laguerre(eps,r0)

    def roots_companion(self):
        return  sortRoots(super(Polynomial, self).roots())

    def roots_companion_laguerre(self, eps=1e-16, r0=None):
        '''Find all roots of the polynomial using a diagonalization of the companion matrix,
        and polish the roots with Laguerre's method
        '''

        roots = self.roots(method = "companion", unique = True)

        # polish using laguerre
        p, dp, ddp = _laguerreInputs(self, eps)
        for i,r in enumerate(roots):
            r = _improveRoot(r, p, dp, ddp, eps)[1]
            roots[i] = r

        return sortRoots(roots)

    def roots_laguerre(self, eps=1e-16, r0=None):
        '''Find all roots of the polynomial using Laguerre's method.
        The algorithm is described in Numerical Recipes
        '''

        # If no initial guess is provided, 
        if r0 == None:
            r0 = self.roots_companion()[0]

        p, dp, ddp = _laguerreInputs(self, eps)

        roots = []
        q = p
        dq = dp
        ddq = ddp

        while q.degree() > 0:
            r = _improveRoot(r0, q, dq, ddq, eps)[1]
            roots.append(r)

            # Deflate the polynomial
            q = q//Polynomial([-r,1])
            dq = q.deriv()
            ddq = dq.deriv()

            # Try the root just found as the next guess.  This rapidly
            # eliminates multiple roots and sets up to find the next closest
            # root.
            r0 = r

        # Polish the roots since the reduction of the polynomial can lead to
        # accumulation of errors.
        for i,r in enumerate(roots):
            r = _improveRoot(r, p, dp, ddp, eps)[1]
            
            if r.imag == 0:
                r = r.real
            roots[i] = r

        return sortRoots(roots)

def sortRoots(roots):
    '''Return a new list of roots sorted first on the real part then the imag.
    '''
    lst = []
    for i,r in enumerate(roots):
        try:
            lst.append((r.real, r.imag, i, r))
        except AttributeError:
            lst.append((r, 0, i, r))
    lst.sort()
    return [r[-1] for r in lst]

def improveRoots(roots, p, eps=1e-16):
    return [improveRoot(r, p, eps)[1] for r in roots]

def improveRoot(root, p, eps=1e-16):
    '''Use Laguerre's method to improve the estimate of polynomial root.

    Returns an estimate of the error and the current estimate of the root.  The
    error will be less than eps if the algorithm converges within the iteration
    limit.

    root - initial estimate of the root
    p - polynomial
    eps - error tolerance
    '''
    p, dp, ddp = _laguerreInputs(p, eps)
    
    z = root
    err, z = _improveRoot(z, p, dp, ddp, eps)
    
    if z.imag == 0:
        z = z.real
    return err, z

def _laguerreInputs(p, eps):
    if p.degree() < 1:
        raise ValueError('The degree of the polynomial must be at least one.')

    if not (0 < eps < 0.001):
        raise ValueError('Tolerance out of range')

    dp = p.deriv()
    ddp = dp.deriv()
    return p, dp, ddp

def _improveRoot(z, p, dp, ddp, eps):
    max_iterations = 100
    iterations = 0
    while iterations < max_iterations:
        try:
            a = _laguerreStep(z, p, dp, ddp)
        except ZeroDivisionError:
            # In this case the derivative is zero so take a step away from that
            # point.
            z = z+1
            a = _laguerreStep(z, p, dp, ddp)
        z = z - a
        iterations += 1

        # Compute the estimated distance from the root.
        err = abs(complex(a))
        if err < eps:
            break

    return err, z

def _laguerreStep(z, p, dp, ddp):
    '''Perform one iteration of Laguerre's method.

    z - current estimate of the root
    p - polynomial
    dp - first derivative of p
    ddp - second derivative of p
    '''
    pz = p(z)

    if pz == 0:
        # The result is zero so a zero has already been located.
        return pz

    dpz = dp(z)
    ddpz = ddp(z)

    G = dpz/pz
    H = G*G - ddpz/pz

    n = p.degree()
    s = np.sqrt(np.absolute((n-1)*(n*H - G*G)))
    d1 = G + s
    d2 = G - s

    if abs(complex(d1)) < abs(complex(d2)):
        return n/d2
    else:
        return n/d1

def gcd(u, v):
    '''gcd(u,v) returns the greatest common divisor of u and v.
    Where u and v are integers or polynomials.
    '''
    # If something other than an integer or a polynomial is fed in, the
    # algorithm may fail to converge.  In that case we want an exception
    # instead of an infinite loop.
    iterations = 0
    max_iterations = 1000
    while v != Polynomial(0):
        u, v = v, u % v
        iterations += 1
        if iterations >= max_iterations:
            raise ValueError('gcd(u,v) failed to converge')

    # return a monic polynomial
    return u/u.coef[-1]


class RationalFunction(object):
    '''Rational function data type that can mix with ordinary numbers and
    polynomials.
    
    Supports the usual arithmetic operations (+,-,*,/,**).
    '''

    # Width of the screen in characters.  Setting this to zero will suppress
    # the 3 line string representation.  If you only want it suppressed on
    # particular instances, set it on the instance instead of the class.
    screenWidth = 80

    def __init__(self, numer, denom=[1]):

        # casting
        if isinstance(numer, RationalFunction):
            ratfun = numer
            numer = ratfun.numer.coef
            denom = ratfun.denom.coef
        if isinstance(numer,Polynomial):
            numer = numer.coef
        if isinstance(denom,Polynomial):
            denom = denom.coef

        numer = Polynomial(numer)
        if numer == Polynomial(0):
            # The zero rational function is always represented as 0/1
            denom = 1
        denom = Polynomial(denom)

        if denom == Polynomial(0):
            raise ZeroDivisionError('Rational function denominator is zero')

        self.numer = numer
        self.denom = denom

        # self.simplify()
        # self.monic()


    def monic(self):
        c = self.denom.coef[-1]
        if c != 1:
            self.numer /= c
            self.denom /= c


    def simplify(self):
        # Remove any common factors.  Avoid numeric only common factors since
        # the rational function needs to be normalied below.
        if self.numer.degree() > 0 and self.denom.degree() > 0:
            d = gcd(self.numer, self.denom)
            if d.degree() > 0:
                self.numer //= d
                self.denom //= d

    def __eq__(self, other):
        return (self.numer == other.numer and self.denom == other.denom)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self

    def __neg__(self):
        return RationalFunction(-self.numer, self.denom)

    def __add__(self, other):
        other = RationalFunction(other)

        # Use algorithm for fractions from Knuth.
        un, ud = self.numer, self.denom
        vn, vd = other.numer, other.denom

        # If either polynomial is a constant, don't use gcd.
        if ud.degree() < 1 or vd.degree() < 1:
            d1 = 1
        else:
            d1 = gcd(ud, vd)
        if d1 == 1 or d1.degree() < 1:
            return RationalFunction(un*vd + ud*vn, ud*vd)

        t = un*(vd//d1) + vn*(ud//d1)
        d2 = gcd(t, d1)
        if d2.degree() < 1:
            return RationalFunction(t, (ud//d1)*vd)

        return RationalFunction(t//d2, (ud//d1)*(vd//d2))

    def __radd__(self, other):
        return RationalFunction(other).__add__(self)

    def __sub__(self, other):
        return self.__add__(-RationalFunction(other))

    def __rsub__(self, other):
        return RationalFunction(other).__add__(-self)

    def __mul__(self, other):
        other = RationalFunction(other)

        un, ud = self.numer, self.denom
        vn, vd = other.numer, other.denom

        # If any constant polynomials are involved, skip the gcd calculation.
        if un.degree() < 1 or ud.degree() < 1 or vn.degree() < 1 or vd.degree() < 1:
            return RationalFunction(un*vn, ud*vd)

        # Use algorithm for fractions from Knuth.
        d1 = gcd(un, vd)
        d2 = gcd(ud, vn)

        if d1.degree() < 1 and d2.degree() < 1:
            return RationalFunction(un*vn, ud*vd)

        if d1.degree() > 0 and d2.degree() < 1:
            return RationalFunction((un//d1)*vn, ud*(vd//d1))

        if d1.degree() < 1 and d2.degree() > 0:
            return RationalFunction(un*(vn//d2), (ud//d2)*vd)

        return RationalFunction((un//d1)*(vn//d2), (ud//d2)*(vd//d1))

    def __rmul__(self, other):
        return RationalFunction(other).__mul__(self)

    def __truediv__(self, other):
        other = RationalFunction(other)
        return self.__mul__(RationalFunction(other.denom, other.numer))

    def __rtruediv__(self, other):
        other = RationalFunction(other)
        return other.__mul__(RationalFunction(self.denom, self.numer))

    def __pow__(self, n):
        if not isinstance(n, int):
            return NotImplemented
        if n == 0:
            return RationalFunction(1)

        # Use the right-to-left binary method from Knuth.
        y = RationalFunction(1)
        if n < 0:
            z = 1//self
            n = -n
        else:
            z = self

        while True:
            odd = n & 1
            n = n >> 1
            if odd:
                y *= z
                if not n:
                    break
            z *= z
        return y

    def __call__(self, x):
        '''Evaluate the function at the point x.  Where x can also be a numpy
        array.  In that case the function is sampled at the points in x.
        '''
        a = self.numer(x)
        b = self.denom(x)

        # Otherwise the divide should work correctly.
        return a / b

    def __repr__(self):
        return 'RationalFunction(%r, %r)' % (self.numer, self.denom)

    def __str__(self):
        if self.denom.degree() == 0:
            # The rational function is a polynomial divided by one so just
            # display the polynomial.
            return str(self.numer)

        n = str(self.numer)
        d = str(self.denom)
        m = max(len(n), len(d), 3)
        if m < self.screenWidth:
            return '%s\n%s\n%s' % (n.center(m).rstrip(), '-'*m,
                                   d.center(m).rstrip())
        else:
            return '(%s) / (%s)' % (n, d)

    def deriv(self):
        '''Return the derivative of the rational function.
        '''
        a = RationalFunction(self.numer.deriv(), self.denom)
        b = RationalFunction(self.denom.deriv(), self.denom)
        return a - self*b