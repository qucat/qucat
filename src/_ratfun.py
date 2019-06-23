# Original copyright:

#-----------------------------------------------------------------------------
# Copyright (c) 2006  Raymond L. Buvel
#
# The ratfun module implements classes for handling polynomials and rational
# functions.
#
# The ratfun module is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
#
# The ratfun nodule is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# ratfun; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA
#-----------------------------------------------------------------------------


__all__ = '''Polynomial RationalFunction gcd polyFromRoots sortRoots
improveRoot improveRoots
'''.split()

import operator, math
import numpy as np
# ------ clnum features disabled for now
# from clnum import mpf, mpq, cmpf, sqrt

#-----------------------------------------------------------------------------
class Polynomial(object):
    '''Polynomial data type that can mix with ordinary numbers.
    
    Supports the usual arithmetic operations (+,-,*,/,%,**).

    Three class level variables can be used to change the format of a polynomial
    when it is converted to a string.  If you want to control the output of a
    single polynomial instance, simply set these variables on the instance.

    varName - this string holds the variable name (defaults to x).
    mulSymbol - symbol used to represent the multiplication of the coefficent
        and the variable.  Defaults to '*'.  Use an empty string if you like
        the operator implicit.
    coefFormat - format specifier used to display floating point values.  Set
        to None if the full string representation is desired.
    '''

    varName = 'x'
    mulSymbol = '*'
    coefFormat = '%.3g'

    def __new__(cls, *args):
        # Handle the zero polynomial.
        if not args:
            return Polynomial._zero

        # If the input is a polynomial, just return it.
        if len(args) == 1 and isinstance(args[0], cls):
            return args[0]

        # Eliminate strings so that generic sequence handling can be used.
        if isinstance(args[0], basestring):
            raise ValueError('Cannot create a polynomial from a string')

        if len(args) == 1:
            # If the input is a sequence, create the polynomial from it.
            try:
                # Note: this handles the case where the input is an iterator.
                seq = tuple(args[0])
            except TypeError:
                seq = None

            if seq is not None:
                coef = _coefFromSeq(seq)
            else:
                # Could be a constant polynomial. If conversion to complex
                # works, accept the value.  Otherwise, complex() will raise an
                # exception.
                val = args[0]
                complex(val)
                if not val:
                    return Polynomial._zero

                # ------ clnum features disabled for now
                # # Floats and complex do not mix with rationals so convert them
                # # to the corresponding clnum type.
                # if isinstance(val, float):
                #     val = mpf(val)
                # elif isinstance(val, complex):
                #     val = cmpf(val)

                coef = (val,)
        else:
            # The args tuple must now be a valid sequence of coefficients or
            # the construction fails with an exception.
            coef = _coefFromSeq(args)

        # Construct the new polynomial from the validated coefficient tuple.
        if not coef:
            return Polynomial._zero

        poly = object.__new__(cls)
        poly._coef = coef
        return poly


    def __eq__(self, other):
        return self._coef == Polynomial(other)._coef

    def __ne__(self, other):
        return self._coef != Polynomial(other)._coef

    # Note: the comparison is only provided to allow sorting lists of
    # polynomials.  Don't know what it would mean in other contexts.
    def __cmp__(self, other):
        other = Polynomial(other)
        n = len(self._coef)
        m = len(other._coef)

        # First, use the degree of the polynomials for comparison if they are
        # different.
        if n < m: return -1
        if n > m: return 1

        # The degree is the same so compute the difference polynomial
        p = self-other

        # The polynomials are equal if the difference is the zero polynomial.
        if not p._coef: return 0

        # The degree of the polynomials is the same but they are not equal.
        # Use the first non-zero coefficient of the highest power to do the
        # comparison.  Note that polynomial normalization makes sure this is
        # the last element in the coefficient tuple.
        c = p._coef[-1]

        # Note: this will raise an exception if complex numbers are invloved
        # since there is no ordering on complex numbers.
        try:
            return cmp(c, 0)
        except TypeError:
            # Make an arbitrary choice that allows a usable sort.  If the
            # coefficient is not complex, an exception will be raised and we
            # have to give up.
            if abs(c.real) >= abs(c.imag):
                return cmp(c.real, 0)
            else:
                return cmp(c.imag, 0)


    def __nonzero__(self):
        return bool(self._coef)


    def __pos__(self):
        return self

    def __neg__(self):
        return Polynomial(map(operator.neg, self._coef))


    def __add__(self, other):
        p0 = self._coef
        # Make polynomials interoperable with rational functions.
        try:
            p1 = Polynomial(other)._coef
        except TypeError:
            return NotImplemented
        n = len(p0)
        m = len(p1)
        if m != n:
            mx = max(n,m)
            if n != mx:
                p0 += (0,)*(mx-n)
            else:
                p1 += (0,)*(mx-m)
        return Polynomial(map(operator.add, p0, p1))

    def __radd__(self, other):
        return Polynomial(other).__add__(self)


    def __sub__(self, other):
        # Make polynomials interoperable with rational functions.
        try:
            other = Polynomial(other)
        except TypeError:
            return NotImplemented
        return self.__add__(-other)

    def __rsub__(self, other):
        return Polynomial(other).__add__(-self)


    def __mul__(self, other):
        p0 = self._coef
        # Make polynomials interoperable with rational functions.
        try:
            p1 = Polynomial(other)._coef
        except TypeError:
            return NotImplemented
        n = len(p0)
        m = len(p1)

        # Since polynomials are normalized, this handles all cases where one of
        # the polynomials is zero.
        if n == 0 or m == 0: return Polynomial._zero

        N = m+n-1 # Length of the result
        r = [0]*N
        # Iterate over the polynomial with the smallest degree.
        if n > m:
            n, m = m, n
            p0, p1 = p1, p0

        for i,a in enumerate(p0):
            if a == 0:
                # This short circuit improves performance for sparse
                # polynomials of high degree.
                continue
            r[i:i+m] = map(operator.add, r[i:i+m], [x*a for x in p1])
        return Polynomial(r)

    def __rmul__(self, other):
        return Polynomial(other).__mul__(self)


    def __divmod__(self, other):
        u = list(self._coef)
        v = Polynomial(other)._coef

        m = len(u)
        n = len(v)

        if n == 0: raise ZeroDivisionError('polynomial division by zero')

        # The following can only happen if u is the zero polynomial.
        if m == 0: return self, self

        # The degree of u is smaller than degree of v.
        if m < n: return Polynomial._zero, self

        # If the degree of v is 0, then it is a constant and must be handled as
        # a special case.
        if n == 1:
            a = _invert(v[0])
            return Polynomial([x*a for x in u]), Polynomial._zero

        # At this point both polynomials have degree at least one and
        # deg(u) >= deg(v).  Use the division algorithm from Knuth.
        m -= 1 # deg(u)
        n -= 1 # deg(v)
        k = m-n
        q = [None]*(k+1)
        vn = _invert(v[n])
        while k >= 0:
            qk = u[n+k]*vn
            q[k] = qk
            j = n+k-1
            while j >= k:
                u[j] -= qk*v[j-k]
                j -= 1
            k -= 1
        return Polynomial(q), Polynomial(u[:n])

    def __rdivmod__(self, other):
        return Polynomial(other).__divmod__(self)

    def __div__(self, other):
        # Make polynomials interoperable with rational functions.
        try:
            other = Polynomial(other)
        except TypeError:
            return NotImplemented
        return self.__divmod__(other)[0]

    def __rdiv__(self, other):
        return Polynomial(other).__divmod__(self)[0]

    def __mod__(self, other):
        return self.__divmod__(other)[1]

    def __rmod__(self, other):
        return Polynomial(other).__divmod__(self)[1]

    # Make true division work the same as regular division.  This will allow
    # the / operator to work with true division turned on.
    def __truediv__(self, other):
        # Make polynomials interoperable with rational functions.
        try:
            other = Polynomial(other)
        except TypeError:
            return NotImplemented
        return self.__divmod__(other)[0]

    def __rtruediv__(self, other):
        return Polynomial(other).__divmod__(self)[0]


    def __pow__(self, n):
        if not isinstance(n, (int,long)):
            return NotImplemented
        if n < 0:
            return NotImplemented
        if n == 0:
            return Polynomial(1)

        # Use the right-to-left binary method from Knuth.
        y = Polynomial(1)
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
        '''Evaluate the polynomial at the point x.  Where x can also be a numpy
        array.  In that case the polynomial is sampled at the points in x.
        '''
        p = self._coef
        n = len(p)-1

        y = 0*x  # produces an array of zeros if x is a numpy array.
        while n >= 0:
            y = y*x + p[n]
            n -= 1
        return y


    def __repr__(self):
        return 'Polynomial%r' % (self._coef,)

    def __str__(self):
        p = self._coef
        if len(p) == 0: return '0'

        lst = []
        for exp,coef in enumerate(p):
            if coef:
                lst.append(self._term_str(coef,exp))
        lst.reverse()
        if lst[0][0] == '+': lst[0] = lst[0][1:]
        return ''.join(lst)


    def _term_str(self, coef, exp):
        '''Return the string form of a single polynomial term.
        '''
        # Handle it this way to cover complex numbers.
        if coef == 1 or coef == -1:
            if coef == 1:
                sign = '+'
            else:
                sign = '-'

            if exp == 0:
                return '%s1' % sign
            elif exp == 1:
                return '%s%s' % (sign, self.varName)
            else:
                return '%s%s^%d' % (sign, self.varName, exp)

        if self.coefFormat:
            # If there is a user-supplied floating point format, use it.
            if isinstance(coef, (mpf,float)):
                s = self.coefFormat % coef

            elif isinstance(coef, (cmpf,complex)):
                s1 = self.coefFormat % coef.real
                s2 = self.coefFormat % coef.imag
                if s2.startswith('-'):
                    s = '(%s%sj)' % (s1,s2)
                else:
                    s = '(%s+%sj)' % (s1,s2)
            else:
                s = str(coef)
        else:
            s = str(coef)

        if not s.startswith('-'):
            s = '+' + s

        if exp == 0:
            return s
        elif exp == 1:
            return '%s%s%s' % (s, self.mulSymbol, self.varName)
        else:
            return '%s%s%s^%d' % (s, self.mulSymbol, self.varName, exp)


    def coef(self, n):
        '''Return the coefficient of x^n.
        '''
        if n < 0:
            raise ValueError('Coefficient of x^%d does not exist' % n)
        if n >= len(self._coef):
            return 0
        return self._coef[n]


    def coefAsPairs(self):
        '''Return an iterator of non-zero coefficients as (coef,exponent) pairs.
        '''
        for exp,coef in enumerate(self._coef):
            if coef:
                yield coef,exp


    def deg(self):
        '''return the degree of the polynomial'''
        # The zero polynomial is represented as an empty tuple so this defines
        # the degree of zero as -1 (good as any other undefined value).
        return len(self._coef) - 1
    deg = property(deg)


    def deriv(self):
        '''Return the derivative of the polynomial.
        '''
        # Handle all the constant forms
        if self.deg < 1: return Polynomial._zero

        lst = []
        for coef,exp in self.coefAsPairs():
            if exp:
                lst.append((coef*exp, exp-1))

        return Polynomial(lst)
    deriv = property(deriv)


    def sample(self, a, b, n):
        '''Evaluate the polynomial at n sample points.

        The samples are evenly spaced along a closed line between the points
        a and b.  The distinction between a line and an interval only matters
        for complex a or b.
        '''
        if n < 2 or a == b:
            raise ValueError('Call the polynomial if you only want one sample')

        if isinstance(a, (int,long)) and isinstance(b, (int,long)):
            # Force float only in cases where both endpoints are integers.
            # Other data types are expected to handle fractional values.
            b = float(b)

        d = (b-a)/(n-1)
        for i in xrange(n):
            x = a + i*d
            yield x, self(x)


    def float(self):
        '''Return a new polynomial with all of the coefficients forced to
        float.  Note: an exception is raised if a coefficient is complex.
        '''
        poly = object.__new__(Polynomial)
        poly._coef = tuple([float(c) for c in self._coef])
        return poly


    def complex(self):
        '''Return a new polynomial with all of the coefficients forced to
        complex.
        '''
        poly = object.__new__(Polynomial)
        poly._coef = tuple([complex(c) for c in self._coef])
        return poly


    def mpf(self, prec=0):
        '''Return a new polynomial with all of the coefficients forced to
        extended precision floating point.  The optional precision parameter
        allows selection of the precision of all the coefficients.
        '''
        poly = object.__new__(Polynomial)
        poly._coef = tuple([mpf(c, prec) for c in self._coef])
        return poly


    def cmpf(self, prec=0):
        '''Return a new polynomial with all of the coefficients forced to
        extended precision complex.  The optional precision parameter allows
        selection of the precision of all the coefficients.
        '''
        poly = object.__new__(Polynomial)
        poly._coef = tuple([cmpf(c, prec=prec) for c in self._coef])
        return poly


    def roots(self, eps=1e-16, start=0):
        '''Find all roots of the polynomial using Laguerre's method.
        '''
        p, dp, ddp, prec, prec2 = _laguerreInputs(self, eps)

        roots = []
        q = p
        dq = dp
        ddq = ddp

        # Allow the user to specify a starting point since guessing the first
        # root can improve the convergence of the process.
        r0 = cmpf(start, prec=prec2)

        while q.deg > 0:
            r = _improveRoot(r0, q, dq, ddq, eps)[1]
            roots.append(r)

            # Deflate the polynomial
            q = q/Polynomial(-r,1)
            dq = q.deriv
            ddq = dq.deriv

            # Try the root just found as the next guess.  This rapidly
            # eliminates multiple roots and sets up to find the next closest
            # root.
            r0 = r

        # Make sure there was no degredation of the precision in the
        # calculations.
        assert r0.prec >= prec2

        # The final precision should reflect the requested tolerance.
        prec = max(prec, 16)

        # Polish the roots since the reduction of the polynomial can lead to
        # accumulation of errors.
        for i,r in enumerate(roots):
            r = _improveRoot(r, p, dp, ddp, eps)[1]
            r = cmpf(r, prec=prec)
            if r.imag == 0:
                r = r.real
            roots[i] = r

        return sortRoots(roots)


    def uniqueRoots(self, eps=1e-16, start=0):
        '''Find the unique roots of the polynomial using Laguerre's method.
        '''
        g = gcd(self, self.deriv)

        if g.deg >= 1:
            p = self / g
        else:
            # There are no multiple roots so use the original polynomial.
            p = self

        return p.roots(eps, start)


    def codeGen(self, varName='x', resultName='y'):
        '''Generate a list of C code to evaluate the polynomial at the point
        varName.
        
        Note: the generated code is also valid Python.
        '''
        p = self._coef
        n = len(p)-1

        if n < 0: return ['%s = 0;' % resultName]

        subs = {'x':varName, 'y':resultName, 'pn':float(p[n])}
        lst = ['%(y)s = %(pn)r;' % subs]
        n -= 1
        while n >= 0:
            if p[n]:
                subs['pn'] = float(p[n])
                lst.append('%(y)s = %(y)s*%(x)s + %(pn)r;' % subs)
            else:
                lst.append('%(y)s = %(y)s*%(x)s;' % subs)
            n -= 1

        return lst


# Create the zero polynomial.
poly = object.__new__(Polynomial)
poly._coef = ()
Polynomial._zero = poly
del poly

#-----------------------------------------------------------------------------
class RationalFunction(object):
    '''Rational function data type that can mix with ordinary numbers and
    polynomials.
    
    Supports the usual arithmetic operations (+,-,*,/,**).
    '''

    # Width of the screen in characters.  Setting this to zero will suppress
    # the 3 line string representation.  If you only want it suppressed on
    # particular instances, set it on the instance instead of the class.
    screenWidth = 80

    def __new__(cls, numer, denom=None):
        if denom is None:
            # In this case the call is being used to perform a type cast to a
            # rational function.  If the input is already a rational function,
            # just return it.
            if isinstance(numer, cls):
                return numer
            denom = 1

        numer = Polynomial(numer)
        if not numer:
            # The zero rational function is always represented as 0/1
            denom = 1
        denom = Polynomial(denom)

        if not denom:
            raise ZeroDivisionError('Rational function denominator is zero')

        # Remove any common factors.  Avoid numeric only common factors since
        # the rational function needs to be normalied below.
        if numer.deg > 0 and denom.deg > 0:
            d = gcd(numer, denom)
            if d.deg > 0:
                numer /= d
                denom /= d

        # Force the denominator to be monic so that equality of rational
        # functions can be evaluated.
        c = denom._coef[-1]
        if c != 1:
            numer /= c
            denom /= c

        # Create and return the instance.
        fun = object.__new__(cls)
        fun._numer = numer
        fun._denom = denom

        return fun


    def __eq__(self, other):
        other = RationalFunction(other)
        return (self._numer == other._numer and self._denom == other._denom)

    def __ne__(self, other):
        return not self.__eq__(other)


    def __nonzero__(self):
        return bool(self._numer)


    def __pos__(self):
        return self

    def __neg__(self):
        return RationalFunction(-self._numer, self._denom)


    def __add__(self, other):
        other = RationalFunction(other)

        # Use algorithm for fractions from Knuth.
        un, ud = self._numer, self._denom
        vn, vd = other._numer, other._denom

        # If either polynomial is a constant, don't use gcd.
        if ud.deg < 1 or vd.deg < 1:
            d1 = 1
        else:
            d1 = gcd(ud, vd)
        if d1 == 1 or d1.deg < 1:
            return RationalFunction(un*vd + ud*vn, ud*vd)

        t = un*(vd/d1) + vn*(ud/d1)
        d2 = gcd(t, d1)
        if d2.deg < 1:
            return RationalFunction(t, (ud/d1)*vd)

        return RationalFunction(t/d2, (ud/d1)*(vd/d2))

    def __radd__(self, other):
        return RationalFunction(other).__add__(self)


    def __sub__(self, other):
        return self.__add__(-RationalFunction(other))

    def __rsub__(self, other):
        return RationalFunction(other).__add__(-self)


    def __mul__(self, other):
        other = RationalFunction(other)

        un, ud = self._numer, self._denom
        vn, vd = other._numer, other._denom

        # If any constant polynomials are involved, skip the gcd calculation.
        if un.deg < 1 or ud.deg < 1 or vn.deg < 1 or vd.deg < 1:
            return RationalFunction(un*vn, ud*vd)

        # Use algorithm for fractions from Knuth.
        d1 = gcd(un, vd)
        d2 = gcd(ud, vn)

        if d1.deg < 1 and d2.deg < 1:
            return RationalFunction(un*vn, ud*vd)

        if d1.deg > 0 and d2.deg < 1:
            return RationalFunction((un/d1)*vn, ud*(vd/d1))

        if d1.deg < 1 and d2.deg > 0:
            return RationalFunction(un*(vn/d2), (ud/d2)*vd)

        return RationalFunction((un/d1)*(vn/d2), (ud/d2)*(vd/d1))

    def __rmul__(self, other):
        return RationalFunction(other).__mul__(self)


    def __div__(self, other):
        other = RationalFunction(other)
        return self.__mul__(RationalFunction(other._denom, other._numer))

    def __rdiv__(self, other):
        other = RationalFunction(other)
        return other.__mul__(RationalFunction(self._denom, self._numer))

    # Make true division work the same as regular division.  This will allow
    # the / operator to work with true division turned on.
    def __truediv__(self, other):
        other = RationalFunction(other)
        return self.__mul__(RationalFunction(other._denom, other._numer))

    def __rtruediv__(self, other):
        other = RationalFunction(other)
        return other.__mul__(RationalFunction(self._denom, self._numer))


    def __pow__(self, n):
        if not isinstance(n, (int,long)):
            return NotImplemented
        if n == 0:
            return RationalFunction(1)

        # Use the right-to-left binary method from Knuth.
        y = RationalFunction(1)
        if n < 0:
            z = 1/self
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
        a = self._numer(x)
        b = self._denom(x)

        # If the results are both integers, return a rational.
        if isinstance(a, (int, long)) and isinstance(b, (int, long)):
            return mpq(a,b)

        # Otherwise the divide should work correctly.
        return a / b


    def __repr__(self):
        return 'RationalFunction(%r, %r)' % (self._numer, self._denom)


    def __str__(self):
        if self._denom.deg == 0:
            # The rational function is a polynomial divided by one so just
            # display the polynomial.
            return str(self._numer)

        n = str(self._numer)
        d = str(self._denom)
        m = max(len(n), len(d), 3)
        if m < self.screenWidth:
            return '%s\n%s\n%s' % (n.center(m).rstrip(), '-'*m,
                                   d.center(m).rstrip())
        else:
            return '(%s) / (%s)' % (n, d)


    # Note: These operations perform the same as the corresponding attributes
    # of mpq.
    def numer(self):
        return self._numer
    numer = property(numer)

    def denom(self):
        return self._denom
    denom = property(denom)


    def deriv(self):
        '''Return the derivative of the rational function.
        '''
        a = RationalFunction(self._numer.deriv, self._denom)
        b = RationalFunction(self._denom.deriv, self._denom)
        return a - self*b
    deriv = property(deriv)


    # Note: integration of a rational function gives a rational function only
    # in some special cases.  Consequently, the integration method is not
    # implemented.


    def sample(self, a, b, n):
        '''Evaluate the function at n sample points.

        The samples are evenly spaced along a closed line between the points
        a and b.  The distinction between a line and an interval only matters
        for complex a or b.
        '''
        if n < 2 or a == b:
            raise ValueError('Call the function if you only want one sample')

        if isinstance(a, (int,long)) and isinstance(b, (int,long)):
            # Force float only in cases where both endpoints are integers.
            # Other data types are expected to handle fractional values.
            b = float(b)

        d = (b-a)/(n-1)
        for i in xrange(n):
            x = a + i*d
            yield x, self(x)


    def float(self):
        '''Return a new function with all of the coefficients forced to
        float.
        '''
        return RationalFunction(self._numer.float(), self._denom.float())


    def complex(self):
        '''Return a new function with all of the coefficients forced to
        complex.
        '''
        return RationalFunction(self._numer.complex(), self._denom.complex())


    def mpf(self, prec=0):
        '''Return a new function with all of the coefficients forced to
        extended precision float.
        '''
        return RationalFunction(self._numer.mpf(prec), self._denom.mpf(prec))


    def cmpf(self, prec=0):
        '''Return a new function with all of the coefficients forced to
        extended precision complex.
        '''
        return RationalFunction(self._numer.cmpf(prec), self._denom.cmpf(prec))

#-----------------------------------------------------------------------------
def gcd(u, v):
    '''gcd(u,v) returns the greatest common divisor of u and v.
    Where u and v are integers or polynomials.
    '''
    # If something other than an integer or a polynomial is fed in, the
    # algorithm may fail to converge.  In that case we want an exception
    # instead of an infinite loop.
    iterations = 1000
    while v:
        u, v = v, u % v
        iterations -= 1
        if iterations < 1:
            raise ValueError('gcd(u,v) failed to converge')
    return u

#-----------------------------------------------------------------------------
def polyFromRoots(roots):
    '''Generate a polynomial from a sequence of roots.

    If the sequence of roots is empty, the constant polynomial 1 is returned.
    '''
    p = Polynomial(1)
    for r in roots:
        p *= Polynomial(-r,1)
    return p

#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
def improveRoots(roots, p, eps=1e-16):
    return [improveRoot(r, p, eps)[1] for r in roots]

#-----------------------------------------------------------------------------
def improveRoot(root, p, eps=1e-16):
    '''Use Laguerre's method to improve the estimate of polynomial root.

    Returns an estimate of the error and the current estimate of the root.  The
    error will be less than eps if the algorithm converges within the iteration
    limit.

    root - initial estimate of the root
    p - polynomial
    eps - error tolerance
    '''
    p, dp, ddp, prec, prec2 = _laguerreInputs(p, eps)
    z = cmpf(root, prec=prec2)
    err, z = _improveRoot(z, p, dp, ddp, eps)

    # The precision should reflect the requested tolerance.
    prec = max(prec, 16)
    z = cmpf(z, prec=prec)
    if z.imag == 0:
        z = z.real
    return err, z

#-----------------------------------------------------------------------------
def _laguerreInputs(p, eps):
    if p.deg < 1:
        raise ValueError('The degree of the polynomial must be at least one.')

    if not (0 < eps < 0.001):
        raise ValueError('Tolerance out of range')

    # The precision is set from the base 10 log of the tolerance
    prec = int(math.ceil(-math.log10(eps)))
    prec2 = prec*prec  # Needed to ensure convergence

    # Force all coefficients to high precision complex numbers.
    p = p.cmpf(prec2)
    dp = p.deriv
    ddp = dp.deriv
    return p, dp, ddp, prec, prec2

#-----------------------------------------------------------------------------
def _improveRoot(z, p, dp, ddp, eps):
    loop = 100
    while loop > 0:
        try:
            a = _laguerreStep(z, p, dp, ddp)
        except ZeroDivisionError:
            # In this case the derivative is zero so take a step away from that
            # point.
            z = z+1
            a = _laguerreStep(z, p, dp, ddp)
        z = z - a
        loop -= 1

        # Compute the estimated distance from the root.
        err = abs(complex(a))
        if err < eps:
            break

    return err, z

#-----------------------------------------------------------------------------
# Note: The algorithm is described in Numerical Recipes

def _laguerreStep(z, p, dp, ddp):
    '''Perform one iteration of Laguerre's method.

    z - current estimate of the root
    p - polynomial
    dp - first derivative of p
    ddp - second derivative of p
    '''
    pz = p(z)

    if not pz:
        # The result is zero so a zero has already been located.
        return pz

    dpz = dp(z)
    ddpz = ddp(z)

    G = dpz/pz
    H = G*G - ddpz/pz

    n = p.deg
    s = np.sqrt((n-1)*(n*H - G*G))
    d1 = G + s
    d2 = G - s

    # No need to do this comparison in high precision.
    if abs(complex(d1)) < abs(complex(d2)):
        return n/d2
    else:
        return n/d1


#-----------------------------------------------------------------------------
# This section contains some helper functions that are not part of the user
# interface.
#-----------------------------------------------------------------------------
# Choose the appropriate type for computing 1/x.

def _invert(x):
    # Integers are converted to rationals
    if isinstance(x, (int, long)):
        return mpq(1,x)

    # All other types should handle the calculation correctly.
    return 1/x

#-----------------------------------------------------------------------------
def _coefFromSeq(seq):
    if not seq:
        return ()

    try:
        n = len(seq[0])
    except TypeError:
        n = 0

    if n == 2:
        return _coefFromPairs(seq)
    elif n > 0:
        raise ValueError('Invalid coefficient sequence')

    # Validate the input sequence.
    lst = []
    for x in seq:
        if isinstance(x, basestring):
            raise ValueError('Cannot create coefficients from a string')
        # Just check that the value can be converted to complex.  This makes it
        # numeric enough to be an acceptable coefficient.
        complex(x)

        # Floats and complex do not mix with rationals so convert them
        # to the corresponding clnum type.
        if isinstance(x, float):
            x = mpf(x)
        elif isinstance(x, complex):
            x = cmpf(x)
        lst.append(x)
    seq = lst

    # Normalize the polynomial by removing zero coefficients of higher powers.
    i = len(seq)-1
    while i >= 0:
        if seq[i]:
            break
        i -= 1

    return tuple(seq[:i+1])

#-----------------------------------------------------------------------------
def _coefFromPairs(seq):
    pairs = {}
    for coef, exp in seq:
        # Check the input for validity.  The coefficient can be any numeric
        # type.  Conversion to complex is the criterian for it to be numeric.
        if not isinstance(exp, (int, long)) or exp < 0:
            raise ValueError('Invalid exponent')
        if isinstance(coef, basestring):
            raise ValueError('Cannot create coefficients from a string')
        complex(coef)

        if exp in pairs:
            raise ValueError('Duplicate exponent')

        # Floats and complex do not mix with rationals so convert them
        # to the corresponding clnum type.
        if isinstance(coef, float):
            coef = mpf(coef)
        elif isinstance(coef, complex):
            coef = cmpf(coef)

        # Accumulate the pairs in a dictionary
        pairs[exp] = coef

    # Find the length of the coefficient list from the largest exponent.
    n = max(pairs) + 1

    # Use the coefficient of the highest power to determine the type of zero to
    # use.  This produces a uniform data type for those cases where the user
    # has been careful to construct the polynomial that way in the first place.
    # In some cases, this can speed up some of the calculations by avoiding
    # type conversions.
    junk,zero = coerce(pairs[n-1],0)

    # Insert the coefficients in the proper places in the list.
    lst = [zero]*n
    for exp, coef in pairs.iteritems():
        lst[exp] = coef

    # Normalize the polynomial by removing zero coefficients of higher powers.
    i = n-1
    while i >= 0:
        if lst[i]:
            break
        i -= 1

    return tuple(lst[:i+1])
