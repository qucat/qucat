from math import floor
import numpy as np
import functools
from warnings import warn
import traceback
np.seterr(divide='raise', invalid='raise')


exponent_to_letter = {
    -18: 'a',
    -15: 'f',
    -12: 'p',
    -9: 'n',
    -6: 'u',
    -3: 'm',
    0: '',
    3: 'k',
    6: 'M',
    9: 'G',
    12: 'T'
}
exponent_to_letter_math = {
    -18: 'a',
    -15: 'f',
    -12: 'p',
    -9: 'n',
    -6: r'$\mu$',
    -3: 'm',
    0: '',
    3: 'k',
    6: 'M',
    9: 'G',
    12: 'T'
}
exponent_to_letter_unicode = {
    -18: 'a',
    -15: 'f',
    -12: 'p',
    -9: 'n',
    -6: u'\u03bc',
    -3: 'm',
    0: '',
    3: 'k',
    6: 'M',
    9: 'G',
    12: 'T'
}

def safely_evaluate(func_to_evaluate):
    @functools.wraps(func_to_evaluate)
    def wrapper_safely_evaluate(self,w, **kwargs):
        try:
            return func_to_evaluate(self,w, **kwargs)
        except FloatingPointError as e:
            warn(str(e)+f"\n\tPerturbing f = {w/2./np.pi} to find finite value")
            perturbation = 1e-14
            while perturbation<0.01:
                try:
                    to_return = func_to_evaluate(self,w*(1.+perturbation), **kwargs)
                    warn(f"Perturbed f = {w/2./np.pi}*(1+{perturbation:.1e}) to find finite value")
                    return to_return
                except FloatingPointError as e:
                    perturbation *= 10
            raise FloatingPointError("Even perturbing the frequency by a percent failed to produce finite value.")
    return wrapper_safely_evaluate


def vectorize(func_to_evaluate):
    @functools.wraps(func_to_evaluate)
    def wrapper_vectorize(self,w, **kwargs):
        try:
            iter(w)
        except TypeError:
            # single
            return func_to_evaluate(self,w,**kwargs)
        else:
            # iterable = True
            return np.vectorize(lambda w_single: func_to_evaluate(self,w_single,**kwargs))(w)
    return wrapper_vectorize


def pretty_value(v, use_power_10=False, use_math=True, use_unicode=False, maximum_info = False):
    if v == 0:
        return '0'
    elif v < 0:
        sign = '-'
        v *= -1
    elif v > 0:
        sign = ''
    exponent = floor(np.log10(v))
    exponent_3 = exponent-(exponent % 3)
    float_part = v/(10**exponent_3)
    if use_power_10 or v >= 1e15 or v < 1e-18:
        if exponent_3 == 0:
            exponent_part = ''
        else:
            if use_math:
                exponent_part = r'$\times 10^{%d}$' % exponent_3
            else:
                exponent_part = r'e%d' % exponent_3
    else:
        if use_unicode:
            exponent_part = ' '+exponent_to_letter_unicode[exponent_3]
        elif use_math:
            exponent_part = ' '+exponent_to_letter_math[exponent_3]
        else:
            exponent_part = ' '+exponent_to_letter[exponent_3]

    if maximum_info == False:
        if float_part >= 10.:
            pretty = "%.0f%s" % (float_part, exponent_part)
        else:
            pretty = "%.1f%s" % (float_part, exponent_part)
        return sign+pretty
    else:
        if ("%.7f"%float_part)[-1] != '0':
            float_part ="%.6f.."%float_part
        else:
            float_part ="%.6f"%float_part
            while float_part[-1]=="0":
                float_part = float_part[:-1]
            if float_part[-1]==".":
                float_part = float_part[:-1]

        return sign + float_part + exponent_part

def shift(to_shift,shift):
    for i,_ in enumerate(to_shift):
        to_shift[i]+= shift
    return to_shift

def to_string(unit,label,value, use_math=True, use_unicode=False, maximum_info = False):

    if unit is not None:
        if use_unicode:
            unit = unit.replace(r'$\Omega$', u"\u03A9")
        if use_math == False:
            unit = unit.replace(r'$\Omega$', 'Ohm')

    if label is None:
        s = ''
    else:
        if use_math:
            s = "$%s$" % (label)
        else:
            s = label

        if value is not None:
            s+='='

    if value is not None:
        s+= pretty_value(
            value, use_math=use_math, use_unicode=use_unicode,maximum_info=maximum_info)

        if unit is not None:
            s+=unit
    return s