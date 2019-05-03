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


def pretty_value(v,is_complex = True, use_power_10=False, use_unicode=True, maximum_info = False):

    if v == 0:
        return '0'

    def get_exponent_3(value):
        value = np.absolute(value)
        exponent = floor(np.log10(value))
        exponent_3 = exponent-(exponent % 3)
        return exponent_3
    
    def exponent_3_to_string(value,exponent_3):
        value = np.absolute(value)
        if use_power_10 or value >= 1e15 or value < 1e-18:
            if exponent_3 == 0:
                exponent_part = ''
            else:
                exponent_part = r'e%d' % exponent_3
        else:
            if use_unicode:
                exponent_part = ' '+exponent_to_letter_unicode[exponent_3]
            else:
                exponent_part = ' '+exponent_to_letter[exponent_3]
        return exponent_part
    
    def get_float_part(v,exponent_3):
        if v == 0:
            return '','0'
        if v < 0:
            sign = '-'
            v *= -1
        elif v > 0:
            sign = ''
        float_part = v/(10**exponent_3)

        if maximum_info == False:
            if float_part >= 100.:
                float_part = "%.0f" % (float_part)
            elif float_part >= 10.:
                float_part = "%.1f" % (float_part)
            else:
                float_part = "%.2f" % (float_part)
            # remove trailing 0s or .
            while float_part[-1]=="0":
                float_part = float_part[:-1]
            if float_part[-1]==".":
                float_part = float_part[:-1]

        else:
            if float_part-float("%.6f"%float_part) != 0:
                # if there is a digit beyond digit 6
                float_part ="%.6f.."%float_part
            else:
                float_part ="%.6f"%float_part
                while float_part[-1]=="0":
                    float_part = float_part[:-1]
                if float_part[-1]==".":
                    float_part = float_part[:-1]
        return [sign,float_part]

    if is_complex:
        exp3 = get_exponent_3(max(np.absolute(np.real(v)),np.absolute(np.imag(v))))
        exponent = exponent_3_to_string(max(np.absolute(np.real(v)),np.absolute(np.imag(v))),exp3)

        sign_r, numbers_r = get_float_part(np.real(v),exp3)
        sign_i, numbers_i = get_float_part(np.imag(v),exp3)

        to_return = ''
        if numbers_r.replace('0','').replace('.','')!='':
            to_return += sign_r+numbers_r 
            if numbers_i.replace('0','').replace('.','')!='':
                if sign_i == '':
                    to_return += '+'
        if numbers_i.replace('0','').replace('.','')!='':
            to_return += sign_i+numbers_i +'i'
        
        return to_return + exponent
    else:
        v = np.real(v)
        exp3 = get_exponent_3(v)
        exponent = exponent_3_to_string(v,exp3)
        sign, numbers = get_float_part(v,exp3)
        return sign+numbers+exponent

def shift(to_shift,shift):
    for i,_ in enumerate(to_shift):
        to_shift[i]+= shift
    return to_shift

def to_string(unit,label,value, use_unicode=True, maximum_info = False):

    if unit is not None:
        if not use_unicode:
            unit = unit.replace(u"\u03A9", 'Ohm')

    if label is None:
        s = ''
    else:
        s = label

        if value is not None:
            s+='='

    if value is not None:
        s+= pretty_value(
            value, use_unicode=use_unicode,maximum_info=maximum_info)

        if unit is not None:
            s+=unit
    return s

if __name__ == "__main__":
    print(pretty_value((1.1110564436955427e-09-3.4692726334931547e-10j),is_complex = True))