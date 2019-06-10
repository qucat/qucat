from math import floor
import numpy as np
import functools
from warnings import warn
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



def refuse_vectorize_kwargs(func_to_evaluate = None,*,exclude = []):
    # Only works for functions which return a list

    def _decorate(func):
        @functools.wraps(func)
        def wrapper_vectorize(self, *args,**kwargs):
            non_iterables = {}
            iterables = {}
            for kw, arg in kwargs.items():
                if kw not in exclude:
                    try:
                        iter(arg)
                        raise ValueError("No iterables are allowed, use a single value for %s"%kw)
                    except TypeError:
                        # not an iterable
                        pass
            
            return func(self,*args,**kwargs)
        return wrapper_vectorize

    if func_to_evaluate:
        return _decorate(func_to_evaluate)
    return _decorate

def vectorize_kwargs(func_to_evaluate = None,*,exclude = []):
    # Only works for functions which return a list

    def _decorate(func):
        @functools.wraps(func)
        def wrapper_vectorize(self, *args,**kwargs):
            non_iterables = {}
            iterables = {}
            for kw, arg in kwargs.items():
                if kw in exclude:
                    non_iterables[kw] = arg
                else:
                    try:
                        iter(arg)
                    except TypeError:
                        # not an iterable
                        non_iterables[kw] = arg
                    else:
                        # is an iterable
                        # Make sure it has the same shape as other iterables
                        if len(iterables)>0:
                            first_iterable = iterables[list(iterables)[0]]
                            if np.array(arg).shape != first_iterable.shape:
                                raise ValueError("Keyword arguments have incompatible shapes: %s %s and %s %s"%(
                                    list(iterables)[0],
                                    first_iterable.shape,
                                    kw,
                                    np.array(arg).shape
                                ))
                        iterables[kw] = np.array(arg)
            
            if len(iterables)==0:
                return func(self,*args,**kwargs)
            else:
                first_iterable = iterables[list(iterables)[0]]
                kwargs_single = non_iterables
                i = 0
                for index,_ in np.ndenumerate(first_iterable):

                    for kw, arg in iterables.items():
                        kwargs_single[kw] = arg[index]

                    to_return_single = func(self,*args,**kwargs_single)

                    if i == 0:
                        try:
                            iter(to_return_single)
                            to_return = np.empty((*first_iterable.shape,*to_return_single.shape), dtype=np.complex128)
                        except TypeError:
                            # not an iterable
                            to_return = np.empty(first_iterable.shape, dtype=np.complex128)
                        i+=1

                    to_return[index] = to_return_single
                for i in range(len(first_iterable.shape)):
                    to_return = np.moveaxis(to_return,0,-1)
                return to_return
        return wrapper_vectorize

    if func_to_evaluate:
        return _decorate(func_to_evaluate)
    return _decorate

def get_exponent_3(value):
    value = np.absolute(value)
    exponent = floor(np.log10(value))
    exponent_3 = exponent-(exponent % 3)
    return exponent_3

def exponent_3_to_string(value,exponent_3,use_power_10,use_unicode):
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

def get_float_part(v,exponent_3,maximum_info):
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
            float_part = "%2.0f." % (float_part)
        elif float_part >= 10.:
            float_part = "%1.1f" % (float_part)
        else:
            float_part = "%0.2f" % (float_part)
        # remove trailing 0s or .
        while float_part[-1]=="0":
            float_part = float_part[:-1]
        if float_part[-1]==".":
            float_part = float_part[:-1]

    else:
        if (v-float("%2.6e"%v) != 0 and float_part >= 100.)\
            or (v-float("%1.6e"%v) != 0 and 100 > float_part >= 10.)\
            or (v-float("%.6e"%v) != 0 and float_part < 10.) :
            # if there is a digit beyond digit 6
            float_part ="%2.6f.."%float_part
        else:
            float_part ="%.6f"%float_part
            while float_part[-1]=="0":
                float_part = float_part[:-1]
            if float_part[-1]==".":
                float_part = float_part[:-1]
    return [sign,float_part]


def pretty_value(v,is_complex = True, use_power_10=False, use_unicode=True, maximum_info = False):

    if v == 0:
        return '0'


    if is_complex:
        exp3 = get_exponent_3(max(np.absolute(np.real(v)),np.absolute(np.imag(v))))
        exponent = exponent_3_to_string(max(np.absolute(np.real(v)),np.absolute(np.imag(v))),exp3,use_power_10,use_unicode)

        sign_r, numbers_r = get_float_part(np.real(v),exp3,maximum_info)
        sign_i, numbers_i = get_float_part(np.imag(v),exp3,maximum_info)

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
        exponent = exponent_3_to_string(v,exp3,use_power_10,use_unicode)
        sign, numbers = get_float_part(v,exp3,maximum_info)
        return sign+numbers+exponent

def shift(to_shift,shift):
    for i,_ in enumerate(to_shift):
        to_shift[i]+= shift
    return to_shift

def to_string(unit,label,value, use_unicode=True, maximum_info = False):

    if unit is not None:
        if not use_unicode:
            unit = unit.replace(u"\u03A9", 'Ohm')
            unit = unit.replace(u'\u03bc', 'u')
                
    if label is None:
        s = ''
    else:
        s = label
        if value is not None:
            s+='='

    if value is not None:
        s+= pretty_value(value, use_unicode=use_unicode,maximum_info=maximum_info)

        if unit is not None:
            s+=unit
    return s
