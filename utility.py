from math import floor
import numpy as np

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

def pretty_value(v, use_power_10=False, use_math=True, use_unicode=False):
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
    if float_part >= 10.:
        pretty = "%.0f%s" % (float_part, exponent_part)
    else:
        pretty = "%.1f%s" % (float_part, exponent_part)
    return sign+pretty

def check_there_are_no_iterables_in_kwarg(**kwargs):
    for el, value in kwargs.items():
        try:
            iter(value)
        except TypeError:
            pass
        else:
            raise ValueError(
                "This function accepts no lists or iterables as input.")

def shift(to_shift,shift):
    for i,_ in enumerate(to_shift):
        to_shift[i]+= shift
    return to_shift

def to_string(unit,label,value, use_math=True, use_unicode=False):

    unit = unit
    if use_unicode:
        unit = unit.replace(r'$\Omega$', u"\u03A9")
    if use_math == False:
        unit = unit.replace(r'$\Omega$', 'Ohm')

    label = label
    if use_math:
        label = "$%s$" % (label)

    if value is not None:
        pvalue = pretty_value(
            value, use_math=use_math, use_unicode=use_unicode)

    if label is None:
        return pvalue+unit
    elif label == '' and value is None:
        return ''
    elif value is None:
        return label
    elif label == '' and value is not None:
        return pvalue+unit
    else:
        return label+'='+pvalue+unit