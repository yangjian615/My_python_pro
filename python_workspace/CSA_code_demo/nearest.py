#!/usr/bin/env python

def nearest(array, value):
    residuals = array - value
    val, idx = min((val, idx) for (idx, val) in enumerate(abs(residuals)))
    return [array[idx], idx]