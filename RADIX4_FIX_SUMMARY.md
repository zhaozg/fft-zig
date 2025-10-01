# Radix-4 FFT Bug Fix Summary

## Problem
The Radix-4 FFT implementation produced incorrect results with magnitudes at approximately 30.7% of expected values for all non-DC frequency components.

## Root Cause
The Radix-4 butterfly algorithm had **swapped imaginary unit multiplication signs** in the final stage:

### Incorrect Implementation (Before Fix)
```zig
data[idx0] = Complex{ .re = temp0.re + temp1.re, .im = temp0.im + temp1.im };
data[idx1] = Complex{ .re = temp2.re - temp3.im, .im = temp2.im + temp3.re };  // Wrong: -i
data[idx2] = Complex{ .re = temp0.re - temp1.re, .im = temp0.im - temp1.im };
data[idx3] = Complex{ .re = temp2.re + temp3.im, .im = temp2.im - temp3.re };  // Wrong: +i
```

### Correct Implementation (After Fix)
```zig
data[idx0] = Complex{ .re = temp0.re + temp1.re, .im = temp0.im + temp1.im };
data[idx1] = Complex{ .re = temp2.re + temp3.im, .im = temp2.im - temp3.re };  // Correct: +i
data[idx2] = Complex{ .re = temp0.re - temp1.re, .im = temp0.im - temp1.im };
data[idx3] = Complex{ .re = temp2.re - temp3.im, .im = temp2.im + temp3.re };  // Correct: -i
```

## Mathematical Explanation
In the Radix-4 DIT FFT butterfly, the final outputs should be:
- `X[k] = a0 + a2`
- `X[k + N/4] = a1 + i·a3`  (multiply a3 by +i)
- `X[k + N/2] = a0 - a2`
- `X[k + 3N/4] = a1 - i·a3`  (multiply a3 by -i)

Where multiplying by +i means: `(re, im) → (-im, re)`
And multiplying by -i means: `(re, im) → (im, -re)`

The bug had these reversed, causing incorrect phase relationships and the characteristic 0.307 magnitude scaling.

## Test Results

### Before Fix
```
Size      16: mag=    5.2263 expected=    8.0000 ratio=0.653281
Size      64: mag=   10.2428 expected=   32.0000 ratio=0.320088
Size     256: mag=   39.4445 expected=  128.0000 ratio=0.308160
Size    1024: mag=  157.4069 expected=  512.0000 ratio=0.307435
Size 1048576: mag=161159.4296 expected=524288.0000 ratio=0.307387
```

### After Fix
```
Size      16: mag=    8.0000 expected=    8.0000 ratio=1.000000
Size      64: mag=   32.0000 expected=   32.0000 ratio=1.000000
Size     256: mag=  128.0000 expected=  128.0000 ratio=1.000000
Size    1024: mag=  512.0000 expected=  512.0000 ratio=1.000000
Size 1048576: mag=524288.0000 expected=524288.0000 ratio=1.000000
```

## Performance Impact
With Radix-4 re-enabled and fixed:
- **4M samples**: 2318ms → 1448ms (**37% faster**)
- **1M samples**: 518ms → 295ms (**43% faster**)  
- **All tests pass**: 25/25 ✅
- **Accuracy**: 100% (ratio = 1.000000 for all sizes)

## Files Modified
1. `src/fft/fft_radix4.zig` - Fixed butterfly algorithm (2 lines changed)
2. `src/fft.zig` - Re-enabled Radix-4 selection
3. `src/fft/fft_parallel.zig` - Re-enabled Radix-4 in parallel paths
4. `src/fft/base.zig` - Re-enabled Radix-4 in base FFT

## Verification
- ✅ DC signals: Perfect (magnitude = N)
- ✅ Single frequency sine waves: Perfect (magnitude = N/2)
- ✅ Multiple frequency signals: Perfect (Parseval's theorem satisfied)
- ✅ All unit tests pass
- ✅ Large data tests (1M-5M samples) pass
- ✅ Performance improved significantly

## Conclusion
The fix was minimal (swapping two sign operations) but critical. The bug was a subtle error in the implementation of the Radix-4 butterfly's imaginary unit rotations. With this fix, Radix-4 FFT now works correctly and provides significant performance improvements for power-of-4 sized inputs.
