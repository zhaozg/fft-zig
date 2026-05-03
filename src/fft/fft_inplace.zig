const std = @import("std");
const math = std.math;
const fft_twiddle = @import("twiddle.zig");
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_mixed = @import("fft_mixed.zig");
const fft_parallel = @import("fft_parallel.zig");
const fft_utils = @import("utils.zig");

// Performance thresholds for algorithm selection
pub const PARALLEL_THRESHOLD = 16384;
pub const HUGE_DATA_THRESHOLD = 1000000; // 1M threshold for huge data optimizations
pub const SIMD_THRESHOLD = 64;
pub const RADIX4_THRESHOLD = 256;
pub const SMALL_FFT_THRESHOLD = 256;

/// 泛型就地 FFT 函数
pub fn fftInPlace(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) error{ InvalidSize, OutOfMemory }!void {
    // 自动选择算法，委托给各模块
    const n = data.len;
    if (n >= HUGE_DATA_THRESHOLD) {
        try fft_parallel.fftHugeDataParallel(T, allocator, data);
    } else if (n >= PARALLEL_THRESHOLD and fft_utils.isPowerOfTwo(n)) {
        try fft_parallel.fftParallelSIMD(T, allocator, data);
    } else if (n >= RADIX4_THRESHOLD and fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4SIMD(T, data);
    } else if (n >= SIMD_THRESHOLD and fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2SIMD(T, data);
    } else if (fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4(T, data);
    } else if (fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2(T, data);
    } else if (n % 4 == 0 and fft_utils.isPowerOfTwo(n / 4)) {
        try fft_radix4.fftRadix4(T, data);
    } else {
        try fft_mixed.fftMixedRadix(T, allocator, data);
    }
}
