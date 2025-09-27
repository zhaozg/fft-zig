//! 基础就地FFT实现（无自动算法选择，仅供内部模块调用）
const std = @import("std");
const Complex = @import("types.zig").Complex;
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_utils = @import("utils.zig");

/// 最基础的就地FFT实现（仅支持power-of-2，内部使用）
pub fn fftInPlaceBase(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    if (n <= 1) return;
    if (fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4SIMD(data);
    } else if (fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2SIMD(data);
    } else {
        return error.InvalidSize;
    }
    _ = allocator; // 保持接口一致
}
