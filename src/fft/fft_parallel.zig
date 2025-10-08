//! 并行与大数据FFT相关实现

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_mixed = @import("fft_mixed.zig");
const fftRadix2SIMD = fft_radix2.fftRadix2SIMD;
const fftRadix4SIMD = fft_radix4.fftRadix4SIMD;
const fftMixedRadix = fft_mixed.fftMixedRadix;
const fft_utils = @import("utils.zig");
pub const PARALLEL_THRESHOLD = 16384;
pub const HUGE_DATA_THRESHOLD = 262144; // Lower threshold for better algorithm selection

pub fn fftParallelSIMD(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    if (n < 16384) {
        if (fft_utils.isPowerOfFour(n)) {
            return fftRadix4SIMD(data);
        } else {
            return fftRadix2SIMD(data);
        }
    }
    if (fft_utils.isPowerOfFour(n)) {
        try fftRadix4SIMD(data);
    } else {
        try fftRadix2SIMD(data);
    }
    _ = allocator;
}

pub fn fftHugeDataParallel(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    // For huge data, use optimized algorithms directly instead of chunking
    // Chunking would break FFT correctness - FFT must be computed on entire dataset
    if (fft_utils.isPowerOfFour(n)) {
        try fftRadix4SIMD(data);
    } else if (fft_utils.isPowerOfTwo(n)) {
        try fftRadix2SIMD(data);
    } else {
        // For non-power-of-2, use mixed radix which has better performance
        try fftMixedRadix(allocator, data);
    }
}

fn fftMixedRadixHuge(allocator: std.mem.Allocator, data: []Complex) !void {
    // Removed - now using direct fftMixedRadix for all sizes
    try fftMixedRadix(allocator, data);
}

// Note: Chunked processing removed - FFT cannot be correctly computed in chunks
// FFT requires processing the entire dataset as a whole transform

const expect = std.testing.expect;

test "Parallel FFT threshold logic" {
    // 这里只做简单阈值分支测试
    try expect(PARALLEL_THRESHOLD == 16384);
    try expect(HUGE_DATA_THRESHOLD == 262144);
}

test "Parallel FFT edge cases" {
    var empty: [0]Complex = undefined;
    try fftParallelSIMD(std.heap.page_allocator, empty[0..]); // 应不报错

    var one = [_]Complex{Complex{ .re = 7.0, .im = 0.0 }};
    try fftParallelSIMD(std.heap.page_allocator, one[0..]);
    try expect(one[0].re == 7.0);

    var not_pow2 = [_]Complex{ Complex{ .re = 1.0, .im = 0.0 }, Complex{ .re = 2.0, .im = 0.0 }, Complex{ .re = 3.0, .im = 0.0 } };
    const result = fftParallelSIMD(std.heap.page_allocator, not_pow2[0..]) catch |err| err;
    try expect(result == error.InvalidSize);
}
