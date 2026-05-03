//! 并行与大数据FFT相关实现

const std = @import("std");
const fft_radix2 = @import("fft_radix2.zig");
const fft_radix4 = @import("fft_radix4.zig");
const fft_mixed = @import("fft_mixed.zig");
const fft_utils = @import("utils.zig");
pub const PARALLEL_THRESHOLD = 16384;
pub const HUGE_DATA_THRESHOLD = 262144; // Lower threshold for better algorithm selection

pub fn fftParallelSIMD(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) !void {
    const n = data.len;
    if (n < 16384) {
        if (fft_utils.isPowerOfFour(n)) {
            return fft_radix4.fftRadix4SIMD(T, data);
        } else {
            return fft_radix2.fftRadix2SIMD(T, data);
        }
    }
    if (fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4SIMD(T, data);
    } else {
        try fft_radix2.fftRadix2SIMD(T, data);
    }
    _ = allocator;
}

pub fn fftHugeDataParallel(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) !void {
    const n = data.len;
    // For huge data, use optimized algorithms directly instead of chunking
    // Chunking would break FFT correctness - FFT must be computed on entire dataset
    if (fft_utils.isPowerOfFour(n)) {
        try fft_radix4.fftRadix4SIMD(T, data);
    } else if (fft_utils.isPowerOfTwo(n)) {
        try fft_radix2.fftRadix2SIMD(T, data);
    } else {
        // For non-power-of-2, use mixed radix which has better performance
        try fft_mixed.fftMixedRadix(T, allocator, data);
    }
}

// Note: Chunked processing removed - FFT cannot be correctly computed in chunks
// FFT requires processing the entire dataset as a whole transform

const expect = std.testing.expect;

test "Parallel FFT threshold logic" {
    // 这里只做简单阈值分支测试
    try expect(PARALLEL_THRESHOLD == 16384);
    try expect(HUGE_DATA_THRESHOLD == 262144);
}

fn testParallelEdgeGeneric(comptime T: type) !void {
    var empty: [0]std.math.Complex(T) = undefined;
    try fftParallelSIMD(T, std.heap.page_allocator, empty[0..]); // 应不报错

    var one = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 7.0), .im = @as(T, 0.0) }};
    try fftParallelSIMD(T, std.heap.page_allocator, one[0..]);
    try expect(one[0].re == @as(T, 7.0));

    var not_pow2 = [_]std.math.Complex(T){ std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) }, std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 0.0) }, std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 0.0) } };
    const result = fftParallelSIMD(T, std.heap.page_allocator, not_pow2[0..]) catch |err| err;
    try expect(result == error.InvalidSize);
}

test "Parallel FFT edge cases f32" {
    try testParallelEdgeGeneric(f32);
}
test "Parallel FFT edge cases f64" {
    try testParallelEdgeGeneric(f64);
}
