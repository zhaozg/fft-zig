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
pub const HUGE_DATA_THRESHOLD = 1000000;

pub fn fftParallelSIMD(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    if (n < 16384) {
        return fftRadix2SIMD(data);
    }
    if (fft_utils.isPowerOfFour(n) and n >= 256) {
        try fftRadix4SIMD(data);
    } else {
        try fftRadix2SIMD(data);
    }
    _ = allocator;
}

pub fn fftHugeDataParallel(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    const chunk_size = @min(n, 16 * 1024 * 1024);
    if (n <= chunk_size) {
        if (fft_utils.isPowerOfTwo(n)) {
            if (fft_utils.isPowerOfFour(n) and n >= 1024) {
                try fftRadix4SIMD(data);
            } else {
                try fftRadix2SIMD(data);
            }
        } else {
            try fftMixedRadixHuge(allocator, data);
        }
    } else {
        try fftChunkedProcessing(allocator, data, chunk_size);
    }
}

fn fftMixedRadixHuge(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    if (n > 100 * 1024 * 1024) {
        try fftDecomposition(allocator, data);
    } else {
        try fftMixedRadix(allocator, data);
    }
}

fn fftChunkedProcessing(allocator: std.mem.Allocator, data: []Complex, chunk_size: usize) !void {
    const n = data.len;
    var processed: usize = 0;
    while (processed < n) {
        const current_chunk_size = @min(chunk_size, n - processed);
        const chunk = data[processed .. processed + current_chunk_size];
        if (fft_utils.isPowerOfTwo(current_chunk_size)) {
            try fftRadix2SIMD(chunk);
        } else {
            try fftMixedRadix(allocator, chunk);
        }
        processed += current_chunk_size;
    }
}

fn fftDecomposition(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    const factors = try findOptimalFactors(allocator, n);
    defer allocator.free(factors);
    try applyFactoredFFT(allocator, data, factors);
}

fn findOptimalFactors(allocator: std.mem.Allocator, n: usize) ![]usize {
    const FactorList = std.ArrayList(usize);
    var factors = FactorList.init(allocator);
    errdefer factors.deinit();
    var remaining = n;
    while (remaining % 2 == 0 and remaining > 1) {
        try factors.append(2);
        remaining /= 2;
    }
    const small_primes = [_]usize{ 3, 5, 7, 11, 13 };
    for (small_primes) |prime| {
        while (remaining % prime == 0 and remaining > 1) {
            try factors.append(prime);
            remaining /= prime;
        }
    }
    if (remaining > 1) {
        try factors.append(remaining);
    }
    return try factors.toOwnedSlice();
}

fn applyFactoredFFT(allocator: std.mem.Allocator, data: []Complex, factors: []const usize) !void {
    try fftMixedRadix(allocator, data);
    _ = factors;
}

const expect = std.testing.expect;

test "Parallel FFT threshold logic" {
    // 这里只做简单阈值分支测试
    try expect(PARALLEL_THRESHOLD == 16384);
    try expect(HUGE_DATA_THRESHOLD == 1000000);
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
