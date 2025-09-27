//! 逆FFT实现

const std = @import("std");
const Complex = @import("types.zig").Complex;
const fft_types = @import("types.zig");

const fftInPlace = @import("../fft.zig").fftInPlace;

pub fn ifftInPlace(allocator: std.mem.Allocator, data: []Complex) !void {
    const n = data.len;
    for (data) |*v| {
        v.im = -v.im;
    }
    try fftInPlace(allocator, data);
    for (data) |*v| {
        v.re = v.re / @as(f64, @floatFromInt(n));
        v.im = -v.im / @as(f64, @floatFromInt(n));
    }
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const TEST_TOLERANCE = 1e-12;

test "IFFT correctness" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var data = [_]Complex{
        Complex{ .re = 1.0, .im = 0.0 },
        Complex{ .re = 0.0, .im = 0.0 },
        Complex{ .re = -1.0, .im = 0.0 },
        Complex{ .re = 0.0, .im = 0.0 },
    };
    try ifftInPlace(allocator, &data);
    try expectApproxEqRel(data[0].re, 0.0, TEST_TOLERANCE);
}

test "IFFT edge cases" {
    var empty: [0]Complex = undefined;
    try ifftInPlace(std.heap.page_allocator, empty[0..]); // 应不报错

    var one = [_]Complex{Complex{ .re = 7.0, .im = 0.0 }};
    try ifftInPlace(std.heap.page_allocator, one[0..]);
    try expectApproxEqRel(one[0].re, 7.0, TEST_TOLERANCE);
}
