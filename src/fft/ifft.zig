//! 逆FFT实现

const std = @import("std");

// 循环依赖修复：fftInPlace 由主接口调度，不在此直接 import

pub fn ifftInPlace(comptime T: type, allocator: std.mem.Allocator, data: []std.math.Complex(T)) !void {
    const n = data.len;
    for (data) |*v| {
        v.im = -v.im;
    }
    try @import("base.zig").fftInPlaceBase(T, allocator, data);
    for (data) |*v| {
        v.re = v.re / @as(T, @floatFromInt(n));
        v.im = -v.im / @as(T, @floatFromInt(n));
    }
}

const expectApproxEqRel = std.testing.expectApproxEqRel;

fn testIFFTCorrectnessGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);

    var data = [_]std.math.Complex(T){
        std.math.Complex(T){ .re = @as(T, 1.0), .im = @as(T, 0.0) },
        std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) },
        std.math.Complex(T){ .re = @as(T, -1.0), .im = @as(T, 0.0) },
        std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) },
    };
    try ifftInPlace(T, allocator, &data);
    try expectApproxEqRel(data[0].re, @as(T, 0.0), tolerance);
}

test "IFFT correctness f32" { try testIFFTCorrectnessGeneric(f32); }
test "IFFT correctness f64" { try testIFFTCorrectnessGeneric(f64); }

fn testIFFTEdgeGeneric(comptime T: type) !void {
    const tolerance: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);
    var empty: [0]std.math.Complex(T) = undefined;
    try ifftInPlace(T, std.heap.page_allocator, empty[0..]); // 应不报错

    var one = [_]std.math.Complex(T){std.math.Complex(T){ .re = @as(T, 7.0), .im = @as(T, 0.0) }};
    try ifftInPlace(T, std.heap.page_allocator, one[0..]);
    try expectApproxEqRel(one[0].re, @as(T, 7.0), tolerance);
}

test "IFFT edge cases f32" { try testIFFTEdgeGeneric(f32); }
test "IFFT edge cases f64" { try testIFFTEdgeGeneric(f64); }
