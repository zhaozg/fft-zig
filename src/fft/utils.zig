//! 常用工具函数
pub fn isPowerOfTwo(n: usize) bool {
    return n > 0 and (n & (n - 1)) == 0;
}

pub fn isPowerOfFour(n: usize) bool {
    if (!isPowerOfTwo(n)) return false;
    return if (@sizeOf(usize) == 8)
        (n & 0x5555555555555555) != 0
    else
        return (n & 0x55555555) != 0;
}

pub fn nextPowerOfTwo(n: usize) usize {
    if (isPowerOfTwo(n)) return n;
    var power: usize = 1;
    while (power < n) power <<= 1;
    return power;
}

const std = @import("std");
const expect = std.testing.expect;
const expectApproxEqRel = std.testing.expectApproxEqRel;

test "Utility functions" {
    try expect(isPowerOfTwo(128));
    try expect(!isPowerOfTwo(100));
    try expect(isPowerOfFour(16));
    try expect(!isPowerOfFour(12));
    try expect(nextPowerOfTwo(9) == 16);
    try expect(nextPowerOfTwo(16) == 16);
}

fn testNormalizeGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    var data = try allocator.alloc(std.math.Complex(T), 4);
    defer allocator.free(data);
    data[0] = std.math.Complex(T){ .re = @as(T, 4.0), .im = @as(T, 8.0) };
    data[1] = std.math.Complex(T){ .re = @as(T, 2.0), .im = @as(T, 4.0) };
    data[2] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
    data[3] = std.math.Complex(T){ .re = @as(T, 6.0), .im = @as(T, 2.0) };
    normalize(T, data, 2);
    try expectApproxEqRel(data[0].re, @as(T, 2.0), @as(T, 1e-6));
    try expectApproxEqRel(data[0].im, @as(T, 4.0), @as(T, 1e-6));
}

test "Normalize f32" { try testNormalizeGeneric(f32); }
test "Normalize f64" { try testNormalizeGeneric(f64); }

fn testCalcMagnitudeGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    var data = try allocator.alloc(std.math.Complex(T), 2);
    defer allocator.free(data);
    const out = try allocator.alloc(T, 2);
    defer allocator.free(out);
    data[0] = std.math.Complex(T){ .re = @as(T, 3.0), .im = @as(T, 4.0) };
    data[1] = std.math.Complex(T){ .re = @as(T, 0.0), .im = @as(T, 0.0) };
    calcMagnitude(T, data, out);
    try expectApproxEqRel(out[0], @as(T, 5.0), @as(T, 1e-6));
    try expectApproxEqRel(out[1], @as(T, 0.0), @as(T, 1e-6));
}

test "CalcMagnitude f32" { try testCalcMagnitudeGeneric(f32); }
test "CalcMagnitude f64" { try testCalcMagnitudeGeneric(f64); }

/// FFT 输入参数校验，n为长度，radix为基数（2或4）
pub fn checkFftInput(n: usize, radix: usize) !void {
    if (n <= 1) return error.InvalidSize;
    if (radix == 2 and !isPowerOfTwo(n)) return error.InvalidSize;
    if (radix == 4 and !isPowerOfFour(n)) return error.InvalidSize;
}

/// 归一化复数数组
pub fn normalize(comptime T: type, data: []std.math.Complex(T), n: usize) void {
    for (data) |*v| {
        v.re /= @as(T, @floatFromInt(n));
        v.im /= @as(T, @floatFromInt(n));
    }
}

/// 计算幅值谱
pub fn calcMagnitude(comptime T: type, data: []const std.math.Complex(T), out: []T) void {
    const math = std.math;
    for (data, 0..) |v, i| {
        out[i] = math.sqrt(v.re * v.re + v.im * v.im);
    }
}

/// 通用bit反转排列，支持任意基数（如2、4）
pub fn bitReversePermuteGeneric(comptime T: type, data: []T, radix: usize) void {
    const n = data.len;
    if (n <= 1) return;
    for (0..n) |i| {
        var j: usize = 0;
        var temp_i = i;
        var temp_n = n;
        while (temp_n > 1) {
            j = j * radix + (temp_i % radix);
            temp_i /= radix;
            temp_n /= radix;
        }
        if (i < j) {
            const temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
    }
}
