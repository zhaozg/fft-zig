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

test "Utility functions" {
    try expect(isPowerOfTwo(128));
    try expect(!isPowerOfTwo(100));
    try expect(isPowerOfFour(16));
    try expect(!isPowerOfFour(12));
    try expect(nextPowerOfTwo(9) == 16);
    try expect(nextPowerOfTwo(16) == 16);
}

/// FFT 输入参数校验，n为长度，radix为基数（2或4）
pub fn checkFftInput(n: usize, radix: usize) !void {
    if (n <= 1) return error.InvalidSize;
    if (radix == 2 and !isPowerOfTwo(n)) return error.InvalidSize;
    if (radix == 4 and !isPowerOfFour(n)) return error.InvalidSize;
}

/// 归一化复数数组
pub fn normalize(data: anytype, n: usize) void {
    for (data) |*v| {
        v.re /= @as(f64, @floatFromInt(n));
        v.im /= @as(f64, @floatFromInt(n));
    }
}

/// 计算幅值谱
pub fn calcMagnitude(data: anytype, out: []f64) void {
    const std = @import("std");
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
