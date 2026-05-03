//! Twiddle Factor 表与相关工具

const std = @import("std");
const math = std.math;

/// 泛型 Twiddle Factor 表
pub fn TwiddleFactorTable(comptime T: type, comptime N: usize) type {
    return struct {
        const Self = @This();
        pub const twiddle_factors: [N / 2]std.math.Complex(T) = init: {
            var factors: [N / 2]std.math.Complex(T) = undefined;
            for (&factors, 0..) |*factor, k| {
                const angle = -2.0 * math.pi * @as(T, @floatFromInt(k)) / @as(T, @floatFromInt(N));
                factor.* = std.math.Complex(T){
                    .re = math.cos(angle),
                    .im = math.sin(angle),
                };
            }
            break :init factors;
        };
        pub const bit_reverse_table: [N]usize = init: {
            var table: [N]usize = undefined;
            for (&table, 0..) |*entry, i| {
                entry.* = bitReverse(i, log2Int(N));
            }
            break :init table;
        };
        fn bitReverse(x: usize, bits: usize) usize {
            var result: usize = 0;
            var temp = x;
            for (0..bits) |_| {
                result = (result << 1) | (temp & 1);
                temp >>= 1;
            }
            return result;
        }
        fn log2Int(n: usize) usize {
            return @ctz(@as(u64, @intCast(n)));
        }
        pub fn getTwiddle(k: usize) std.math.Complex(T) {
            return twiddle_factors[k];
        }
        pub fn getBitReverse(i: usize) usize {
            return bit_reverse_table[i];
        }
    };
}

/// 生成单个twiddle factor
pub fn genTwiddleFactor(comptime T: type, theta: T, k: usize) std.math.Complex(T) {
    return std.math.Complex(T){
        .re = math.cos(theta * @as(T, @floatFromInt(k))),
        .im = math.sin(theta * @as(T, @floatFromInt(k))),
    };
}

/// 批量生成twiddle table
pub fn genTwiddleTable(comptime T: type, theta: T, count: usize, out: []std.math.Complex(T)) void {
    for (0..count) |k| {
        out[k] = genTwiddleFactor(T, theta, k);
    }
}

/// 运行时 Twiddle 表
pub fn TwiddleTable(comptime T: type) type {
    return struct {
        cos_table: []T,
        sin_table: []T,
        size: usize,
        pub fn init(allocator: std.mem.Allocator, table_size: usize) !@This() {
            const cos_table = try allocator.alloc(T, table_size);
            const sin_table = try allocator.alloc(T, table_size);
            for (0..table_size) |k| {
                const angle = -2.0 * math.pi * @as(T, @floatFromInt(k)) / @as(T, @floatFromInt(table_size * 2));
                cos_table[k] = math.cos(angle);
                sin_table[k] = math.sin(angle);
            }
            return @This(){
                .cos_table = cos_table,
                .sin_table = sin_table,
                .size = table_size,
            };
        }
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            allocator.free(self.cos_table);
            allocator.free(self.sin_table);
        }
        pub fn get(self: *const @This(), k: usize, n: usize) std.math.Complex(T) {
            const index = (k * self.size) / (n / 2);
            return std.math.Complex(T){
                .re = self.cos_table[index],
                .im = self.sin_table[index],
            };
        }
    };
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const expectApproxEqAbs = std.testing.expectApproxEqAbs;
const expect = std.testing.expect;

/// 测试泛型 Twiddle Factor 表
fn testTwiddleFactorTableGeneric(comptime T: type) !void {
    const F32_TOLERANCE: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);

    // N=2, only twiddle_factors[0] exists, angle=0
    const TwiddleTable2 = TwiddleFactorTable(T, 2);
    try expectApproxEqRel(@as(T, @cos(0.0)), TwiddleTable2.twiddle_factors[0].re, F32_TOLERANCE);

    // N=4, twiddle_factors[0] (angle=0), [1] (angle=-π/2)
    const TwiddleTable4 = TwiddleFactorTable(T, 4);
    try expectApproxEqRel(@as(T, @cos(0.0)), TwiddleTable4.twiddle_factors[0].re, F32_TOLERANCE);
    try expectApproxEqAbs(@as(T, @cos(-2.0 * std.math.pi * 1 / 4.0)), TwiddleTable4.twiddle_factors[1].re, F32_TOLERANCE);

    // N=8, test angle=0, -π/2, -π
    const TwiddleTable8 = TwiddleFactorTable(T, 8);
    try expectApproxEqRel(@as(T, @cos(0.0)), TwiddleTable8.twiddle_factors[0].re, F32_TOLERANCE);
    try expectApproxEqAbs(@as(T, @cos(-2.0 * std.math.pi * 2 / 8.0)), TwiddleTable8.twiddle_factors[2].re, F32_TOLERANCE); // -π/2
}

test "Twiddle factor table f32" {
    try testTwiddleFactorTableGeneric(f32);
}

test "Twiddle factor table f64" {
    try testTwiddleFactorTableGeneric(f64);
}

test "TwiddleFactorTable edge cases" {
    // 测试 f32 和 f64
    try testTwiddleFactorTableGeneric(f32);
    try testTwiddleFactorTableGeneric(f64);
}

/// 测试泛型 TwiddleTable
fn testTwiddleTableGeneric(comptime T: type) !void {
    const allocator = std.testing.allocator;
    const F32_TOLERANCE: T = if (T == f32) @as(T, 1e-6) else @as(T, 1e-12);

    const table_size = 16;
    var table = try TwiddleTable(T).init(allocator, table_size);
    defer table.deinit(allocator);

    // 测试获取 twiddle factor
    const twiddle = table.get(1, 8);
    const expected_angle = -2.0 * std.math.pi * @as(T, @floatFromInt(1)) / @as(T, @floatFromInt(8));
    try expectApproxEqRel(@as(T, @cos(expected_angle)), twiddle.re, F32_TOLERANCE);
    try expectApproxEqRel(@as(T, @sin(expected_angle)), twiddle.im, F32_TOLERANCE);
}

test "TwiddleTable f32" {
    try testTwiddleTableGeneric(f32);
}

test "TwiddleTable f64" {
    try testTwiddleTableGeneric(f64);
}
