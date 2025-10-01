//! Twiddle Factor 表与相关工具

const std = @import("std");
const math = std.math;
const Complex = @import("types.zig").Complex;

pub fn TwiddleFactorTable(comptime N: usize) type {
    return struct {
        const Self = @This();
        pub const twiddle_factors: [N / 2]Complex = init: {
            var factors: [N / 2]Complex = undefined;
            for (&factors, 0..) |*factor, k| {
                const angle = -2.0 * math.pi * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(N));
                factor.* = Complex{
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
        pub fn getTwiddle(k: usize) Complex {
            return twiddle_factors[k];
        }
        pub fn getBitReverse(i: usize) usize {
            return bit_reverse_table[i];
        }
    };
}

/// 生成单个twiddle factor
pub fn genTwiddleFactor(theta: f64, k: usize) Complex {
    return Complex{
        .re = math.cos(theta * @as(f64, @floatFromInt(k))),
        .im = math.sin(theta * @as(f64, @floatFromInt(k))),
    };
}

/// 批量生成twiddle table
pub fn genTwiddleTable(theta: f64, count: usize, out: []Complex) void {
    for (0..count) |k| {
        out[k] = genTwiddleFactor(theta, k);
    }
}

pub const TwiddleTable = struct {
    cos_table: []f64,
    sin_table: []f64,
    size: usize,
    pub fn init(allocator: std.mem.Allocator, table_size: usize) !TwiddleTable {
        const cos_table = try allocator.alloc(f64, table_size);
        const sin_table = try allocator.alloc(f64, table_size);
        for (0..table_size) |k| {
            const angle = -2.0 * math.pi * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(table_size * 2));
            cos_table[k] = math.cos(angle);
            sin_table[k] = math.sin(angle);
        }
        return TwiddleTable{
            .cos_table = cos_table,
            .sin_table = sin_table,
            .size = table_size,
        };
    }
    pub fn deinit(self: *TwiddleTable, allocator: std.mem.Allocator) void {
        allocator.free(self.cos_table);
        allocator.free(self.sin_table);
    }
    pub fn get(self: *const TwiddleTable, k: usize, n: usize) Complex {
        const index = (k * self.size) / (n / 2);
        return Complex{
            .re = self.cos_table[index],
            .im = self.sin_table[index],
        };
    }
};

test "Twiddle factor table" {
    const TwiddleTable16 = TwiddleFactorTable(16);
    const twiddle_0 = TwiddleTable16.twiddle_factors[1];
    try expectApproxEqRel(@cos(-2.0 * std.math.pi / 16.0), twiddle_0.re, TEST_TOLERANCE);
}

const expectApproxEqRel = std.testing.expectApproxEqRel;
const TEST_TOLERANCE = 1e-12;

test "TwiddleFactorTable edge cases" {
    // N=2, only twiddle_factors[0] exists, angle=0
    const TwiddleTable2 = TwiddleFactorTable(2);
    try std.testing.expectApproxEqRel(@cos(0.0), TwiddleTable2.twiddle_factors[0].re, 1e-12);

    // N=4, twiddle_factors[0] (angle=0), [1] (angle=-π/2)
    const TwiddleTable4 = TwiddleFactorTable(4);
    try std.testing.expectApproxEqRel(@cos(0.0), TwiddleTable4.twiddle_factors[0].re, 1e-12);
    try std.testing.expectApproxEqRel(@cos(-2.0 * std.math.pi * 1 / 4.0), TwiddleTable4.twiddle_factors[1].re, 1e-12);

    // N=8, test angle=0, -π/2, -π
    const TwiddleTable8 = TwiddleFactorTable(8);
    try std.testing.expectApproxEqRel(@cos(0.0), TwiddleTable8.twiddle_factors[0].re, 1e-12);
    try std.testing.expectApproxEqRel(@cos(-2.0 * std.math.pi * 2 / 8.0), TwiddleTable8.twiddle_factors[2].re, 1e-12); // -π/2
}
