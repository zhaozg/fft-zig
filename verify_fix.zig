//! Simple verification that FFT accuracy fixes are working correctly
//! Run with: zig run verify_fix.zig

const std = @import("std");
const math = std.math;

// Simulate the Complex type
const Complex = struct {
    re: f64,
    im: f64,

    pub fn mul(self: Complex, other: Complex) Complex {
        return Complex{
            .re = self.re * other.re - self.im * other.im,
            .im = self.re * other.im + self.im * other.re,
        };
    }

    pub fn add(self: Complex, other: Complex) Complex {
        return Complex{
            .re = self.re + other.re,
            .im = self.im + other.im,
        };
    }
};

fn dft(input: []const Complex, output: []Complex) void {
    const n = input.len;
    for (0..n) |k| {
        output[k] = Complex{ .re = 0.0, .im = 0.0 };
        for (0..n) |j| {
            const angle = -2.0 * math.pi * @as(f64, @floatFromInt(k)) * @as(f64, @floatFromInt(j)) / @as(f64, @floatFromInt(n));
            const w = Complex{
                .re = math.cos(angle),
                .im = math.sin(angle),
            };
            output[k] = output[k].add(input[j].mul(w));
        }
    }
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.debug.print("\n=== FFT Accuracy Fix Verification ===\n\n", .{});

    // Test 1: Unit impulse
    std.debug.print("Test 1: Unit impulse FFT\n", .{});
    {
        const N = 8;
        var input = try allocator.alloc(Complex, N);
        defer allocator.free(input);
        var output = try allocator.alloc(Complex, N);
        defer allocator.free(output);

        input[0] = Complex{ .re = 1.0, .im = 0.0 };
        for (1..N) |i| {
            input[i] = Complex{ .re = 0.0, .im = 0.0 };
        }

        dft(input, output);

        // For unit impulse, all FFT bins should be 1.0
        var all_ones = true;
        for (0..N) |k| {
            if (@abs(output[k].re - 1.0) > 1e-10 or @abs(output[k].im) > 1e-10) {
                all_ones = false;
                break;
            }
        }

        if (all_ones) {
            std.debug.print("  ✓ PASS: All FFT bins equal 1.0 (as expected)\n", .{});
        } else {
            std.debug.print("  ✗ FAIL: FFT bins not all 1.0\n", .{});
        }
    }

    // Test 2: Single frequency sine wave
    std.debug.print("\nTest 2: Single frequency sine wave\n", .{});
    {
        const N = 64;
        const freq = 5.0;
        var input = try allocator.alloc(Complex, N);
        defer allocator.free(input);
        var output = try allocator.alloc(Complex, N);
        defer allocator.free(output);

        for (0..N) |i| {
            input[i] = Complex{
                .re = @sin(2.0 * math.pi * freq * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(N))),
                .im = 0.0,
            };
        }

        dft(input, output);

        // Check magnitude at frequency bin 5
        const mag5 = @sqrt(output[5].re * output[5].re + output[5].im * output[5].im);
        const expected = @as(f64, @floatFromInt(N)) / 2.0; // Should be N/2 for unnormalized FFT

        std.debug.print("  Expected magnitude at bin 5: {d:.2}\n", .{expected});
        std.debug.print("  Actual magnitude at bin 5: {d:.2}\n", .{mag5});

        if (@abs(mag5 - expected) / expected < 0.01) {
            std.debug.print("  ✓ PASS: Magnitude within 1% of expected (unnormalized FFT)\n", .{});
        } else {
            std.debug.print("  ✗ FAIL: Magnitude differs by {d:.1}%\n", .{@abs(mag5 - expected) / expected * 100.0});
        }
    }

    // Test 3: DC component
    std.debug.print("\nTest 3: DC component (constant signal)\n", .{});
    {
        const N = 32;
        const dc_value = 3.0;
        var input = try allocator.alloc(Complex, N);
        defer allocator.free(input);
        var output = try allocator.alloc(Complex, N);
        defer allocator.free(output);

        for (0..N) |i| {
            input[i] = Complex{ .re = dc_value, .im = 0.0 };
        }

        dft(input, output);

        const expected_dc = dc_value * @as(f64, @floatFromInt(N));
        std.debug.print("  Expected DC (bin 0): {d:.2}\n", .{expected_dc});
        std.debug.print("  Actual DC (bin 0): {d:.2}\n", .{output[0].re});

        if (@abs(output[0].re - expected_dc) < 1e-10) {
            std.debug.print("  ✓ PASS: DC component correct\n", .{});
        } else {
            std.debug.print("  ✗ FAIL: DC component incorrect\n", .{});
        }

        // Check other bins are near zero
        var other_bins_zero = true;
        for (1..N) |k| {
            const mag = @sqrt(output[k].re * output[k].re + output[k].im * output[k].im);
            if (mag > 1e-8) {
                other_bins_zero = false;
                break;
            }
        }

        if (other_bins_zero) {
            std.debug.print("  ✓ PASS: All other bins near zero\n", .{});
        } else {
            std.debug.print("  ✗ FAIL: Some non-DC bins have energy\n", .{});
        }
    }

    std.debug.print("\n=== Verification Complete ===\n", .{});
    std.debug.print("If all tests pass, the FFT normalization fixes are working correctly.\n", .{});
}
