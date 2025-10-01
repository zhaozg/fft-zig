const std = @import("std");
const fft_radix4 = @import("src/fft/fft_radix4.zig");
const Complex = @import("src/fft/types.zig").Complex;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    
    const sizes = [_]usize{ 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576 };
    
    for (sizes) |size| {
        const freq: usize = 5;
        
        const data = try allocator.alloc(Complex, size);
        defer allocator.free(data);
        
        for (0..size) |i| {
            const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(size));
            data[i] = Complex{ 
                .re = @sin(2.0 * std.math.pi * @as(f64, @floatFromInt(freq)) * t),
                .im = 0.0
            };
        }
        
        try fft_radix4.fftRadix4SIMD(data);
        
        const mag = @sqrt(data[freq].re * data[freq].re + data[freq].im * data[freq].im);
        const expected = @as(f64, @floatFromInt(size)) / 2.0;
        const ratio = mag / expected;
        
        std.debug.print("Size {d:7}: mag={d:10.4} expected={d:10.4} ratio={d:.6}\n", .{size, mag, expected, ratio});
    }
}
