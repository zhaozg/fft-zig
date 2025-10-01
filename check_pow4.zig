const std = @import("std");

fn isPowerOfTwo(n: usize) bool {
    return n > 0 and (n & (n - 1)) == 0;
}

fn isPowerOfFour(n: usize) bool {
    if (!isPowerOfTwo(n)) return false;
    return if (@sizeOf(usize) == 8)
        (n & 0x5555555555555555) != 0
    else
        return (n & 0x55555555) != 0;
}

pub fn main() !void {
    const sizes = [_]usize{ 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576 };
    for (sizes) |size| {
        std.debug.print("{d}: pow2={} pow4={} (binary: {b})\n", .{size, isPowerOfTwo(size), isPowerOfFour(size), size});
    }
}
