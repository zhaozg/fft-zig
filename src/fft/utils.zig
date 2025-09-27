//! 常用工具函数
pub fn isPowerOfTwo(n: usize) bool {
    return n > 0 and (n & (n - 1)) == 0;
}

pub fn isPowerOfFour(n: usize) bool {
    if (!isPowerOfTwo(n)) return false;
    return if(@sizeOf(usize) == 8)
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
