#!/bin/bash
# Helper script to format and test FFT code after accuracy fixes
# 
# This script should be run after the FFT accuracy fixes to:
# 1. Format the modified source files
# 2. Build the project
# 3. Run all tests
#
# Requirements: Zig 0.14.1 or 0.15.1

set -e  # Exit on error

echo "=== Formatting modified FFT source files ==="
zig fmt src/fft/fft_radix2.zig
zig fmt src/fft/fft_r2c.zig
zig fmt src/fft/fft_mixed.zig

echo ""
echo "=== Building project ==="
zig build

echo ""
echo "=== Running tests ==="
zig build test

echo ""
echo "=== All checks passed! ==="
echo "The FFT accuracy issues have been fixed:"
echo "  - Forward FFT no longer incorrectly normalizes"
echo "  - fft_r2c now properly converts output for all sizes"
echo "  - Bluestein algorithm chirp mirroring fixed"
echo "  - Tests now validate correct FFT magnitude"
