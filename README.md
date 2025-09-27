# fft-zig

## 依赖关系图

```
types.zig
   ↑
   │
 ┌───────────────┬───────────────┬───────────────┬───────────────┐
 │               │               │               │               │
fft_radix2.zig  fft_radix4.zig  utils.zig      base.zig
   ↑               ↑               ↑               ↑
   │               │               │               │
   └───────┬───────┴───────┬───────┴───────┬───────┘
           │               │               │
      fft_parallel.zig     │           ifft.zig
           ↑               │               ↑
           │               │               │
      fft_mixed.zig────────┘               │
           ↑                               │
           └───────────────┬───────────────┘
                           │
                      fft_r2c.zig
```

fft.zig → (types, utils, fft_radix2, fft_radix4, fft_parallel, fft_mixed, ifft, fft_r2c, base)
