# BlueSphere
Blue Noise Sampling of Sphere Surface.

# Implementation

This code computes using 16 cores, and uses 16-wide SIMD (AVX512) on each core, for conceptually computing 256 candidates "simultaneously."


[![Watch the video](https://img.youtube.com/vi/2ZWt0uVzpAQ/maxresdefault.jpg)](https://youtu.be/2ZWt0uVzpAQ)

# Requirements

Requires AVX512F to run. Not every CPU has this. E.g. AMD CPUs at time of writing, cannot execute this code.

# Building

```
$ make
$ ./bluesphere
```

