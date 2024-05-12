# Wanderer

> Wanderer problem, also known as *Random Walk* or *Drunkard's Walk*, is a random process that describes a path that consists of a succession of random steps on some mathematical space.
> -- <cite>[Wikipedia][1]</cite>

[1]: https://en.wikipedia.org/wiki/Random_walk

This repo is one of the projects for the *Numerical Algorithms* course I'm taking at my University. Besides solving the wanderer problem it also implements the following algorithms for solving systems of linear equations:

| Algorithm \ Representation           | 2D Vector | HashMap Sparse |
|--------------------------------------|-----------|----------------|
| Jacobi                               | ✅         | ✅              |
| Gauss-Seidel                         | ✅         | ✅              |
| Gauss Elimination (without pivoting) | ✅         | ✅              |
| Gauss Elimination (partial pivot)    | ✅         | ✅              |

# Requirements
 - [Rust](https://www.rust-lang.org/)
 - [Cargo](https://doc.rust-lang.org/cargo/)
 - [Python 3+](https://www.python.org/)
 - [matplotlib](https://matplotlib.org/)
 - Linux OS (preferably)

# Quick Setup

Download this repository using:
```bash
$ git clone https://github.com/wedkarz02/wanderer.git
```
or use the *Download ZIP* option from the GitHub repository page.

To compile the project, use *cargo build*. I highly recommend compiling in *--release* mode due to better execution speed.

```bash
$ cargo build --release
```
Run the executable directly or with *cargo run --release*. See [src/main.rs](https://github.com/wedkarz02/wanderer/blob/main/src/main.rs) for avaliable command line arguments.

# License

If not directly stated otherwise, everything in this project is under the MIT License. See the [LICENSE](https://github.com/wedkarz02/wanderer/blob/main/LICENSE) file for more info.

# Cool links and references
 - Numerical Algorithms lecture presentations
 - https://en.wikipedia.org/wiki/Random_walk
 - https://math.dartmouth.edu/~doyle/docs/walks/walks.pdf
 - https://en.wikipedia.org/wiki/Jacobi_method
 - https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
 - https://en.wikipedia.org/wiki/Gaussian_elimination
