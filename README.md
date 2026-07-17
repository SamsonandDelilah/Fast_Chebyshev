## 📖 Documentation / Dokumentation

For a deep dive into the mathematical formulations and implementation theory, check out our full documentation:
 - **[English Comprehensive Guide (PDF-ready)](https://github.com)**
 - **[Deutsche Kurzerklärung (PDF-bereit)](https://github.com)**
 
# Fast_Chebyshev

A blazing-fast, zero-overhead Chebyshev polynomial approximation toolchain written in pure Rust.

Designed for high-performance computing, numerical simulations, and statistical modeling (such as ultra-fast CDF evaluations), `Fast_Chebyshev` features a compile-time optimized variant that completely eliminates heap allocations and fully leverages CPU instruction pipelines.

## 🚀 Key Features

* **Zero-Allocation Stack Turbo:** Use `OptimizedChebyshev<N>` with Const Generics to evaluate functions with **0 bytes of heap memory**.
* **Microsecond-Level Performance:** Evaluates polynomial expansions, analytical derivatives, and integrals in under **20 nanoseconds**.
* **Analytical Derivatives & Integrals:** Instantly generate new Chebyshev structs representing the exact derivative $f'(x)$ or the definite integral $F(x)$ without numerical instability.
* **Clenshaw's Algorithm:** Uses the highly stable Clenshaw recurrence method for all polynomial evaluations, preventing precision loss at high frequencies.
* **100% Safe & Idiomatic Rust:** Fully compatible with stable Rust, zero warning output under strict `cargo clippy` audits.

---

## 📊 Benchmarks

*Benchmarks were executed using `Criterion` on an x86_64 CPU with native optimizations (`-C target-cpu=native`) evaluating a 16-degree polynomial function.*

| Implementation               | Evaluation Method | Allocation Type     | Execution Time |
| :--------------------------- | :---------------- | :------------------ | :------------- |
| **`OriginalChebyshev`**      | Dynamic Range     | Heap (`Vec<f64>`)   | **~26.3 ns**   |
| **`OptimizedChebyshev<16>`** | Const Generics    | Stack (`[f64; 16]`) | **~19.4 ns**   |

### Why is it so fast?
By baking the interval scaling factor directly into the Clenshaw coefficients during generation, `Fast_Chebyshev` avoids secondary array traversals. Combined with Rust's compile-time loop unrolling via Const Generics, data flows directly through the CPU's AVX registers without cache thrashing.

---

## 🛠️ Usage

Add `Fast_Chebyshev` to your local environment or include the source files into your architecture.

### 1. Function Approximation & Evaluation

```rust
use std::f64::consts::PI;
use cheby::OptimizedChebyshev;

fn main() {
    let a = 0.0;
    let b = PI / 2.0;
    const N: usize = 16;

    // Approximate the sinus function over [0, PI/2]
    let approx = OptimizedChebyshev::<N>::new(a, b, |x| x.sin());

    // Evaluate the polynomial at any point in under 20ns
    let test_x = PI / 4.0;
    let result = approx.eval(test_x);
    
    println!("Approximated sin(pi/4): {}", result);
    println!("Exact CPU sin(pi/4):        {}", test_x.sin());
}
```

### 2. Generating the Analytical Derivative $f'(x)$

```rust
// Generate a new, standalone Chebyshev struct representing the exact derivative
let derivative = approx.derive();

// Evaluates directly to cos(x) in ~19.4ns!
let slope = derivative.eval(test_x); 
println!("Approximated slope (cos): {}", slope);
```

### 3. Generating the Definite Integral $F(x)$ (Perfect for CDFs)

```rust
// Generate the analytical integral curve F(x) where F(a) = 0
let integral = approx.integrate();

// Perfect for cumulative distribution functions (CDFs)
let area = integral.eval(test_x); 
println!("Area under the curve: {}", area);
```

---

## 📐 Mathematical Background

Chebyshev approximation works by mapping an arbitrary interval $[a, b]$ onto the canonical domain $[-1, 1]$ via:

$$y = \frac{2x - a - b}{b - a}$$

The function is sampled at the **Chebyshev nodes (roots)** to minimize the Runge phenomenon (maximum error bound oscillation):

$$x_k = \cos\left(\frac{\pi (k + 0.5)}{N}\right)$$

Evaluation is performed via **Clenshaw's Recurrence Formula**, a backwards recurrence structure that computes the polynomial sequence without explicitly constructing powers of $x$:

$$d_j = 2y \cdot d_{j+1} - d_{j+2} + c_j$$

---

## 📜 License

This project is licensed under the MIT License - feel free to use it in commercial, scientific, or private software stacks.
