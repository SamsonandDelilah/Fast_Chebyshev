Here is the complete, high-density documentation in English. It maintains the exact same rigorous mathematical depth, structural clarity, and tailored guidance for beginners.
You can easily copy this Markdown text and convert it into a PDF using any converter (such as Pandoc, VS Code Markdown PDF, or an online tool) to create a professional reference manual for your repository.

------------------------------
## PDF Documentation: Chebyshev Approximation & Clenshaw Toolchain
A Comprehensive Guide for Numerical Beginners, Statistical Modeling, and High-Performance Computing
------------------------------
## 1. Why Chebyshev Approximation? (Motivation)
Function approximation is a cornerstone of numerical mathematics. In many scientific fields, computationally expensive or analytically complex functions (such as $e^x$, $\sin(x)$, or complex probability density functions used in statistical CDFs) must be evaluated millions of times within simulations or real-time architectures.
Chebyshev polynomials solve this bottleneck. They allow us to substitute a continuous, complex function f(x) over a bounded interval $[a,b]$ with an N-degree polynomial that can be evaluated with incredibly low latency while maintaining global precision.
## When is Chebyshev Advantageous?

* Minimizing the Maximum Error (Minimax Principle): Unlike Taylor polynomials, which are highly accurate around a single expansion point but explode exponentially near the borders, Chebyshev approximations distribute the error evenly across the entire interval. Mathematically, they come remarkably close to the absolute optimal "Minimax Polynomial".
* Eliminating Runge's Phenomenon: In classical polynomial interpolation using equally spaced points, high-degree polynomials tend to oscillate wildly near the boundaries of the interval. Chebyshev nodes naturally bunch up at the borders, completely neutralizing this instability.
* Analytical Toolchain Integration: Once the Chebyshev coefficients are generated, you can analytically derive the exact coefficients for the derivative or integral of the function. This entirely bypasses the need for numerical differentiation or integration, which are notoriously prone to noise and truncation errors.

## When is it NOT Advantageous?

* Discontinuities and Sharp Kinks: If a function possesses jumps, singularities, or non-differentiable sharp bends (like a ReLU function or $\vert{}x\vert{}$), Chebyshev series suffer from slow convergence and ringing artifacts near the kink (Gibbs phenomenon).
* Infinite Intervals: The standard framework is strictly bounded to a compact, finite interval $[a,b]$. To approximate functions over infinite or semi-infinite domains (e.g., $[0, \infty)$), one must deploy modified, rational Chebyshev methods.

------------------------------
## 2. Mathematical Formulations of the Architecture
The entire framework relies on mapping the real-world physical domain $x \in [a,b]$ onto the canonical Chebyshev domain $y \in [-1, 1]$.
## 2.1 Domain Transformation
The bijective mapping from the physical coordinate x to the normalized coordinate y (and vice versa) is defined as:
$$y = \frac{2x - a - b}{b - a} \quad \Longleftrightarrow \quad x = \frac{a+b}{2} + \frac{b-a}{2} \cdot y$$ 
## 2.2 Coefficient Generation via Node Sampling
To compute the polynomial coefficients $c_j$, the function is sampled at the roots (zeros) of the Chebyshev polynomials. These Chebyshev nodes are geometrically clustered towards the boundaries of the domain:
$$y_k = \cos\left( \frac{\pi (k + 0.5)}{N} \right) \quad \text{for } k = 0, 1, \dots, N-1$$ 
The corresponding function samples are computed as $f_k = f(x(y_k))$. The array of coefficients $c_j$ is then extracted via a Discrete Cosine Transformation (DCT):
$$c_j = \frac{2}{N} \sum_{k=0}^{N-1} f_k \cdot \cos\left( \frac{\pi j (k + 0.5)}{N} \right) \quad \text{for } j = 0, 1, \dots, N-1$$ 
Note on Normalization: In standard polynomial expansions, the zero-th coefficient c₀ is conventionally weighted by a factor of 0.5 (hence the coeffs *= 0.5; scaling step in the source code).
## 2.3 Clenshaw's Backward Recurrence Evaluation (eval)
To evaluate the polynomial safely without explicitly constructing unstable, rounding-error-prone powers of $x^j$, we implement Clenshaw's recurrence method. We initialize $d_N = 0$ and $d_{N+1} = 0$, then iterate backward from j = N-1 down to j = 1:
$$d_j = 2y \cdot d_{j+1} - d_{j+2} + c_j$$ 
The final function value at physical coordinate x (mapped to y) is synthesized in the very last step using the accumulated state variables:
$$f(x) \approx y \cdot d_1 - d_2 + c_0$$ 
## 2.4 Analytical Derivative Generation (derive)
The coefficients of the exact analytical derivative $c'_j$ can be synthesized directly from the original coefficients $c_j$ using a backward recurrence relationship:
$$c'_{j} = c'_{j+2} + 2 \cdot (j+1) \cdot c_{j+1} \quad \text{for } j = N-2, N-3, \dots, 1$$ 
The boundary condition for the zero-th derivative coefficient is defined as: c'₀ = 0.5 ⋅ c'₂ + c₁. Because this derivative is calculated in the normalized $[-1,1]$ space, it must be scaled by the inner derivative of our physical domain transformation:
$$\text{scale} = \frac{2}{b-a} \quad \Longrightarrow \quad c'_j \leftarrow c'_j \cdot \text{scale}$$ 
## 2.5 Analytical Integral / Antiderivative Generation (integrate)
The coefficients of the exact analytical integral $C_j$ are generated via a forward relationship. They scale inversely proportional to the frequency step j:
$$C_j = \frac{c_{j-1} - c_{j+1}}{2j} \quad \text{for } j = 1, 2, \dots, N-2$$ 
Boundary handling for the highest frequency is given by: $C_{N-1} = \frac{c_{N-2}}{2(N-1)}$. The constant of integration C₀ is dynamically shifted to guarantee the boundary condition F(a) = 0 (clamping the integral to zero at the left bound, which is mandatory for building precise Cumulative Distribution Functions). To achieve this, the integral polynomial is evaluated at x = a (where y = -1) via Clenshaw's recurrence, and the resulting offset is subtracted from C₀. Finally, the scaling factor of the transformation is applied:
$$\text{scale} = \frac{b-a}{2} \quad \Longrightarrow \quad C_j \leftarrow C_j \cdot \text{scale}$$ 
------------------------------
## 3. Recommended Textbooks & Academic Resources
For newcomers looking to master spectral methods, orthogonal polynomials, and numerical approximation theory, these textbooks represent the international gold standard:
## Core English Literature (The Global Benchmarks)

   1. "Trefethen, L. N. (2013). Approximation Practice and Chebfun Guide." SIAM.
   * The definitive, modern masterpiece on the subject. Lloyd N. Trefethen is the pioneer of modern Chebyshev-driven numerics. The book is incredibly intuitive, avoids dry, dense lemma-dumping, and masterfully explains why Chebyshev methods dominate high-end numerical tools.
   2. "Boyd, J. P. (2001). Chebyshev and Fourier Spectral Methods." Dover Publications.
   * The undisputed classic. Perfect for engineers, computer scientists, and software developers. It focuses heavily on how to convert these equations into actual computer source code. (Bonus: The author has made this textbook officially available as a free PDF download through his university portal).
   3. "Press, W. H., et al. (2007). Numerical Recipes: The Art of Scientific Computing (3rd Edition)." Cambridge University Press.
   * The programmer’s bible. Chapter 5.8 is exclusively dedicated to Chebyshev approximation and Clenshaw recurrence. Excellent for understanding how our Rust code mirrors traditional scientific libraries while achieving safe, bare-metal acceleration.
   
## Complementary Specialized Reading

   1. "Rivlin, T. J. (1990). Chebyshev Polynomials: From Approximation Theory to Algebra and Number Theory." Wiley-Interscience.
   * Ideal for users who want to dive deep into the underlying pure algebraic and structural properties of Chebyshev polynomials beyond engineering.
   2. "Mason, J. C., & Handscomb, D. C. (2002). Chebyshev Polynomials." CRC Press.
   * A rigorous, comprehensive text that covers both the continuous and discrete theories of these algorithms, providing highly detailed proofs for derivative and integration matrices.
   