use std::f64::consts::PI;

// =========================================================================
// VERSION 1: Deine originale Implementierung (Dynamisch mit Vec)
// =========================================================================
pub struct OriginalChebyshev {
    a: f64,
    b: f64,
    pub coeffs: Vec<f64>,
}

impl OriginalChebyshev {
    #[allow(clippy::needless_range_loop)]
    pub fn new<F>(a: f64, b: f64, n: usize, func: F) -> Self
    where
        F: Fn(f64) -> f64,
    {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);
        let samples: Vec<f64> = (0..n).map(|k| {
            let theta = PI * (k as f64 + 0.5) / n as f64;
            mid + half * theta.cos()
        }).map(func).collect();

        let mut coeffs = vec![0.0; n];
        for j in 0..n {
            let j_f64 = j as f64;
            let mut s = 0.0;
            for k in 0..n {
                let theta = PI * j_f64 * (k as f64 + 0.5) / n as f64;
                s += samples[k] * theta.cos();
            }
            coeffs[j] = 2.0 * s / n as f64;
        }
        
        if !coeffs.is_empty() {
            // Nur das erste Element wird halbiert
            coeffs[0] *= 0.5;
        }
        
        Self { a, b, coeffs }
    }

    pub fn eval(&self, x: f64) -> f64 {
        let y = (2.0 * x - self.a - self.b) / (self.b - self.a);
        let y2 = 2.0 * y;
        let mut d = 0.0;
        let mut dd = 0.0;
        for &cj in self.coeffs.iter().skip(1).rev() {
            let tmp = d;
            d = y2 * d - dd + cj;
            dd = tmp;
        }
        y * d - dd + self.coeffs[0] // Korrigiert: Nur das erste Element addieren
    }

    #[allow(clippy::needless_range_loop)]
    pub fn derive(&self) -> Self {
        let n = self.coeffs.len();
        let mut d_coeffs = vec![0.0; n];
        
        if n < 2 {
            return Self { a: self.a, b: self.b, coeffs: d_coeffs };
        }

        let scale = 2.0 / (self.b - self.a);

        // Die Intervallskalierung wird direkt bei der Berechnung angewendet!
        for j in (1..n - 1).rev() {
            let next_plus_two = if j + 2 < n { d_coeffs[j + 2] } else { 0.0 };
            d_coeffs[j] = next_plus_two + scale * 2.0 * ((j + 1) as f64) * self.coeffs[j + 1];
        }
        
        let nxt_0 = if 2 < n { d_coeffs[2] } else { 0.0 };
        d_coeffs[0] = 0.5 * nxt_0 + scale * self.coeffs[1];

        Self { a: self.a, b: self.b, coeffs: d_coeffs }
    }

    #[allow(clippy::needless_range_loop)]
    pub fn integrate(&self) -> Self {
        let n = self.coeffs.len();
        let mut i_coeffs = vec![0.0; n];

        if n < 2 {
            return Self { a: self.a, b: self.b, coeffs: i_coeffs };
        }

        let scale = 0.5 * (self.b - self.a);

        // Vorwärts-Integration inklusive direkter Skalierung
        let c2_val = if 2 < n { self.coeffs[2] } else { 0.0 };
        i_coeffs[1] = scale * (self.coeffs[0] * 2.0 - 0.5 * c2_val);

        for j in 2..n - 1 {
            i_coeffs[j] = scale * (self.coeffs[j - 1] - self.coeffs[j + 1]) / (2.0 * j as f64);
        }

        i_coeffs[n - 1] = scale * self.coeffs[n - 2] / (2.0 * (n - 1) as f64);

        // Integrationskonstante bestimmen
        let mut d = 0.0;
        let mut dd = 0.0;
        for &cj in i_coeffs.iter().skip(1).rev() {
            let tmp = d;
            d = -2.0 * d - dd + cj;
            dd = tmp;
        }
        i_coeffs[0] = d + dd;

        Self { a: self.a, b: self.b, coeffs: i_coeffs }
    }
}

// =========================================================================
// VERSION 2: Die optimierte Implementierung (Flexibel für jedes N)
// =========================================================================
pub struct OptimizedChebyshev<const N: usize> {
    a: f64,
    b: f64,
    pub coeffs: [f64; N],
}

impl<const N: usize> OptimizedChebyshev<N> {
    #[allow(clippy::needless_range_loop)]
    pub fn new<F>(a: f64, b: f64, func: F) -> Self
    where
        F: Fn(f64) -> f64,
    {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);
        let n_f64 = N as f64;
        let inv_n = 1.0 / n_f64;

        let mut samples = [0.0; N];
        let mut theta_k = [0.0; N];

        for k in 0..N {
            let theta = PI * (k as f64 + 0.5) * inv_n;
            theta_k[k] = theta;
            let x = mid + half * theta.cos();
            samples[k] = func(x);
        }

        let mut coeffs = [0.0; N];
        for j in 0..N {
            let j_f64 = j as f64;
            let mut s = 0.0;
            for k in 0..N {
                let theta = j_f64 * theta_k[k];
                s += samples[k] * theta.cos();
            }
            coeffs[j] = 2.0 * s * inv_n;
        }
        
        if N > 0 {
            // Nur das erste Element wird halbiert
            coeffs[0] *= 0.5;
        }

        Self { a, b, coeffs }
    }

    pub fn eval(&self, x: f64) -> f64 {
        let y = (2.0 * x - self.a - self.b) / (self.b - self.a);
        let y2 = 2.0 * y;

        let mut d = 0.0;
        let mut dd = 0.0;

        for j in (1..N).rev() {
            let tmp = d;
            d = y2 * d - dd + self.coeffs[j];
            dd = tmp;
        }

        y * d - dd + self.coeffs[0] // Korrigiert: Nur das erste Element addieren
    }

    #[allow(clippy::needless_range_loop)]
    pub fn derive(&self) -> Self {
        let mut d_coeffs = [0.0; N];
        
        if N < 2 {
            return Self { a: self.a, b: self.b, coeffs: d_coeffs };
        }

        let scale = 2.0 / (self.b - self.a);

        // Skalierung direkt einrechnen spart eine komplette Schleife!
        for j in (1..N - 1).rev() {
            let next_plus_two = if j + 2 < N { d_coeffs[j + 2] } else { 0.0 };
            d_coeffs[j] = next_plus_two + scale * 2.0 * ((j + 1) as f64) * self.coeffs[j + 1];
        }
        
        let nxt_0 = if 2 < N { d_coeffs[2] } else { 0.0 };
        d_coeffs[0] = 0.5 * nxt_0 + scale * self.coeffs[1];

        Self { a: self.a, b: self.b, coeffs: d_coeffs }
    }

    #[allow(clippy::needless_range_loop)]
    pub fn integrate(&self) -> Self {
        let mut i_coeffs = [0.0; N];

        if N < 2 {
            return Self { a: self.a, b: self.b, coeffs: i_coeffs };
        }

        let scale = 0.5 * (self.b - self.a);

        // Skalierung direkt einrechnen spart eine komplette Schleife!
        let c2_val = if 2 < N { self.coeffs[2] } else { 0.0 };
        i_coeffs[1] = scale * (self.coeffs[0] * 2.0 - 0.5 * c2_val);

        for j in 2..N - 1 {
            i_coeffs[j] = scale * (self.coeffs[j - 1] - self.coeffs[j + 1]) / (2.0 * j as f64);
        }

        if N > 1 {
            i_coeffs[N - 1] = scale * self.coeffs[N - 2] / (2.0 * (N - 1) as f64);
        }

        // Integrationskonstante bestimmen
        let mut d = 0.0;
        let mut dd = 0.0;
        for j in (1..N).rev() {
            let tmp = d;
            d = -2.0 * d - dd + i_coeffs[j];
            dd = tmp;
        }
        i_coeffs[0] = d + dd;

        Self { a: self.a, b: self.b, coeffs: i_coeffs }
    }
}
