Da du die Dokumentation als PDF-Kurzerklärung ausgelegt haben möchtest, ist Markdown die perfekte Quelle. Du kannst den folgenden Inhalt ganz einfach kopieren und in einen beliebigen Markdown-zu-PDF-Konverter (wie Pandoc, VS Code Markdown PDF oder Online-Tools) einfügen.
Der Text ist mathematisch exakt, modular aufgebaut und bietet Neueinsteigern den perfekten theoretischen Unterbau zu deinem 19-Nanosekunden-Code.

------------------------------
## PDF-Dokumentation: Chebyshev-Approximation & Clenshaw-Toolchain
Ein Leitfaden für Numerik-Einsteiger, statistische Modellierung und Hochleistungs-Computing
------------------------------
## 1. Wozu Chebyshev-Approximation? (Motivation)
Die mathematische Approximation von Funktionen ist ein Grundpfeiler der numerischen Mathematik. Häufig müssen analytisch komplexe oder rechenintensive Funktionen ($e^x$, $\sin(x)$, Wahrscheinlichkeitsdichten für statistische CDFs) millionenfach in Simulationen oder Echtzeitsystemen ausgewertet werden.
Hier kommen Chebyshev-Polynome ins Spiel: Sie erlauben es, eine kontinuierliche Funktion $f(x)$ auf einem Intervall $[a,b]$ durch ein Polynom $N$-ten Grades zu ersetzen, das extrem schnell berechnet werden kann, während der Approximationsfehler minimal bleibt.
## Wann ist Chebyshev von Vorteil?

* Minimierung des maximalen Fehlers (Minimax-Prinzip): Im Gegensatz zu Taylor-Polynomen, die nur um einen einzigen Punkt herum extrem exakt sind und zu den Rändern hin explodieren, verteilen Chebyshev-Polynome den Fehler gleichmäßig über das gesamte Intervall. Sie kommen dem mathematisch optimalen „Minimax-Polynom“ verblüffend nahe.
* Vermeidung des Runge-Phänomens: Bei einer klassischen Polynom-Interpolation mit gleichabständigen Punkten neigen Polynome höheren Grades an den Intervallrändern zu extremen Oszillationen (Runge-Phänomen). Chebyshev-Knoten verhindern dies vollständig.
* Analytische Toolchain: Aus den berechneten Koeffizienten lassen sich die exakten Koeffizienten für die mathematische Ableitung oder das Integral berechnen, ohne dass man numerisch differenzieren oder integrieren muss (was oft instabil ist).

## Wann ist es eher NICHT von Vorteil?

* Unstetigkeiten und Knicke: Besitzt eine Funktion Sprungstellen, Singularitäten oder nicht-differenzierbare Knicke (z. B. eine ReLU-Funktion oder $\vert{}x\vert{}$), konvergiert die Chebyshev-Approximation nur sehr langsam (Gibbssches Phänomen).
* Unendliche Intervalle: Die Standard-Methodik ist strikt auf ein kompaktes, endliches Intervall $[a,b]$ begrenzt. Für unendliche Intervalle (wie $[0, \infty)$) muss man auf abgewandelte, rationale Chebyshev-Funktionen ausweichen.

------------------------------
## 2. Die mathematischen Formeln unserer Implementierung
Die Transformation basiert darauf, dass das reale Intervall $x \in [a,b]$ auf das kanonische Chebyshev-Intervall $y \in [-1, 1]$ abgebildet wird.
## 2.1 Intervall-Transformation
Die Abbildung eines Punktes $x$ auf den normierten Bereich $y$ (und umgekehrt) erfolgt über:
$$y = \frac{2x - a - b}{b - a} \quad \Longleftrightarrow \quad x = \frac{a+b}{2} + \frac{b-a}{2} \cdot y$$ 
## 2.2 Bestimmung der Koeffizienten (Knotenpunkte)
Um die Koeffizienten $c_j$ für das Polynom zu finden, wird die Funktion an den Nullstellen (Wurzeln) der Chebyshev-Polynome gesampelt. Diese Chebyshev-Knoten verdichten sich mathematisch zu den Intervallrändern hin:
$$y_k = \cos\left( \frac{\pi (k + 0.5)}{N} \right) \quad \text{für } k = 0, 1, \dots, N-1$$ 
Die zugehörigen Funktionswerte sind $f_k = f(x(y_k))$. Die Koeffizienten $c_j$ der Reihe berechnen sich anschließend über die Diskrete Kosinustransformation (DCT):
$$c_j = \frac{2}{N} \sum_{k=0}^{N-1} f_k \cdot \cos\left( \frac{\pi j (k + 0.5)}{N} \right) \quad \text{für } j = 0, 1, \dots, N-1$$ 
Hinweis zur Normierung: Der Koeffizient $c_0$ wird in der Summe der Polynomreihe traditionell mit dem Faktor $0.5$ gewichtet (daher die Zeile coeffs[0] *= 0.5; im Code).
## 2.3 Die Clenshaw-Rückwärts-Auswertung (eval)
Um das Polynom stabil und ohne das ungenaue Aufsummieren von $x^j$-Potenzen auszuwerten, nutzen wir die Clenshaw-Rekursion. Wir starten bei $j = N-1$ und rechnen rückwärts bis $j = 1$:
$$d_j = 2y \cdot d_{j+1} - d_{j+2} + c_j \quad \text{mit } d_N = 0, \; d_{N+1} = 0$$ 
Der finale Funktionswert am Punkt $x$ (bzw. transformiert $y$) ergibt sich aus dem letzten Schritt:
$$f(x) \approx y \cdot d_1 - d_2 + c_0$$ 
## 2.4 Analytische Ableitung (derive)
Die Koeffizienten der Ableitung $c'_j$ lassen sich direkt aus den Koeffizienten $c_j$ der Originalfunktion über eine Rückwärts-Verknüpfung bestimmen:
$$c'_{j} = c'_{j+2} + 2 \cdot (j+1) \cdot c_{j+1} \quad \text{für } j = N-2, N-3, \dots, 1$$ 
Für den nullten Koeffizienten gilt: $c'_0 = 0.5 \cdot c'_2 + c_1$. Da die Ableitung auf dem Intervall $[-1,1]$ berechnet wurde, muss sie mit der inneren Ableitung der Intervall-Transformation skaliert werden:
$$\text{scale} = \frac{2}{b-a} \quad \Longrightarrow \quad c'_j \leftarrow c'_j \cdot \text{scale}$$ 
## 2.5 Analytisches Integral / Stammfunktion (integrate)
Die Koeffizienten des Integrals $C_j$ werden über eine Vorwärts-Beziehung bestimmt. Sie verhalten sich antiproportional zur Frequenz $j$:
$$C_j = \frac{c_{j-1} - c_{j+1}}{2j} \quad \text{für } j = 1, 2, \dots, N-2$$ 
Spezialfälle für die Ränder: $C_{N-1} = \frac{c_{N-2}}{2(N-1)}$. Der Koeffizient $C_0$ wird so gewählt, dass die Bedingung $F(a) = 0$ (Nullpunkt am linken Intervallrand) erfüllt ist. Hierzu wird das Polynom mittels Clenshaw an der Stelle $x = a$ (entspricht $y = -1$) ausgewertet. Die Abweichung wird von $C_0$ subtrahiert.
Zusätzlich gilt auch hier die Skalierung der äußeren Transformation:
$$\text{scale} = \frac{b-a}{2} \quad \Longrightarrow \quad C_j \leftarrow C_j \cdot \text{scale}$$ 
------------------------------
## 3. Literaturempfehlungen & Aktuelle Textbooks
Für Neueinsteiger, die tiefer in die faszinierende Welt der orthogonalen Polynome und der spektralen Numerik eintauchen möchten, sind folgende Werke der Goldstandard:
## Englische Fachliteratur (Internationaler Standard)

   1. "Trefethen, L. N. (2013). Approximation Practice and Chebfun Guide." SIAM.
   * Das absolute Kultbuch zum Thema. Lloyd N. Trefethen ist der Urvater der modernen Chebyshev-Numerik. Das Buch ist extrem anschaulich geschrieben, verzichtet auf unnötig trockenes Lemma-Wälzen und erklärt perfekt, warum Chebyshev-Methoden die numerische Welt dominieren.
   2. "Boyd, J. P. (2001). Chebyshev and Fourier Spectral Methods." Dover Publications.
   * Der unangefochtene Klassiker. Perfekt für Entwickler und Ingenieure. Es erklärt detailliert, wie man die Algorithmen in Computercode gießt und analysiert. (Bonus: Der Autor hat das Buch als kostenloses PDF über seine Universitätsseite freigegeben).
   3. "Press, W. H., et al. (2007). Numerical Recipes: The Art of Scientific Computing (3rd Edition)." Cambridge University Press.
   * Das Handbuch für Programmierer. Kapitel 5.8 widmet sich ausschließlich der Chebyshev-Approximation und Clenshaw-Auswertung. Ideal, um die Parallelen zwischen unserem Rust-Code und klassischen C++/Fortran-Bibliotheken zu verstehen.
   
## Deutsche Fachliteratur

   1. "Schwarz, H. R., Köckler, J. (2011). Numerische Mathematik." Vieweg+Teubner.
   * Ein hervorragendes Standardlehrwerk an deutschen Universitäten. Es erklärt die Intervall-Transformationen und die diskrete Kosinustransformation (DCT) mathematisch präzise auf Deutsch.
   2. "Deuflhard, P., Hohmann, A. (2019). Numerische Mathematik 1: Eine algorithmisch orientierte Einführung." De Gruyter.
   * Sehr strukturierter Aufbau, der den Fokus stark darauf legt, wie mathematische Sätze effizient in Algorithmen übersetzt werden. Hervorragend geeignet, um das Prinzip der Fehlerminimierung zu verinnerlichen.
   
