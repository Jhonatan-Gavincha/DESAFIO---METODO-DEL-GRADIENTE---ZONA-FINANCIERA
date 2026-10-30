# Planificación financiera con Gradientes Conjugados Precondicionados (GCP)

**Universidad Mayor de San Andrés · Facultad de Ciencias Puras y Naturales · Carrera de Informática**
**Desafío académico de Métodos Numéricos · Mención: Ingeniería de Sistemas · Modalidad individual**
**Estudiante:** Jhonatan Pongo Gavincha · **Tolerancia:** 1e-6

Aplicación web interactiva que modela la distribución de presupuesto entre seis áreas de una empresa como un sistema lineal `Ax = b` con `A` simétrica definida positiva (SPD). Lo resuelve con el método de Gradientes Conjugados Precondicionado (Suñagua, 2020), lo contrasta con Gradientes Conjugados simple y con un método directo de referencia (Gauss con pivoteo parcial), y lo analiza en tres escenarios: Ideal, Bajo estrés y Mal condicionado.

## 1. Estructura

```
├── index.html   # 10 secciones del desafío
├── style.css    # Estilos y diseño responsivo
├── script.js    # Lógica numérica, validaciones, interfaz y pruebas QA
└── README.md
```

Sin dependencias ni compilación: HTML, CSS y JavaScript puros.

## 2. Ejecución local

Abrir `index.html` en un navegador, o servir la carpeta:

```
python3 -m http.server 8000   # luego abrir http://localhost:8000
```

## 3. Publicación (GitHub Pages)

Subir los archivos a un repositorio público, ir a **Settings → Pages**, elegir la rama `main` y la carpeta `/ (root)`. El sitio queda en `https://<usuario>.github.io/<repositorio>/`.

## 4. Secciones de la página

| N.º | Sección | Contenido |
|---|---|---|
| 1 | Caso de estudio | Contexto, justificación financiera del modelo y procedencia de los datos |
| 2 | Formulación | Variables (MM Bs.), Ax = b, interpretación de coeficientes, prueba teórica y verificación numérica de SPD (autovalores y Sylvester) en los 3 escenarios |
| 3 | Modelado organizacional (Ingeniería de Sistemas) | Procesos y flujo interdepartamental, flujo de caja, valor organizacional f(x*) y reglas de negocio validadas |
| 4 | Entrada de datos | Edición de A, b y x0, tolerancia, Nmáx, precondicionador, matriz M, opción de simetría automática |
| 5 | Tres escenarios | Carga rápida y análisis del estrés frente al ideal |
| 6 | Solución de referencia | Método directo (Gauss con pivoteo parcial) |
| 7 | GCP vs GC | Vectores, iteraciones, residuos, κ, gráfico de convergencia y tabla comparativa de los tres escenarios |
| 8 | Análisis de errores | Ea, Er, Ep por variable y residuo relativo global |
| 9 | QA | Casos de prueba con Entrada, Esperado, Obtenido y Conclusión, ejecutados al cargar y con un botón |
| 10 | Conclusiones y referencias | Interpretación técnica y financiera, Suñagua (2020) y fuente de datos |

## 5. Modelo financiero

Se minimiza `f(x) = 1/2 x^T A x - b^T x` (x en millones de Bs.). Su gradiente `Ax - b` se anula en el equilibrio, es decir en la solución de `Ax = b`.

| Variable | Área |
|---|---|
| x1 | Tecnología |
| x2 | Operaciones |
| x3 | Marketing |
| x4 | Logística |
| x5 | Personal |
| x6 | Innovación |

- `a_ii > 0`: costo propio del área.
- `a_ij = a_ji < 0`: complementariedad entre áreas consecutivas de la cadena de valor (Tecnología–Operaciones, Operaciones–Marketing, Marketing–Logística, Logística–Personal, Personal–Innovación).
- `b = A · x_meta`: metas ponderadas que reproducen la asignación objetivo de la gerencia.
- Los datos son sintéticos, con fines académicos.

### Escenarios

Los tres usan `A = D·C·D`, con `C` tridiagonal y `D` diagonal positiva.

| Escenario | C (diagonal / vecinas) | D | Observación |
|---|---|---|---|
| Ideal | 4 / −0,7 | cercana a 1 | κ(A) ≈ 2,3 |
| Bajo estrés | 4,4 / −1,25 | cercana a 1 | Inflación e interdependencia; metas b recortadas (Logística sube). κ(A) ≈ 3,3: solo se desplaza el equilibrio |
| Mal condicionado | 4 / −1 | diag(1, 3, 10, 30, 100, 300) | Mal escalamiento diagonal; κ(A) ≈ 9,8·10⁴ y Jacobi lo reduce a ≈ 2,6 |

**Prueba de SPD.** `C` es simétrica y estrictamente diagonal dominante (4 > 2·0,7; 4,4 > 2·1,25; 4 > 2·1), por tanto definida positiva. Como `D` es invertible, `vᵀAv = (Dv)ᵀC(Dv) > 0` para todo `v ≠ 0`. El programa lo confirma con autovalores, menores principales (Sylvester) y Cholesky.

## 6. Algoritmos (`script.js`)

| Función | Descripción |
|---|---|
| `solveDirect` | Eliminación de Gauss con pivoteo parcial (solución de referencia) |
| `cholesky`, `minorsOf` | Verificación de definición positiva y menores principales |
| `eigenSym`, `condSPD` | Autovalores (rotaciones de Jacobi) y κ₂ = λmax/λmin |
| `transformedJacobi` | Matriz D^{-1/2} A D^{-1/2}, semejante a M⁻¹A, para κ(M⁻¹A) |
| `cg` | GCP; con `precond = false` equivale al GC simple |
| `validateMatrix` | Simetría, valores finitos y definición positiva, con mensajes al usuario |
| `renderExtras`, `renderConclusions` | Verificación SPD, flujo de caja, valor organizacional, reglas de negocio y análisis del estrés |
| `runQA` | Pruebas de software |

Algoritmo GCP (Suñagua, 2020), con `M = diag(A)`:

```
r0 = b - A x0;  resolver M z0 = r0;  p0 = z0
para k = 0, 1, ...
    alpha_k = (r_k·z_k) / (p_k·A p_k)
    x_{k+1} = x_k + alpha_k p_k
    r_{k+1} = r_k - alpha_k A p_k
    si ||r_{k+1}||2 / ||b||2 < tol -> parar
    resolver M z_{k+1} = r_{k+1}
    beta_k = (r_{k+1}·z_{k+1}) / (r_k·z_k)
    p_{k+1} = z_{k+1} + beta_k p_k
```

## 7. Resultados de referencia (tol = 1e-6, x0 = 0)

| Escenario | κ(A) | κ(M⁻¹A) | Iter. GC | Iter. GCP | Error abs. máx. |
|---|---|---|---|---|---|
| Ideal | 2,33 | 1,92 | 6 | 6 | ~4e-16 |
| Bajo estrés | 3,35 | 3,10 | 6 | 6 | ~7e-16 |
| Mal condicionado | 98 087 | 2,64 | 7 | 6 | ~4e-16 |

Con n = 6, GC termina en a lo sumo 6 pasos en aritmética exacta; por eso la ventaja de Jacobi se aprecia sobre todo en el número de condición y no tanto en las iteraciones.

## 8. Pruebas de software (QA)

1. **Carga de escenarios:** Ideal, Bajo estrés y Mal condicionado válidos y convergentes.
2. **Matriz no válida:** no simétrica, simétrica no definida positiva y coeficiente no numérico; se notifica sin colapsar. En la interfaz se puede probar desmarcando "Mantener simetría automática".
3. **Edición dinámica:** modificar `a11` y `b1` recalcula referencia, GCP y errores; además, convergencia desde un x0 distinto de cero.
4. **Precondicionamiento:** en el escenario mal condicionado, iteraciones GCP ≤ iteraciones GC.

## 9. Limitaciones

- n = 6 fijo (cumple n ≥ 4).
- Precondicionadores disponibles: Jacobi y ninguno.
- Datos sintéticos, sin recomendación de inversión.

## 10. Referencias

- Suñagua S., P. (2020). *Método de Gradientes Conjugados Precondicionado.* Revista Boliviana de Matemática – UMSA, 04, 2–7.
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.). SIAM.
