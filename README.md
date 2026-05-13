# 🛡️ MalwareODE Simulator

> **Simulador interactivo de propagación y supresión de malware en redes LAN mediante EDO de 2do Orden Amortiguada**

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![PyQt5](https://img.shields.io/badge/PyQt5-5.15%2B-green?logo=qt)](https://pypi.org/project/PyQt5/)
[![SymPy](https://img.shields.io/badge/SymPy-1.12%2B-orange)](https://www.sympy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.7%2B-red)](https://matplotlib.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 📋 Descripción general

**MalwareODE Simulator** es una aplicación de escritorio desarrollada en Python que modela matemáticamente la dinámica de propagación y neutralización de malware en una red de área local (LAN), utilizando como base teórica la **ecuación diferencial ordinaria (EDO) de segundo orden amortiguada forzada**:

```
m·y''(t) + b·y'(t) + k·y(t) = f(t)
```

El proyecto surge de la analogía física entre el **oscilador masa-resorte-amortiguador** y el comportamiento de una red bajo ataque informático. Cada parámetro de la EDO representa un elemento concreto del ecosistema de ciberseguridad:

| Parámetro | Significado matemático | Significado en ciberseguridad |
|-----------|------------------------|-------------------------------|
| `y(t)` | Función incógnita | Proporción de nodos comprometidos en el instante `t` (0 = red sana, 1 = 100% infectada) |
| `y'(t)` | Primera derivada | Velocidad de propagación del malware |
| `y''(t)` | Segunda derivada | Aceleración del contagio |
| `m` | Inercia / masa | Complejidad y latencia de la red (topología, cantidad de nodos) |
| `b` | Coeficiente de amortiguamiento | Efectividad del firewall y respuesta del equipo SOC |
| `k` | Rigidez del sistema | Capacidad de auto-recuperación (antivirus, sandboxing, snapshots) |
| `f(t)` | Fuerza externa forzante | Intensidad del ataque externo a lo largo del tiempo |

El simulador resuelve esta EDO tanto **analíticamente** (con SymPy, paso a paso, mostrando cada derivación algebraica) como **numéricamente** (con integración de Euler), y presenta el resultado en tiempo real sobre un gráfico interactivo.

---

## ✨ Características principales

### Simulación interactiva en tiempo real
- Tres sliders independientes para ajustar `m`, `b` y `k` con realimentación visual inmediata
- Identificación automática del régimen de amortiguamiento según el discriminante `Δ = b² − 4mk`
- Estadísticas dinámicas: pico máximo de infección, instante del pico y tiempo de estabilización

### Seis tipos de fuerza externa f(t)
- **Exponencial** `A·e^(-αt)`: modela un ataque DDoS que se mitiga con el tiempo
- **Coseno** `A·cos(ωt)`: campañas de phishing periódicas y cíclicas
- **Seno** `A·sin(ωt)`: oleadas graduales de intentos de acceso no autorizado
- **Polinómica** `A·t^n·e^(-0.25t)`: malware tipo gusano que escala progresivamente
- **Impulso de Dirac** `P·δ(t−t₀)`: exploit de día cero ejecutado en un instante exacto, aproximado numéricamente como pulso rectangular de altura `P/ε` y duración `ε`
- **Personalizada**: el usuario ingresa cualquier expresión arbitraria en `t`

### Solver simbólico automático (Coeficientes Indeterminados)
Para funciones personalizadas, el simulador analiza automáticamente la forma de `f(t)` y aplica el Método de Coeficientes Indeterminados (Zill / Ibarra Escutia):
- Detecta la familia funcional (exponencial, trigonométrica, polinómica, exp·trig)
- Detecta resonancia y aplica la forma modificada `t^n · yp`
- Deriva `yp`, `yp'`, `yp''` simbólicamente
- Sustituye en la EDO y resuelve el sistema algebraico para los coeficientes `A`, `B`, `C`...
- Convierte la expresión exacta en función numérica evaluable

Para funciones que no admiten Coeficientes Indeterminados (Dirac, formas mixtas), aplica automáticamente **Variación de Parámetros** (Lagrange) con integración numérica acumulativa.

### Desarrollo analítico paso a paso (actualizado en tiempo real)
El panel derecho muestra la derivación matemática completa:
1. **Paso 1** — Planteamiento y clasificación de la EDO
2. **Paso 2** — Ecuación característica, discriminante y solución complementaria `yc`
3. **Paso 3** — Propuesta de `yp`, derivadas, sustitución y resolución de coeficientes
4. **Paso 4** — Solución general `y(t) = yc(t) + yp(t)`
5. **Paso 5** — Aplicación de condiciones iniciales, valores de `C₁` y `C₂`, solución única

### Condiciones iniciales contextualizadas
- `y(0)`: porcentaje de nodos ya comprometidos al inicio del análisis
- `y'(0)`: velocidad de propagación inicial (positivo = expansión activa, negativo = ya contenido)

---

## 🧮 Fundamento matemático

### Clasificación de la EDO

La ecuación `m·y''(t) + b·y'(t) + k·y(t) = f(t)` pertenece a la categoría:

- **Ordinaria**: variable independiente escalar `t`
- **Lineal**: los coeficientes `m`, `b`, `k` no dependen de `y`
- **Orden 2**: involucra hasta la segunda derivada
- **Coeficientes constantes**: `m`, `b`, `k` ∈ ℝ, no dependen de `t`
- **No homogénea**: `f(t) ≠ 0`

### Los tres regímenes de amortiguamiento

La solución depende del discriminante `Δ = b² − 4mk`:

```
Δ < 0  →  SUBAMORTIGUADO   →  Raíces complejas r = α ± iω
           La red oscila antes de estabilizarse.
           El malware puede reinfectar nodos ya limpios.
           yc = e^(αt)·[C₁·sin(ωt) + C₂·cos(ωt)]

Δ > 0  →  SOBREAMORTIGUADO  →  Raíces reales distintas r₁, r₂
           Recuperación monótona pero lenta.
           yc = C₁·e^(r₁t) + C₂·e^(r₂t)

Δ = 0  →  CRÍTICO           →  Raíz real doble r
           Neutralización óptima: máxima velocidad sin oscilaciones.
           Fórmula del firewall ideal: b = 2·√(m·k)
           yc = (C₁ + C₂·t)·e^(rt)
```

### Método de Coeficientes Indeterminados

Basado en los textos de **Zill (2017)** e **Ibarra Escutia (2010)**, el método es válido cuando `f(t)` pertenece a una familia funcional cerrada bajo diferenciación:

| Familia de f(t) | Propuesta para yp |
|-----------------|-------------------|
| `A·e^(αt)` | `B·e^(αt)` |
| `A·sin(ωt)` o `A·cos(ωt)` | `B·sin(ωt) + C·cos(ωt)` |
| `Aₙtⁿ + ... + A₀` | `Bₙtⁿ + ... + B₀` |
| `e^(αt)·[sin/cos]` | `e^(αt)·[B·sin + C·cos]` |
| Resonancia (α = raíz) | Propuesta × `t^(multiplicidad)` |

### Variación de Parámetros (Lagrange)

Para funciones arbitrarias se aplica la fórmula de Lagrange:

```
yp(t) = y₁(t)·u₁(t) + y₂(t)·u₂(t)

u₁(t) = ∫ [−y₂(τ)·f(τ)] / [m·W(τ)] dτ
u₂(t) = ∫ [y₁(τ)·f(τ)]  / [m·W(τ)] dτ

W(t) = y₁·y₂' − y₁'·y₂  (Wronskiano)
```

donde `y₁`, `y₂` son las soluciones independientes de la ecuación homogénea.

---

## 🖥️ Requisitos del sistema

| Requisito | Versión mínima |
|-----------|---------------|
| Python | 3.9 |
| PyQt5 | 5.15 |
| matplotlib | 3.7 |
| numpy | 1.24 |
| sympy | 1.12 |
| scipy | 1.10 |

> **Sistemas operativos soportados:** Windows 10/11, macOS 12+, Ubuntu 20.04+

---

## 🚀 Instalación y uso

### 1. Clonar el repositorio

```bash
git clone https://github.com/tu-usuario/malware-ode-simulator.git
cd malware-ode-simulator
```

### 2. Crear entorno virtual (recomendado)

```bash
python -m venv venv

# Windows
venv\Scripts\activate

# macOS / Linux
source venv/bin/activate
```

### 3. Instalar dependencias

```bash
pip install -r requirements.txt
```

### 4. Ejecutar el simulador

```bash
python src/malware_ode_simulator.py
```

---

## 📦 Empaquetar como .exe (Windows)

El simulador puede distribuirse como ejecutable independiente usando **PyInstaller**:

```bash
pip install pyinstaller

pyinstaller --onefile --windowed \
            --name "MalwareODE_Simulator" \
            src/malware_ode_simulator.py
```

El ejecutable se genera en `dist/MalwareODE_Simulator.exe`.

> **Nota:** el archivo `.exe` tendrá un tamaño de aproximadamente 40–60 MB ya que incluye Python, PyQt5 y todas las librerías científicas empaquetadas.

---

## 🎮 Guía de uso rápido

### Panel izquierdo — Controles

| Control | Descripción |
|---------|-------------|
| Slider **m** | Inercia de la red. Arrastrar hacia la derecha = red más compleja |
| Slider **b** | Intensidad del firewall. Valor óptimo: `b = 2·√(m·k)` |
| Slider **k** | Recuperación automática. Alto = antivirus más agresivo |
| **y(0)** | Porcentaje de nodos ya infectados al inicio |
| **y'(0)** | Velocidad de propagación inicial. Negativo = ya se está conteniendo |
| **Chips f(t)** | Tipo de ataque externo |
| **Campo texto** | Expresión personalizada de `f(t)` en Python/SymPy |
| **Presets** | Carga automáticamente los 3 escenarios clásicos |

### Presets incluidos

| Preset | Parámetros | Situación |
|--------|-----------|-----------|
| Firewall débil | m=1, b=1, k=2 | Red sin protección adecuada, infección oscila |
| Respuesta lenta | m=1, b=3, k=2 | Firewall demasiado conservador, recuperación tardía |
| Configuración óptima | m=1, b=2, k=1 | Amortiguamiento crítico: neutralización más rápida posible |

### Expresiones válidas para f(t) personalizada

```python
6*exp(-2*t)                    # exponencial simple
3*sin(t) + 2*cos(2*t)          # trigonométrica mixta
t**2 * exp(-0.5*t)             # polinómica con decaimiento
exp(-t) * cos(3*t)             # exponencial × trigonométrica
0.5 * sin(t) * exp(-0.2*t)     # con amortiguamiento exponencial
t**3 - 2*t + 1                 # polinomio puro
```

---

## 📁 Estructura del repositorio

```
malware-ode-simulator/
│
├── 📄 README.md                        ← Este archivo
├── 📄 requirements.txt                 ← Dependencias Python
├── 📄 .gitignore                       ← Archivos ignorados por Git
│
├── 📂 src/                             ← Código fuente principal
│   └── 📄 ciberseguridad_EDO.py        ← Aplicación completa (901 líneas)
│
├── 📂 docs/                            ← Documentación
│   ├── 📄 fundamento_matematico.md     ← Teoría de EDOs y analogía
│   └── 📂 img/                         ← Capturas de pantalla
│       ├── 🖼️ screenshot_main.png
│       ├── 🖼️ screenshot_subamortiguado.png
│       ├── 🖼️ screenshot_sobreamortiguado.png
│       ├── 🖼️ screenshot_critico.png
│       └── 🖼️ screenshot_impulso.png
│
│
├── 📂 script/                         ← Utilidades de empaquetado
    ├── 📄 build_exe_windows.bat        ← Script para generar .exe en Windows

```

---

## 🏗️ Arquitectura del código

El código de `malware_ode_simulator.py` está organizado en cinco capas funcionales:

```
┌─────────────────────────────────────────────────────┐
│                   CAPA 5: ENTRY POINT               │
│   main() → QApplication + QPalette + MainWindow     │
├─────────────────────────────────────────────────────┤
│               CAPA 4: VENTANA PRINCIPAL             │
│   MainWindow → _build_ui(), _left_panel(),          │
│                _right_panel(), _do_update()         │
├─────────────────────────────────────────────────────┤
│                CAPA 3: HELPERS DE UI                │
│   make_label(), make_slider(), make_groupbox(),     │
│   make_step_frame(), make_math_box(), make_chip()   │
├─────────────────────────────────────────────────────┤
│              CAPA 2: CANVAS MATPLOTLIB              │
│   ODECanvas (FigureCanvas) → update_plot()          │
├─────────────────────────────────────────────────────┤
│              CAPA 1: MATEMÁTICA PURA                │
│   _parse_expr()           → string → SymPy          │
│   _extract_alpha/omega()  → detectar parámetros     │
│   solve_yp_symbolic()     → Coef. Indeterminados    │
│   variation_of_params()   → Método de Lagrange      │
│   eval_f_numeric()        → evaluar f(t)            │
│   solve_ode_numeric()     → Integración Euler       │
│   build_all_steps()       → texto de los 5 pasos    │
└─────────────────────────────────────────────────────┘
```

### Flujo de datos principal

```
Usuario mueve slider
        │
        ▼
   _sched()  ──[debounce 80ms]──▶  _do_update()
                                        │
                          ┌─────────────┼─────────────┐
                          ▼             ▼             ▼
                   eval_f_numeric()  _f_expr_str()    │
                   (array f(t))      (string SymPy)   │
                          │             │             │
                          │             ▼             │
                          │    solve_yp_symbolic()    │
                          │    ┌──────────────────┐   │
                          │    │ Detectar familia │   │
                          │    │ Proponer yp      │   │
                          │    │ Derivar yp, yp'' │   │
                          │    │ Sustituir en EDO │   │
                          │    │ Resolver A,B,C   │   │
                          │    └──────────────────┘   │
                          │             │             │
                          └──────┬──────┘             │
                                 ▼                    │
                        solve_ode_numeric()           │
                        (y(t) con CI)                 │
                                 │                    │
                    ┌────────────┼───────────┐        │
                    ▼            ▼           ▼        │
              ODECanvas    _update_stats  build_all_steps()
              .update_plot()             → math_boxes[0..4]
```

---

## 🔬 Descripción técnica de cada función

### `_parse_expr(expr_str)` — Línea 63
Convierte el string ingresado por el usuario (`"6*exp(-2*t)"`) en un objeto matemático de SymPy. Normaliza prefijos de numpy (`np.exp` → `exp`) y mapea el texto `t` a la variable simbólica global `_t_sym`.

### `_extract_alpha(expr)` — Línea 78
Recorre el árbol de expresión SymPy buscando nodos `exp(...)`. Extrae el coeficiente de `t` en el argumento del exponencial. Para `e^(-2t)` retorna `α = -2`.

### `_extract_omega(expr)` — Línea 85
Similar al anterior pero para funciones trigonométricas. Para `sin(3t)` retorna `ω = 3`.

### `solve_yp_symbolic(expr_str, m, b, k)` — Línea 93
Motor principal del solver analítico. Implementa el Método de Coeficientes Indeterminados:
1. Clasifica `f(t)` en su familia funcional
2. Verifica resonancia comparando α contra las raíces de `m·r² + b·r + k = 0`
3. Construye la propuesta `yp` con coeficientes simbólicos A, B, C, D
4. Calcula `yp'` y `yp''` con `sp.diff()`
5. Forma la ecuación `m·yp'' + b·yp' + k·yp = f(t)`
6. Resuelve con `sp.solve()` para A, B, C, D
7. Retorna la expresión exacta + texto del desarrollo + función numérica via `lambdify`

### `variation_of_params(m, b, k, f_arr, t_arr)` — Línea 269
Método de Lagrange para `f(t)` arbitraria. Construye `y₁`, `y₂` (soluciones independientes de la homogénea) según el régimen, calcula el Wronskiano `W = y₁·y₂' − y₁'·y₂` y obtiene `yp` integrando numéricamente con `np.cumsum` (suma de Riemann con paso `dt`).

### `eval_f_numeric(t_val, ftype, ...)` — Línea 307
Evalúa `f(t)` en un único punto escalar. La aproximación del Impulso de Dirac `δ(t−t₀)` se implementa como pulso rectangular de altura `A/ε` y ancho `ε`: cuando `|t − t₀| < ε/2` retorna `A/ε`, de lo contrario `0`.

### `solve_ode_numeric(m, b, k, y0, dy0, ...)` — integración completa
Reduce la EDO de 2do orden al sistema:
- `u₁' = u₂` (donde `u₁ = y`, `u₂ = y'`)
- `u₂' = [f(t) − b·u₂ − k·u₁] / m`

Integra con **método de Euler hacia adelante** con `dt = 25/800 ≈ 0.031` unidades de tiempo.

### `build_all_steps(...)` — Línea 347
Toma todos los resultados calculados y construye el texto de cada uno de los 5 pasos analíticos para mostrar en la interfaz. Calcula `C₁` y `C₂` resolviendo el sistema de condiciones iniciales `y(0) = y₀`, `y'(0) = dy₀`.

### `ODECanvas` — Línea 489
Clase que hereda de `FigureCanvasQTAgg`: el gráfico matplotlib integrado como widget Qt. Las tres curvas (`y(t)`, `yp(t)`, `f(t)`) se crean con datos vacíos y se actualizan eficientemente con `set_data()` + `draw()`, evitando recrear el gráfico completo en cada render.

### `MainWindow._sched()` — Línea 783
Implementa **debounce** de 80ms. Cada cambio de slider reinicia el temporizador. Cuando el usuario deja de interactuar, el timer dispara `_do_update()` una sola vez, evitando recálculos innecesarios mientras se arrastra el slider.

### `MainWindow._do_update()` — Línea 785
Orquestador central. En cada actualización: evalúa `f(t)`, calcula `yp` (simbólico o numérico según la familia), integra la EDO completa, actualiza el gráfico, las estadísticas y el desarrollo analítico.

---

## 📊 Interpretación de los resultados

### Badge de régimen
El indicador en la parte superior identifica automáticamente el comportamiento del sistema:
- 🔴 **SUBAMORTIGUADO**: el firewall es insuficiente, habrá brotes recurrentes
- 🟡 **SOBREAMORTIGUADO**: la contención es excesiva o lenta, mayor ventana de vulnerabilidad
- 🟢 **AMORTIGUAMIENTO CRÍTICO**: configuración óptima, neutralización más rápida sin rebotes

### Métricas en tiempo real
- **Δ = b²−4mk**: discriminante numérico. Monitorear su signo para ajustar `b`
- **Pico de infección**: porcentaje máximo de red comprometida y cuándo ocurre
- **Tiempo estable**: instante en que la infección cae por debajo del 1% (red efectivamente sana)

### Curvas del gráfico
- 🔴 **y(t)**: solución completa — comportamiento real de la infección con las condiciones iniciales
- 🟣 **yp(t)**: solución particular — comportamiento de largo plazo forzado por el ataque
- 🔵 **f(t)** (escalada al 10%): intensidad del ataque externo a lo largo del tiempo
- 🟢 **umbral 1%**: línea de referencia que define cuándo la red se considera estable

---

## 🎓 Contexto académico

Este proyecto es una implementación práctica de los modelos matemáticos presentados en:

- **Zill, D.G.** (2017). *Ecuaciones Diferenciales con Problemas de Valores en la Frontera*. Cengage Learning. Capítulos 4 y 5.
- **Ibarra Escutia, M.** (2010). *Ecuaciones Diferenciales Ordinarias*. Método de Coeficientes Indeterminados y Variación de Parámetros.

La analogía entre el oscilador amortiguado y la dinámica de propagación de malware ha sido estudiada en la literatura de modelado matemático de ciberseguridad, donde los modelos SIR (Susceptible-Infected-Recovered) pueden aproximarse por EDOs de segundo orden bajo ciertas condiciones de linealización.

---

## 🤝 Contribuciones

Las contribuciones son bienvenidas. Para contribuir:

1. Hacer fork del repositorio
2. Crear una rama: `git checkout -b feature/nueva-funcion`
3. Hacer commit: `git commit -m "Agrega soporte para f(t) = log(t)"`
4. Push a la rama: `git push origin feature/nueva-funcion`
5. Abrir un Pull Request

## 👤 Autor

Desarrollado como proyecto académico de **Ingeniería en Sistemas de Información** — aplicación de Matemática Aplicada y Ecuaciones Diferenciales al modelado de sistemas de ciberseguridad.

---

<div align="center">
  <sub>Construido con Python · PyQt5 · SymPy · Matplotlib · NumPy · SciPy</sub>
</div>
