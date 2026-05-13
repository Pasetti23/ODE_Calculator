import sys
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QSlider, QGroupBox, QScrollArea, QFrame, QSizePolicy,
    QPushButton, QLineEdit, QTextEdit
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QColor, QPalette


# COLORES

CLR_BG       = "#0f1117"
CLR_PANEL    = "#1a1d27"
CLR_PANEL2   = "#13151f"
CLR_BORDER   = "#2a2d3e"
CLR_TEXT     = "#e2e8f0"
CLR_TEXT_MUT = "#8892a4"
CLR_TEXT_DIM = "#4a5568"
CLR_ACCENT   = "#4f9cf9"
CLR_RED      = "#f05252"
CLR_GREEN    = "#0ea579"
CLR_YELLOW   = "#f5a623"
CLR_PURPLE   = "#9b8afb"
CLR_CYAN     = "#22d3ee"
CLR_ORANGE   = "#fb923c"
PLOT_Y       = "#f05252"
PLOT_F       = "#4f9cf9"
PLOT_YP      = "#a78bfa"
PLOT_STAB    = "#0ea579"
PLOT_BG      = "#0f1117"
PLOT_GRID    = "#1e2535"
PLOT_AXES    = "#2a2d3e"


# CONSTANTES NUMÉRICAS

N_POINTS = 800
T_MAX    = 25.0
T_ARRAY  = np.linspace(0, T_MAX, N_POINTS)


#  CAPA MATEMÁTICA: RESOLVER SIMBÓLICO DE yp


_t_sym = sp.Symbol('t', real=True)

def _parse_expr(expr_str):
    """Convierte el string del usuario en expresión SymPy."""
    s = expr_str.strip()
    for old, new in [
        ('np.exp','exp'), ('np.sin','sin'), ('np.cos','cos'),
        ('np.sqrt','sqrt'), ('np.pi','pi'), ('np.log','log'),
        ('np.abs','Abs'), ('Math.exp','exp'), ('Math.sin','sin'),
        ('Math.cos','cos'),
    ]:
        s = s.replace(old, new)
    return sp.sympify(s, locals={
        't': _t_sym, 'exp': sp.exp, 'sin': sp.sin, 'cos': sp.cos,
        'sqrt': sp.sqrt, 'pi': sp.pi, 'log': sp.log, 'Abs': sp.Abs
    })

def _extract_alpha(expr):
    """Extrae α de la parte exponencial de expr (coeficiente de t en exp(α·t))."""
    for a in sp.preorder_traversal(expr):
        if isinstance(a, sp.exp):
            return a.args[0].coeff(_t_sym)
    return None

def _extract_omega(expr):
    """Extrae ω de sin(ω·t) o cos(ω·t)."""
    for cls in (sp.sin, sp.cos):
        for a in sp.preorder_traversal(expr):
            if isinstance(a, cls):
                return a.args[0].coeff(_t_sym)
    return sp.Integer(1)

def solve_yp_symbolic(expr_str, m_val, b_val, k_val):
    """
    Resuelve yp(t) por el Método de Coeficientes Indeterminados usando SymPy.

    Flujo:
      1. Parsear f(t) desde el string del usuario.
      2. Detectar la familia funcional: exponencial, trigonométrica, polinómica,
         exp·trig, etc.
      3. Proponer la forma de yp con coeficientes indeterminados (A, B, C…).
      4. Calcular yp', yp'' simbólicamente.
      5. Sustituir en m·yp'' + b·yp' + k·yp = f(t).
      6. Resolver el sistema algebraico para A, B, C…
      7. Retornar yp final + texto del desarrollo paso a paso.

    Retorna: (yp_sympy_expr, steps_text, yp_numeric_fn)
      yp_sympy_expr : expresión SymPy de yp (None si falló)
      steps_text    : string con el desarrollo completo para mostrar en la UI
      yp_numeric_fn : función Python yp(t_array) evaluable numéricamente
    """
    A, B, C, D = sp.symbols('A B C D', real=True)
    m, b, k    = sp.Float(m_val), sp.Float(b_val), sp.Float(k_val)

    # ── 1. Parsear ────────────────────────────────────────────────────────────
    try:
        expr = _parse_expr(expr_str)
        expr = sp.expand(expr)
    except Exception as e:
        return None, f"[Error al parsear f(t)]: {e}\nVerificá la sintaxis.", None

    lines = []   # acumulador del desarrollo textual
    lines.append(f"f(t) = {expr}")
    lines.append("")

    has_exp = expr.has(sp.exp)
    has_sin = expr.has(sp.sin)
    has_cos = expr.has(sp.cos)
    try:
        sp.Poly(expr, _t_sym); is_poly = True
        poly_deg = sp.degree(sp.Poly(expr, _t_sym), _t_sym)
    except Exception:
        is_poly = False; poly_deg = 0

    # ── 2. Detectar familia y construir propuesta ─────────────────────────────
    if has_exp and not has_sin and not has_cos:
        # ── Caso A: Exponencial (pura o con polinomio) ────────────────────────
        alpha = _extract_alpha(expr)
        poly_part = sp.expand(expr / sp.exp(alpha * _t_sym))
        try:
            deg = sp.degree(sp.Poly(sp.expand(poly_part), _t_sym), _t_sym)
        except Exception:
            deg = 0

        lines.append(f"Familia: EXPONENCIAL-POLINÓMICA")
        lines.append(f"  f(t) = P(t)·e^({alpha}t),  P(t) = {poly_part},  grado = {deg}")
        lines.append(f"  Según Coef. Indeterminados → yp = Q(t)·e^({alpha}t)")
        lines.append(f"  donde Q(t) tiene el mismo grado que P(t).")

        # Verificar resonancia: ¿α es raíz de mr²+br+k=0?
        r_s      = sp.Symbol('r')
        char_eq  = m * r_s**2 + b * r_s + k
        ch_roots = sp.solve(char_eq, r_s)
        mult     = sum(1 for r in ch_roots if sp.simplify(alpha - r) == 0)

        if mult > 0:
            lines.append(f"\n  ⚠ RESONANCIA: α={alpha} es raíz de la ec. característica")
            lines.append(f"  Multiplicidad = {mult} → propuesta × t^{mult}")
            yp_prop  = (A if deg == 0 else A*_t_sym + B) * _t_sym**mult * sp.exp(alpha*_t_sym)
            unknowns = [A] if deg == 0 else [A, B]
        else:
            if deg == 0:
                yp_prop = A * sp.exp(alpha*_t_sym); unknowns = [A]
            elif deg == 1:
                yp_prop = (A*_t_sym + B) * sp.exp(alpha*_t_sym); unknowns = [A, B]
            else:
                yp_prop = (A*_t_sym**2 + B*_t_sym + C) * sp.exp(alpha*_t_sym); unknowns = [A, B, C]

    elif (has_sin or has_cos) and not has_exp:
        # ── Caso B: Trigonométrica pura ───────────────────────────────────────
        omega = _extract_omega(expr)
        fam   = ("SENO" if has_sin and not has_cos
                 else "COSENO" if has_cos and not has_sin else "TRIG MIXTA")
        lines.append(f"Familia: {fam}  ω = {omega}")
        lines.append(f"  Propuesta siempre incluye sin Y cos:")
        lines.append(f"  yp = A·sin({omega}t) + B·cos({omega}t)")
        yp_prop  = A*sp.sin(omega*_t_sym) + B*sp.cos(omega*_t_sym)
        unknowns = [A, B]

    elif has_exp and (has_sin or has_cos):
        # ── Caso C: Exponencial × Trigonométrica ─────────────────────────────
        alpha = _extract_alpha(expr)
        omega = _extract_omega(expr)
        lines.append(f"Familia: EXPONENCIAL·TRIGONOMÉTRICA")
        lines.append(f"  α = {alpha},  ω = {omega}")
        lines.append(f"  Propuesta: yp = e^({alpha}t)·[A·sin({omega}t) + B·cos({omega}t)]")
        yp_prop  = sp.exp(alpha*_t_sym) * (A*sp.sin(omega*_t_sym) + B*sp.cos(omega*_t_sym))
        unknowns = [A, B]

    elif is_poly:
        # ── Caso D: Polinómica ────────────────────────────────────────────────
        lines.append(f"Familia: POLINÓMICA de grado {poly_deg}")
        if poly_deg == 0:
            yp_prop = A; unknowns = [A]
            lines.append("  Propuesta: yp = A")
        elif poly_deg == 1:
            yp_prop = A*_t_sym + B; unknowns = [A, B]
            lines.append("  Propuesta: yp = At + B")
        elif poly_deg == 2:
            yp_prop = A*_t_sym**2 + B*_t_sym + C; unknowns = [A, B, C]
            lines.append("  Propuesta: yp = At² + Bt + C")
        else:
            yp_prop = A*_t_sym**3 + B*_t_sym**2 + C*_t_sym + D; unknowns = [A, B, C, D]
            lines.append("  Propuesta: yp = At³ + Bt² + Ct + D")
    else:
        # ── Caso E: No reconocida → variación de parámetros ──────────────────
        lines.append("Familia: NO RECONOCIDA por Coef. Indeterminados.")
        lines.append("→ Se aplica Variación de Parámetros (Lagrange) numéricamente.")
        return None, "\n".join(lines), None

    # ── 3. Derivadas de yp ────────────────────────────────────────────────────
    dyp  = sp.diff(yp_prop, _t_sym)
    ddyp = sp.diff(yp_prop, _t_sym, 2)

    lines.append("")
    lines.append("Derivadas de la propuesta:")
    lines.append(f"  yp   = {yp_prop}")
    lines.append(f"  yp'  = {sp.expand(dyp)}")
    lines.append(f"  yp'' = {sp.expand(ddyp)}")

    # ── 4. Sustituir en la EDO ────────────────────────────────────────────────
    lhs = sp.expand(m*ddyp + b*dyp + k*yp_prop)
    lines.append("")
    lines.append(f"Sustituyendo en {m_val}·yp'' + {b_val}·yp' + {k_val}·yp = f(t):")
    lines.append(f"  LI = {lhs}")
    lines.append(f"  LD = {sp.expand(expr)}")

    # ── 5. Resolver el sistema ────────────────────────────────────────────────
    try:
        equation = sp.Eq(sp.expand(lhs), sp.expand(expr))
        sol      = sp.solve(equation, unknowns)

        if isinstance(sol, dict):
            sol_dict = sol
        elif isinstance(sol, list) and len(sol) > 0:
            first = sol[0]
            if isinstance(first, (list, tuple)):
                sol_dict = dict(zip(unknowns, first))
            else:
                sol_dict = {unknowns[0]: first}
        else:
            sol_dict = {}

        if sol_dict:
            lines.append("")
            lines.append("Igualando coeficientes y despejando:")
            for var, val in sol_dict.items():
                lines.append(f"  {var} = {sp.nsimplify(val, rational=True)}")

            yp_final = sp.simplify(yp_prop.subs(sol_dict))
            lines.append("")
            lines.append(f"▶  yp(t) = {yp_final}")

            # Función numérica a partir de la expresión simbólica
            try:
                yp_fn = sp.lambdify(_t_sym, yp_final, modules=['numpy', 'sympy'])
                return yp_final, "\n".join(lines), yp_fn
            except Exception:
                return yp_final, "\n".join(lines), None
        else:
            lines.append("No se encontró solución única → posible resonancia no detectada.")
            return None, "\n".join(lines), None

    except Exception as e:
        lines.append(f"[Error al resolver el sistema]: {e}")
        return None, "\n".join(lines), None


def variation_of_params(m, b, k, f_arr, t_arr):
    """
    Calcula yp(t) numéricamente por Variación de Parámetros (Lagrange).
    Se usa como fallback cuando la familia de f(t) no es reconocida por
    el método de Coeficientes Indeterminados.

    Fórmula:
        yp(t) = y1(t)·u1(t) + y2(t)·u2(t)
        u1(t) = ∫ [−y2(τ)·f(τ) / (m·W(τ))] dτ
        u2(t) = ∫ [y1(τ)·f(τ)  / (m·W(τ))] dτ
        W = y1·y2' − y1'·y2  (Wronskiano)
    """
    disc = b**2 - 4*m*k
    r    = -b / (2*m)
    dt   = t_arr[1] - t_arr[0]

    if disc < -1e-6:
        omega = np.sqrt(-disc) / (2*m)
        y1  = np.exp(r*t_arr) * np.sin(omega*t_arr)
        y2  = np.exp(r*t_arr) * np.cos(omega*t_arr)
        dy1 = r*y1 + omega*np.exp(r*t_arr)*np.cos(omega*t_arr)
        dy2 = r*y2 - omega*np.exp(r*t_arr)*np.sin(omega*t_arr)
    elif disc > 1e-6:
        sq  = np.sqrt(disc)
        r1, r2 = (-b+sq)/(2*m), (-b-sq)/(2*m)
        y1 = np.exp(r1*t_arr); y2 = np.exp(r2*t_arr)
        dy1 = r1*y1; dy2 = r2*y2
    else:
        y1  = np.exp(r*t_arr); y2 = t_arr*np.exp(r*t_arr)
        dy1 = r*y1;            dy2 = np.exp(r*t_arr)*(1 + r*t_arr)

    W  = y1*dy2 - dy1*y2
    W  = np.where(np.abs(W) < 1e-14, 1e-14, W)
    u1 = np.cumsum(-y2*f_arr/(m*W)) * dt
    u2 = np.cumsum( y1*f_arr/(m*W)) * dt
    return u1*y1 + u2*y2


def eval_f_numeric(t_val, ftype, A, p2, t0, eps, custom_expr):
    """Evalúa f(t) numéricamente en un punto escalar."""
    try:
        if ftype == "Exponencial":
            return A * np.exp(-p2 * t_val)
        elif ftype == "Coseno":
            return A * np.cos(p2 * t_val)
        elif ftype == "Seno":
            return A * np.sin(p2 * t_val)
        elif ftype == "Polinómica":
            return A * np.power(t_val, max(1, int(round(p2)))) * np.exp(-0.25*t_val)
        elif ftype == "Impulso Dirac":
            return (A/eps) if abs(t_val - t0) < eps/2 else 0.0
        elif ftype == "Personalizada":
            return float(eval(custom_expr.replace('np.','').replace('sin(','np.sin(')
                              .replace('cos(','np.cos(').replace('exp(','np.exp('),
                              {"t": t_val, "np": np, "sin": np.sin, "cos": np.cos,
                               "exp": np.exp, "sqrt": np.sqrt, "log": np.log,
                               "pi": np.pi, "__builtins__": {}}))
    except Exception:
        return 0.0


def _f_expr_str(ftype, A, p2, t0, eps, custom_expr):
    """Devuelve el string de f(t) para ser parseado por SymPy."""
    if ftype == "Exponencial":
        return f"{A}*exp(-{p2}*t)"
    elif ftype == "Coseno":
        return f"{A}*cos({p2}*t)"
    elif ftype == "Seno":
        return f"{A}*sin({p2}*t)"
    elif ftype == "Polinómica":
        n = max(1, int(round(p2)))
        return f"{A}*t**{n}*exp(-0.25*t)"
    elif ftype == "Personalizada":
        return custom_expr
    else:
        return None   # Dirac no tiene forma simbólica cerrada


def solve_ode_numeric(m, b, k, y0_val, dy0_val, ftype, A, p2, t0, eps, custom_expr):
    """Integra la EDO completa con método de Euler (incluye condiciones iniciales)."""
    dt = T_MAX / (N_POINTS - 1)
    y  = np.zeros(N_POINTS); v = np.zeros(N_POINTS)
    y[0] = y0_val; v[0] = dy0_val
    for i in range(1, N_POINTS):
        f  = eval_f_numeric(T_ARRAY[i-1], ftype, A, p2, t0, eps, custom_expr)
        ac = (f - b*v[i-1] - k*y[i-1]) / m
        v[i] = v[i-1] + ac*dt
        y[i] = y[i-1] + v[i]*dt
    return y


def build_all_steps(m, b, k, y0_val, dy0_val,
                    ftype, A, p2, t0, eps, custom_expr,
                    yp_sympy, yp_steps_txt, yp_fn, yp_arr):
    """
    Genera los textos de los 5 pasos + régimen para la UI.
    Recibe los resultados ya calculados por solve_yp_symbolic.
    """
    disc = b**2 - 4*m*k
    r    = -b / (2*m)

    # ── PASO 1 ────────────────────────────────────────────────────────────────
    p1 = (f"Modelo de ciberseguridad — oscilador amortiguado forzado:\n\n"
          f"  {m:.2f}·y''(t) + {b:.2f}·y'(t) + {k:.2f}·y(t) = f(t)\n\n"
          f"Clasificación:\n"
          f"  • Ordinaria (ODE)\n"
          f"  • Orden 2 (involucra y'', y', y)\n"
          f"  • Lineal — coeficientes no dependen de y\n"
          f"  • Coeficientes constantes: m={m:.2f}, b={b:.2f}, k={k:.2f}\n"
          f"  • No homogénea: f(t) ≠ 0\n\n"
          f"Estrategia de solución:\n"
          f"  y(t) = yc(t)  [complementaria]  +  yp(t)  [particular]")

    # ── PASO 2: yc ────────────────────────────────────────────────────────────
    p2_txt  = f"Ecuación característica:\n  {m:.2f}·r² + {b:.2f}·r + {k:.2f} = 0\n\n"
    p2_txt += f"Discriminante:\n  Δ = ({b:.2f})² − 4·({m:.2f})·({k:.2f}) = {disc:.4f}\n\n"

    if disc < -1e-6:
        omega = np.sqrt(-disc)/(2*m)
        p2_txt += f"Δ < 0 → Raíces COMPLEJAS CONJUGADAS:\n  r = {r:.4f} ± {omega:.4f}i\n\n"
        p2_txt += f"yc(t) = e^({r:.4f}t)·[C₁·sin({omega:.4f}t) + C₂·cos({omega:.4f}t)]"
        yc_str  = f"e^({r:.4f}t)·[C₁·sin({omega:.4f}t)+C₂·cos({omega:.4f}t)]"
        regime  = "SUBAMORTIGUADO"
        regime_detail = (f"Δ < 0 → raíces complejas.\nLa red OSCILA antes de estabilizarse.\n"
                         f"El malware reinfecta nodos ya limpios.\n\n"
                         f"b² = {b**2:.3f}  <  4mk = {4*m*k:.3f}")
    elif disc > 1e-6:
        sq = np.sqrt(disc); r1,r2 = (-b+sq)/(2*m), (-b-sq)/(2*m)
        p2_txt += f"Δ > 0 → Raíces REALES DISTINTAS:\n  r₁={r1:.4f},  r₂={r2:.4f}\n\n"
        p2_txt += f"yc(t) = C₁·e^({r1:.4f}t) + C₂·e^({r2:.4f}t)"
        yc_str  = f"C₁·e^({r1:.4f}t)+C₂·e^({r2:.4f}t)"
        regime  = "SOBREAMORTIGUADO"
        regime_detail = (f"Δ > 0 → raíces reales distintas.\nRecuperación MONÓTONA y lenta.\n"
                         f"Sin rebotes, pero mayor ventana de vulnerabilidad.\n\n"
                         f"b² = {b**2:.3f}  >  4mk = {4*m*k:.3f}")
    else:
        p2_txt += f"Δ = 0 → Raíz REAL DOBLE:\n  r = {r:.4f}\n\n"
        p2_txt += f"yc(t) = (C₁ + C₂·t)·e^({r:.4f}t)"
        yc_str  = f"(C₁+C₂·t)·e^({r:.4f}t)"
        regime  = "AMORTIGUAMIENTO CRÍTICO"
        regime_detail = (f"Δ = 0 → raíz real doble.\nNeutralización ÓPTIMA: máxima velocidad\n"
                         f"sin oscilaciones. Firewall perfectamente calibrado.\n\n"
                         f"b_óptimo = 2√(mk) = {2*np.sqrt(m*k):.4f}\n"
                         f"b actual = {b:.4f}")

    # ── PASO 3: yp ────────────────────────────────────────────────────────────
    p3_txt = yp_steps_txt if yp_steps_txt else "No se pudo calcular yp."

    # Valores numéricos de yp para mostrar
    if yp_arr is not None and len(yp_arr) > 1:
        peak_yp = float(np.max(np.abs(yp_arr)))
        t_peak  = T_ARRAY[int(np.argmax(np.abs(yp_arr)))]
        v5  = yp_arr[min(int(5/T_MAX*N_POINTS),  N_POINTS-1)]
        v10 = yp_arr[min(int(10/T_MAX*N_POINTS), N_POINTS-1)]
        p3_txt += (f"\n\nValores numéricos de yp(t):\n"
                   f"  yp(5)    = {v5:.5f}  ({v5*100:.2f}% de red)\n"
                   f"  yp(10)   = {v10:.5f}  ({v10*100:.2f}% de red)\n"
                   f"  |yp|máx  = {peak_yp:.5f}  en t = {t_peak:.2f}s")

    # ── PASO 4: Solución general ──────────────────────────────────────────────
    yp_short = str(yp_sympy) if yp_sympy is not None else "yp(t) [numérica]"
    p4_txt = (f"Solución general (superposición lineal):\n\n"
              f"  y(t) = yc(t) + yp(t)\n\n"
              f"  yc(t) = {yc_str}\n\n"
              f"  yp(t) = {yp_short}\n\n"
              f"Las constantes C₁ y C₂ se determinan a continuación\n"
              f"aplicando las condiciones iniciales de la red.")

    # ── PASO 5: Condiciones iniciales ─────────────────────────────────────────
    yp0  = float(yp_arr[0])  if yp_arr is not None and len(yp_arr) > 0 else 0.0
    dyp0 = float((yp_arr[1]-yp_arr[0])/(T_MAX/(N_POINTS-1))) if yp_arr is not None and len(yp_arr) > 1 else 0.0

    C1, C2 = 0.0, 0.0
    if disc < -1e-6:
        om = np.sqrt(-disc)/(2*m)
        C1 = y0_val - yp0
        C2 = (dy0_val - dyp0 - r*C1) / (om or 1e-12)
        yf = (f"y(t) = e^({r:.4f}t)·[{C1:.5f}·sin({om:.4f}t) + {C2:.5f}·cos({om:.4f}t)]\n"
              f"       + {yp_short}")
    elif disc > 1e-6:
        sq = np.sqrt(disc); r1,r2 = (-b+sq)/(2*m), (-b-sq)/(2*m)
        den = (r2-r1) or 1e-12
        C2  = (dy0_val - dyp0 - r1*(y0_val-yp0)) / den
        C1  = y0_val - yp0 - C2
        yf  = (f"y(t) = {C1:.5f}·e^({r1:.4f}t) + {C2:.5f}·e^({r2:.4f}t)\n"
               f"       + {yp_short}")
    else:
        C1 = y0_val - yp0
        C2 = dy0_val - dyp0 - r*C1
        yf = (f"y(t) = ({C1:.5f} + {C2:.5f}·t)·e^({r:.4f}t)\n"
              f"       + {yp_short}")

    net  = ("Alta exposición — respuesta urgente." if y0_val > 0.4
            else "Infección moderada — firewall puede contenerla." if y0_val > 0.15
            else "Red mayormente sana — ataque en etapa temprana.")
    rate = ("Propagación ACTIVA." if dy0_val > 0.08
            else "Contención ya iniciada." if dy0_val < 0 else "Inicio lento del brote.")

    p5_txt = (f"Condiciones iniciales:\n\n"
              f"  y(0)  = {y0_val*100:.1f}% nodos comprometidos  →  y₀ = {y0_val:.4f}\n"
              f"  y'(0) = {dy0_val:.3f}  →  {rate}\n"
              f"  yp(0) = {yp0:.5f}\n\n"
              f"Sistema de ecuaciones en C₁, C₂:\n"
              f"  y(0)  = [yc(0) evaluada] + yp(0) = y₀\n"
              f"  y'(0) = [yc'(0) evaluada] + yp'(0) = dy₀\n\n"
              f"  C₁ = {C1:.6f}\n"
              f"  C₂ = {C2:.6f}\n\n"
              f"▶ Solución única:\n{yf}\n\n"
              f"Interpretación: {net}")

    return {
        "p1": p1, "p2": p2_txt, "p3": p3_txt, "p4": p4_txt, "p5": p5_txt,
        "regime": regime, "regime_detail": regime_detail, "disc": disc
    }



#  CANVAS MATPLOTLIB


class ODECanvas(FigureCanvas):
    def __init__(self):
        self.fig = Figure(figsize=(9, 3.0), facecolor=PLOT_BG, tight_layout=True)
        self.ax  = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        ax = self.ax
        ax.set_facecolor(PLOT_BG)
        ax.tick_params(colors=CLR_TEXT_MUT, labelsize=8)
        for sp_ in ax.spines.values(): sp_.set_color(PLOT_AXES)
        ax.grid(True, color=PLOT_GRID, lw=0.5, alpha=0.7)
        ax.set_xlim(0, T_MAX); ax.set_ylim(-0.08, 0.6)
        ax.set_xlabel("tiempo t", color=CLR_TEXT_MUT, fontsize=9)
        ax.set_ylabel("y(t) — red comprometida", color=CLR_TEXT_MUT, fontsize=9)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,_: f"{x*100:.0f}%"))

        self.line_y,  = ax.plot([], [], color=PLOT_Y,  lw=2.2, label="y(t) — infección", zorder=5)
        self.line_yp, = ax.plot([], [], color=PLOT_YP, lw=1.6, ls="-.", alpha=0.75, label="yp(t) — particular", zorder=4)
        self.line_f,  = ax.plot([], [], color=PLOT_F,  lw=1.2, ls="--", alpha=0.55, label="f(t) × 0.1 — ataque", zorder=3)
        ax.plot(T_ARRAY, [0.01]*N_POINTS, color=PLOT_STAB, lw=0.8, ls=":", alpha=0.5, label="umbral 1%", zorder=2)
        ax.fill_between(T_ARRAY, 0, 0.01, alpha=0.05, color=PLOT_STAB)
        ax.legend(loc="upper right", fontsize=7.5, framealpha=0.3,
                  facecolor=CLR_PANEL, edgecolor=CLR_BORDER,
                  labelcolor=CLR_TEXT, ncol=2)
        self.pk_dot = ax.scatter([], [], color=PLOT_Y, s=55, zorder=10, edgecolors="white", linewidths=0.8)
        self.pk_ann = ax.annotate("", xy=(0,0), xytext=(1,0.1),
                        color=PLOT_Y, fontsize=7.5,
                        arrowprops=dict(arrowstyle="->", color=PLOT_Y, lw=0.8),
                        bbox=dict(boxstyle="round,pad=0.3", fc=CLR_PANEL, ec=PLOT_Y, alpha=0.85))

    def update_plot(self, y_arr, f_arr, yp_arr):
        self.line_y.set_data(T_ARRAY, y_arr)
        self.line_yp.set_data(T_ARRAY, np.clip(yp_arr, -0.5, 1.1))
        self.line_f.set_data(T_ARRAY, np.clip(f_arr*0.10, -0.2, 0.5))
        peak_v = float(np.max(y_arr)); peak_i = int(np.argmax(y_arr))
        self.pk_dot.set_offsets([[T_ARRAY[peak_i], peak_v]])
        self.pk_ann.xy = (T_ARRAY[peak_i], peak_v)
        tx = min(T_ARRAY[peak_i]+1.5, T_MAX-4)
        self.pk_ann.set_position((tx, peak_v+0.04))
        self.pk_ann.set_text(f"Pico\n{peak_v*100:.1f}%\nt={T_ARRAY[peak_i]:.1f}s")
        self.ax.set_ylim(-0.06, min(max(0.15, peak_v*1.3), 1.15))
        self.draw()



#  HELPERS DE UI


def make_label(text, size=12, color=CLR_TEXT, bold=False):
    lbl = QLabel(text)
    f = QFont(); f.setPointSize(size); f.setBold(bold)
    lbl.setFont(f)
    lbl.setStyleSheet(f"color:{color};background:transparent;")
    return lbl

def make_slider(mn, mx, val, decimals=1):
    s = QSlider(Qt.Horizontal); factor = 10**decimals
    s.setMinimum(int(mn*factor)); s.setMaximum(int(mx*factor)); s.setValue(int(val*factor))
    s.setStyleSheet("""QSlider::groove:horizontal{background:#2a2d3e;height:4px;border-radius:2px}
        QSlider::handle:horizontal{background:#4f9cf9;width:14px;height:14px;border-radius:7px;margin:-5px 0}
        QSlider::sub-page:horizontal{background:#4f9cf9;border-radius:2px}""")
    return s, factor

def make_val_lbl(text, color=CLR_ACCENT):
    lbl = QLabel(text)
    lbl.setFont(QFont("Consolas", 11))
    lbl.setStyleSheet(f"color:{color};background:transparent;font-weight:bold;")
    lbl.setMinimumWidth(52); lbl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
    return lbl

def make_groupbox(title):
    gb = QGroupBox(title)
    gb.setStyleSheet(f"""QGroupBox{{color:{CLR_TEXT_MUT};border:0.5px solid {CLR_BORDER};
        border-radius:10px;margin-top:10px;padding-top:10px;font-size:10px;
        font-weight:bold;letter-spacing:1px}}
        QGroupBox::title{{subcontrol-origin:margin;left:12px;padding:0 6px;background:{CLR_PANEL}}}""")
    return gb

def make_step_frame(num, title, color=CLR_BORDER):
    frame = QFrame()
    frame.setStyleSheet(f"QFrame{{background:{CLR_PANEL2};border-left:2px solid {color};border-radius:0 8px 8px 0}}")
    lay = QVBoxLayout(frame); lay.setSpacing(4); lay.setContentsMargins(12,10,12,10)
    lay.addWidget(make_label(num, 8, CLR_TEXT_DIM))
    lay.addWidget(make_label(title, 11, bold=True))
    return frame, lay

def make_math_box(color=CLR_CYAN):
    box = QTextEdit(); box.setReadOnly(True)
    box.setFont(QFont("Consolas", 10))
    box.setStyleSheet(f"QTextEdit{{background:{CLR_PANEL};color:{color};border:0.5px solid {CLR_BORDER};border-radius:6px;padding:8px}}")
    box.setMaximumHeight(180); box.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
    return box

def make_chip(text, active=False):
    btn = QPushButton(text); btn.setCheckable(True); btn.setChecked(active)
    btn.setStyleSheet(f"""QPushButton{{background:{CLR_PANEL2};color:{CLR_TEXT_MUT};
        border:0.5px solid {CLR_BORDER};border-radius:12px;padding:4px 11px;font-size:10px}}
        QPushButton:checked{{background:#1a2e4a;color:{CLR_ACCENT};border-color:{CLR_ACCENT}}}
        QPushButton:hover{{background:#1f2336}}""")
    return btn



#  VENTANA PRINCIPAL

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulador de Malware — EDO 2do Orden Amortiguada")
        self.setMinimumSize(1100, 740); self.resize(1240, 900)
        # Estado
        self.m_val=1.0; self.b_val=1.0; self.k_val=2.0
        self.A_val=0.5; self.p2_val=1.0
        self.y0_val=0.10; self.dy0_val=0.05
        self.t0_val=3.0;  self.eps_val=0.2
        self.ftype="Exponencial"
        self.custom_expr="0.5*exp(-t)"
        # Debounce
        self._timer = QTimer(); self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._do_update)
        self._build_ui()
        self._do_update()

    # ─── UI ──────────────────────────────────────────────────────────────────
    def _build_ui(self):
        cw = QWidget(); self.setCentralWidget(cw)
        cw.setStyleSheet(f"background:{CLR_BG};")
        root = QVBoxLayout(cw); root.setSpacing(8); root.setContentsMargins(12,12,12,12)

        # Título
        tr = QHBoxLayout()
        tr.addWidget(make_label("Simulador de propagación de malware", 15, bold=True))
        self.regime_badge = make_label("SUBAMORTIGUADO", 9, CLR_RED, bold=True)
        self.regime_badge.setStyleSheet(
            f"color:{CLR_RED};background:#2d1a1a;border:0.5px solid {CLR_RED};"
            "border-radius:10px;padding:2px 10px;font-size:9px;font-weight:bold;")
        self.eq_lbl = make_label("", 9, CLR_TEXT_MUT)
        self.eq_lbl.setFont(QFont("Consolas", 9))
        tr.addWidget(self.regime_badge); tr.addStretch(); tr.addWidget(self.eq_lbl)
        root.addLayout(tr)

        self.canvas = ODECanvas(); self.canvas.setFixedHeight(225)
        root.addWidget(self.canvas)

        # Stats
        sr = QHBoxLayout(); sr.setSpacing(8)
        self.s_delta = self._stat("—","Δ = b²−4mk")
        self.s_peak  = self._stat("—","Pico de infección")
        self.s_stab  = self._stat("—","Tiempo estable (<1%)")
        for w,_ in [self.s_delta,self.s_peak,self.s_stab]: sr.addWidget(w)
        root.addLayout(sr)

        body = QHBoxLayout(); body.setSpacing(10)
        body.addLayout(self._left_panel(), stretch=5)
        body.addWidget(self._right_panel(), stretch=6)
        root.addLayout(body)

    def _stat(self, val, lbl):
        f = QFrame()
        f.setStyleSheet(f"background:{CLR_PANEL};border:0.5px solid {CLR_BORDER};border-radius:8px;")
        ly = QVBoxLayout(f); ly.setContentsMargins(10,7,10,7)
        v = make_label(val, 13, bold=True); v.setAlignment(Qt.AlignCenter)
        l = make_label(lbl, 8, CLR_TEXT_MUT); l.setAlignment(Qt.AlignCenter)
        ly.addWidget(v); ly.addWidget(l)
        return f, v

    def _left_panel(self):
        col = QVBoxLayout(); col.setSpacing(8)

        # ── m, b, k ──────────────────────────────────────────────────────────
        gb = make_groupbox("PARÁMETROS DEL SISTEMA")
        ly = QVBoxLayout(gb); ly.setSpacing(5); ly.setContentsMargins(10,14,10,10)
        params = [
            ("m — inercia de red", 0.1,5.0,1.0,"m_val","sl_m","vl_m",CLR_ACCENT,
             "Latencia y complejidad de la red.\nAlto = red grande/compleja.\nMalware tarda más en propagarse y en eliminarse."),
            ("b — firewall / defensa", 0.1,10.0,1.0,"b_val","sl_b","vl_b",CLR_ACCENT,
             "Efectividad del firewall + respuesta humana.\nAlto = contención rápida.\nb² = 4mk → configuración ÓPTIMA."),
            ("k — recuperación automática", 0.1,8.0,2.0,"k_val","sl_k","vl_k",CLR_ACCENT,
             "Auto-limpieza: antivirus, sandboxing, snapshots.\nAlto = la red se cura sola rápido."),
        ]
        for lbl_t,mn,mx,dflt,attr,sln,vln,clr,tip in params:
            row = QHBoxLayout()
            lbl = make_label(lbl_t,10,CLR_TEXT_MUT); lbl.setToolTip(tip)
            sl,f = make_slider(mn,mx,dflt); vl = make_val_lbl(f"{dflt:.1f}",clr)
            setattr(self,sln,sl); setattr(self,vln,vl)
            def on(v,f=f,vl=vl,a=attr): val=v/f; vl.setText(f"{val:.1f}"); setattr(self,a,val); self._sched()
            sl.valueChanged.connect(on)
            row.addWidget(lbl,4); row.addWidget(sl,5); row.addWidget(vl,1); ly.addLayout(row)
        # Presets
        pr = QHBoxLayout(); pr.addWidget(make_label("Preset:",9,CLR_TEXT_DIM))
        for lt,vals in [("Firewall débil",(1,1,2)),("Respuesta lenta",(1,3,2)),("Óptimo",(1,2,1))]:
            btn = QPushButton(lt)
            btn.setStyleSheet(f"font-size:10px;padding:3px 8px;background:{CLR_PANEL2};color:{CLR_TEXT_MUT};border:0.5px solid {CLR_BORDER};border-radius:6px;")
            btn.clicked.connect(lambda _,v=vals: self._preset(*v)); pr.addWidget(btn)
        ly.addLayout(pr); col.addWidget(gb)

        # ── Condiciones iniciales ─────────────────────────────────────────────
        gb2 = make_groupbox("CONDICIONES INICIALES DE LA RED")
        ly2 = QVBoxLayout(gb2); ly2.setSpacing(5); ly2.setContentsMargins(10,14,10,10)
        ics = [
            ("y(0) — nodos comprometidos",0.0,0.95,0.10,"y0_val","sl_y0","vl_y0",
             lambda v:f"{v*100:.0f}%",CLR_PURPLE,2,
             "% de equipos ya infectados al inicio.\n0%=sana, 95%=comprometida."),
            ("y'(0) — velocidad de propagación",-0.5,0.5,0.05,"dy0_val","sl_dy0","vl_dy0",
             lambda v:f"{'+' if v>=0 else ''}{v:.2f}",CLR_PURPLE,2,
             "Tasa de propagación en t=0.\n+: expansión activa.  −: ya contenido."),
        ]
        for lbl_t,mn,mx,dflt,attr,sln,vln,fmt,clr,dec,tip in ics:
            row = QHBoxLayout()
            lbl = make_label(lbl_t,10,CLR_TEXT_MUT); lbl.setToolTip(tip)
            sl,f = make_slider(mn,mx,dflt,dec); vl = make_val_lbl(fmt(dflt),clr)
            setattr(self,sln,sl); setattr(self,vln,vl)
            def on_ic(v,f=f,vl=vl,a=attr,fn=fmt): val=v/f; vl.setText(fn(val)); setattr(self,a,val); self._sched()
            sl.valueChanged.connect(on_ic)
            row.addWidget(lbl,4); row.addWidget(sl,5); row.addWidget(vl,1); ly2.addLayout(row)
        col.addWidget(gb2)

        # ── f(t) ──────────────────────────────────────────────────────────────
        gb3 = make_groupbox("FUERZA DEL ATAQUE f(t)")
        ly3 = QVBoxLayout(gb3); ly3.setSpacing(5); ly3.setContentsMargins(10,14,10,10)
        cr = QHBoxLayout(); cr.setSpacing(3)
        self.chips = {}
        for ft in ["Exponencial","Coseno","Seno","Polinómica","Impulso Dirac","Personalizada"]:
            btn = make_chip(ft, ft=="Exponencial")
            btn.clicked.connect(lambda _,t=ft: self._set_ftype(t))
            self.chips[ft]=btn; cr.addWidget(btn)
        ly3.addLayout(cr)
        self.fx_hint = make_label("",8,CLR_TEXT_DIM); self.fx_hint.setWordWrap(True); ly3.addWidget(self.fx_hint)
        self.fx_inp = QLineEdit("0.5*exp(-t)")
        self.fx_inp.setFont(QFont("Consolas",10))
        self.fx_inp.setStyleSheet(f"background:{CLR_PANEL2};color:{CLR_CYAN};border:0.5px solid {CLR_BORDER};border-radius:6px;padding:5px 8px;")
        self.fx_inp.setPlaceholderText("f(t) en t — ej: 3*sin(2*t), t**2+1, exp(-t)*cos(t)")
        self.fx_inp.textChanged.connect(self._sched)
        ly3.addWidget(self.fx_inp)
        for lbl_t,mn,mx,dflt,attr,sln,vln,clr in [
            ("Amplitud A",0.1,3.0,0.5,"A_val","sl_A","vl_A",CLR_ORANGE),
            ("Parámetro α/ω/n",0.1,5.0,1.0,"p2_val","sl_p2","vl_p2",CLR_YELLOW),
        ]:
            row = QHBoxLayout()
            sl,f = make_slider(mn,mx,dflt); vl = make_val_lbl(f"{dflt:.1f}",clr)
            setattr(self,sln,sl); setattr(self,vln,vl)
            def on_f(v,f=f,vl=vl,a=attr): val=v/f; vl.setText(f"{val:.1f}"); setattr(self,a,val); self._sched()
            sl.valueChanged.connect(on_f)
            row.addWidget(make_label(lbl_t,10,CLR_TEXT_MUT),3); row.addWidget(sl,5); row.addWidget(vl,1); ly3.addLayout(row)
        # Dirac extras
        self.dirac_frame = QFrame()
        self.dirac_frame.setStyleSheet(f"background:{CLR_PANEL2};border-radius:6px;border:0.5px solid {CLR_BORDER};")
        df = QVBoxLayout(self.dirac_frame); df.setContentsMargins(8,8,8,8); df.setSpacing(4)
        for lbl_t,mn,mx,dflt,attr,sln,vln in [
            ("t₀ — momento del golpe",0.5,20.0,3.0,"t0_val","sl_t0","vl_t0"),
            ("ε — ancho del pulso",0.05,1.0,0.2,"eps_val","sl_eps","vl_eps"),
        ]:
            row = QHBoxLayout()
            sl,f = make_slider(mn,mx,dflt,2); vl = make_val_lbl(f"{dflt:.2f}",CLR_RED)
            setattr(self,sln,sl); setattr(self,vln,vl)
            def on_d(v,f=f,vl=vl,a=attr): val=v/f; vl.setText(f"{val:.2f}"); setattr(self,a,val); self._sched()
            sl.valueChanged.connect(on_d)
            row.addWidget(make_label(lbl_t,10,CLR_TEXT_MUT),3); row.addWidget(sl,5); row.addWidget(vl,1); df.addLayout(row)
        note = make_label("δ(t−t₀) ≈ pulso rect. altura=P/ε, duración=ε → 0\n'Martillazo' informático: exploit día cero en t₀.",8,CLR_TEXT_DIM)
        note.setWordWrap(True); df.addWidget(note)
        self.dirac_frame.setVisible(False); ly3.addWidget(self.dirac_frame)
        col.addWidget(gb3); col.addStretch()
        return col

    def _right_panel(self):
        scroll = QScrollArea(); scroll.setWidgetResizable(True)
        scroll.setStyleSheet(f"QScrollArea{{border:none;background:{CLR_BG};}}"
                             f"QScrollBar:vertical{{background:{CLR_PANEL};width:6px;border-radius:3px;}}"
                             f"QScrollBar::handle:vertical{{background:{CLR_BORDER};border-radius:3px;}}")
        content = QWidget(); content.setStyleSheet(f"background:{CLR_BG};")
        lay = QVBoxLayout(content); lay.setSpacing(8); lay.setContentsMargins(0,0,8,0)
        hdr = QHBoxLayout()
        hdr.addWidget(make_label("Desarrollo analítico paso a paso",12,bold=True))
        hdr.addWidget(make_label("— actualizado en tiempo real",9,CLR_TEXT_DIM))
        hdr.addStretch(); lay.addLayout(hdr)
        steps_def = [
            ("PASO 1","Planteamiento de la EDO",CLR_BORDER,CLR_CYAN),
            ("PASO 2","Solución complementaria yc — ecuación característica",CLR_ACCENT,CLR_CYAN),
            ("PASO 3","Solución particular yp — Coeficientes Indeterminados",CLR_GREEN,CLR_GREEN),
            ("PASO 4","Solución general  y(t) = yc + yp",CLR_PURPLE,CLR_CYAN),
            ("PASO 5","Condiciones iniciales → solución única",CLR_YELLOW,CLR_YELLOW),
        ]
        self.math_boxes = []
        for num,title,color,tclr in steps_def:
            frame,fl = make_step_frame(num,title,color)
            box = make_math_box(tclr); self.math_boxes.append(box)
            fl.addWidget(box); lay.addWidget(frame)
        reg_frame,rl = make_step_frame("RÉGIMEN","Clasificación y significado en ciberseguridad",CLR_RED)
        self.regime_box = make_math_box(CLR_YELLOW); rl.addWidget(self.regime_box)
        lay.addWidget(reg_frame); lay.addStretch()
        scroll.setWidget(content); return scroll

    # ─── Lógica ──────────────────────────────────────────────────────────────
    def _sched(self): self._timer.start(80)

    def _do_update(self):
        custom = self.fx_inp.text().strip()

        # 1. Array numérico de f(t)
        f_arr = np.array([eval_f_numeric(t, self.ftype, self.A_val, self.p2_val,
                                          self.t0_val, self.eps_val, custom)
                          for t in T_ARRAY])

        # 2. Calcular yp simbólico + numérico
        expr_str = _f_expr_str(self.ftype, self.A_val, self.p2_val,
                                self.t0_val, self.eps_val, custom)

        if expr_str and self.ftype != "Impulso Dirac":
            yp_sym, yp_steps, yp_fn = solve_yp_symbolic(
                expr_str, self.m_val, self.b_val, self.k_val)
            if yp_fn is not None:
                try:
                    yp_arr = np.array([float(yp_fn(t)) for t in T_ARRAY])
                    yp_arr = np.nan_to_num(yp_arr, nan=0.0, posinf=0.0, neginf=0.0)
                except Exception:
                    yp_arr = variation_of_params(self.m_val, self.b_val, self.k_val, f_arr, T_ARRAY)
            else:
                yp_arr = variation_of_params(self.m_val, self.b_val, self.k_val, f_arr, T_ARRAY)
        else:
            # Dirac u otras sin forma simbólica → variación de parámetros
            yp_sym   = None
            yp_steps = (f"f(t) = P·δ(t−{self.t0_val:.1f})\n\n"
                        "La función Delta de Dirac no admite Coef. Indeterminados.\n"
                        "Se usa Variación de Parámetros numéricamente:\n\n"
                        "  yp(t) = y₁·u₁ + y₂·u₂\n"
                        "  u₁ = ∫[−y₂·f/(m·W)]dt,  u₂ = ∫[y₁·f/(m·W)]dt")
            yp_arr = variation_of_params(self.m_val, self.b_val, self.k_val, f_arr, T_ARRAY)

        # 3. EDO completa
        y_arr = solve_ode_numeric(self.m_val, self.b_val, self.k_val,
                                   self.y0_val, self.dy0_val,
                                   self.ftype, self.A_val, self.p2_val,
                                   self.t0_val, self.eps_val, custom)

        # 4. Gráfico
        self.canvas.update_plot(y_arr, f_arr, yp_arr)

        # 5. Stats
        self._update_stats(y_arr)

        # 6. Desarrollo analítico
        steps = build_all_steps(self.m_val, self.b_val, self.k_val,
                                 self.y0_val, self.dy0_val,
                                 self.ftype, self.A_val, self.p2_val,
                                 self.t0_val, self.eps_val, custom,
                                 yp_sym, yp_steps, None, yp_arr)
        for i, key in enumerate(["p1","p2","p3","p4","p5"]):
            self.math_boxes[i].setText(steps[key])
        self.regime_box.setText(
            f"{steps['regime']}\n\n{steps['regime_detail']}\n\nΔ = {steps['disc']:.5f}")

    def _update_stats(self, y_arr):
        disc = self.b_val**2 - 4*self.m_val*self.k_val
        self.s_delta[1].setText(f"{disc:.3f}")
        pv = float(np.max(y_arr)); pi = int(np.argmax(y_arr))
        self.s_peak[1].setText(f"{pv*100:.1f}%  t={T_ARRAY[pi]:.1f}s")
        si = next((i for i in range(pi,N_POINTS) if abs(y_arr[i])<0.01), None)
        self.s_stab[1].setText(f"t≈{T_ARRAY[si]:.1f}s" if si else f">{T_MAX}s")
        self.eq_lbl.setText(f"{self.m_val:.1f}·y''+{self.b_val:.1f}·y'+{self.k_val:.1f}·y=f(t)")
        if disc < -0.01:
            self.regime_badge.setText("SUBAMORTIGUADO")
            self.regime_badge.setStyleSheet(f"color:{CLR_RED};background:#2d1a1a;border:0.5px solid {CLR_RED};border-radius:10px;padding:2px 10px;font-size:9px;font-weight:bold;")
        elif disc > 0.01:
            self.regime_badge.setText("SOBREAMORTIGUADO")
            self.regime_badge.setStyleSheet(f"color:{CLR_YELLOW};background:#2d2510;border:0.5px solid {CLR_YELLOW};border-radius:10px;padding:2px 10px;font-size:9px;font-weight:bold;")
        else:
            self.regime_badge.setText("AMORTIGUAMIENTO CRÍTICO")
            self.regime_badge.setStyleSheet(f"color:{CLR_GREEN};background:#0d2d1f;border:0.5px solid {CLR_GREEN};border-radius:10px;padding:2px 10px;font-size:9px;font-weight:bold;")

    def _preset(self, m, b, k):
        self.m_val,self.b_val,self.k_val = float(m),float(b),float(k)
        for a,sln,vln in [("m_val","sl_m","vl_m"),("b_val","sl_b","vl_b"),("k_val","sl_k","vl_k")]:
            getattr(self,sln).setValue(int(getattr(self,a)*10))
            getattr(self,vln).setText(f"{getattr(self,a):.1f}")
        self._sched()

    def _set_ftype(self, ftype):
        self.ftype = ftype
        for ft,btn in self.chips.items(): btn.setChecked(ft==ftype)
        hints = {
            "Exponencial":   "f(t)=A·e^(-αt) — DDoS que se mitiga solo.",
            "Coseno":        "f(t)=A·cos(ωt) — ataque periódico cíclico.",
            "Seno":          "f(t)=A·sin(ωt) — oleadas graduales de acceso.",
            "Polinómica":    "f(t)=A·t^n·e^(-0.25t) — malware en escalada.",
            "Impulso Dirac": "f(t)=P·δ(t−t₀) — exploit instantáneo (día cero).",
            "Personalizada": "Escribí cualquier f(t): 6*exp(-2*t), 3*sin(t), t**2+1, sin(t)*exp(-t)…",
        }
        exprs = {
            "Exponencial":   "0.5*exp(-t)",
            "Coseno":        "0.5*cos(1.0*t)",
            "Seno":          "0.5*sin(1.0*t)",
            "Polinómica":    "0.1*t**2*exp(-0.25*t)",
            "Impulso Dirac": "",
            "Personalizada": "6*exp(-2*t)",
        }
        self.fx_hint.setText(hints.get(ftype,""))
        self.fx_inp.setText(exprs.get(ftype,""))
        self.fx_inp.setEnabled(ftype != "Impulso Dirac")
        self.dirac_frame.setVisible(ftype=="Impulso Dirac")
        self._sched()




def main():
    app = QApplication(sys.argv); app.setStyle("Fusion")
    palette = QPalette()
    for role, hx in [
        (QPalette.Window, CLR_BG),(QPalette.WindowText, CLR_TEXT),
        (QPalette.Base, CLR_PANEL),(QPalette.AlternateBase, CLR_PANEL2),
        (QPalette.Text, CLR_TEXT),(QPalette.Button, CLR_PANEL),
        (QPalette.ButtonText, CLR_TEXT),(QPalette.Highlight, CLR_ACCENT),
        (QPalette.HighlightedText,"#ffffff"),
        (QPalette.ToolTipBase, CLR_PANEL2),(QPalette.ToolTipText, CLR_TEXT),
    ]: palette.setColor(role, QColor(hx))
    app.setPalette(palette)
    win = MainWindow(); win.show(); sys.exit(app.exec_())

if __name__ == "__main__":
    main()
    