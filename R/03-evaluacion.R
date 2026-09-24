# 03-evaluacion.R
# Evaluación: medidas(), ljung_box(), jarque_bera(), durbin_watson(), validar_errores()
#
# Todas las pruebas devuelven una lista con los mismos campos, para que
# .reportar_prueba() pueda escribir los seis elementos exigidos por el enunciado:
#   nombre, H0, H1, formula, estadistico, gl, valor_critico, valor_p, region, decision
# (más n y los datos propios de cada prueba). H0, H1, formula y region vienen en LaTeX
# ($...$) porque el consumidor es informe.qmd.

# ---- Formato y auxiliares internas ----------------------------------------------------

# Coma decimal (convención de las notas del curso). Estadísticos con 4 cifras.
.fmt_num <- function(x, digitos = 4) {
  formatC(x, format = "f", digits = digitos, decimal.mark = ",")
}

# Valores p con 3 cifras; «< 0,001» cuando corresponda.
.fmt_p <- function(p) {
  if (is.na(p)) {
    "no disponible"
  } else if (p < 0.001) {
    "< 0,001"
  } else {
    .fmt_num(p, 3)
  }
}

# Decisión al 5 % a partir de un lógico «rechaza».
.decision <- function(rechaza) {
  if (rechaza) "Se rechaza $H_0$" else "No se rechaza $H_0$"
}

# ---- ljung_box ------------------------------------------------------------------------

#' ljung_box(r, T, m, p)
#'
#' Descripción: prueba conjunta de Ljung–Box de que las primeras m autocorrelaciones
#'   son nulas (Clase 3, Parte II).
#'
#' Ecuaciones:
#'   Q_m = T (T + 2) sum_{h=1}^{m} r_h^2 / (T - h),  Q_m ~ chi^2_{m-p} bajo H0.
#'   Grados de libertad: m - p, con p el número de parámetros que el método estimó sobre
#'   la muestra. Sobre la serie original p = 0 (gl = m). Sobre los errores de un método,
#'   cada parámetro estimado vuelve a los errores artificialmente menos autocorrelacionados
#'   que las perturbaciones; referir Q_m a chi^2_m en vez de chi^2_{m-p} deja la prueba
#'   conservadora (Clase 3).
#'
#' Sobre T: es el número de valores de la sucesión de la que salen las r_h, no el T de la
#'   serie original. Con errores de un método es el número de errores sin NA. Es el mismo
#'   n que usa Box.test de stats (la longitud de x).
#'
#' @param r  vector numérico con las autocorrelaciones r_1, r_2, ...; se usan las m primeras.
#' @param T  número de observaciones de la sucesión que generó las r_h.
#' @param m  número de rezagos de la prueba (entero, m > p).
#' @param p  número de parámetros estimados por el método (entero >= 0; 0 sobre la serie).
#' @return lista con los campos comunes de las pruebas (ver cabecera del archivo) y
#'   además n (= T), m y p.
#'
#' Referencia: enunciado, sección 2(d); Clase 3, Parte II (Ljung y Box, 1978).
ljung_box <- function(r, T, m, p) {
  stopifnot(
    "r debe ser un vector numérico sin NA" = is.numeric(r) && length(r) > 0L && !anyNA(r),
    "T debe ser un entero positivo" =
      is.numeric(T) && length(T) == 1L && !is.na(T) && T >= 1 && T == floor(T),
    "m debe ser un entero positivo" =
      is.numeric(m) && length(m) == 1L && !is.na(m) && m >= 1 && m == floor(m),
    "p debe ser un entero no negativo" =
      is.numeric(p) && length(p) == 1L && !is.na(p) && p >= 0 && p == floor(p)
  )
  stopifnot(
    "m debe ser mayor que p: si no, Ljung-Box no tiene grados de libertad (m - p > 0)" = m > p,
    "r debe tener al menos m autocorrelaciones" = length(r) >= m,
    "T debe ser mayor que m" = T > m
  )

  h <- seq_len(m)
  Q <- T * (T + 2) * sum(r[h]^2 / (T - h))
  gl <- m - p
  valor_critico <- stats::qchisq(0.95, df = gl)
  valor_p <- stats::pchisq(Q, df = gl, lower.tail = FALSE)

  list(
    nombre = "Prueba de Ljung–Box",
    H0 = r"($H_0:\ \rho_1 = \rho_2 = \cdots = \rho_m = 0$)",
    H1 = r"($H_1:\ \exists\, h \le m$ con $\rho_h \neq 0$)",
    formula = r"($Q_m = T(T+2)\sum_{h=1}^{m} \dfrac{r_h^2}{T-h}\ \overset{H_0}{\sim}\ \chi^2_{m-p}$)",
    estadistico = Q,
    gl = gl,
    valor_critico = valor_critico,
    valor_p = valor_p,
    region = sprintf(r"($Q_m > \chi^2_{0{,}95;\,%d} = %s$)", as.integer(gl), .fmt_num(valor_critico)),
    decision = .decision(Q > valor_critico),
    n = T, m = m, p = p
  )
}

# ---- jarque_bera ----------------------------------------------------------------------

#' jarque_bera(e)
#'
#' Descripción: prueba de normalidad de Jarque–Bera sobre los errores de un método.
#'
#' Ecuaciones (momentos centrales con divisor N, no N - 1):
#'   m_k = (1/N) sum (e_t - ebar)^k,  A = m_3 / m_2^(3/2),  K = m_4 / m_2^2,
#'   JB = N/6 [ A^2 + (K - 3)^2 / 4 ]  ~  chi^2_2 bajo H0 (asintótica).
#'   Bajo normalidad A = 0 y K = 3.
#'
#' Advertencia: es una prueba asintótica. Con N < 20 la función devuelve el campo
#'   `advertencia` (cadena; NULL si N >= 20) y la decisión debe reportarse con esa
#'   salvedad. La normalidad no mejora el pronóstico puntual: solo habilita intervalos con
#'   cuantiles normales (Clase 3).
#'
#' @param e  vector numérico sin NA (errores de un paso ya sin el calentamiento).
#' @return lista con los campos comunes de las pruebas y además n (= N), A, K y
#'   advertencia.
#'
#' Referencia: enunciado, sección 2(d); Clase 3, Parte III (Jarque y Bera, 1987).
jarque_bera <- function(e) {
  stopifnot(
    "e debe ser un vector numérico sin NA" = is.numeric(e) && length(e) > 0L && !anyNA(e),
    "se necesitan al menos 3 errores" = length(e) >= 3L
  )
  N <- length(e)
  desv <- e - mean(e)
  m2 <- mean(desv^2)
  stopifnot("los errores son todos iguales: asimetría y curtosis no están definidas" = m2 > 0)
  m3 <- mean(desv^3)
  m4 <- mean(desv^4)
  A <- m3 / m2^(3 / 2)
  K <- m4 / m2^2
  JB <- N / 6 * (A^2 + (K - 3)^2 / 4)
  valor_critico <- stats::qchisq(0.95, df = 2)

  list(
    nombre = "Prueba de Jarque–Bera",
    H0 = r"($H_0:\ A = 0\ \text{y}\ K = 3$ (errores normales))",
    H1 = r"($H_1:\ A \neq 0\ \text{o}\ K \neq 3$)",
    formula = r"($JB = \dfrac{N}{6}\left[A^2 + \dfrac{(K-3)^2}{4}\right]\ \overset{H_0}{\sim}\ \chi^2_2$, con $A = m_3/m_2^{3/2}$ y $K = m_4/m_2^2$)",
    estadistico = JB,
    gl = 2L,
    valor_critico = valor_critico,
    valor_p = stats::pchisq(JB, df = 2, lower.tail = FALSE),
    region = sprintf(r"($JB > \chi^2_{0{,}95;\,2} = %s$)", .fmt_num(valor_critico)),
    decision = .decision(JB > valor_critico),
    n = N, A = A, K = K,
    advertencia = if (N < 20L) {
      sprintf("N = %d < 20: Jarque–Bera es una prueba asintótica y su decisión no es confiable con tan pocos errores.", N)
    } else {
      NULL
    }
  )
}
