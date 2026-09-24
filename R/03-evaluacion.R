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
#'   n que usa Box.test() (la longitud de x).
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
