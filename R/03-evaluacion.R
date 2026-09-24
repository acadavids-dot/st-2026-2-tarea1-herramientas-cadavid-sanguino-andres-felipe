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

# ---- durbin_watson --------------------------------------------------------------------

#' durbin_watson(e, dL = NULL, dU = NULL)
#'
#' Descripción: estadístico de Durbin–Watson sobre los errores (o residuos) e_1, ..., e_N
#'   y decisión por regiones con las cotas de la tabla.
#'
#' Ecuaciones:
#'   d = sum_{t=2}^{N} (e_t - e_{t-1})^2 / sum_{t=1}^{N} e_t^2,   d ~ 2 (1 - r_1).
#'   d cerca de 2: sin autocorrelación de orden 1; cerca de 0: positiva; cerca de 4: negativa.
#'   El denominador es la suma de cuadrados sin centrar, como en la definición del
#'   estadístico (en residuos de una regresión con intercepto la media ya es cero).
#'
#' Decisión: d no tiene una distribución exacta libre de la matriz de diseño, así que no
#'   hay valor p; se compara con las cotas d_L y d_U de la tabla de Savin y White (1977) para
#'   N y k' (número de regresores sin contar el intercepto), al 5 % de una cola. Las cotas
#'   NO están en R base y no se aproximan: se dan como argumentos.
#'   - Si d < 2 se contrasta H1: rho > 0 con d: se rechaza H0 si d < d_L, no se rechaza si
#'     d > d_U, y la prueba es inconclusa si d_L <= d <= d_U.
#'   - Si d >= 2 se contrasta H1: rho < 0 con 4 - d y las mismas cotas.
#'   Sin dL y dU la función devuelve d y la decisión queda «No evaluada».
#'
#' @param e   vector numérico sin NA.
#' @param dL,dU  cotas inferior y superior (0 < dL < dU < 2), ambas o ninguna.
#' @return lista con los campos comunes de las pruebas (valor_critico = c(dL, dU) y
#'   valor_p = NA) y además n, lado ("positiva" o "negativa": la autocorrelación contrastada)
#'   y rho_aprox = 1 - d/2 (la lectura d ~ 2(1 - r_1)).
#'
#' Referencia: Clase 4, Parte VI (DW = 0,7063 en la tendencia lineal); Savin y White (1977).
durbin_watson <- function(e, dL = NULL, dU = NULL) {
  stopifnot(
    "e debe ser un vector numérico sin NA" = is.numeric(e) && length(e) > 0L && !anyNA(e),
    "se necesitan al menos 3 errores" = length(e) >= 3L,
    "los errores no pueden ser todos cero" = sum(e^2) > 0,
    "dL y dU se dan juntos o ninguno" = is.null(dL) == is.null(dU)
  )
  if (!is.null(dL)) {
    stopifnot("dL y dU deben cumplir 0 < dL < dU < 2" =
                is.numeric(dL) && is.numeric(dU) && length(dL) == 1L && length(dU) == 1L &&
                !anyNA(c(dL, dU)) && dL > 0 && dL < dU && dU < 2)
  }

  N <- length(e)
  d <- sum(diff(e)^2) / sum(e^2)
  lado <- if (d < 2) "positiva" else "negativa"
  d_eval <- if (d < 2) d else 4 - d  # con d >= 2 se contrasta la cola negativa con 4 - d
  simbolo <- if (d < 2) "d" else "4-d"

  if (is.null(dL)) {
    decision <- "No evaluada: faltan las cotas $d_L$ y $d_U$ de la tabla"
    region <- r"(Se rechaza $H_0$ si $d < d_L$; no se rechaza si $d > d_U$; inconclusa si $d_L \le d \le d_U$, con $d_L$ y $d_U$ de la tabla de Savin y White (1977) para $N$ y $k'$.)"
    valor_critico <- c(dL = NA_real_, dU = NA_real_)
  } else {
    decision <- if (d_eval < dL) {
      "Se rechaza $H_0$"
    } else if (d_eval > dU) {
      "No se rechaza $H_0$"
    } else {
      "Prueba inconclusa ($d_L \\le d \\le d_U$)"
    }
    region <- sprintf(
      "Se rechaza $H_0$ si $%s < d_L = %s$; no se rechaza si $%s > d_U = %s$; inconclusa si $d_L \\le %s \\le d_U$.",
      simbolo, .fmt_num(dL), simbolo, .fmt_num(dU), simbolo
    )
    valor_critico <- c(dL = dL, dU = dU)
  }

  list(
    nombre = "Prueba de Durbin–Watson",
    H0 = r"($H_0:\ \rho = 0$ (sin autocorrelación de orden 1))",
    H1 = if (lado == "positiva") r"($H_1:\ \rho > 0$)" else r"($H_1:\ \rho < 0$)",
    formula = r"($d = \dfrac{\sum_{t=2}^{N}(e_t - e_{t-1})^2}{\sum_{t=1}^{N} e_t^2} \approx 2(1 - r_1)$; sin distribución exacta bajo $H_0$, se decide con las cotas $d_L$ y $d_U$)",
    estadistico = d,
    gl = NA_real_,
    valor_critico = valor_critico,
    valor_p = NA_real_,
    region = region,
    decision = decision,
    n = N, lado = lado, rho_aprox = 1 - d / 2
  )
}

# ---- medidas --------------------------------------------------------------------------

#' medidas(y, yhat, y_entrenamiento, s = 1)
#'
#' Descripción: MSE, MAD, MAPE y MASE de un pronóstico, sobre los pares (y_t, yhat_t) sin NA.
#'
#' Ecuaciones (e_t = y_t - yhat_t, N = número de pares sin NA):
#'   MSE  = (1/N) sum e_t^2
#'   MAD  = (1/N) sum |e_t|
#'   MAPE = (100/N) sum |e_t / y_t|                  (en %)
#'   MASE = MAD / [ (1/(T-s)) sum_{t=s+1}^{T} |Y_t - Y_{t-s}| ]
#'   donde el denominador se calcula sobre `y_entrenamiento` (el tramo de estimación, con
#'   T = su longitud), no sobre `y`: es el MAD dentro de muestra del pronóstico ingenuo
#'   (s = 1) o ingenuo estacional (s = frecuencia). MASE < 1 significa que el método gana
#'   al ingenuo (Clase 3, Parte III; Hyndman y Koehler, 2006).
#'
#' Tratamiento de NA: se descartan los pares con yhat = NA (el calentamiento del método) y N
#'   cuenta los que quedan. y debe venir completo: un NA en y es un error, no calentamiento.
#'
#' MAPE con ceros: si algún y_t de los pares usados vale 0, el MAPE no está definido; se
#'   devuelve NA con una advertencia (el dato no se omite en silencio).
#'
#' @param y   observaciones del tramo evaluado (vector numérico sin NA).
#' @param yhat  pronósticos de esas observaciones (misma longitud; NA en el calentamiento).
#' @param y_entrenamiento  serie del tramo de estimación, para la escala del MASE.
#' @param s   período del ingenuo de referencia: 1 (series no estacionales) o la frecuencia.
#' @return tibble de una fila con N, MSE, MAD, MAPE y MASE; el atributo `escala_mase` guarda
#'   el denominador del MASE.
#'
#' Referencia: Clase 3, Parte III (tres medidas usuales y la que falta).
medidas <- function(y, yhat, y_entrenamiento, s = 1) {
  stopifnot(
    "y debe ser un vector numérico sin NA" = is.numeric(y) && length(y) > 0L && !anyNA(y),
    "yhat debe ser numérico y de la misma longitud que y" =
      is.numeric(yhat) && length(yhat) == length(y),
    "y_entrenamiento debe ser un vector numérico sin NA" =
      is.numeric(y_entrenamiento) && length(y_entrenamiento) > 0L && !anyNA(y_entrenamiento),
    "s debe ser un entero positivo" =
      is.numeric(s) && length(s) == 1L && !is.na(s) && s >= 1 && s == floor(s)
  )
  stopifnot("y_entrenamiento debe tener más de s observaciones" = length(y_entrenamiento) > s)

  usar <- !is.na(yhat)
  N <- sum(usar)
  stopifnot("no hay ningún par (y, yhat) sin NA para evaluar" = N >= 1L)
  y_u <- y[usar]
  e <- y_u - yhat[usar]

  if (any(y_u == 0)) {
    warning("MAPE no definido: ", sum(y_u == 0), " de los ", N,
            " valores observados son 0; se devuelve NA", call. = FALSE)
    mape <- NA_real_
  } else {
    mape <- 100 * mean(abs(e / y_u))
  }

  mad <- mean(abs(e))
  escala <- mean(abs(diff(y_entrenamiento, lag = s)))
  stopifnot("la escala del MASE es cero: y_entrenamiento es constante" = escala > 0)

  res <- tibble::tibble(N = as.integer(N), MSE = mean(e^2), MAD = mad, MAPE = mape,
                        MASE = mad / escala)
  attr(res, "escala_mase") <- escala
  res
}
