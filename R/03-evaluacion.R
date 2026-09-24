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
    formula = r"($JB = \dfrac{N}{6}\left[A^2 + \dfrac{(K-3)^2}{4}\right]\ \overset{H_0}{\sim}\ \chi^2_2$, con $A = m_3/m_2^{3/2}$ la asimetría y $K = m_4/m_2^2$ la curtosis muestrales)",
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
    region <- r"(Se rechaza $H_0$ si $d < d_L$; no se rechaza si $d > d_U$; inconclusa si $d_L \le d \le d_U$, con $d_L$ y $d_U$ de la tabla de Savin y White (1977) para $N$ y $k'$)"
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
      "Se rechaza $H_0$ si $%s < d_L = %s$; no se rechaza si $%s > d_U = %s$; inconclusa si $d_L \\le %s \\le d_U$",
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

# ---- Pruebas auxiliares (internas) ----------------------------------------------------

# Prueba t de media cero sobre los errores: H0: mu_e = 0 contra H1: mu_e != 0 (dos colas).
#   t = ebar / (s_e / sqrt(N)) ~ t_{N-1} bajo H0, con s_e la desviación estándar muestral
#   (divisor N - 1). Región de rechazo al 5 %: |t| > qt(0.975, N - 1).
.t_media_cero <- function(e) {
  stopifnot(
    "e debe ser un vector numérico sin NA" = is.numeric(e) && length(e) > 0L && !anyNA(e),
    "se necesitan al menos 2 errores" = length(e) >= 2L
  )
  N <- length(e)
  s <- stats::sd(e)
  stopifnot("los errores son todos iguales: la prueba t no está definida" = s > 0)
  t <- mean(e) / (s / sqrt(N))
  gl <- N - 1L
  valor_critico <- stats::qt(0.975, df = gl)

  list(
    nombre = "Prueba t de media cero de los errores",
    H0 = r"($H_0:\ \mu_e = 0$)",
    H1 = r"($H_1:\ \mu_e \neq 0$)",
    formula = r"($t = \dfrac{\bar{e}}{s_e/\sqrt{N}}\ \overset{H_0}{\sim}\ t_{N-1}$)",
    estadistico = t,
    gl = gl,
    valor_critico = valor_critico,
    valor_p = 2 * stats::pt(-abs(t), df = gl),
    region = sprintf(r"($|t| > t_{0{,}975;\,%d} = %s$)", gl, .fmt_num(valor_critico)),
    decision = .decision(abs(t) > valor_critico),
    n = N, media = mean(e), desv = s
  )
}

# Contraste individual de la autocorrelación en el rezago h con la banda de ruido blanco
# (Clase 3, «Prueba individual y sus dos bandas»): H0: rho_h = 0 contra H1: rho_h != 0,
#   T* = r_h / desv(r_h), con desv(r_h) = 1 / sqrt(n) bajo ruido blanco, T* ~ N(0, 1)
#   aproximadamente. Rechaza si |T*| > qnorm(0.975), es decir |r_h| > 1,96 / sqrt(n): es la
#   misma banda que dibuja correlograma(). `n` es el de la sucesión de la que sale `r`.
.contraste_rh <- function(r, n, h) {
  stopifnot(
    "r debe ser un vector numérico sin NA" = is.numeric(r) && length(r) > 0L && !anyNA(r),
    "n debe ser un entero positivo" = is.numeric(n) && length(n) == 1L && !is.na(n) && n >= 1 && n == floor(n),
    "h debe ser un entero entre 1 y length(r)" =
      is.numeric(h) && length(h) == 1L && !is.na(h) && h >= 1 && h <= length(r) && h == floor(h)
  )
  z <- r[h] * sqrt(n)
  valor_critico <- stats::qnorm(0.975)

  list(
    nombre = sprintf("Contraste individual de la autocorrelación en el rezago h = %d", h),
    H0 = sprintf(r"($H_0:\ \rho_{%d} = 0$)", h),
    H1 = sprintf(r"($H_1:\ \rho_{%d} \neq 0$)", h),
    formula = r"($T^* = \dfrac{r_h}{\sqrt{1/n}} = \sqrt{n}\, r_h\ \overset{H_0}{\approx}\ N(0,1)$)",
    estadistico = z,
    gl = NA_real_,
    valor_critico = valor_critico,
    valor_p = 2 * stats::pnorm(-abs(z)),
    region = sprintf(r"($|T^*| > z_{0{,}975} = %s$, es decir $|r_{%d}| > %s$)",
                     .fmt_num(valor_critico), h, .fmt_num(valor_critico / sqrt(n))),
    decision = .decision(abs(z) > valor_critico),
    n = n, h = h, r_h = r[h]
  )
}

# ---- .reportar_prueba -----------------------------------------------------------------

# Escribe una prueba de hipótesis en los seis elementos exigidos, en este orden, como
# markdown listo para `cat()` en un chunk con results = "asis":
#   1. H0 y H1 en términos del parámetro.
#   2. Estadístico con su fórmula y su distribución bajo H0, con los grados de libertad.
#   3. Región de rechazo al 5 % con el valor crítico calculado.
#   4. Valor observado y valor p.
#   5. Decisión.
#   6. Lectura en términos de la serie o del método.
# Los elementos 1 a 5 salen del objeto que devuelve la prueba; el sexto, la lectura, SIEMPRE
# se escribe a mano para cada caso y entra como argumento: la función se niega a producir el
# bloque sin ella, para que una prueba no se reporte solo con la salida de una función.
#
# obj:     lista devuelta por ljung_box(), jarque_bera(), durbin_watson(), .t_media_cero()
#          o .contraste_rh() (mismos campos).
# lectura: cadena no vacía, escrita por quien analiza.
.reportar_prueba <- function(obj, lectura) {
  campos <- c("nombre", "H0", "H1", "formula", "estadistico", "gl", "valor_p", "region", "decision")
  stopifnot(
    "obj no es el resultado de una prueba: faltan campos" =
      is.list(obj) && all(campos %in% names(obj)),
    "lectura debe ser una cadena no vacía escrita a mano para este caso" =
      is.character(lectura) && length(lectura) == 1L && !is.na(lectura) && nzchar(trimws(lectura))
  )

  titulo <- if (is.null(obj$n)) {
    sprintf("**%s**", obj$nombre)
  } else {
    sprintf("**%s** (n = %d)", obj$nombre, as.integer(obj$n))
  }

  gl_txt <- if (is.na(obj$gl)) {
    ""
  } else if (!is.null(obj$m) && !is.null(obj$p)) {
    sprintf("; $m = %d$, $p = %d$, luego los grados de libertad son $m - p = %d$",
            as.integer(obj$m), as.integer(obj$p), as.integer(obj$gl))
  } else {
    sprintf("; grados de libertad: $%d$", as.integer(obj$gl))
  }

  p_txt <- if (is.na(obj$valor_p)) {
    "sin valor p exacto (se decide con las cotas tabuladas)"
  } else if (obj$valor_p < 0.001) {
    "valor p < 0,001"
  } else {
    paste0("valor p = ", .fmt_p(obj$valor_p))
  }

  decision <- paste0(obj$decision, " (nivel de significancia del 5 %).")
  if (!is.null(obj$advertencia)) {
    decision <- paste0(decision, " *Advertencia: ", obj$advertencia, "*")
  }

  paste(
    titulo,
    "",
    paste0("1. **Hipótesis.** ", obj$H0, "; ", obj$H1, "."),
    paste0("2. **Estadístico y distribución bajo $H_0$.** ", obj$formula, gl_txt, "."),
    paste0("3. **Región de rechazo al 5 %.** ", obj$region, "."),
    paste0("4. **Valor observado.** estadístico observado = ", .fmt_num(obj$estadistico),
           "; ", p_txt, "."),
    paste0("5. **Decisión.** ", decision),
    paste0("6. **Lectura.** ", trimws(lectura)),
    "",
    sep = "\n"
  )
}

# ---- validar_errores ------------------------------------------------------------------

#' validar_errores(e, p, T_serie = NULL, dL = NULL, dU = NULL)
#'
#' Descripción: validación de los errores de un paso de un método (Clase 3, Parte III:
#'   los errores deben ser ruido blanco). Sobre los errores sin NA calcula y devuelve:
#'   - gráfico de los errores en el tiempo, con la línea en cero;
#'   - correlograma (ACF a mano, PACF) con la banda calculada sobre N, el número de errores;
#'   - prueba t de media cero (t = ebar / (s_e / sqrt(N)), t_{N-1});
#'   - Ljung–Box con m - p grados de libertad, m = min(floor(N/4), 24) y T = N;
#'   - Jarque–Bera (con advertencia si N < 20);
#'   - Durbin–Watson (por regiones si se dan dL y dU).
#'   N se declara explícitamente en la salida: es el número de errores, no el T de la serie.
#'
#' Los NA solo pueden estar al inicio (calentamiento del método); se descartan y las
#'   pruebas se hacen sobre los N errores restantes.
#'
#' @param e   errores de un paso e_t = y_t - yhat_t (vector; NA en el calentamiento).
#' @param p   número de parámetros que el método estimó (convención del README): fija los
#'            grados de libertad m - p de Ljung–Box. Debe ser menor que m.
#' @param T_serie  tamaño de la serie, opcional y solo informativo: se devuelve junto a N
#'            para declarar «N errores de una serie de T observaciones».
#' @param dL,dU  cotas de Durbin–Watson (ver durbin_watson()); opcionales.
#' @return lista con n, T_serie, e (errores sin NA), correlograma (lista de correlograma()),
#'   t_media, ljung_box, jarque_bera, durbin_watson, grafico_errores y grafico (errores y
#'   correlograma en un solo panel). El gráfico solo se imprime si la sesión es interactiva;
#'   guardarlo en disco es tarea de ejemplos.R.
#'
#' Referencia: enunciado, sección 4(a); Clase 3, Parte III.
validar_errores <- function(e, p, T_serie = NULL, dL = NULL, dU = NULL) {
  stopifnot(
    "e debe ser un vector numérico" = is.numeric(e) && length(e) > 0L,
    "p debe ser un entero no negativo" =
      is.numeric(p) && length(p) == 1L && !is.na(p) && p >= 0 && p == floor(p),
    "T_serie debe ser NULL o un entero positivo" =
      is.null(T_serie) || (is.numeric(T_serie) && length(T_serie) == 1L && !is.na(T_serie) &&
                             T_serie >= 1 && T_serie == floor(T_serie))
  )

  # correlograma() valida que los NA sean solo de calentamiento y descarta esos NA:
  # su n es el número de errores, N, sobre el que se calcula todo lo demás.
  cg <- correlograma(e)
  N <- cg$n
  stopifnot("T_serie no puede ser menor que el número de errores" = is.null(T_serie) || T_serie >= N)
  if (cg$m <= p) {
    stop(sprintf(paste0("con N = %d errores, m = min(floor(N/4), 24) = %d no supera p = %d: ",
                        "Ljung-Box no tiene grados de libertad (m - p > 0)"),
                 N, cg$m, as.integer(p)), call. = FALSE)
  }
  tiempo <- (length(e) - N + 1L):length(e)
  e_ok <- e[tiempo]

  lb <- ljung_box(cg$acf, N, cg$m, p)
  lb$nombre <- "Prueba de Ljung–Box sobre los errores"

  g_err <- ggplot2::ggplot(data.frame(t = tiempo, e = e_ok), ggplot2::aes(x = t, y = e)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_line(colour = "#1F4E79", linewidth = 0.5) +
    ggplot2::geom_point(colour = "#1F4E79", size = 1.3) +
    ggplot2::labs(title = "Errores de un paso", x = "t (índice de la observación)",
                  y = "e_t = y_t - yhat_t") +
    .tema_tarea()

  # Al anidar el panel de correlograma() se pierde su subtítulo, así que se vuelve a escribir
  # aquí: N (y T si se dio) y m quedan declarados en la propia figura.
  grafico <- patchwork::wrap_plots(g_err, cg$grafico, ncol = 1, heights = c(1, 2.2)) +
    patchwork::plot_annotation(
      title = "Validación de los errores de un paso",
      subtitle = sprintf(
        "N = %d errores%s; m = %d rezagos.\nLíneas discontinuas: banda de ruido blanco al 95 %% (±%s)",
        N, if (is.null(T_serie)) "" else sprintf(" (serie de T = %d observaciones)", as.integer(T_serie)),
        cg$m, formatC(cg$banda, format = "f", digits = 4, decimal.mark = ",")
      ),
      theme = .tema_tarea()
    )
  if (interactive()) {
    print(grafico)
  }

  list(
    n = N, T_serie = T_serie, e = e_ok, correlograma = cg,
    t_media = .t_media_cero(e_ok),
    ljung_box = lb,
    jarque_bera = jarque_bera(e_ok),
    durbin_watson = durbin_watson(e_ok, dL, dU),
    grafico_errores = g_err, grafico = grafico
  )
}

# ---- Auxiliares de los ejemplos y del informe ------------------------------------------

# medidas() con el aviso de MAPE indefinido explicado: cuando la serie tiene ceros (es el
# caso de discoveries) medidas() devuelve MAPE = NA y avisa; aquí ese aviso, y solo ese, se
# convierte en un mensaje. Cualquier otra advertencia sigue saliendo.
.medir <- function(y, yhat, y_entrenamiento, s = 1) {
  withCallingHandlers(
    medidas(y, yhat, y_entrenamiento, s),
    warning = function(w) {
      if (grepl("MAPE no definido", conditionMessage(w), fixed = TRUE)) {
        message("  [esperado] MAPE no definido: hay ceros entre los datos evaluados; se reporta NA.")
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Tabla de medidas de un ejemplo: método dentro de la muestra, método en validación y
# referente en validación (mismo tramo y mismo h).
.tabla_medidas <- function(m_in, m_val, m_ref, nombre_ref) {
  as.data.frame(rbind(
    cbind(caso = "Método, un paso (estimación)", m_in),
    cbind(caso = "Método, h pronósticos (validación)", m_val),
    cbind(caso = paste0(nombre_ref, ", h pronósticos (validación)"), m_ref)
  ))
}

# Cotas d_L y d_U de Durbin–Watson al 5 % (una cola) para los N errores de cada ejemplo y
# k' regresores sin contar el intercepto. Fuente: tabla de Savin y White (1977). La tabla no
# cubre todos los n (después de 40 solo trae múltiplos de 5), así que se calcularon de forma
# exacta con la definición de Durbin y Watson (1951) y la integral de Imhof (1961); ese
# cálculo reproduce 12 filas de la tabla con diferencias de a lo sumo 0,005 (ver README).
# Se redondean a 2 decimales, como la tabla. Para los métodos de suavizamiento se usa
# k' = 1 (aproximación declarada: no hay matriz de diseño).
# La tabla vive dentro de la función (ninguna función lee objetos globales). Devuelve la
# fila del ejemplo `id` (N, k, dL, dU) y se detiene si el número de errores cambió: las
# cotas solo valen para ese N.
.cotas <- function(id, N) {
  cotas_dw <- list(
    ej01 = list(N = 87, k = 1, dL = 1.63, dU = 1.67),
    ej02 = list(N = 86, k = 1, dL = 1.63, dU = 1.67),
    ej03 = list(N = 87, k = 1, dL = 1.63, dU = 1.67),
    ej04 = list(N = 74, k = 1, dL = 1.60, dU = 1.65),
    ej05 = list(N = 86, k = 1, dL = 1.63, dU = 1.67),
    ej06 = list(N = 20, k = 2, dL = 1.10, dU = 1.54),
    ej07 = list(N = 72, k = 1, dL = 1.59, dU = 1.65),
    ej08 = list(N = 76, k = 1, dL = 1.60, dU = 1.65),
    ej09 = list(N = 131, k = 1, dL = 1.70, dU = 1.73)
  )
  fila <- cotas_dw[[id]]
  stopifnot("no hay cotas de Durbin-Watson para este ejemplo" = !is.null(fila),
            "el número de errores cambió: recalcular las cotas de Durbin-Watson" = N == fila$N)
  fila
}


# Prueba t de un coeficiente de una tendencia con error estándar robusto (HAC), en el mismo
# formato que las demás pruebas para que .reportar_prueba() la escriba en seis elementos.
#   simbolo     el coeficiente en LaTeX sin signos de dólar, p. ej. r"(\beta_1)", "a" o r"(\theta)"
#   estimacion  estimación del coeficiente
#   ee          error estándar robusto (columna ee_robusto de la tabla de ajustar_tendencia())
#   gl          grados de libertad, T - p (parametros$gl)
#   n           tamaño de la muestra de la regresión, T (se muestra en el título de la prueba)
# Bajo H0 el estadístico se refiere a la t de Student con T - p grados de libertad (elección
# documentada en ajustar_tendencia(): con T pequeño es más conservadora que la normal). La
# aproximación es asintótica: con errores autocorrelacionados el HAC corrige la varianza pero
# no vuelve exacta la distribución.
.t_coeficiente <- function(simbolo, estimacion, ee, gl, n) {
  stopifnot(
    "simbolo debe ser una cadena" = is.character(simbolo) && length(simbolo) == 1L,
    "estimacion y ee deben ser números" = is.numeric(estimacion) && is.numeric(ee) &&
      length(estimacion) == 1L && length(ee) == 1L && !anyNA(c(estimacion, ee)),
    "ee debe ser positivo" = ee > 0,
    "gl debe ser un entero positivo" = is.numeric(gl) && length(gl) == 1L && gl >= 1 && gl == floor(gl),
    "n debe ser un entero mayor que gl" = is.numeric(n) && length(n) == 1L && n > gl && n == floor(n)
  )
  t <- estimacion / ee
  valor_critico <- stats::qt(0.975, df = gl)
  list(
    nombre = sprintf("Prueba t del coeficiente $%s$ con error estándar robusto", simbolo),
    H0 = sprintf(r"($H_0:\ %s = 0$)", simbolo),
    H1 = sprintf(r"($H_1:\ %s \neq 0$)", simbolo),
    formula = sprintf(r"($t = \dfrac{\hat{%s}}{\widehat{EE}_{HAC}(\hat{%s})}\ \overset{H_0}{\approx}\ t_{T-p}$)", simbolo, simbolo),
    estadistico = t,
    gl = gl,
    valor_critico = valor_critico,
    valor_p = 2 * stats::pt(-abs(t), df = gl),
    region = sprintf(r"($|t| > t_{0{,}975;\,%d} = %s$)", as.integer(gl), .fmt_num(valor_critico)),
    decision = .decision(abs(t) > valor_critico),
    n = n
  )
}
