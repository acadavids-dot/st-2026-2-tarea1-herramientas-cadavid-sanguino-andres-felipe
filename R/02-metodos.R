# 02-metodos.R
# Métodos de pronóstico: ajustar_media(), ajustar_mm(), ajustar_ses(), ajustar_dmm(),
# ajustar_tendencia(), ajustar_holt() y optimizar()
#
# Contrato común (enunciado, sección 2): cada ajustar_*() recibe un vector numérico y
# devuelve una lista de clase "metodo_pronostico" con
#   yhat        vector de largo length(y); yhat[t] pronostica Y_t con información hasta t - 1
#               (calentamiento en NA)
#   pronosticar function(h) con los pronósticos Yhat_{T+1}, ..., Yhat_{T+h}
#   parametros  lista con los parámetros usados y los estados finales
# más metodo (nombre), errores (y - yhat) y n_param (el p de Ljung–Box; ver el README).

# ---- Validaciones y constructores internos --------------------------------------------

# y: numérico, sin NA ni infinitos y con al menos dos datos. Cada método usa lo que
# necesita de esto y añade sus propias restricciones (ventana, constantes, positividad).
# Valida el objeto ORIGINAL (si se convirtiera antes, un vector de texto pasaría a NA o,
# peor, c("1", "2") pasaría como número) y devuelve y como vector numérico simple, sin
# los atributos de ts / tibble: los métodos solo usan los valores.
.validar_serie <- function(y) {
  stopifnot(
    "y debe ser un vector numérico" = is.numeric(y) && is.null(dim(y)),
    "y no puede contener NA" = !anyNA(y),
    "y no puede contener valores infinitos" = all(is.finite(y)),
    "y debe tener al menos 2 observaciones" = length(y) >= 2L
  )
  as.numeric(y)
}

# Lógicos auxiliares para los stopifnot() nombrados de cada función: así el mensaje de
# error dice exactamente qué argumento falló.
.es_entero <- function(x, minimo) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) && x >= minimo
}
.es_constante <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x < 1
}

# Horizonte h: entero positivo.
.validar_h <- function(h) {
  stopifnot("h debe ser un entero positivo" = .es_entero(h, 1))
}

# Pronosticador de los métodos de nivel: repite un valor h veces. La fábrica existe para
# que el closure capture solo `valor` (y no todo el entorno de ajustar_*, con y, yhat...).
.pronosticador_constante <- function(valor) {
  force(valor)
  function(h) {
    .validar_h(h)
    rep(valor, h)
  }
}

# Pronosticador de los métodos con pendiente: nivel + pendiente * (1, ..., h).
.pronosticador_lineal <- function(nivel, pendiente) {
  force(nivel)
  force(pendiente)
  function(h) {
    .validar_h(h)
    nivel + pendiente * seq_len(h)
  }
}

# Arma el objeto que devuelven los ocho métodos. errores = y - yhat conserva los NA del
# calentamiento; los consumidores (medidas, validar_errores) los descartan.
.nuevo_metodo <- function(metodo, y, yhat, pronosticar, parametros, n_param) {
  structure(
    list(metodo = metodo, yhat = yhat, errores = y - yhat, pronosticar = pronosticar,
         parametros = parametros, n_param = n_param),
    class = "metodo_pronostico"
  )
}

# ---- ajustar_media --------------------------------------------------------------------

#' ajustar_media(y)
#'
#' Descripción: pronóstico con la media de todos los datos observados hasta t (media
#'   simple; nivel fijo). Es el método para una serie con nivel aproximadamente constante.
#'
#' Ecuaciones:
#'   Yhat_{t+1} = Ybar_t = (1/t) sum_{i<=t} Y_i, actualizada de forma recursiva:
#'   Ybar_t = Ybar_{t-1} + (Y_t - Ybar_{t-1}) / t.
#'   La recursión evita recalcular mean(y[1:t]) en cada paso (T sumas de largo t, O(T^2)),
#'   y da el mismo número salvo por redondeo de punto flotante.
#'
#' Inicialización: Ybar_1 = Y_1.
#' Calentamiento: yhat[1] = NA (no hay información previa); yhat[2] = Y_1.
#' Pronóstico extramuestral: rep(Ybar_T, h); el nivel es fijo, no hay pendiente.
#'
#' @param y  vector numérico ordenado en el tiempo, sin NA (al menos 2 datos).
#' @return objeto de clase "metodo_pronostico" (ver cabecera del archivo). En parametros:
#'   media (Ybar_T) y trayectoria (Ybar_1, ..., Ybar_T). n_param = 0: nada se elige sobre
#'   la muestra, así que Ljung–Box sobre los errores usa chi^2_m.
#'
#' Referencia: Clase 3, media simple.
ajustar_media <- function(y) {
  y <- .validar_serie(y)

  n <- length(y)
  yhat <- rep(NA_real_, n)
  trayectoria <- numeric(n)
  media <- y[1]
  trayectoria[1] <- media
  for (t in 2:n) {
    yhat[t] <- media                       # Yhat_t = Ybar_{t-1}: solo usa datos hasta t - 1
    media <- media + (y[t] - media) / t    # Ybar_t con el dato nuevo
    trayectoria[t] <- media
  }

  .nuevo_metodo(
    metodo = "Media simple", y = y, yhat = yhat,
    pronosticar = .pronosticador_constante(media),
    parametros = list(media = media, trayectoria = trayectoria),
    n_param = 0L
  )
}

# ---- ajustar_mm -----------------------------------------------------------------------

#' ajustar_mm(y, k)
#'
#' Descripción: pronóstico con la media móvil de orden k, el promedio de los k datos más
#'   recientes (nivel local). Es el método para una serie con nivel que cambia despacio.
#'
#' Ecuaciones:
#'   MM_t(k) = (1/k) sum_{i=0}^{k-1} Y_{t-i},   Yhat_{t+1} = MM_t(k)  para t >= k.
#'   Un solo recorrido con la forma recursiva (Proposición de la Clase 3):
#'   MM_t = MM_{t-1} + (Y_t - Y_{t-k}) / k. La primera media se calcula una sola vez.
#'
#' Inicialización: MM_k = (Y_1 + ... + Y_k) / k.
#' Calentamiento: yhat[1:k] = NA; el primer pronóstico es yhat[k + 1] = MM_k.
#' Pronóstico extramuestral: rep(MM_T, h).
#'
#' @param y  vector numérico ordenado, sin NA.
#' @param k  ventana: entero con 2 <= k <= length(y).
#' @return objeto de clase "metodo_pronostico". En parametros: k, mm (trayectoria MM_t,
#'   con NA para t < k) y ultima (MM_T). n_param = 1 (se estima k).
#'
#' Referencia: Clase 3, medias móviles y su forma recursiva.
ajustar_mm <- function(y, k) {
  y <- .validar_serie(y)
  stopifnot("k debe ser un entero mayor o igual a 2" = .es_entero(k, 2))
  n <- length(y)
  stopifnot("k no puede superar length(y)" = k <= n)

  mm <- rep(NA_real_, n)
  mm[k] <- sum(y[1:k]) / k
  for (t in seq_len(n - k) + k) {          # vacío si k == n
    mm[t] <- mm[t - 1L] + (y[t] - y[t - k]) / k
  }
  yhat <- c(NA_real_, mm[-n])              # yhat[t] = MM_{t-1}: NA para t <= k

  .nuevo_metodo(
    metodo = sprintf("Media móvil (k = %d)", as.integer(k)), y = y, yhat = yhat,
    pronosticar = .pronosticador_constante(mm[n]),
    parametros = list(k = k, mm = mm, ultima = mm[n]),
    n_param = 1L
  )
}

# ---- ajustar_ses ----------------------------------------------------------------------

#' ajustar_ses(y, alpha)
#'
#' Descripción: suavizamiento exponencial simple: el pronóstico corrige el anterior con una
#'   fracción alpha del último error (nivel local con pesos geométricos decrecientes).
#'
#' Ecuaciones (forma de corrección de error):
#'   e_t = Y_t - Yhat_t,   Yhat_{t+1} = Yhat_t + alpha e_t.
#'   Equivale al promedio ponderado finito de esta inicialización:
#'   Yhat_{t+1} = alpha sum_{i=0}^{t-2} (1 - alpha)^i Y_{t-i} + (1 - alpha)^{t-1} Y_1,
#'   con pesos que suman 1.
#'
#' Inicialización: Yhat_2 = Y_1.
#' Calentamiento: yhat[1] = NA; yhat[2] = Y_1.
#' Pronóstico extramuestral: rep(Yhat_{T+1}, h), con Yhat_{T+1} = alpha Y_T + (1 - alpha) Yhat_T.
#'
#' @param y      vector numérico ordenado, sin NA.
#' @param alpha  constante de suavizamiento, estrictamente entre 0 y 1.
#' @return objeto de clase "metodo_pronostico". En parametros: alpha y siguiente
#'   (Yhat_{T+1}). n_param = 1 (se estima alpha).
#'
#' Referencia: Clase 3, suavizamiento exponencial simple.
ajustar_ses <- function(y, alpha) {
  y <- .validar_serie(y)
  stopifnot("alpha debe ser un número en el intervalo (0, 1)" = .es_constante(alpha))

  n <- length(y)
  yhat <- rep(NA_real_, n)
  pron <- y[1]                             # Yhat_2 = Y_1
  yhat[2] <- pron
  for (t in 2:n) {
    pron <- pron + alpha * (y[t] - pron)   # Yhat_{t+1} = Yhat_t + alpha e_t
    if (t < n) yhat[t + 1L] <- pron
  }                                        # al salir, pron = Yhat_{T+1}

  .nuevo_metodo(
    metodo = sprintf("SES (alpha = %s)", format(alpha)), y = y, yhat = yhat,
    pronosticar = .pronosticador_constante(pron),
    parametros = list(alpha = alpha, siguiente = pron),
    n_param = 1L
  )
}

# ---- ajustar_dmm ----------------------------------------------------------------------

#' ajustar_dmm(y, k)
#'
#' Descripción: doble media móvil de orden k; corrige el rezago de la media móvil ante una
#'   tendencia lineal local combinando MM y la media de las últimas k medias móviles.
#'
#' Ecuaciones (Definición 3 de la Clase 4):
#'   MM_t(k)  = (1/k) sum_{i=0}^{k-1} Y_{t-i}                para t >= k,
#'   DMM_t(k) = (1/k) sum_{i=0}^{k-1} MM_{t-i}(k)            para t >= 2k - 1,
#'   Ehat_t = 2 MM_t - DMM_t,   beta1hat(t) = 2/(k - 1) (MM_t - DMM_t),
#'   Yhat_{t+1} = Ehat_t + beta1hat(t)                        para t >= 2k - 1.
#'   Con tendencia lineal MM rezaga (k-1)/2 y DMM rezaga k-1 (Proposición 9); 2 MM - DMM
#'   recupera el nivel y la diferencia, la pendiente.
#'   Un solo recorrido con sumas corridas: MM_t = MM_{t-1} + (Y_t - Y_{t-k})/k y
#'   DMM_t = DMM_{t-1} + (MM_t - MM_{t-k})/k.
#'
#' Inicialización: MM_k = mean(Y_1..Y_k); DMM_{2k-1} = mean(MM_k..MM_{2k-1}).
#' Calentamiento: yhat[1:(2k - 1)] = NA; el primer pronóstico es yhat[2k].
#' Pronóstico extramuestral: Yhat_{T+h} = Ehat_T + beta1hat(T) h.
#'
#' @param y  vector numérico ordenado, sin NA.
#' @param k  ventana: entero con k >= 2 (con k = 1 el factor 2/(k-1) no existe) y
#'           2k - 1 <= length(y).
#' @return objeto de clase "metodo_pronostico". En parametros: k, mm, dmm, nivel (Ehat_t) y
#'   pendiente (beta1hat(t)), todos de largo T con NA donde no están definidos.
#'   n_param = 1 (se estima k).
#'
#' Referencia: Clase 4, Parte VIII (dobles medias móviles).
ajustar_dmm <- function(y, k) {
  y <- .validar_serie(y)
  stopifnot("k debe ser un entero mayor o igual a 2" = .es_entero(k, 2))
  n <- length(y)
  stopifnot("2k - 1 no puede superar length(y)" = 2 * k - 1 <= n)

  mm <- rep(NA_real_, n)
  dmm <- rep(NA_real_, n)
  for (t in k:n) {                         # un solo recorrido para MM y DMM
    mm[t] <- if (t == k) sum(y[1:k]) / k else mm[t - 1L] + (y[t] - y[t - k]) / k
    if (t == 2 * k - 1) {
      dmm[t] <- sum(mm[k:t]) / k           # primera DMM: se calcula una sola vez
    } else if (t > 2 * k - 1) {
      dmm[t] <- dmm[t - 1L] + (mm[t] - mm[t - k]) / k
    }
  }
  nivel <- 2 * mm - dmm
  pendiente <- 2 / (k - 1) * (mm - dmm)
  yhat <- c(NA_real_, (nivel + pendiente)[-n])  # yhat[t] = Ehat_{t-1} + beta1hat(t-1)

  .nuevo_metodo(
    metodo = sprintf("Doble media móvil (k = %d)", as.integer(k)), y = y, yhat = yhat,
    pronosticar = .pronosticador_lineal(nivel[n], pendiente[n]),
    parametros = list(k = k, mm = mm, dmm = dmm, nivel = nivel, pendiente = pendiente),
    n_param = 1L
  )
}

# ---- Tendencias por mínimos cuadrados -------------------------------------------------

# Matriz de diseño a mano para t = (t_1, ...): [1, t] en la lineal y en la exponencial
# (que es lineal en logaritmos) y [1, t, t^2] en la cuadrática.
.matriz_diseno <- function(t, tipo) {
  if (tipo == "cuadratica") cbind(1, t, t^2) else cbind(1, t)
}

# Lleva un valor de la escala de la regresión (lineal o logarítmica) a la escala de Y:
# identidad en las lineales; exp() por el factor de retransformación en la exponencial.
.a_escala_y <- function(valor, tipo, factor) {
  if (tipo == "exponencial") exp(valor) * factor else valor
}

# Pronosticador de las tendencias: evalúa la tendencia en t = n + 1, ..., n + h. Fábrica
# aparte para que el closure capture solo beta, n, tipo y factor.
.pronosticador_tendencia <- function(beta, n, tipo, factor) {
  force(beta); force(n); force(tipo); force(factor)
  function(h) {
    .validar_h(h)
    .a_escala_y(drop(.matriz_diseno(n + seq_len(h), tipo) %*% beta), tipo, factor)
  }
}

# Varianza robusta de Newey–West (HAC) con núcleo de Bartlett, escrita a mano:
#   Gamma_l = sum_{t=l+1}^{T} u_t u_{t-l} x_t x_{t-l}',  l = 0, ..., L
#   S = Gamma_0 + sum_{l=1}^{L} (1 - l/(L + 1)) (Gamma_l + Gamma_l')
#   V = (X'X)^{-1} S (X'X)^{-1}
# Sin corrección por grados de libertad (no se multiplica por T/(T - p)). Cada fila de
# XU es u_t x_t, así que Gamma_l = XU[(l+1):T, ]' XU[1:(T-l), ] sin bucles sobre t.
# Los pesos de Bartlett 1 - l/(L + 1) mantienen S semidefinida positiva.
.varianza_hac <- function(X, u, L) {
  n <- nrow(X)
  XU <- X * u                              # u se recicla por columnas: fila t = u_t x_t
  S <- crossprod(XU)                       # Gamma_0
  for (l in seq_len(L)) {
    G <- crossprod(XU[(l + 1L):n, , drop = FALSE], XU[1:(n - l), , drop = FALSE])
    S <- S + (1 - l / (L + 1)) * (G + t(G))
  }
  XtXi <- solve(crossprod(X))
  XtXi %*% S %*% XtXi
}

#' ajustar_tendencia(y, tipo = c("lineal", "cuadratica", "exponencial"), corregir_sesgo = FALSE)
#'
#' Descripción: tendencia determinista en t = 1, ..., T estimada por mínimos cuadrados con
#'   las ecuaciones normales. Es el método para una serie con tendencia global estable.
#'
#' Ecuaciones:
#'   lineal:       Y_t = b0 + b1 t + e_t
#'   cuadrática:   Y_t = b0 + b1 t + b2 t^2 + e_t
#'   exponencial:  ln Y_t = a + theta t + e_t,  con b0 = e^a y b1 = e^theta (Y = b0 b1^t).
#'   bhat = solve(crossprod(X), crossprod(X, z)), con X construida a mano y z = y
#'   (o ln y en la exponencial). No se usa la regresión de stats.
#'
#' Exponencial: se detiene si algún Y_t <= 0. Ojo: e^{ahat + thetahat t} estima la MEDIANA
#'   condicional de Y_t, no la media (Proposición 4 de la Clase 4). Con
#'   corregir_sesgo = TRUE se multiplica por e^{sigma2hat_ln / 2}, con sigma2hat_ln =
#'   SCR/(T - 2) de la regresión en logaritmos (corrección lognormal: supone errores
#'   normales; la de Duan, que no lo supone, no se implementa). Solo aplica a "exponencial".
#'
#' Calentamiento: ninguno. yhat son los VALORES AJUSTADOS de la regresión sobre todo el
#'   tramo de estimación, así que NO son pronósticos con información hasta t - 1: la
#'   regresión usó toda la muestra, incluidos los datos posteriores a t. Por eso el MSE
#'   dentro de muestra de una tendencia es optimista frente al de los métodos de
#'   suavizamiento y no debe compararse con él sin decirlo. errores = y - yhat en la
#'   escala original.
#' Pronóstico extramuestral: la tendencia evaluada en t = T + 1, ..., T + h (multiplicada
#'   por el factor de sesgo en la exponencial).
#'
#' Errores estándar: ordinario, sqrt(sigma2hat [(X'X)^{-1}]_jj) con sigma2hat = SCR/(T - p);
#'   y robusto HAC (Newey–West, núcleo de Bartlett, L = floor(4 (T/100)^(2/9)) rezagos, sin
#'   corrección por grados de libertad; ver .varianza_hac). El estadístico t usa el EE
#'   robusto y su valor p la distribución t con T - p grados de libertad (elección
#'   documentada: con T pequeño es más conservadora que la normal). En la exponencial la
#'   tabla, R^2, sigma2hat, DW y residuos están en la escala logarítmica, que es la de
#'   la regresión y donde valen los supuestos.
#'
#' @param y              vector numérico ordenado, sin NA (y > 0 en la exponencial).
#' @param tipo           "lineal", "cuadratica" o "exponencial".
#' @param corregir_sesgo TRUE/FALSE; ver arriba.
#' @return objeto de clase "metodo_pronostico". En parametros:
#'   tipo, coef (vector de la regresión, comparable con los coeficientes de la regresión de R), coef_exp (b0 y b1 en
#'   la escala original, solo exponencial), tabla (tibble: coeficiente, estimacion,
#'   ee_ordinario, ee_robusto, t_robusto, valor_p), r2, sigma2 (SCR/(T - p)), dw (estadístico
#'   de Durbin–Watson de los residuos), residuos, L (rezagos HAC), gl (T - p), factor_sesgo.
#'   n_param = 2 (lineal, exponencial) o 3 (cuadrática).
#'
#' Referencia: Clase 4, Parte VI (tendencias, HAC y retransformación); enunciado 2(b).
ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial"),
                              corregir_sesgo = FALSE) {
  tipo <- match.arg(tipo)
  y <- .validar_serie(y)
  stopifnot("corregir_sesgo debe ser TRUE o FALSE" =
              is.logical(corregir_sesgo) && length(corregir_sesgo) == 1L && !is.na(corregir_sesgo))
  if (tipo == "exponencial") {
    stopifnot("la tendencia exponencial exige y > 0 en todos los datos" = all(y > 0))
  }

  n <- length(y)
  p <- if (tipo == "cuadratica") 3L else 2L
  stopifnot("length(y) debe superar el número de coeficientes" = n > p)

  t <- seq_len(n)
  z <- if (tipo == "exponencial") log(y) else y   # variable sobre la que se hace la regresión
  X <- .matriz_diseno(t, tipo)
  beta <- drop(solve(crossprod(X), crossprod(X, z)))  # ecuaciones normales
  ajustado <- drop(X %*% beta)
  u <- z - ajustado                                    # residuos de la regresión
  scr <- sum(u^2)
  gl <- n - p
  sigma2 <- scr / gl
  r2 <- 1 - scr / sum((z - mean(z))^2)

  L <- floor(4 * (n / 100)^(2 / 9))                    # regla del enunciado para Bartlett
  ee_ord <- sqrt(sigma2 * diag(solve(crossprod(X))))
  ee_rob <- sqrt(diag(.varianza_hac(X, u, L)))
  t_rob <- beta / ee_rob

  nombres <- if (tipo == "exponencial") c("a", "theta") else paste0("beta", 0:(p - 1L))
  tabla <- tibble::tibble(
    coeficiente = nombres, estimacion = beta, ee_ordinario = ee_ord, ee_robusto = ee_rob,
    t_robusto = t_rob, valor_p = 2 * stats::pt(-abs(t_rob), df = gl)
  )

  # Factor de retransformación: 1 salvo en la exponencial con corrección lognormal.
  factor <- if (tipo == "exponencial" && corregir_sesgo) exp(sigma2 / 2) else 1
  yhat <- .a_escala_y(ajustado, tipo, factor)

  parametros <- list(
    tipo = tipo, coef = stats::setNames(beta, nombres),
    coef_exp = if (tipo == "exponencial") c(beta0 = exp(beta[[1]]), beta1 = exp(beta[[2]])),
    tabla = tabla, r2 = r2, sigma2 = sigma2, dw = durbin_watson(u)$estadistico,
    residuos = u, L = L, gl = gl, factor_sesgo = factor
  )

  .nuevo_metodo(
    metodo = sprintf("Tendencia %s", tipo), y = y, yhat = yhat,
    pronosticar = .pronosticador_tendencia(beta, n, tipo, factor),
    parametros = parametros, n_param = p
  )
}
