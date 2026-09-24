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
