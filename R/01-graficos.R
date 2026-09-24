# 01-graficos.R
# Gráficos: graficar_serie(), correlograma()

# ---- Tema común -----------------------------------------------------------------------

# Tema de todas las figuras del proyecto, definido una sola vez aquí. Partiendo de
# theme_minimal (sin fondo gris) se fija el fondo blanco de forma explícita para que los
# PNG guardados con ggsave() no queden transparentes.
.tema_tarea <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(colour = "grey30"),
      plot.caption     = ggplot2::element_text(colour = "grey40", hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Marcas del eje de fechas (1 de enero de años «redondos»). El paso es el menor de una
# lista de pasos redondos que deja como máximo unas 8 marcas en el rango de la serie,
# para que el eje se lea igual en una serie de 24 años que en una de 289. Las marcas caen
# en múltiplos del paso (1880, 1900, ...) y no a partir del primer dato (1871, 1891, ...),
# que es lo que haría date_breaks = "20 years".
.marcas_anios <- function(fecha) {
  anio <- as.integer(format(range(fecha), "%Y"))
  pasos <- c(1, 2, 5, 10, 20, 25, 50, 100)
  paso <- pasos[which(pasos >= (anio[2] - anio[1]) / 8)[1]]
  marcas <- seq(ceiling(anio[1] / paso) * paso, floor(anio[2] / paso) * paso, by = paso)
  as.Date(sprintf("%04d-01-01", as.integer(marcas)))
}

# ---- graficar_serie -------------------------------------------------------------------

#' graficar_serie(datos, titulo)
#'
#' Descripción: gráfico de la serie en el tiempo, pensado para leerse sin el texto que
#'   lo acompaña: título, eje x con años, eje y rotulado con la unidad y un pie
#'   (caption) con la fuente y el número de observaciones.
#'
#' Decisiones:
#'   - El rótulo del eje y es el atributo `unidad` tal como lo dio el usuario a
#'     leer_serie(), que por eso debe incluir la magnitud y su unidad.
#'   - Eje x: «Año» cuando la serie es anual y «Fecha» en el resto; las marcas se
#'     rotulan con el año, en múltiplos de un paso redondo (ver .marcas_anios()).
#'   - Devuelve el objeto ggplot; guardarlo en disco es tarea de ejemplos.R.
#'
#' @param datos   tibble devuelto por leer_serie() (con sus atributos intactos).
#' @param titulo  cadena con el título del gráfico.
#' @return objeto `ggplot`.
#'
#' Referencia: enunciado, sección 2(b).
graficar_serie <- function(datos, titulo) {
  stopifnot(
    "datos debe ser el tibble que devuelve leer_serie()" =
      is.data.frame(datos) && all(c("fecha", "y") %in% names(datos)),
    "datos no tiene los atributos frecuencia, fuente y unidad de leer_serie()" =
      !is.null(attr(datos, "frecuencia")) && !is.null(attr(datos, "fuente")) &&
        !is.null(attr(datos, "unidad")),
    "titulo debe ser una cadena de texto" =
      is.character(titulo) && length(titulo) == 1L && !is.na(titulo)
  )

  # El pie se parte en líneas para que una fuente larga no se salga del gráfico. Se quita
  # el punto final de la fuente (las citas suelen traerlo) para no escribir «..».
  pie <- paste(
    strwrap(sprintf("Fuente: %s. n = %d observaciones.",
                    sub("[.]\\s*$", "", attr(datos, "fuente")), nrow(datos)),
            width = 110),
    collapse = "\n"
  )

  ggplot2::ggplot(datos, ggplot2::aes(x = fecha, y = y)) +
    ggplot2::geom_line(colour = "#1F4E79", linewidth = 0.6) +
    ggplot2::scale_x_date(breaks = .marcas_anios(datos$fecha), date_labels = "%Y") +
    ggplot2::labs(
      title   = titulo,
      x       = if (attr(datos, "frecuencia") == 1) "Año" else "Fecha",
      y       = paste(strwrap(attr(datos, "unidad"), width = 40), collapse = "
"),  # unidades largas: varias líneas
      caption = pie
    ) +
    .tema_tarea()
}

# ---- Gráfico final de cada ejemplo ------------------------------------------------------

# Gráfico final del protocolo (enunciado 4(a)): la serie completa, el ajuste dentro de la
# muestra (yhat sobre el tramo de estimación) y los h pronósticos sobre el tramo de
# validación, con una línea vertical en el corte. Si se da `referente` (el ingenuo, o el
# ingenuo estacional) se dibuja también sobre el tramo de validación para poder compararlo.
#   datos    tibble de leer_serie() con la serie COMPLETA (estimación + validación)
#   yhat_est yhat del método sobre el tramo de estimación (largo T_est; NA en el calentamiento)
#   pron     los h pronósticos extramuestrales (largo h = nrow(datos) - T_est)
#   referente  pronósticos del referente sobre el tramo de validación (largo h) o NULL
#   nombre_referente  rótulo de la leyenda: «Referente ingenuo» o «Referente ingenuo estacional»
# La línea del corte queda a mitad de camino entre la última fecha de estimación y la
# primera de validación, así no coincide con ningún dato.
.grafico_pronostico <- function(datos, yhat_est, pron, referente = NULL, titulo,
                                  nombre_referente = "Referente ingenuo") {
  n <- nrow(datos)
  h <- length(pron)
  T_est <- n - h
  stopifnot(
    "yhat_est debe tener un valor por cada observación del tramo de estimación" = length(yhat_est) == T_est,
    "referente debe tener el mismo largo que pron" = is.null(referente) || length(referente) == h,
    "h debe ser al menos 1" = h >= 1L
  )
  f <- datos$fecha
  obs <- "Serie observada"; aju <- "Ajuste dentro de la muestra"
  pro <- "Pronóstico"; ref <- nombre_referente
  df <- rbind(
    data.frame(fecha = f, valor = datos$y, tipo = obs),
    data.frame(fecha = f[seq_len(T_est)], valor = yhat_est, tipo = aju),
    data.frame(fecha = f[T_est + seq_len(h)], valor = pron, tipo = pro)
  )
  if (!is.null(referente)) {
    df <- rbind(df, data.frame(fecha = f[T_est + seq_len(h)], valor = referente, tipo = ref))
  }
  df <- df[!is.na(df$valor), ]
  df$tipo <- factor(df$tipo, levels = c(obs, aju, pro, ref))
  corte <- as.Date(mean(as.numeric(f[c(T_est, T_est + 1L)])))

  graficar_serie(datos, titulo) +
    ggplot2::geom_vline(xintercept = corte, linetype = "dotted", colour = "grey30") +
    ggplot2::geom_line(data = df, ggplot2::aes(x = fecha, y = valor, colour = tipo,
                                               linetype = tipo), linewidth = 0.6) +
    ggplot2::geom_point(data = df[df$tipo == pro, ], ggplot2::aes(x = fecha, y = valor,
                                                                  colour = tipo), size = 1.6) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c("#1F4E79", "#E08E0B", "#B03A2E", "grey45"), c(obs, aju, pro, ref)),
      name = NULL, drop = FALSE) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(c("solid", "solid", "solid", "dashed"), c(obs, aju, pro, ref)),
      name = NULL, drop = FALSE) +
    ggplot2::annotate("text", x = corte, y = Inf, label = " validación →", hjust = 0,
                      vjust = 1.6, size = 3, colour = "grey30") +
    ggplot2::theme(legend.position = "bottom")
}

# Residuos contra valores ajustados de una tendencia (enunciado 4(a), supuesto de varianza
# constante): una nube sin forma alrededor de cero apoya la varianza constante; un abanico o
# una curva la contradicen. Ambos ejes están en la escala de la regresión (en la exponencial,
# logaritmos), que es donde valen los supuestos.
.grafico_residuos_ajustados <- function(ajustado, residuos, titulo) {
  stopifnot("ajustado y residuos deben tener el mismo largo" = length(ajustado) == length(residuos))
  ggplot2::ggplot(data.frame(a = ajustado, r = residuos), ggplot2::aes(x = a, y = r)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_point(colour = "#1F4E79", size = 1.8) +
    ggplot2::labs(title = titulo, x = "Valores ajustados (escala de la regresión)",
                  y = "Residuos") +
    .tema_tarea()
}

# ---- correlograma ---------------------------------------------------------------------

# Un panel del correlograma: barras verticales desde cero para los rezagos 1, ..., m
# (el rezago 0 no se dibuja) y la banda como dos líneas discontinuas en ±banda. Se fuerza
# que el eje y incluya la banda aunque todas las barras queden dentro de ella.
.panel_correlograma <- function(valores, banda, titulo, etiqueta_y) {
  df <- data.frame(h = seq_along(valores), v = valores)
  ggplot2::ggplot(df, ggplot2::aes(x = h, y = v)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_segment(ggplot2::aes(xend = h, yend = 0), colour = "#1F4E79",
                          linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(-banda, banda), linetype = "dashed",
                        colour = "#B03A2E") +
    ggplot2::scale_x_continuous(breaks = function(lim) {
      marcas <- pretty(lim)
      marcas[marcas >= 1]  # sin marca en el rezago 0
    }) +
    ggplot2::expand_limits(y = c(-banda, banda)) +
    ggplot2::labs(title = titulo, x = "Rezago h", y = etiqueta_y) +
    .tema_tarea()
}

#' correlograma(datos, m)
#'
#' Descripción: correlograma de una serie, o de los errores de un método, en un panel de
#'   dos gráficos (ACF arriba, PACF abajo).
#'
#' Ecuaciones:
#'   ACF muestral a mano, con divisor único T (Definición 9 de la Clase 3):
#'     r_h = sum_{t=h+1}^{n} (y_t - ybar)(y_{t-h} - ybar) / sum_{t=1}^{n} (y_t - ybar)^2,
#'     h = 1, ..., m.
#'   El numerador tiene n - h sumandos y el denominador n, por lo que r_h queda sesgada
#'   hacia cero; se prefiere porque garantiza una matriz de autocorrelaciones semidefinida
#'   positiva (Clase 3).
#'   PACF: pacf(y, lag.max = m, plot = FALSE)$acf, de stats.
#'   Banda: la misma línea que dibuja plot.acf (método de stats), clim0 <- qnorm((1 + ci)/2)/sqrt(x$n.used)
#'   con ci = 0.95, es decir +/- qnorm(0.975)/sqrt(n) ~ 1.96/sqrt(n). No es un estadístico:
#'   es el intervalo de la distribución asintótica de r_h bajo ruido blanco (Bartlett),
#'   r_h ~ N(0, 1/n). Se dibuja igual en los dos paneles.
#'
#' Sobre n: es el número de valores de la sucesión que se grafica, no el T de la serie
#'   original. Con errores de un método, los NA del calentamiento se eliminan y n es el
#'   número de errores que quedan; la ACF, la PACF y la banda se calculan sobre esos n.
#'   Los NA solo pueden estar al inicio (calentamiento): quitar NA del medio cambiaría
#'   qué observaciones están a h períodos de distancia, por eso se rechazan.
#'
#' @param datos  tibble de leer_serie() (se usa la columna `y`) o vector numérico.
#' @param m      número de rezagos (entero, 1 <= m < n). Por defecto min(floor(n / 4), 24).
#' @return lista con `acf` (vector r_h, h = 1..m), `pacf` (vector, h = 1..m), `banda`
#'   (semiancho), `n`, `m` y `grafico` (panel patchwork, subtítulo con n y m).
#'
#' Referencia: enunciado, sección 2(c); Clase 3, Parte II.
correlograma <- function(datos, m = NULL) {
  y <- if (is.data.frame(datos)) datos$y else datos
  stopifnot("datos debe ser el tibble de leer_serie() o un vector numérico" =
              is.numeric(y) && length(y) > 0L)
  y <- as.numeric(y)

  valido <- !is.na(y)
  stopifnot("y no puede ser todo NA" = any(valido))
  primero <- which(valido)[1]
  stopifnot("los NA solo pueden estar al inicio (calentamiento), no entre los valores" =
              all(valido[primero:length(y)]))
  y <- y[primero:length(y)]
  n <- length(y)

  if (is.null(m)) {
    m <- min(floor(n / 4), 24)
  }
  stopifnot(
    "m debe ser un entero positivo menor que el número de observaciones sin NA" =
      is.numeric(m) && length(m) == 1L && !is.na(m) && m >= 1 && m == floor(m) && m < n,
    "la serie es constante: la ACF no está definida" = sum((y - mean(y))^2) > 0
  )

  # ACF a mano, divisor único: el denominador suma los n cuadrados y cada numerador
  # suma los n - h productos disponibles.
  yc <- y - mean(y)
  den <- sum(yc^2)
  r <- vapply(seq_len(m), function(h) sum(yc[(h + 1):n] * yc[1:(n - h)]) / den, numeric(1))

  # PACF: la única llamada a stats permitida para esto (pacf con plot = FALSE).
  p <- as.numeric(stats::pacf(y, lag.max = m, plot = FALSE)$acf)

  banda <- stats::qnorm((1 + 0.95) / 2) / sqrt(n)

  grafico <- patchwork::wrap_plots(
    .panel_correlograma(r, banda, "Autocorrelación muestral (ACF)", "r_h"),
    .panel_correlograma(p, banda, "Autocorrelación parcial (PACF)", "PACF"),
    ncol = 1
  ) +
    patchwork::plot_annotation(
      title = "Correlograma",
      subtitle = sprintf(
        "n = %d observaciones, m = %d rezagos. Líneas discontinuas: banda de ruido blanco al 95 %% (±%s)",
        n, as.integer(m), formatC(banda, format = "f", digits = 4, decimal.mark = ",")
      ),
      theme = .tema_tarea()
    )

  list(acf = r, pacf = p, banda = banda, n = n, m = as.integer(m), grafico = grafico)
}
