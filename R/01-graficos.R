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
      y       = attr(datos, "unidad"),
      caption = pie
    ) +
    .tema_tarea()
}
