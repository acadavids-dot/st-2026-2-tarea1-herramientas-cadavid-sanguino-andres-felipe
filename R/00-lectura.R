# 00-lectura.R
# Lectura de series: leer_serie()

# ---- Auxiliares internas --------------------------------------------------------------

# Paso de calendario que corresponde a cada frecuencia soportada (1, 4, 12).
# Es el único lugar donde se listan las frecuencias soportadas: si la frecuencia no es
# una de ellas, se detiene con un mensaje que dice cuáles sí lo son.
.paso_calendario <- function(frecuencia) {
  switch(as.character(frecuencia),
    "1"  = "year",
    "4"  = "quarter",
    "12" = "month",
    stop("frecuencia no soportada (", frecuencia, "): se soportan 1 (anual), ",
         "4 (trimestral) y 12 (mensual)", call. = FALSE)
  )
}

# Convierte el índice de un objeto ts en fechas Date (primer día del período).
#   frecuencia 1  -> AAAA-01-01
#   frecuencia 4  -> mes 3 * (trimestre - 1) + 1, día 1
#   frecuencia 12 -> mes = ciclo, día 1
# Año y ciclo salen de floor(time(x) + 1e-8) y cycle(x), no de aritmética de decimales:
# time() guarda 1971.5 para el segundo trimestre de 1971, y trocear ese decimal a mano
# puede caer en el período equivocado por error de punto flotante. El 1e-8 protege el
# caso contrario, un año exacto representado como 1970.9999999.
.fechas_de_ts <- function(x) {
  frecuencia <- stats::frequency(x)
  .paso_calendario(frecuencia)  # aborta si la frecuencia no es soportada
  anio  <- as.integer(floor(as.numeric(stats::time(x)) + 1e-8))
  ciclo <- as.integer(stats::cycle(x))
  mes <- switch(as.character(frecuencia),
    "1"  = 1L,
    "4"  = 3L * (ciclo - 1L) + 1L,
    "12" = ciclo
  )
  as.Date(sprintf("%04d-%02d-01", anio, mes))
}

# Infiere la frecuencia de una columna de fechas a partir de la mediana de las
# diferencias en días: ~365 -> 1, ~91 -> 4, ~30 -> 12. Se usa la mediana y rangos, no
# igualdad, porque los meses miden 28 a 31 días, los trimestres 90 a 92 y los años
# 365 o 366.
.frecuencia_de_fechas <- function(fecha) {
  mediana <- stats::median(as.numeric(diff(fecha)))
  if (mediana >= 365 && mediana <= 366) {
    1
  } else if (mediana >= 89 && mediana <= 92) {
    4
  } else if (mediana >= 28 && mediana <= 31) {
    12
  } else {
    stop("no se pudo inferir la frecuencia: la mediana de diferencias entre fechas es ",
         mediana, " días y se soportan ~365 (anual), ~91 (trimestral) y ~30 (mensual)",
         call. = FALSE)
  }
}

# Verifica que las fechas no tengan NA, sean estrictamente crecientes y, si se da
# `paso` ("year", "quarter" o "month"), que coincidan con seq(fecha[1], by = paso).
# Con paso = NULL solo se hacen las dos primeras comprobaciones (la frecuencia de un
# CSV se infiere de las fechas, y para eso deben ser válidas y crecientes).
# Ante una falla se detiene indicando la primera posición que la incumple.
.verificar_fechas <- function(fecha, paso = NULL) {
  falla <- which(is.na(fecha))
  if (length(falla) > 0L) {
    stop("fecha no válida (NA) en la posición ", falla[1], call. = FALSE)
  }
  falla <- which(diff(as.numeric(fecha)) <= 0) + 1L
  if (length(falla) > 0L) {
    stop("las fechas no son estrictamente crecientes: falla la posición ", falla[1],
         " (", format(fecha[falla[1]]), " no es posterior a ", format(fecha[falla[1] - 1L]),
         ")", call. = FALSE)
  }
  if (!is.null(paso)) {
    esperado <- seq(fecha[1], by = paso, length.out = length(fecha))
    falla <- which(fecha != esperado)
    if (length(falla) > 0L) {
      stop("las fechas no son equiespaciadas (paso = ", paso, "): la primera posición ",
           "que falla es ", falla[1], " (", format(fecha[falla[1]]), "; se esperaba ",
           format(esperado[falla[1]]), ")", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ---- leer_serie -----------------------------------------------------------------------

#' leer_serie(x, fuente, unidad)
#'
#' Descripción: convierte una serie de tiempo en el tibble de trabajo del proyecto, con
#'   una fila por observación, y verifica que las fechas sean crecientes y
#'   equiespaciadas según la frecuencia; si no lo son, se detiene.
#'
#' Entradas admitidas:
#'   - un objeto `ts` univariado de frecuencia 1, 4 o 12. Las fechas salen de
#'     time(x) y cycle(x) (ver .fechas_de_ts()): el primer día de cada año, trimestre o mes.
#'   - la ruta a un archivo .csv con columnas `fecha` (AAAA-MM-DD) y `valor`. La
#'     frecuencia se infiere de la mediana de diferencias en días (~365, ~91, ~30).
#'     Las fechas deben caer siempre en el mismo día del mes: seq() por meses desde un
#'     día 29 a 31 desborda a otro mes y la verificación lo reportaría como falla.
#'
#' Verificación de fechas: sin NA, estrictamente crecientes y iguales a
#'   seq(fecha[1], by = paso, length.out = n), con paso "year", "quarter" o "month".
#'   El mensaje de error indica la primera posición que falla.
#'
#' @param x       objeto `ts` univariado, o ruta (cadena) a un .csv con `fecha` y `valor`.
#'                Una serie multivariada (p. ej. Seatbelts) debe entrar ya reducida a una
#'                columna, como Seatbelts[, "DriversKilled"].
#' @param fuente  cadena: fuente original tal como la documenta help() de la serie.
#' @param unidad  cadena: magnitud y unidad de la serie. graficar_serie() la usa tal cual
#'                como rótulo del eje y, por eso conviene incluir la magnitud, p. ej.
#'                "Caudal anual (10^8 m^3)".
#'
#' @return tibble con columnas
#'   t      entero, 1, ..., n
#'   fecha  Date
#'   y      numérico
#'   y los atributos `frecuencia` (1, 4 o 12), `fuente` y `unidad`. Los atributos se
#'   pierden con algunas operaciones de dplyr; graficar_serie() los exige por eso.
#'
#' Referencia: enunciado, sección 2(a).
leer_serie <- function(x, fuente, unidad) {
  stopifnot(
    "fuente debe ser una cadena de texto" =
      is.character(fuente) && length(fuente) == 1L && !is.na(fuente),
    "unidad debe ser una cadena de texto" =
      is.character(unidad) && length(unidad) == 1L && !is.na(unidad)
  )

  if (stats::is.ts(x)) {
    stopifnot(
      "x debe ser una serie univariada; para una serie multivariada use una sola columna, p. ej. Seatbelts[, \"DriversKilled\"]" =
        NCOL(x) == 1L,
      "x debe ser numérica" = is.numeric(x)
    )
    frecuencia <- stats::frequency(x)
    fecha <- .fechas_de_ts(x)
    y <- as.numeric(x)
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("no existe el archivo: ", x, call. = FALSE)
    }
    crudo <- utils::read.csv(x, stringsAsFactors = FALSE)
    stopifnot(
      "el CSV debe tener las columnas 'fecha' y 'valor'" =
        all(c("fecha", "valor") %in% names(crudo)),
      "el CSV debe tener al menos dos filas" = nrow(crudo) >= 2L,
      "la columna 'valor' del CSV debe ser numérica" = is.numeric(crudo$valor)
    )
    fecha <- as.Date(crudo$fecha)
    .verificar_fechas(fecha)  # NA y orden, antes de inferir la frecuencia
    frecuencia <- .frecuencia_de_fechas(fecha)
    y <- as.numeric(crudo$valor)
  } else {
    stop("x debe ser un objeto ts o la ruta (cadena) a un archivo .csv", call. = FALSE)
  }

  .verificar_fechas(fecha, .paso_calendario(frecuencia))

  datos <- tibble::tibble(t = seq_along(y), fecha = fecha, y = y)
  attr(datos, "frecuencia") <- frecuencia
  attr(datos, "fuente") <- fuente
  attr(datos, "unidad") <- unidad
  datos
}
