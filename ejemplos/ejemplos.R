# ejemplos.R
# Se ejecuta desde la raíz del repositorio: source("ejemplos/ejemplos.R")
#
# Bloque 0: verificaciones de las funciones contra R (se detiene si alguna falla).
# Bloques 1 a 8: un ejemplo por método, con el protocolo del enunciado (sección 4(a)).
# Bloque 9: contraejemplo (SES sobre una serie con tendencia y estacionalidad).
# Las figuras se guardan en figs/ con el nombre ejNN-metodo-serie-tipo.png.

invisible(lapply(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)), source))
stopifnot("ejemplos.R debe ejecutarse desde la raíz del repositorio" = exists("leer_serie"))

dir.create("figs", showWarnings = FALSE)
set.seed(2026)  # ningún paso de este script es aleatorio; se fija por regla del enunciado
t_inicio <- proc.time()

# ---- Auxiliares de este script ---------------------------------------------------------

# Guarda una figura en figs/ con tamaño y resolución fijos.
.guardar_fig <- function(g, nombre, alto = 5) {
  ggplot2::ggsave(file.path("figs", nombre), g, width = 8, height = alto, dpi = 150)
}

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

# Una línea por prueba: estadístico, grados de libertad, valor p y decisión.
.linea <- function(obj) {
  cat(sprintf("  %-40s est = %s  gl = %s  p = %s  -> %s\n", obj$nombre, .fmt_num(obj$estadistico),
              if (is.na(obj$gl)) "-" else format(obj$gl), .fmt_p(obj$valor_p),
              gsub("[$]", "", obj$decision)))
}

# Cotas d_L y d_U de Durbin–Watson al 5 % (una cola) para los N errores de cada ejemplo y
# k' regresores sin contar el intercepto. Fuente: tabla de Savin y White (1977). La tabla no
# cubre todos los n (después de 40 solo trae múltiplos de 5), así que se calcularon de forma
# exacta con la definición de Durbin y Watson (1951) y la integral de Imhof (1961); ese
# cálculo reproduce 12 filas de la tabla con diferencias de a lo sumo 0,005 (ver README).
# Se redondean a 2 decimales, como la tabla. Para los métodos de suavizamiento se usa
# k' = 1 (aproximación declarada: no hay matriz de diseño).
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
# Devuelve dL y dU del ejemplo `id` y se detiene si el número de errores cambió: las cotas
# de arriba solo valen para ese N.
.cotas <- function(id, N) {
  fila <- cotas_dw[[id]]
  stopifnot("no hay cotas de Durbin-Watson para este ejemplo" = !is.null(fila),
            "el número de errores cambió: recalcular las cotas de Durbin-Watson" = N == fila$N)
  fila
}

# Fila del resumen final: MASE del método y del referente en el tramo de validación.
.fila_resumen <- function(ejemplo, serie, fit, m_val, m_ref) {
  data.frame(ejemplo = ejemplo, serie = serie, metodo = fit$metodo,
             MASE_metodo = m_val$MASE, MASE_referente = m_ref$MASE)
}

resumen <- list()

# ======================================================================================
# Bloque 0. Verificaciones contra R
# ======================================================================================
# Es el único lugar (junto al informe) donde aparecen acf(), lm() y Box.test(). Cada fila
# guarda la diferencia máxima observada y su tolerancia; si alguna la supera, el script se
# detiene. Las equivalencias de los métodos se miden como error RELATIVO a max|y|: con datos
# del orden de 1e4 un ulp ya mide ~2e-12, así que una cota absoluta de 1e-12 no la cumple
# ninguna implementación correcta (ver README, Convenciones).
cat("\n== Bloque 0: verificaciones contra R ==\n")
fila_v <- function(prueba, dif, tol) data.frame(prueba = prueba, dif = dif, tol = tol, ok = dif < tol)

series0 <- list(Nile = Nile, LakeHuron = LakeHuron, austres = austres, discoveries = discoveries,
                airmiles = airmiles, JohnsonJohnson = JohnsonJohnson, AirPassengers = AirPassengers)

# ACF a mano, banda y Ljung–Box contra acf() y Box.test()
dif_acf <- dif_banda <- dif_Q <- dif_Qp <- 0
for (x in series0) {
  y <- as.numeric(x); n <- length(y)
  for (m in unique(c(4, min(floor(n / 4), 24)))) {
    cg <- correlograma(y, m)
    ref <- acf(y, lag.max = m, plot = FALSE)
    dif_acf <- max(dif_acf, max(abs(cg$acf - ref$acf[-1])))
    dif_banda <- max(dif_banda, abs(cg$banda - qnorm((1 + 0.95) / 2) / sqrt(ref$n.used)))
    for (p in 0:min(3, m - 1)) {
      lb <- ljung_box(cg$acf, n, m, p)
      bt <- Box.test(y, lag = m, type = "Ljung-Box", fitdf = p)
      dif_Q <- max(dif_Q, abs(lb$estadistico - bt$statistic) / bt$statistic)
      dif_Qp <- max(dif_Qp, abs(lb$valor_p - bt$p.value))
    }
  }
}

# Equivalencias de los métodos (error relativo a max|y|)
dif_mm <- dif_ses <- dif_pesos <- dif_dmm <- dif_holt <- 0
for (x in series0) {
  y <- as.numeric(x); n <- length(y); esc <- max(abs(y))
  for (k in 2:12) {                                     # MM recursiva vs promedio directo por ventana
    directa <- rep(NA_real_, n)
    for (t in (k + 1):n) directa[t] <- mean(y[(t - k):(t - 1)])
    dif_mm <- max(dif_mm, max(abs(ajustar_mm(y, k)$yhat - directa), na.rm = TRUE) / esc)
    if (2 * k - 1 <= n) {                               # DMM (un recorrido) vs definición directa
      mm <- rep(NA_real_, n); for (t in k:n) mm[t] <- mean(y[(t - k + 1):t])
      dm <- rep(NA_real_, n); for (t in (2 * k - 1):n) dm[t] <- mean(mm[(t - k + 1):t])
      yh <- c(NA, (2 * mm - dm + 2 / (k - 1) * (mm - dm))[-n])
      dif_dmm <- max(dif_dmm, max(abs(ajustar_dmm(y, k)$yhat - yh), na.rm = TRUE) / esc)
    }
  }
  for (a in c(0.05, 0.15, 0.5, 0.98)) {                 # SES vs promedio ponderado finito
    f <- ajustar_ses(y, a)
    pond <- vapply(2:n, function(t) {
      i <- 0:(t - 2)
      a * sum((1 - a)^i * y[t - i]) + (1 - a)^(t - 1) * y[1]
    }, 0)
    suma_pesos <- vapply(2:n, function(t) a * sum((1 - a)^(0:(t - 2))) + (1 - a)^(t - 1), 0)
    dif_pesos <- max(dif_pesos, max(abs(suma_pesos - 1)))
    dif_ses <- max(dif_ses, max(abs(c(f$yhat[3:n], f$parametros$siguiente) - pond)) / esc)
  }
  for (ab in list(c(0.1, 0.1), c(0.3, 0.4), c(0.9, 0.05))) {   # Holt vs forma de corrección de error
    f <- ajustar_holt(y, ab[1], ab[2])
    e <- y - f$yhat; L <- f$parametros$nivel; Tt <- f$parametros$pendiente; t <- 2:n
    dif_holt <- max(dif_holt, max(abs(L[t] - (L[t - 1] + Tt[t - 1] + ab[1] * e[t])),
                                  abs(Tt[t] - (Tt[t - 1] + ab[1] * ab[2] * e[t]))) / esc)
  }
}

# Tendencias: coeficientes contra lm()
dif_lm <- 0
for (cs in list(list(LakeHuron, "lineal"), list(airmiles, "cuadratica"), list(airmiles, "exponencial"),
                list(JohnsonJohnson, "exponencial"), list(Nile, "cuadratica"))) {
  y <- as.numeric(cs[[1]]); t <- seq_along(y)
  fit_lm <- switch(cs[[2]], lineal = lm(y ~ t), cuadratica = lm(y ~ t + I(t^2)), exponencial = lm(log(y) ~ t))
  dif_lm <- max(dif_lm, max(abs(unname(ajustar_tendencia(y, cs[[2]])$parametros$coef) - unname(coef(fit_lm)))))
}

# Cifras publicadas en las notas: Clase 3 (precios de cierre, mayo-diciembre) y Clase 4
# (créditos al consumo, once meses). Las notas publican 4 cifras significativas o 2 decimales,
# por eso la tolerancia es media unidad de la última cifra publicada.
precio <- c(18.43, 19.98, 19.51, 20.63, 19.78, 21.25, 21.18, 22.14)
credito <- c(133, 155, 165, 171, 194, 231, 274, 312, 313, 333, 343)
m_mm3 <- medidas(precio, ajustar_mm(precio, 3)$yhat, precio)
ses15 <- ajustar_ses(precio, 0.15); yh_ses <- ses15$yhat; yh_ses[2:3] <- NA   # origen común t = 4
m_ses15 <- medidas(precio, yh_ses, precio)
f_dmm3 <- ajustar_dmm(credito, 3); m_dmm3 <- medidas(credito, f_dmm3$yhat, credito)
f_holt <- ajustar_holt(credito, 0.3, 0.4); m_holt <- medidas(credito, f_holt$yhat, credito)

verif <- do.call(rbind, list(
  fila_v("ACF a mano vs acf() (7 series, 2 valores de m)", dif_acf, 1e-12),
  fila_v("Banda vs qnorm(0,975)/sqrt(n.used) de plot.acf", dif_banda, 1e-12),
  fila_v("Q de Ljung-Box vs Box.test, p = 0..3 (relativa)", dif_Q, 1e-12),
  fila_v("Valor p de Ljung-Box vs Box.test con fitdf = p", dif_Qp, 1e-12),
  fila_v("MM recursiva vs directa, k = 2..12 (relativa)", dif_mm, 1e-12),
  fila_v("SES vs promedio ponderado finito (relativa)", dif_ses, 1e-12),
  fila_v("SES: los pesos suman 1", dif_pesos, 1e-12),
  fila_v("DMM un recorrido vs definición directa (relativa)", dif_dmm, 1e-12),
  fila_v("Holt vs forma de corrección de error (relativa)", dif_holt, 1e-12),
  fila_v("Coeficientes de las tendencias vs lm()", dif_lm, 1e-8),
  fila_v("Clase 3: MM(3) MSE, MAD y MAPE", max(abs(c(m_mm3$MSE - 1.1621, m_mm3$MAD - 0.9780, (m_mm3$MAPE - 4.607) / 10))), 5e-5),
  fila_v("Clase 3: SES(0,15) en origen común, MSE", abs(m_ses15$MSE - 3.3706), 5e-5),
  fila_v("Clase 4: DMM(3) MSE y pronósticos", max(abs(c(m_dmm3$MSE - 841.24, f_dmm3$pronosticar(3) - c(356.56, 370.00, 383.44)))), 5e-3),
  fila_v("Clase 4: Holt(0,3; 0,4) MSE y pronósticos", max(abs(c(m_holt$MSE - 1133.08, f_holt$pronosticar(3) - c(387.15, 414.97, 442.80)))), 5e-3)
))
verif_imp <- verif; verif_imp$dif <- format(verif$dif, digits = 3, scientific = TRUE)
verif_imp$tol <- format(verif$tol, scientific = TRUE)
print(verif_imp, row.names = FALSE, right = FALSE)
stopifnot("alguna verificación del Bloque 0 superó su tolerancia" = all(verif$ok))
cat("Bloque 0: todas las verificaciones pasan.\n")

# ======================================================================================
# Ejemplo 1. Media simple sobre discoveries
# ======================================================================================
cat("\n== Ejemplo 1: media simple sobre discoveries ==\n")
d <- leer_serie(discoveries,
                fuente = "The World Almanac and Book of Facts, 1975 Edition, pages 315-318 (paquete datasets, discoveries)",
                unidad = "Descubrimientos importantes por año (número)")
n_total <- nrow(d); h <- min(12, floor(0.2 * n_total)); T_est <- n_total - h
est <- d[seq_len(T_est), ]; val <- d[T_est + seq_len(h), ]
cat(sprintf("T = %d, h = %d, T_est = %d; corte entre %s y %s\n", n_total, h, T_est, d$fecha[T_est], d$fecha[T_est + 1]))
.guardar_fig(graficar_serie(d, "Ejemplo 1. Descubrimientos importantes por año, 1860-1959"), "ej01-media-discoveries-serie.png")
cg <- correlograma(d)
.guardar_fig(cg$grafico, "ej01-media-discoveries-correlograma.png", alto = 6)
.linea(ljung_box(cg$acf, n_total, cg$m, 0))

fit <- ajustar_media(est$y)
pron <- fit$pronosticar(h); ref <- rep(est$y[T_est], h)
m_in <- .medir(est$y, fit$yhat, est$y); m_val <- .medir(val$y, pron, est$y); m_ref <- .medir(val$y, ref, est$y)
print(.tabla_medidas(m_in, m_val, m_ref, "Ingenuo"), digits = 5, row.names = FALSE)
cot <- .cotas("ej01", sum(!is.na(fit$errores)))
ve <- validar_errores(fit$errores, fit$n_param, T_est, cot$dL, cot$dU)
cat(sprintf("Validación de errores (N = %d):\n", ve$n))
for (p in list(ve$t_media, ve$ljung_box, ve$jarque_bera, ve$durbin_watson)) .linea(p)
.guardar_fig(ve$grafico, "ej01-media-discoveries-errores.png", alto = 7)
.guardar_fig(.grafico_pronostico(d, fit$yhat, pron, ref, "Ejemplo 1. Media simple: ajuste y pronóstico"),
             "ej01-media-discoveries-pronostico.png")
resumen[[1]] <- .fila_resumen(1, "discoveries", fit, m_val, m_ref)

# ======================================================================================
# Ejemplo 2. Media móvil sobre Nile
# ======================================================================================
cat("\n== Ejemplo 2: media móvil sobre Nile ==\n")
d <- leer_serie(Nile,
                fuente = "Durbin J, Koopman SJ (2001). Time Series Analysis by State Space Methods. Oxford University Press (paquete datasets, Nile)",
                unidad = "Caudal anual del Nilo en Asuán (10^8 m^3)")
n_total <- nrow(d); h <- min(12, floor(0.2 * n_total)); T_est <- n_total - h
est <- d[seq_len(T_est), ]; val <- d[T_est + seq_len(h), ]
cat(sprintf("T = %d, h = %d, T_est = %d; corte entre %s y %s\n", n_total, h, T_est, d$fecha[T_est], d$fecha[T_est + 1]))
.guardar_fig(graficar_serie(d, "Ejemplo 2. Caudal anual del Nilo en Asuán, 1871-1970"), "ej02-mm-nile-serie.png")
cg <- correlograma(d)
.guardar_fig(cg$grafico, "ej02-mm-nile-correlograma.png", alto = 6)
.linea(ljung_box(cg$acf, n_total, cg$m, 0))

opt <- optimizar(est$y, "mm")
cat(sprintf("Óptimo: k = %d (MSE = %s); en el borde de la rejilla: %s\n", opt$optimo$k, .fmt_num(opt$optimo$mse), opt$en_borde))
.guardar_fig(opt$grafico, "ej02-mm-nile-optimizacion.png")
fit <- ajustar_mm(est$y, opt$optimo$k)
pron <- fit$pronosticar(h); ref <- rep(est$y[T_est], h)
m_in <- .medir(est$y, fit$yhat, est$y); m_val <- .medir(val$y, pron, est$y); m_ref <- .medir(val$y, ref, est$y)
print(.tabla_medidas(m_in, m_val, m_ref, "Ingenuo"), digits = 5, row.names = FALSE)
cot <- .cotas("ej02", sum(!is.na(fit$errores)))
ve <- validar_errores(fit$errores, fit$n_param, T_est, cot$dL, cot$dU)
cat(sprintf("Validación de errores (N = %d):\n", ve$n))
for (p in list(ve$t_media, ve$ljung_box, ve$jarque_bera, ve$durbin_watson)) .linea(p)
.guardar_fig(ve$grafico, "ej02-mm-nile-errores.png", alto = 7)
.guardar_fig(.grafico_pronostico(d, fit$yhat, pron, ref, "Ejemplo 2. Media móvil: ajuste y pronóstico"),
             "ej02-mm-nile-pronostico.png")
resumen[[2]] <- .fila_resumen(2, "Nile", fit, m_val, m_ref)

# ======================================================================================
# Ejemplo 3. Suavizamiento exponencial simple sobre Nile
# ======================================================================================
cat("\n== Ejemplo 3: suavizamiento exponencial simple sobre Nile ==\n")
# Misma serie del ejemplo 2 (el enunciado permite dos ejemplos por serie): d, est y val
# ya están cargados; se repite el gráfico y el correlograma con el nombre de este ejemplo.
cat(sprintf("T = %d, h = %d, T_est = %d; corte entre %s y %s\n", n_total, h, T_est, d$fecha[T_est], d$fecha[T_est + 1]))
.guardar_fig(graficar_serie(d, "Ejemplo 3. Caudal anual del Nilo en Asuán, 1871-1970"), "ej03-ses-nile-serie.png")
.guardar_fig(cg$grafico, "ej03-ses-nile-correlograma.png", alto = 6)

opt <- optimizar(est$y, "ses")
cat(sprintf("Óptimo: alpha = %s (MSE = %s); en el borde de la rejilla: %s\n", format(opt$optimo$alpha),
            .fmt_num(opt$optimo$mse), opt$en_borde))
.guardar_fig(opt$grafico, "ej03-ses-nile-optimizacion.png")
fit <- ajustar_ses(est$y, opt$optimo$alpha)
pron <- fit$pronosticar(h); ref <- rep(est$y[T_est], h)
m_in <- .medir(est$y, fit$yhat, est$y); m_val <- .medir(val$y, pron, est$y); m_ref <- .medir(val$y, ref, est$y)
print(.tabla_medidas(m_in, m_val, m_ref, "Ingenuo"), digits = 5, row.names = FALSE)
cot <- .cotas("ej03", sum(!is.na(fit$errores)))
ve <- validar_errores(fit$errores, fit$n_param, T_est, cot$dL, cot$dU)
cat(sprintf("Validación de errores (N = %d):\n", ve$n))
for (p in list(ve$t_media, ve$ljung_box, ve$jarque_bera, ve$durbin_watson)) .linea(p)
.guardar_fig(ve$grafico, "ej03-ses-nile-errores.png", alto = 7)
.guardar_fig(.grafico_pronostico(d, fit$yhat, pron, ref, "Ejemplo 3. Suavizamiento exponencial simple: ajuste y pronóstico"),
             "ej03-ses-nile-pronostico.png")
resumen[[3]] <- .fila_resumen(3, "Nile", fit, m_val, m_ref)

# ======================================================================================
# Ejemplo 4. Doble media móvil sobre austres
# ======================================================================================
cat("\n== Ejemplo 4: doble media móvil sobre austres ==\n")
d <- leer_serie(austres,
                fuente = "Brockwell PJ, Davis RA (1996). Introduction to Time Series and Forecasting. Springer (paquete datasets, austres)",
                unidad = "Residentes de Australia (miles)")
n_total <- nrow(d); h <- min(12, floor(0.2 * n_total)); T_est <- n_total - h
est <- d[seq_len(T_est), ]; val <- d[T_est + seq_len(h), ]
cat(sprintf("T = %d, h = %d, T_est = %d; corte entre %s y %s\n", n_total, h, T_est, d$fecha[T_est], d$fecha[T_est + 1]))
.guardar_fig(graficar_serie(d, "Ejemplo 4. Residentes de Australia, trimestral"), "ej04-dmm-austres-serie.png")
cg <- correlograma(d)
.guardar_fig(cg$grafico, "ej04-dmm-austres-correlograma.png", alto = 6)
.linea(ljung_box(cg$acf, n_total, cg$m, 0))

opt <- optimizar(est$y, "dmm")
cat(sprintf("Óptimo: k = %d (MSE = %s); en el borde de la rejilla: %s\n", opt$optimo$k, .fmt_num(opt$optimo$mse), opt$en_borde))
.guardar_fig(opt$grafico, "ej04-dmm-austres-optimizacion.png")
fit <- ajustar_dmm(est$y, opt$optimo$k)
pron <- fit$pronosticar(h); ref <- rep(est$y[T_est], h)
m_in <- .medir(est$y, fit$yhat, est$y); m_val <- .medir(val$y, pron, est$y); m_ref <- .medir(val$y, ref, est$y)
print(.tabla_medidas(m_in, m_val, m_ref, "Ingenuo"), digits = 5, row.names = FALSE)
cot <- .cotas("ej04", sum(!is.na(fit$errores)))
ve <- validar_errores(fit$errores, fit$n_param, T_est, cot$dL, cot$dU)
cat(sprintf("Validación de errores (N = %d):\n", ve$n))
for (p in list(ve$t_media, ve$ljung_box, ve$jarque_bera, ve$durbin_watson)) .linea(p)
.guardar_fig(ve$grafico, "ej04-dmm-austres-errores.png", alto = 7)
.guardar_fig(.grafico_pronostico(d, fit$yhat, pron, ref, "Ejemplo 4. Doble media móvil: ajuste y pronóstico"),
             "ej04-dmm-austres-pronostico.png")
resumen[[4]] <- .fila_resumen(4, "austres", fit, m_val, m_ref)

# ======================================================================================
# Ejemplo 5. Tendencia lineal sobre LakeHuron
# ======================================================================================
cat("\n== Ejemplo 5: tendencia lineal sobre LakeHuron ==\n")
d <- leer_serie(LakeHuron,
                fuente = "Brockwell PJ, Davis RA (1991). Time Series: Theory and Methods, 2nd ed. Springer, page 555 (paquete datasets, LakeHuron)",
                unidad = "Nivel del lago Hurón (pies)")
n_total <- nrow(d); h <- min(12, floor(0.2 * n_total)); T_est <- n_total - h
est <- d[seq_len(T_est), ]; val <- d[T_est + seq_len(h), ]
cat(sprintf("T = %d, h = %d, T_est = %d; corte entre %s y %s\n", n_total, h, T_est, d$fecha[T_est], d$fecha[T_est + 1]))
.guardar_fig(graficar_serie(d, "Ejemplo 5. Nivel del lago Hurón, 1875-1972"), "ej05-tendencia-lineal-lakehuron-serie.png")
cg <- correlograma(d)
.guardar_fig(cg$grafico, "ej05-tendencia-lineal-lakehuron-correlograma.png", alto = 6)
.linea(ljung_box(cg$acf, n_total, cg$m, 0))

fit <- ajustar_tendencia(est$y, "lineal")
print(as.data.frame(fit$parametros$tabla), digits = 5)
cat(sprintf("R2 = %s, sigma2 = %s, DW de los residuos = %s, L (HAC) = %d\n", .fmt_num(fit$parametros$r2),
            .fmt_num(fit$parametros$sigma2), .fmt_num(fit$parametros$dw), fit$parametros$L))
.guardar_fig(.grafico_residuos_ajustados(est$y - fit$parametros$residuos, fit$parametros$residuos,
                                         "Ejemplo 5. Residuos contra valores ajustados"),
             "ej05-tendencia-lineal-lakehuron-residuos.png")
pron <- fit$pronosticar(h); ref <- rep(est$y[T_est], h)
m_in <- .medir(est$y, fit$yhat, est$y); m_val <- .medir(val$y, pron, est$y); m_ref <- .medir(val$y, ref, est$y)
print(.tabla_medidas(m_in, m_val, m_ref, "Ingenuo"), digits = 5, row.names = FALSE)
cot <- .cotas("ej05", sum(!is.na(fit$errores)))
ve <- validar_errores(fit$errores, fit$n_param, T_est, cot$dL, cot$dU)
cat(sprintf("Validación de errores (N = %d):\n", ve$n))
for (p in list(ve$t_media, ve$ljung_box, ve$jarque_bera, ve$durbin_watson)) .linea(p)
.guardar_fig(ve$grafico, "ej05-tendencia-lineal-lakehuron-errores.png", alto = 7)
.guardar_fig(.grafico_pronostico(d, fit$yhat, pron, ref, "Ejemplo 5. Tendencia lineal: ajuste y pronóstico"),
             "ej05-tendencia-lineal-lakehuron-pronostico.png")
resumen[[5]] <- .fila_resumen(5, "LakeHuron", fit, m_val, m_ref)
