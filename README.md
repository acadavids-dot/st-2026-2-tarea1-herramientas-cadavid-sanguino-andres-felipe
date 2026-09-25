https://github.com/acadavids-dot/st-2026-2-tarea1-herramientas-cadavid-sanguino-andres-felipe

# Tarea 1 — Caja de herramientas de pronóstico

Series de Tiempo (3009297), UNAL Medellín, 2026-II. Prof. Juan Pablo Valencia Arango.

## Qué contiene el repositorio

| Archivo / carpeta        | Contenido                                                                                                                                                                                            | Dependencias                                                                                                     |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `README.md`            | este archivo                                                                                                                                                                                         | —                                                                                                               |
| `.gitignore`           | archivos y carpetas que no se versionan (`.Rhistory`, `.RData`, `.Rproj.user/`, `informe/informe_files/`, `informe/.quarto/`, etc.)                                                        | —                                                                                                               |
| `R/00-lectura.R`       | `leer_serie()`                                                                                                                                                                                     | base (`stats`, `utils`), `tibble`                                                                          |
| `R/01-graficos.R`      | `graficar_serie()`, `correlograma()`, tema `.tema_tarea()` y los auxiliares de gráficos (`.grafico_pronostico`, `.grafico_residuos_ajustados`)                                            | `ggplot2`, `patchwork`, `stats` (`pacf`, `qnorm`)                                                      |
| `R/02-metodos.R`       | los ocho`ajustar_*()` y `optimizar()`                                                                                                                                                            | base,`tibble`, `ggplot2` (gráficos de `optimizar()`), `stats` (`pt`)                                  |
| `R/03-evaluacion.R`    | `medidas()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()`, `validar_errores()` y sus auxiliares (`.reportar_prueba`, `.medir`, `.tabla_medidas`, `.cotas`, `.contraste_rh`) | base,`tibble`, `ggplot2`, `patchwork`, `stats`                                                           |
| `ejemplos/ejemplos.R`  | carga`R/`, corre el Bloque 0 de verificaciones contra R base y los 9 ejemplos (8 + contraejemplo), guarda las figuras y `sesion-info.txt`                                                        | los anteriores;`stats::acf()`, `stats::lm()` y `stats::Box.test()` **solo** en el Bloque 0           |
| `informe/informe.qmd`  | fuente del informe: introducción y convenciones, verificaciones contra R, los 9 ejemplos, tabla resumen, referencias y`sessionInfo()`                                                             | los anteriores;`knitr`/`rmarkdown`/Quarto únicamente como herramienta de render, nunca dentro de un método |
| `informe/informe.html` | informe renderizado, autocontenido (`embed-resources`)                                                                                                                                             | —                                                                                                               |
| `figs/`                | 44 figuras que genera`ejemplos.R` (serie, correlograma, optimización, errores y pronóstico de cada ejemplo)                                                                                      | —                                                                                                               |
| `sesion-info.txt`      | `sessionInfo()` de la corrida de cierre de `ejemplos.R`                                                                                                                                          | —                                                                                                               |

De los paquetes permitidos por el enunciado (`dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `patchwork`), este trabajo
solo necesitó `tibble`, `ggplot2` y `patchwork`; `dplyr`, `tidyr` y `purrr` no se usan en ningún archivo.

## Cómo se corre

Desde la raíz del repositorio, con un R que tenga instalados `tibble`, `ggplot2` y `patchwork`:

```bash
Rscript --vanilla -e 'source("ejemplos/ejemplos.R")'
```

Corre de principio a fin sin intervención, imprime el Bloque 0 de verificaciones (se detiene con `stopifnot()` si
alguna supera su tolerancia), los 9 ejemplos y el resumen final, y deja 44 figuras en `figs/` y la corrida de
`sessionInfo()` en `sesion-info.txt`. Tiempo observado: **23,6 segundos** (R 4.6.0, Windows 11, máquina de esta
entrega).

Para volver a renderizar el informe (Quarto necesita un R con esos mismos paquetes visible en el `PATH`):

```bash
quarto render informe/informe.qmd
```

## Cómo se usan las funciones

Ejemplo mínimo: leer una serie, ajustar un método y pedir un pronóstico extramuestral.

```r
invisible(lapply(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)), source))

d   <- leer_serie(Nile, fuente = "Durbin J, Koopman SJ (2001). Time Series Analysis by State Space Methods. Oxford University Press",
                  unidad = "Caudal anual del Nilo en Asuán (10^8 m^3)")
fit <- ajustar_ses(d$y, alpha = 0.24)
fit$pronosticar(5)   # Yhat_{T+1}, ..., Yhat_{T+5}
```

`d` es un `tibble` con `t`, `fecha` e `y`, y los atributos `frecuencia`, `fuente` y `unidad`. `fit` es una lista de
clase `"metodo_pronostico"` con `yhat` (el ajuste dentro de muestra, con el calentamiento en `NA`), `errores`,
`parametros` y la función `pronosticar(h)`; ninguna de las funciones de `R/` lee objetos globales. Para elegir la
constante o la ventana en vez de fijarla a mano: `opt <- optimizar(d$y, "ses")` y `opt$optimo$alpha`.

## Convenciones que fijan los números

Cuando el enunciado y las notas de clase no fijaban una decisión, se tomó una y se documentó en el encabezado de la
función correspondiente y aquí. Estas son las que determinan las cifras de `ejemplos.R` y del informe.

- **Tolerancia de las equivalencias numéricas contra R:** error **relativo** a $\max|y|$, con cota $10^{-12}$ (no
  absoluto). Con datos del orden de $10^4$ (p. ej. `austres`) un ulp ya mide cerca de $1,8\times10^{-12}$, así que una
  cota absoluta de $10^{-12}$ no la cumpliría ninguna implementación correcta (diagnóstico en `scratch/diag-tol.R`,
  fuera del repo).
- **ACF a mano:** divisor único $T$ (no $T-h$); verificada contra `acf()` con diferencia relativa máxima del orden de
  $10^{-15}$. Banda de ruido blanco $\pm z_{0{,}975}/\sqrt{n}$, con $n$ el número de valores graficados (los errores,
  no la serie, cuando se grafican errores). $m$ por defecto $\min\{\lfloor n/4\rfloor, 24\}$.
- **Ljung–Box:** $Q_m$ se refiere a $\chi^2_{m-p}$, con $p$ el número de parámetros estimados por el método (0 en
  media simple; 1 en media móvil, SES y doble media móvil; 2 en tendencia lineal y exponencial y en Holt; 3 en
  tendencia cuadrática); sobre la serie original siempre $p=0$. Verificada contra `Box.test(..., fitdf = p)`.
- **Durbin–Watson:** sin valor p exacto; se decide con las cotas $d_L$ y $d_U$ al 5 % de Savin y White (1977). Como
  la tabla no cubre todos los $N$ de este trabajo, las cotas se calcularon de forma exacta (definición de Durbin y
  Watson 1951 más la integral de Imhof 1961) y ese cálculo reproduce 12 filas de la tabla con diferencias de a lo
  sumo $0{,}005$; se redondean a dos decimales, como la tabla. Para los métodos de suavizamiento se usa $k'=1$
  (aproximación declarada: no hay matriz de diseño).
- **MASE:** la escala es el error absoluto medio del referente ingenuo de un paso, calculada **en el tramo de
  estimación** (nunca en validación); $s=1$ en series no estacionales y $s=$ frecuencia en las estacionales
  (`JohnsonJohnson`, $s=4$; `AirPassengers`, $s=12$).
- **MAPE:** si algún valor evaluado es cero, `medidas()` devuelve `NA` con una advertencia explícita en vez de omitir
  el dato en silencio (ocurre en `discoveries`).
- **`optimizar()`:** minimiza el MSE de un paso **dentro del tramo de estimación**, con origen común: el MSE de
  todas las filas de la rejilla se calcula desde el primer período en que la fila más exigente ya produce
  pronóstico ($t=13$ en media móvil, $t=24$ en doble media móvil, $t=2$ en SES y Holt). En caso de empate gana la
  primera fila de la rejilla.
- **Tendencias:** `yhat` son los valores ajustados de una regresión sobre toda la muestra de estimación, no
  pronósticos con información hasta $t-1$; por eso su MSE dentro de muestra es optimista frente al de los métodos
  de suavizamiento. El error estándar robusto es HAC de Newey–West a mano, núcleo de Bartlett,
  $L=\lfloor 4\,(T/100)^{2/9}\rfloor$ rezagos, **sin** corrección por grados de libertad; los valores p usan la $t$
  de Student con $T-p$ grados de libertad. Coeficientes verificados contra `coef(lm(...))` (diferencia menor a
  $10^{-8}$). La tendencia exponencial se estima sobre $\ln Y_t$: $e^{\hat a+\hat\theta t}$ estima la mediana
  condicional, no la media; `corregir_sesgo = TRUE` multiplica por $e^{\hat\sigma^2_{\ln}/2}$ (no la corrección de
  Duan de la Clase 4).
- **Holt:** verificado contra su forma de corrección de error ($L_t = L_{t-1}+\hat T_{t-1}+\alpha e_t$,
  $\hat T_t = \hat T_{t-1}+\alpha\beta e_t$) con diferencia relativa menor a $10^{-12}$.
- **Referente:** el pronóstico ingenuo $\hat Y_{T+h}=Y_T$, o el ingenuo estacional (el último ciclo repetido) en las
  series estacionales, calculado sobre el mismo tramo y el mismo horizonte $h$ que el método.

## Resumen de resultados

MASE de validación de cada método frente a su referente, sobre el mismo tramo y el mismo horizonte $h=\min\{12,
\lfloor 0{,}2\,T\rfloor\}$ (ver el detalle de cada ejemplo en el informe):

| # | Serie              | Método                    | Parámetros                                | MASE método | MASE referente |
| - | ------------------ | -------------------------- | ------------------------------------------ | ------------ | -------------- |
| 1 | `discoveries`    | Media simple               | —                                         | 0,9192       | 1,1365         |
| 2 | `Nile`           | Media móvil               | $k=2$ (borde)                            | 0,8467       | 0,8355         |
| 3 | `Nile`           | SES                        | $\alpha=0{,}24$                          | 0,8062       | 0,8355         |
| 4 | `austres`        | Doble media móvil         | $k=2$ (borde)                            | 1,8632       | 6,0912         |
| 5 | `LakeHuron`      | Tendencia lineal           | —                                         | 2,0667       | 2,1642         |
| 6 | `airmiles`       | Tendencia cuadrática      | —                                         | 1,2055       | 4,4959         |
| 7 | `JohnsonJohnson` | Tendencia exponencial | —                                         | 3,5142       | 6,5286         |
| 8 | `austres`        | Holt                       | $\alpha=0{,}95$ (borde), $\beta=0{,}8$ | 1,4842       | 6,0912         |
| C | `AirPassengers`  | SES (contraejemplo)        | $\alpha=0{,}98$ (borde)                  | 2,5143       | 1,5709         |

Siete de los nueve casos tienen un MASE del método menor que el de su referente. Los dos que no: la media móvil del
ejemplo 2, por un margen mínimo (0,8467 contra 0,8355; el patrón de `Nile` igual justifica el método, como se
argumenta en el informe), y el contraejemplo, a propósito, por el margen más amplio (2,5143 contra 1,5709): es la
firma numérica de aplicar un suavizamiento sin tendencia ni estacionalidad a una serie que tiene las dos cosas.

## Declaración de uso de IA

Este trabajo se desarrolló de forma conjunta entre el estudiante y Claude (Anthropic, modelo Sonnet 5), utilizado como
asistente de programación en sesiones de Claude Code sobre este repositorio.

**Punto de partida.** El proyecto partió de una propuesta inicial aportada por el estudiante: borradores de las
funciones de `R/`, de `ejemplos.R`, del informe y del README, y la especificación del proyecto (`CLAUDE.md`).

**Reparto del trabajo.** El estudiante dirigió el proceso: aprobó el plan de cada fase antes de la implementación,
tomó las decisiones de entrega (repositorio público, fecha confirmada con el profesor), revisó y modificó partes del
código, solicitó una auditoría completa del código y redactó las interpretaciones de los ejemplos 1 a 6 del informe.
Claude revisó, corrigió, reescribió y amplió el código de la propuesta inicial, ejecutó las verificaciones contra las
funciones de R, integró las interpretaciones del estudiante al informe convirtiendo sus cifras en código en línea,
redactó los demás textos del informe (las interpretaciones de los ejemplos 7 y 8 y del contraejemplo, y los textos
descriptivos de todos los ejemplos) y explicó el código al estudiante para que pueda entenderlo y sustentarlo.

| Fase | Trabajo realizado | Qué se verificó |
| ---- | ----------------- | --------------- |
| 0 | Estructura del repositorio, `.gitignore` y README esqueleto | Comparación contra la estructura exigida |
| 1 | `leer_serie()`, `graficar_serie()` y `correlograma()` con sus auxiliares | Fechas de las 15 series contra aritmética de `start()`/`end()`; ACF a mano contra `acf()`; banda idéntica a `plot.acf()`; casos de error |
| 2 | `ljung_box()`, `jarque_bera()`, `durbin_watson()`, `medidas()` y `validar_errores()` con sus auxiliares | $Q_m$ contra `Box.test()`; JB y DW contra `tseries`/`lmtest` (solo en `scratch/`); reproducción de las cifras de la Clase 3 |
| 3 | Los ocho métodos y `optimizar()` | Cifras de las Clases 3 y 4; equivalencias de los métodos con formulaciones de referencia; tendencias contra `lm()`; HAC contra `sandwich::NeweyWest` (solo en `scratch/`) |
| 4 | Revisión de las series candidatas con gráficos y correlogramas; cálculo de las cotas $d_L$/$d_U$ | Gráficos y correlogramas de las series candidatas; selección final y cotas contra la tabla de Savin y White |
| 5 | `ejemplos/ejemplos.R` con el Bloque 0, los 9 ejemplos y las figuras | `source("ejemplos/ejemplos.R")` corre limpio; 14 verificaciones del Bloque 0; figuras revisadas |
| 6 | Informe completo; interpretaciones de los ejemplos 1 a 6 redactadas por el estudiante | Renderizado completo; contraste del texto con las cifras reales de la corrida; toda cifra del informe sale de código |

Las verificaciones se ejecutaron en las sesiones de trabajo con Claude y el estudiante revisó sus resultados.

**Responsabilidad.** El estudiante revisó el trabajo, es responsable de su contenido y puede explicar cada función y
cada ejemplo. La responsabilidad final sobre las interpretaciones presentadas, su comprensión y su explicación en la
sustentación corresponde al estudiante.
