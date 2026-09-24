# ejemplos.R
# Se ejecuta desde la raíz del repositorio: source("ejemplos/ejemplos.R")

invisible(lapply(sort(list.files("R", pattern = "\\.R$", full.names = TRUE)), source))
