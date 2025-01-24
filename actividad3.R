#Actividad 3 - Equipo 9 Lote 7

library(gtsummary)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)



setwd('~/Documents/Màster Bioinformàtica/Assignatures/Estadística y R/taller grupal/mubio02_act3/')
getwd()

datos <- read.csv("Dataset expresi¢n genes.csv")

head(datos)
str(datos)
summary(datos)

is.na(datos)
colSums(is.na((datos)))
rowSums(is.na(datos))

# 2. Realizar una pca
# Cargar librerías
library(stats)
install.packages("factoextra")
library(factoextra)

#Eliminar las columnas X e id del dataset:
datos <- datos %>%
  select(-X, -id)

# Seleccionar únicamente las columnas de expresión génica
genes <- datos %>%
  select(starts_with("AQ_"))

# Estandarizar los datos
genes_scaled <- scale(genes)

#################### 1 -  Realizar el PCA ###############

pca <- prcomp(genes_scaled, center = TRUE, scale. = TRUE)

# Resumen del PCA
summary(pca)

# Extraer la varianza explicada: Scree plot
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100))

# El factor 1 explica el 52% de la varianza y el factor 2 ya solo nos explica el 6.5%.
# Cogeríamos el factor 1 y 2

# Gráfico de las variables (cargas) con colores según su contribución
library(factoextra)
fviz_pca_var(pca, 
             col.var = "contrib",  # Colorear las variables según su contribución
             gradient.cols = c("blue", "yellow", "red"),
             title = "Contribución de las Variables a los Componentes Principales")

#Agrupación por clustering de variables:

# Aplicar k-means con k=3 (puedes ajustar el número de clústeres)
set.seed(123)  # Fijar semilla para reproducibilidad
kmeans_result <- kmeans(pca$x[, 1:2], centers = 3)  # Usamos los primeros dos componentes para clustering

# Gráfico de las observaciones (puntuaciones) con los clústeres
fviz_pca_ind(pca, 
             col.ind = factor(kmeans_result$cluster),  # Convertir a factor para colores discretos
             palette = "jco",  # Paleta de colores
             addEllipses = TRUE,  # Agregar elipses de confianza para los clústeres
             title = "Observaciones Agrupadas por Clústeres en el PCA")





#-----------tablas (OPTATIVAS)----------------------

# TABLA - Varianza explicada por cada componente
var_exp <- data.frame(
  Componente = paste0("PC", 1:length(pca$sdev)),
  R2 = round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 2)
)

# Crear la tabla con kable (puedes ajustar el formato)
library(knitr)
kable(var_exp, caption = "Varianza explicada por cada componente principal")



# TABLA - Cargas de las variables
cargas <- as.data.frame(pca$rotation)

# Seleccionar las cargas para los primeros componentes (ajustar si necesitas más componentes)
cargas_tabla <- cargas[, 1:2]  # Componente 1 y 2 como ejemplo
cargas_tabla$Variable <- rownames(cargas)

# Reorganizar columnas para la tabla
cargas_tabla <- cargas_tabla %>%
  select(Variable, PC1, PC2)  # Renombra PC1, PC2 según corresponda

# Crear tabla
kable(cargas_tabla, caption = "Cargas de las variables en los componentes principales")



# ------------------------ 2. Gráficos descriptivos PCA ---------------------------

#Gráficos descriptivos de los PCA

# Cargar las librerías necesarias
library(dplyr)
library(gtsummary)

# Crear terciles para los componentes principales
datos <- datos %>%
  mutate(
    PC1_tercile = cut(pca$x[, 1], breaks = quantile(pca$x[, 1], probs = c(0, 0.33, 0.66, 1)),
                      labels = c("t1", "t2", "t3"), include.lowest = TRUE),
    PC2_tercile = cut(pca$x[, 2], breaks = quantile(pca$x[, 2], probs = c(0, 0.33, 0.66, 1)),
                      labels = c("t1", "t2", "t3"), include.lowest = TRUE)
  )

# Seleccionar solo las variables de genes y los terciles
genes <- datos %>% select(starts_with("AQ_"), PC1_tercile, PC2_tercile)

# Generar tabla descriptiva para el Componente 1
tbl_PC1 <- genes %>%
  tbl_summary(
    by = PC1_tercile,  # Agrupar por terciles del Componente 1
    statistic = all_continuous() ~ "{mean} ({sd})",  # Media y desviación estándar
    digits = all_continuous() ~ 1,
    missing = "no"  # Excluir valores faltantes
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test",  # Prueba de Kruskal-Wallis
    pvalue_fun = function(x) sprintf("%.3f", x)  # Formato decimal estándar para valores p
  )

# Generar tabla descriptiva para el Componente 2
tbl_PC2 <- genes %>%
  tbl_summary(
    by = PC2_tercile,  # Agrupar por terciles del Componente 2
    statistic = all_continuous() ~ "{mean} ({sd})",  # Media y desviación estándar
    digits = all_continuous() ~ 1,
    missing = "no"  # Excluir valores faltantes
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test",  # Prueba de Kruskal-Wallis
    pvalue_fun = function(x) sprintf("%.3f", x)  # Mostrar valores p con 3 decimales
  )

# Combinar ambas tablas
tbl_combinadas <- tbl_merge(
  tbls = list(tbl_PC1, tbl_PC2),
  tab_spanner = c("**Component 1**", "**Component 2**")
)

# Añadir el título a la tabla
tbl_combinadas <- tbl_combinadas %>%
  modify_caption("**Tabla descriptiva de genes agrupados por terciles**")

tbl_combinadas

# EN EL RMARKDOWN las tablas de percentiles estan expresadas en notación científica (que me parece más descriptivo)



######## 3 - Regresión logística ##########

# Observando el dataset vemos que no existe la columna metastasis pero en su lugar existe extensión:
colnames(datos)
datos$extension

# Crear variable dependiente (metástasis: sí/no) en función de la variable extensión:
datos <- datos %>%
  mutate(metastasis = ifelse(extension == "metastasico", 1, 0))

# Verificar la nueva columna
table(datos$metastasis)

# Modelo de regresión logística
log_model <- glm(metastasis ~ PC1_tercile + edad + trat + tumor, 
                 data = datos, family = binomial)

# Resumen del modelo
summary(log_model)

# Crear tabla con gtsummary
tbl_reg <- tbl_regression(log_model, 
                          exponentiate = TRUE,
                          conf.level = 0.95)

tbl_reg






# Evaluar calidad del modelo
#library(pROC)
#roc_curve <- roc(data$metastasis, predict(log_model, type = "response"))
#auc(roc_curve)  # Área bajo la curva
#plot(roc_curve)




