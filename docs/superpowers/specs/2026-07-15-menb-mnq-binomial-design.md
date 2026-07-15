# Diseño: MaxSPRT binomial MENB vs MNQ para Pyrexia

## Objetivo

Agregar a la app Shiny una pestaña fija que compare MENB con MNQ para Pyrexia mediante un MaxSPRT binomial condicional. El estimando es la proporción relativa de reportes VAERS con Pyrexia, no la incidencia clínica entre personas vacunadas.

También se corregirá la presentación para dejar claro que, en el diseño binomial del paper, la probabilidad nula de que un evento caiga del lado expuesto es `1 / (z + 1)` y solo vale `1/2` cuando `z = 1`. El horizonte binomial se denotará por `N`, y la frontera se describirá como dependiente de `alpha`, `N` y `z`.

## Modelo estadístico

Para cada mes `i` se observan:

- `A_i`: total de reportes MENB;
- `B_i`: total de reportes MNQ;
- `X_i`: reportes MENB con Pyrexia;
- `Y_i`: reportes MNQ con Pyrexia;
- `n_i = X_i + Y_i`: total mensual de reportes de Pyrexia entre ambas vacunas.

Bajo igualdad de proporciones de reporte:

```text
X_i | n_i ~ Binomial(n_i, p0_i)
p0_i = A_i / (A_i + B_i)
```

La alternativa usa un único riesgo relativo de reporte `RR > 1` durante la vigilancia:

```text
p_i(RR) = RR A_i / (RR A_i + B_i)
```

En cada mes se acumula la log-verosimilitud condicional de todos los estratos observados hasta ese momento y se maximiza numéricamente sobre `RR >= 1`. El LLR vale cero cuando el máximo restringido ocurre en `RR = 1`.

## Calibración secuencial

La frontera se calibrará mediante Monte Carlo bajo `H0`, manteniendo fijos los totales mensuales `n_i` y las probabilidades `p0_i`. Cada réplica simulará `X_i ~ Binomial(n_i, p0_i)`, reconstruirá toda la trayectoria mensual del GLR y registrará su máximo.

La frontera será el menor valor empírico cuya probabilidad simulada de cruce sea a lo sumo el alfa nominal, siguiendo la convención conservadora ya usada por la app. La pestaña reutilizará el control general de alfa y tendrá una semilla fija para reproducibilidad. El texto de la interfaz describirá el control de alfa como condicional a los márgenes y al calendario mensual observados.

## Interfaz

La nueva pestaña se llamará `Binomial MENB vs MNQ` y no tendrá selectores de vacuna o síntoma. Mostrará:

1. una explicación del estimando y de la hipótesis nula;
2. métricas acumuladas de reportes totales, reportes de Pyrexia, proporciones y RR de reporte;
3. un gráfico de proporciones mensuales de Pyrexia para MENB y MNQ;
4. un gráfico de la trayectoria del LLR binomial, la frontera Monte Carlo y el primer cruce;
5. un resumen de la decisión secuencial;
6. una tabla mensual con conteos, proporciones y evidencia acumulada;
7. una advertencia de que VAERS no contiene denominadores de dosis ni permite concluir incidencia o causalidad.

El análisis usará 2016--2025 y será independiente del preset y del mes inicial elegidos en la pestaña Poisson.

## Componentes

- `R/sequential.R`: funciones puras para calcular la log-verosimilitud estratificada, maximizar el RR y calibrar la frontera.
- `R/data.R`: preparación explícita de las cuatro series MENB/MNQ necesarias.
- `app.R`: reactivos, gráficos, métricas, decisión y tabla de la pestaña.
- `tests/testthat/`: pruebas numéricas del modelo, calibración reproducible, datos reales e interfaz.
- `presentation/maxsprt.typ`: aclaraciones de `z`, del caso `z = 1` y del horizonte `N`.
- `MEMORY.md` y README de la app: registro de la decisión metodológica y de su interpretación.

## Validación

Las pruebas cubrirán:

- LLR nulo cuando los conteos coinciden con la asignación esperada;
- LLR positivo ante exceso MENB;
- equivalencia con la fórmula binomial simple cuando `p0_i` es constante;
- invariancia ante meses sin Pyrexia;
- frontera reproducible y conservadora en la simulación usada para calibrarla;
- integridad de los conteos MENB/MNQ cargados desde Parquet;
- presencia y funcionamiento básico de la pestaña Shiny;
- compilación de la presentación Typst y ejecución completa de `testthat`.

## Fuera de alcance

- estimar incidencia clínica o riesgo causal;
- generalizar la pestaña a otras vacunas o síntomas;
- usar `p0 = 1/2` sin considerar los distintos volúmenes de reportes;
- usar directamente las tablas binomiales del paper, que suponen un `z` fijo.
