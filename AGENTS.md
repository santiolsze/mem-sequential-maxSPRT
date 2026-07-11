# Objetivo del repo

Este repositorio corresponde a un trabajo universitario.

El objetivo es presentar el paper `paper.pdf` en una presentacion hecha con Typst y realizar experimentacion en R para comprobar resultados por cuenta propia sobre un dataset publico. La experimentacion debe servir tambien como material demostrable para mostrar en clase.

# Memoria del proyecto

Usar `MEMORY.md` como memoria persistente del proyecto para transferir contexto entre conversaciones de agente.

Si `MEMORY.md` existe, leerlo antes de tomar decisiones sobre el estado del trabajo, experimentos, presentacion, pendientes o proximos pasos. Actualizarlo cuando haya decisiones relevantes, hallazgos importantes o cambios de rumbo que deban sobrevivir entre conversaciones.

## Resumen del paper

Kulldorff y colaboradores proponen el *Maximized Sequential Probability Ratio Test* (MaxSPRT) para la vigilancia posterior a la comercialización de medicamentos y vacunas. El problema surge porque los ensayos clínicos previos a la aprobación suelen tener muestras pequeñas y poblaciones seleccionadas: por ello pueden no detectar eventos adversos raros, graves o concentrados en subgrupos. El desafío estadístico es alertar pronto sobre un posible exceso de riesgo sin generar muchas falsas alarmas por examinar repetidamente los datos acumulados.

El punto de partida es el SPRT clásico de Wald. Para datos Poisson, bajo la hipótesis nula el número acumulado de eventos tiene una media esperada que incorpora población expuesta y riesgo basal; bajo la alternativa, esa media se multiplica por un riesgo relativo (RR) fijado de antemano. En cada momento se calcula la razón de verosimilitudes y se compara con fronteras superior e inferior. Sin embargo, su desempeño depende de forma crítica del RR elegido para la alternativa. Si el riesgo real no coincide con el especificado, el método puede tardar demasiado en detectar una señal o incluso aceptar la nula aunque exista un incremento relevante.

Los autores ilustran esa dificultad con datos históricos de la Vaccine Safety Datalink, usados para simular vigilancia semanal después de la vacuna Pediarix. Para fiebre dentro de los 28 días posteriores, un SPRT clásico con alternativa RR = 2 puede concluir rápidamente que no hay evidencia de ese aumento, mientras que con RR = 1.2 detecta una señal. Pero son poco satisfactorias para farmacovigilancia, donde el tamaño real del exceso de riesgo se desconoce antes de iniciar la vigilancia.

El MaxSPRT reemplaza esa alternativa puntual por una alternativa compuesta, RR > 1. En el modelo Poisson, maximiza la verosimilitud sobre todos los riesgos relativos mayores que uno. Cuando los eventos observados superan a los esperados, el estimador es RR = observado/esperado y el estadístico logarítmico resulta de comparar ese máximo con la nula; de lo contrario, vale cero. Se emite una señal cuando alcanza un valor crítico único antes de un límite preespecificado de vigilancia, expresado como número esperado de eventos bajo la nula. A diferencia del SPRT clásico, no hay una frontera inferior: si no aparece señal, el seguimiento termina al llegar a ese límite.

Los valores críticos del MaxSPRT se calibran para mantener el nivel de significación deseado pese a las múltiples revisiones. Para Poisson, los autores los obtienen mediante cálculos numéricos iterativos exactos, aprovechando que una señal solo puede surgir cuando ocurre un evento. Presentan tablas para distintos niveles de alfa y horizontes de vigilancia. Un horizonte mayor aumenta la potencia para RR pequeños, pero también exige un umbral más alto y puede prolongar el seguimiento. El trabajo reporta además potencia, tiempo esperado hasta la señal y duración media de la vigilancia.

También se desarrolla una versión binomial para situaciones sin una estimación fiable del número esperado de eventos. Allí se comparan períodos o personas expuestas con controles no expuestos emparejados, y se condiciona al total de eventos. El límite se expresa como el total de eventos y los valores críticos se derivan mediante una cadena de Markov. Esta versión resulta útil cuando faltan estimaciones previas del riesgo basal en la población expuesta. En Pediarix, el MaxSPRT detectó fiebre a las 13 semanas y síntomas neurológicos a las 42, evitando elegir un RR alternativo arbitrario. Los autores enfatizan que una señal requiere investigación epidemiológica posterior: puede reflejar riesgo real, cambios de codificación, diagnóstico o sesgos de diseño. La elección de controles históricos, contemporáneos o autocontrolados introduce distintos confusores; por ello, MaxSPRT es una herramienta de alerta temprana, no una prueba definitiva de daño clínico.
