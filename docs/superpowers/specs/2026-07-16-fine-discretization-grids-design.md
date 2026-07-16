# Grillas finas para el experimento de discretización

## Objetivo

Extender el selector `Cantidad de looks K` con las opciones 25.000 y 50.000 para observar mejor la convergencia hacia la vigilancia continua.

## Diseño

- Agregar `25000` y `50000` como opciones no seleccionadas por defecto.
- Mantener el comportamiento actual para grillas con máximo menor que 25.000.
- Cuando alguna grilla seleccionada tenga `K >= 25000`, limitar automáticamente las repeticiones efectivas a 5.000.
- Informar en la interfaz cuando el límite reduzca la cantidad solicitada, mencionando la grilla fina seleccionada y las 5.000 repeticiones efectivas.
- Actualizar la documentación para reflejar el nuevo máximo y la regla de carga.

## Verificación

- Probar que ambas opciones aparecen en la interfaz y no están seleccionadas inicialmente.
- Probar que `K = 25000` y `K = 50000` aplican el límite de 5.000 repeticiones.
- Probar que las grillas menores conservan la cantidad solicitada y que la regla existente para `K >= 10000` sigue siendo compatible.
- Ejecutar la suite de pruebas de la aplicación Shiny.
