# Estructura del proyecto

Ver [AGENTS.md](AGENTS.md) para el objetivo del repo y el resumen del paper, y
[MEMORY.md](MEMORY.md) como memoria persistente de decisiones y hallazgos.

## Presentacion (Typst)

`presentation/`

- `maxsprt.typ`: fuente de la presentacion (deck de clase sobre el paper).
- `deck_setup.typ`: configuracion/tema del deck.
- `maxsprt.pdf`: PDF compilado de la presentacion.
- `make_figs.py`: genera las figuras en `figures/` a partir de los resultados
  de `experimentacion/`.
- `figures/`: imagenes usadas por el deck (`paper_figure_1_signals.png`,
  `sprt_fiebre.png`, `sprt_neuro.png`, `sprt_trayectoria.png`).

El paper original esta en la raiz: `paper.pdf`.

## App Shiny (demo interactiva)

`experimentacion/shiny/`

- `app.R`: UI y server de la app. Se lanza desde la raiz del repo con
  `shiny::runApp("experimentacion/shiny", host = "127.0.0.1", port = 3838)`.
- `R/data.R`: carga y prepara los datasets Parquet locales para los paneles de
  datos reales.
- `R/sequential.R`: implementacion de MaxSPRT/SPRT (Poisson y binomial) usada
  por los paneles de simulacion.
- `www/app.css`: estilos.
- `tests/`: suite `testthat`, se corre con
  `Rscript experimentacion/shiny/tests/run_tests.R`.
- `README.md`: detalle de los 4 paneles (Datos reales, Diagnostico Poisson,
  Potencia simulada, Binomial simulado) y sus limitaciones de interpretacion.

## Generacion de datos procesados

Datos crudos en `experimentacion/data/raw/` (VAERS: ZIPs anuales 2016-2025;
openFDA: metadata + descargas diarias), datos procesados en
`experimentacion/data/processed/` (Parquet, particionado). Ver
`experimentacion/data/README.md` para el detalle completo de fuentes, esquema
de join de VAERS y particiones.

Pipeline, en orden:

1. `experimentacion/scripts/fetch_openfda_eliquis_comparator.R`: descarga /
   refresca el agregado mensual openFDA ELIQUIS/Haemorrhage (unico agregado
   local para ese par producto/evento) a partir de la API `api.fda.gov`.
2. `experimentacion/scripts/build_parquet_datasets.py`: convierte los datos
   crudos locales (openFDA + VAERS) a Parquet comprimido (`zstd`), genera las
   particiones (`year`, `month`, `vax_type_slug`, `medicine_slug`, etc.) y
   escribe el manifiesto `experimentacion/data/processed/parquet_manifest.json`
   (row counts, columnas, particiones, compresion).
3. `experimentacion/scripts/replicate_paper_tables_1_3.R`: replicacion exacta
   y deterministica en base R de las Tablas 1-3 del paper (valores criticos,
   potencia, tiempos de senal/vigilancia). No depende de los datos VAERS/openFDA;
   escribe resultados en `experimentacion/results/paper_tables/`.

La app Shiny (`R/data.R`) solo lee los Parquet ya generados en
`experimentacion/data/processed/`; no descarga ni procesa datos en runtime.
