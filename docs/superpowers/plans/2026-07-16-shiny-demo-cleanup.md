# Shiny Demo Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Simplificar la demo quitando el diagnóstico Poisson y extender el experimento de discretización hasta 50.000 looks con una carga acotada.

**Architecture:** La interfaz y sus salidas se mantienen en `app.R`; la política de repeticiones efectivas permanece encapsulada en `discretization_effective_reps()`. Las pruebas verifican tanto el HTML de la interfaz como la regla computacional sin ejecutar grillas gigantes.

**Tech Stack:** R, Shiny, testthat.

## Global Constraints

- `K = 25000` y `K = 50000` son opcionales y no quedan seleccionados por defecto.
- Toda selección con `K >= 25000` usa como máximo 5.000 repeticiones.
- Se conserva `poisson_gof()` y sus pruebas aunque se retire su pestaña.

---

### Task 1: Grillas finas

**Files:**
- Modify: `experimentacion/shiny/tests/testthat/test-discretization.R`
- Modify: `experimentacion/shiny/tests/testthat/test-app.R`
- Modify: `experimentacion/shiny/R/discretization.R`
- Modify: `experimentacion/shiny/app.R`
- Modify: `experimentacion/shiny/README.md`

**Interfaces:**
- Consumes: `discretization_effective_reps(looks, requested_reps)`.
- Produces: opciones UI 25.000/50.000 y máximo efectivo de 5.000 repeticiones.

- [ ] Agregar expectativas para las dos opciones nuevas y el límite de 5.000.
- [ ] Ejecutar las pruebas focalizadas y confirmar que fallan por faltar el comportamiento.
- [ ] Actualizar la función, el selector, el mensaje de carga y el README.
- [ ] Ejecutar nuevamente las pruebas focalizadas y confirmar que pasan.

### Task 2: Retirar diagnóstico Poisson

**Files:**
- Modify: `experimentacion/shiny/tests/testthat/test-app.R`
- Modify: `experimentacion/shiny/app.R`
- Modify: `experimentacion/shiny/README.md`

**Interfaces:**
- Consumes: estructura `navset_tab()` de la aplicación.
- Produces: app sin pestaña, guía ni salidas de diagnóstico Poisson.

- [ ] Agregar expectativas de ausencia para la pestaña y sus identificadores de salida.
- [ ] Ejecutar la prueba de interfaz y confirmar que falla por el contenido todavía presente.
- [ ] Eliminar el panel, su mención en la guía y los tres `render*` asociados.
- [ ] Actualizar el README sin borrar `poisson_gof()` ni sus pruebas unitarias.
- [ ] Ejecutar la suite completa con `Rscript experimentacion/shiny/tests/run_tests.R`.
