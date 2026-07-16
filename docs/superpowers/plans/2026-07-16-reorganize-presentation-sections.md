# Reorganize Presentation Sections Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Organize the Typst presentation into seven clearly signposted narrative sections.

**Architecture:** Extend the existing section-divider helper with a compact section counter, insert four missing dividers, and replace the three existing dividers with audience-facing titles. Keep every technical content slide in its current order.

**Tech Stack:** Typst, Touying, project-local deck theme.

## Global Constraints

- Use seven narrative blocks: context, classical SPRT, Poisson MaxSPRT, binomial variant, threshold calculation, discrete looks, and real-data demo.
- Preserve all existing technical slide content and unrelated working-tree changes.
- Avoid generic divider titles such as `Sección` and transitional language such as `Volvemos al caso Poisson`.
- Verify every new divider visually after compiling.

---

### Task 1: Add seven narrative dividers

**Files:**
- Modify: `presentation/maxsprt.typ`

**Interfaces:**
- Consumes: `separador_seccion`, existing slide order, and the deck theme.
- Produces: seven section slides numbered `01 / 07` through `07 / 07`.

- [ ] **Step 1:** Add an `etapa` argument and a compact counter treatment to `separador_seccion`.
- [ ] **Step 2:** Insert or replace dividers at the seven approved narrative boundaries.
- [ ] **Step 3:** Compile with `typst compile presentation/maxsprt.typ presentation/maxsprt.pdf`; expect exit status 0.
- [ ] **Step 4:** Render and inspect all seven divider slides for wrapping, overflow, and hierarchy.
- [ ] **Step 5:** Recompile after any visual corrections and confirm the final page count.
