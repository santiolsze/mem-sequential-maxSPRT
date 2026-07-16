# Improve Presentation Introduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the generic problem slide with a concise two-slide opening that establishes the CDC's institutional importance and motivates post-market pharmacovigilance.

**Architecture:** Modify only the opening section of the existing Typst deck. Preserve the established theme and flow into the relative-risk slide, then compile and visually inspect the affected slides.

**Tech Stack:** Typst, Touying, project-local presentation sources.

## Global Constraints

- Describe CDC accurately as a major U.S. federal public-health agency.
- Do not imply that the paper itself is a CDC publication.
- State that MaxSPRT arose from a concrete surveillance need in the CDC-sponsored Vaccine Safety Datalink.
- Keep source links out of the audience-facing narrative.
- Preserve all unrelated user changes in the working tree.

---

### Task 1: Rewrite the opening

**Files:**
- Modify: `presentation/maxsprt.typ:36-41`

**Interfaces:**
- Consumes: the existing title slide and Touying theme.
- Produces: two introductory slides leading into `¿Qué es el riesgo relativo?`.

- [ ] **Step 1:** Replace the current generic problem slide with one CDC-context slide and one pharmacovigilance-motivation slide.
- [ ] **Step 2:** Compile with `typst compile presentation/maxsprt.typ presentation/maxsprt.pdf` and expect exit status 0.
- [ ] **Step 3:** Render the opening pages and inspect them for overflow, wrapping, and narrative continuity.
- [ ] **Step 4:** Correct any visual defects and compile again.
