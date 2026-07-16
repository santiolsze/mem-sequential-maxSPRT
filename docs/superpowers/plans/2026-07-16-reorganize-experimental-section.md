# Reorganize Experimental Section Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the two mixed dataset slides with three focused slides that prepare the audience before switching to the Shiny application.

**Architecture:** Present sources and targets first, comparator construction second, and interpretation limits third. Move the Shiny separator after these slides so it becomes the handoff from the deck to the live application.

**Tech Stack:** Typst, Touying, current Shiny data definitions.

## Global Constraints

- Use the current 2016--2025 FAERS/openFDA and VAERS series.
- Keep observed and expected quantities explicitly at report level.
- Explain the non-ELIQUIS and MNQ comparators separately.
- Mention the fixed binomial MENB--MNQ panel.
- Do not claim incidence, clinical risk, or causality.

---

### Task 1: Replace and reorder the experimental slides

**Files:**
- Modify: `presentation/maxsprt.typ`

- [ ] **Step 1:** Create a source/target overview slide without accumulated result totals.
- [ ] **Step 2:** Create a comparator-construction slide with the two expected-count formulas and the binomial-panel distinction.
- [ ] **Step 3:** Create an interpretation slide contrasting spontaneous reports with the VSD data in the paper.
- [ ] **Step 4:** Move the Shiny separator after the three slides.
- [ ] **Step 5:** Compile and inspect all three slides and both neighboring separators.
