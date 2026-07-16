# Replace Table Summary Slide Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the slide summarizing Tables 1--3 with a legible selection of Poisson MaxSPRT critical values from Table 1 only.

**Architecture:** Remove the conceptual table and Table 3 formulas. Insert one four-column native Typst table with representative horizons and the three alpha levels reproduced in the local results.

**Tech Stack:** Typst, locally replicated CSV results.

## Global Constraints

- Show only Table 1 critical values.
- Use values from `experimentacion/results/paper_tables/table_1_critical_values.csv`.
- Select rows that cover small, medium, and large horizons without overcrowding the slide.
- Preserve all unrelated content.

---

### Task 1: Replace the summary slide

**Files:**
- Modify: `presentation/maxsprt.typ`

- [ ] **Step 1:** Replace the slide title with an audience-facing description of the table.
- [ ] **Step 2:** Replace all body content with ten representative Table 1 rows and three alpha columns.
- [ ] **Step 3:** Compile and visually inspect the complete table at full size.
