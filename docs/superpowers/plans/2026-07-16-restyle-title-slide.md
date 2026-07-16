# Restyle Title Slide Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Improve the title slide's visual hierarchy without removing or rewriting any of its visible content.

**Architecture:** Keep the existing Touying title bar and replace the plain body paragraphs with an intentional Typst composition: prominent subtitle, accent rule, emphasized authors, and a visually subordinate full citation.

**Tech Stack:** Typst, Touying, existing deck palette and typography.

## Global Constraints

- Preserve the exact title, subtitle, author names, and full paper citation.
- Do not alter the visual language of the rest of the deck.
- Keep the citation legible but clearly secondary.
- Preserve all unrelated working-tree changes.

---

### Task 1: Restyle the title slide

**Files:**
- Modify: `presentation/maxsprt.typ`

**Interfaces:**
- Consumes: the existing title frame and project palette.
- Produces: one redesigned title slide with unchanged content.

- [ ] **Step 1:** Replace the plain title-slide body with explicit spacing, typography, an accent rule, and a citation block.
- [ ] **Step 2:** Compile `presentation/maxsprt.pdf` and expect exit status 0.
- [ ] **Step 3:** Render the title slide at full size and inspect hierarchy, wrapping, and footer clearance.
- [ ] **Step 4:** Correct any visual issues and recompile.
