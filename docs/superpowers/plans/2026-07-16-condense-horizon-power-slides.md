# Condense Horizon and Power Slides Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove the repetitive power/trade-off slide and preserve its only new idea on the preceding threshold table slide.

**Architecture:** Add one concise power implication below the existing threshold explanation, then delete the standalone slide. Preserve all following content and section boundaries.

**Tech Stack:** Typst, Touying.

## Global Constraints

- Do not repeat the `larger T -> larger V` chain.
- Do not claim that extending surveillance costs only computation.
- Preserve the distinction between threshold calibration and power against moderate risks.

---

### Task 1: Merge the power implication into the threshold slide

**Files:**
- Modify: `presentation/maxsprt.typ`

- [ ] **Step 1:** Add one concise sentence explaining that a larger horizon helps detect more moderate risk increases but requires accumulating evidence for longer.
- [ ] **Step 2:** Delete the standalone `Potencia y el trade-off de siempre` slide.
- [ ] **Step 3:** Compile and visually inspect the modified threshold slide.
