# Binomial MaxSPRT Alpha Replication Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce exactly in base R the binomial MaxSPRT type-I error calibration and explain every Markov-chain step for the paper's `z = 1`, `N = 10`, `alpha = 0.05` example.

**Architecture:** A self-contained sourceable script will define the one-sided binomial LLR, exact absorbing-chain propagation, candidate-boundary calibration, and a printable walkthrough. A focused `testthat` file will verify probability conservation, direct path enumeration for `N = 10`, conservativeness, and the Table 4 critical value `2.77259`.

**Tech Stack:** Base R, `testthat` only for automated verification.

## Global Constraints

- Do not depend on Shiny or modify active Shiny functions.
- Calculate alpha exactly; do not use Monte Carlo simulation.
- Treat `LLR >= V` as rejection and stop propagating rejected paths.
- Preserve the paper's discreteness: choose a reachable LLR threshold with actual alpha at most nominal alpha.

---

### Task 1: Exact binomial Markov chain

**Files:**
- Create: `experimentacion/scripts/replicate_binomial_alpha.R`
- Test: `experimentacion/tests/testthat/test-replicate-binomial-alpha.R`

**Interfaces:**
- Produces: `binomial_maxsprt_llr_exact(c, n, z)`, `binomial_alpha_markov(V, N, z, keep_states)`, and `calibrate_binomial_boundary(alpha, N, z)`.

- [x] Write tests for the LLR, probability conservation, and agreement between the Markov recursion and all `2^10` event paths.
- [x] Run the focused test and confirm it fails because the script does not yet exist.
- [x] Implement the one-sided LLR and exact two-branch transition from `(n,c)` to `(n+1,c)` or `(n+1,c+1)`.
- [x] Run the focused test and confirm it passes.

### Task 2: Critical value and pedagogical walkthrough

**Files:**
- Modify: `experimentacion/scripts/replicate_binomial_alpha.R`
- Modify: `experimentacion/tests/testthat/test-replicate-binomial-alpha.R`

**Interfaces:**
- Produces: a calibration result containing `critical_value`, `actual_alpha`, candidate table, and step history; direct execution prints the derivation and writes CSV tables under `experimentacion/results/binomial_alpha/`.

- [x] Add a failing test requiring critical value `2.77259` for `z = 1`, `N = 10`, nominal alpha `0.05`.
- [x] Enumerate all reachable positive LLR values, evaluate exact crossing probability immediately above each, and select the smallest threshold with actual alpha no greater than nominal alpha.
- [x] Print inputs, formulas, transitions, per-event absorbed probability, survival probability, accumulated alpha, and the Table 4 comparison.
- [ ] Run the focused test, full experiment test suite, and the script; confirm all pass and inspect the generated CSV values.
