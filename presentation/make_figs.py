from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(7)

FIGURES_DIR = Path(__file__).parent / "figures"

BLUE = "#2f6f9f"
DARK = "#1f2933"
GRAY = "#9aa5b1"
LIGHT = "#c9d9e6"

plt.rcParams.update({
    "font.size": 13,
    "axes.edgecolor": DARK,
    "axes.labelcolor": DARK,
    "text.color": DARK,
    "xtick.color": DARK,
    "ytick.color": DARK,
})

# --- Figure 1: two SPRT trajectories with different RR alternatives ---
fig, ax = plt.subplots(figsize=(6.5, 4))

t = np.arange(0, 40)
# simulated-looking LLR path that trends up with noise
walk = np.cumsum(np.random.normal(0.28, 0.55, size=len(t)))
walk = walk - walk[0]

ax.plot(t, walk, color=BLUE, linewidth=2.2, label="log-verosimilitud (LLR)")
ax.axhline(2.77, color=DARK, linestyle="--", linewidth=1.4)
ax.text(1, 2.95, "límite crítico (rechazo de H0)", fontsize=10, color=DARK)

# mark crossing point
cross_idx = np.argmax(walk > 2.77) if np.any(walk > 2.77) else len(t) - 1
ax.scatter([t[cross_idx]], [walk[cross_idx]], color="#b3541e", zorder=5, s=45)
ax.annotate("señal", xy=(t[cross_idx], walk[cross_idx]), xytext=(t[cross_idx]+2, walk[cross_idx]+1.3),
            fontsize=10, color="#b3541e", arrowprops=dict(arrowstyle="->", color="#b3541e"))

ax.set_xlabel("semana de vigilancia")
ax.set_ylabel("LLR acumulado")
ax.set_ylim(-2, 8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "sprt_trayectoria.png", dpi=170)
plt.close(fig)

# --- Figure 2: sensitivity to the choice of RR in classical SPRT ---
fig, ax = plt.subplots(figsize=(6.8, 4.2))

t = np.arange(0, 60)
base = np.cumsum(np.random.normal(0.10, 0.5, size=len(t)))
base -= base[0]

llr_low = base * 1.55   # behaves like a "sensitive" small-RR alternative
llr_high = base * 0.35  # behaves like a "conservative" large-RR alternative

ax.plot(t, llr_low, color=BLUE, linewidth=2.2, label="alternativa RR=1.2 (sensible)")
ax.plot(t, llr_high, color="#b3541e", linewidth=2.2, label="alternativa RR=2.0 (conservadora)")
ax.axhline(2.77, color=DARK, linestyle="--", linewidth=1.2)
ax.set_xlabel("semana de vigilancia")
ax.set_ylabel("LLR acumulado")
ax.legend(frameon=False, fontsize=10, loc="upper left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "sprt_sensibilidad.png", dpi=170)
plt.close(fig)

print("done")
