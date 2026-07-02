#import "@preview/touying:0.7.3": *
#import themes.metropolis: *

#let deck(body) = {
  set text(font: "New Computer Modern", size: 18pt)
  show math.equation: set text(font: "New Computer Modern Math")
  set par(justify: false, leading: 0.65em)
  show raw: set text(size: 11pt)
  show: metropolis-theme.with(
    aspect-ratio: "16-9",
    config-info(
      title: [MaxSPRT: vigilancia secuencial de seguridad de vacunas],
      subtitle: [Kulldorff et al. (2011), Sequential Analysis],
      author: [Santiago Olszevicki],
      institution: [],
      logo: [],
    ),
    config-colors(
      primary: rgb("#2f6f9f"),
      primary-light: rgb("#c9d9e6"),
      secondary: rgb("#1f2933"),
      neutral-lightest: rgb("#fbfcfe"),
      neutral-dark: rgb("#1f2933"),
      neutral-darkest: rgb("#111827"),
    ),
    footer: [MaxSPRT · Vigilancia de fármacos y vacunas],
  )
  body
}
