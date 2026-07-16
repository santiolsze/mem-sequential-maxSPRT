#import "deck_setup.typ": deck
#import "@preview/cetz:0.4.1"
#show: deck

// colores para distinguir el RR de la alternativa/estimado y los esperados mu
#let aRR = rgb("#c0392b")   // RR (alternativa / a estimar)
#let aMU = rgb("#1f8a70")   // mu (esperados bajo H0)
#let cRR = text(fill: aRR)[RR]
#let cMt = text(fill: aMU)[$mu_t$]
#let cMi = text(fill: aMU)[$mu_i$]

#let separador_seccion(titulo, detalle) = {
  v(1fr)
  align(center)[
    #block(width: 78%)[
      #line(length: 40%, stroke: 1.4pt + rgb("#2f6f9f"))
      #v(14pt)
      #text(size: 30pt, weight: "bold", fill: rgb("#111827"))[#titulo]
      #v(8pt)
      #text(size: 16pt, fill: rgb("#4b5563"))[#detalle]
      #v(14pt)
      #line(length: 40%, stroke: 1.4pt + rgb("#2f6f9f"))
    ]
  ]
  v(1fr)
}

= MaxSPRT: vigilancia secuencial de seguridad de vacunas

#v(0.45fr)

#block(width: 100%)[
  #text(size: 25pt, weight: "medium", fill: rgb("#111827"))[Un test de razón de verosimilitud secuencial]
  #v(3pt)
  #text(size: 25pt, weight: "medium", fill: rgb("#111827"))[con alternativa compuesta]

  #v(18pt)
  #line(length: 18%, stroke: 2pt + rgb("#2f6f9f"))
  #v(18pt)

  #text(size: 20pt, weight: "bold", fill: rgb("#2f6f9f"))[
    Santiago Dandois y Santiago Olszevicki
  ]
]

#v(0.55fr)

#block(
  width: 92%,
  inset: (left: 14pt, top: 8pt, bottom: 8pt),
  stroke: (left: 2pt + rgb("#c9d9e6")),
)[
  #text(size: 13pt, fill: rgb("#4b5563"))[
    Basado en Kulldorff, Davis, Kolczak, Lewis, Lieu y Platt (2011),
    _A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance_,
    Sequential Analysis 30(1).
  ]
]

#v(0.25fr)

== Introducción

#separador_seccion[Un problema real de salud pública][El paper, sus autores, el CDC y por qué hace falta vigilancia después de la aprobación]

== Esto no nace de un ejercicio académico

El *CDC* —Centers for Disease Control and Prevention— es una de las principales agencias federales de salud pública de Estados Unidos. Vigila enfermedades y convierte datos clínicos en decisiones de alcance nacional.

- Presupuesto solicitado para 2025: *USD 9.683 millones*.
- Su Vaccine Safety Datalink reúne historias clínicas de *más de 10 millones de personas por año* y registros de *más de 180 millones de vacunaciones*.
- MaxSPRT nació para una necesidad concreta de ese sistema: *revisar semanalmente la seguridad de vacunas sin multiplicar las falsas alarmas*.

== ¿Por qué hace falta farmacovigilancia?

Antes de la aprobación, los ensayos clínicos permiten detectar efectos frecuentes. Pero su tamaño y su población seleccionada pueden ocultar eventos que aparecen:

- una vez cada decenas o cientos de miles de aplicaciones;
- muchos meses después;
- solamente en embarazadas, niños, adultos mayores u otros subgrupos.

Cuando el producto empieza a utilizarse masivamente, la pregunta pasa a ser:

#align(center)[
  #text(size: 22pt, weight: "bold", fill: rgb("#2f6f9f"))[
    ¿Podemos detectar un riesgo raro lo antes posible,
    sin generar una alarma cada vez que miramos los datos?
  ]
]


== SPRT clásico

#separador_seccion[El punto de partida: SPRT clásico][Cómo funciona y por qué fijar el riesgo relativo puede retrasar o perder una señal]

== ¿Qué es el riesgo relativo?

$ "RR" = ("casos observados")/("casos esperados sin efecto de la vacuna") $

- $"RR" = 1$: no hay exceso, todo como se esperaba.
- $"RR" = 2$: el doble de casos que los esperados.
- $"RR" = 1.2$: un 20% más de casos que lo esperado.

== SPRT clásico (Wald, 1945-47)

$C_t$: eventos adversos observados hasta $t$ (en $D$ días); $#cMt$ los esperados bajo $H_0$. El SPRT contrasta dos hipótesis *simples*:

$ H_0: C_t tilde "Poisson"(#cMt) quad (#cRR=1) quad "vs." quad H_A: C_t tilde "Poisson"(#cRR #cMt) quad (#cRR" fija, a priori") $

== La tasa basal puede cambiar en el tiempo

En general, bajo $H_0$ puede plantearse un *proceso de Poisson no homogéneo*:

$ C_t tilde "Poisson"(#cMt), quad #cMt = integral_0^t lambda_0(s) dif s $

La intensidad basal $lambda_0(s)$ puede variar por estacionalidad, volumen de vacunación, campañas, días hábiles, etc. No se supone una tasa constante: $#cMt$ acumula la cantidad esperada de eventos hasta $t$.

Lo importante es que, para $s<t$, los incrementos en intervalos disjuntos sean independientes y satisfagan

$ C_t - C_s tilde "Poisson"(#cMt - mu_s) quad "bajo " H_0. $

Aquí $C_t$ es el conteo aleatorio y $c_t$ es el valor acumulado efectivamente observado.

== SPRT clásico: construcción de la verosimilitud

El estadístico es la razón de verosimilitudes $"LR"_t = P("datos"|H_A) \/ P("datos"|H_0)$. Para armarla, partimos $[0,t]$ en intervalos $i=1,...,t$ con $n_i$ eventos y $#cMi$ esperados, independientes:

$ n_i tilde "Poisson"(#cMi) " bajo " H_0, quad n_i tilde "Poisson"(#cRR #cMi) " bajo " H_A $

La verosimilitud conjunta es la productoria de las densidades Poisson:

$ L_t (#cRR) = product_(i=1)^t (e^(-#cRR #cMi) (#cRR #cMi)^(n_i))/(n_i !) $

== SPRT clásico: de la verosimilitud al estadístico

Reagrupando la productoria (exponenciales se suman, potencias de $#cRR$ se acumulan):

$ L_t (#cRR) = e^(-#cRR #cMt) #cRR^(c_t) product_(i=1)^t #cMi^(n_i)/n_i! $

con $#cMt = sum_i #cMi$ y $c_t = sum_i n_i$. El último factor no depende de $#cRR$ y se cancela en el cociente:

$ "LR"_t = (L_t (#cRR))/(L_t (1)) = e^((1-#cRR)#cMt) #cRR^(c_t) quad => quad "LLR"_t = (1-#cRR)#cMt + c_t ln(#cRR) $

Regla: señal si $"LLR"_t >= ln[(1-beta)\/alpha]=2.77$; se acepta $H_0$ si $"LLR"_t <= ln[beta\/(1-alpha)]=-1.56$ ($alpha=0.05, beta=0.20$).

Hay que fijar el $#cRR$ de antemano — ese "detalle" es el problema.

== El problema del SPRT clásico (I): fiebre

#grid(columns: (58%, 42%), gutter: 1em,
  image("figures/sprt_fiebre.png", width: 100%),
  [
    Pediarix, ≈650.000 chicos (Vaccine Safety Datalink). El RR verdadero termina en *≈ 1.16*.

    - con $H_A: "RR"=1.2$ (cerca del real) → *señal a las 13 sem*.
    - con $H_A: "RR"=2.0$ (lejos del real) → cruza el límite inferior y *acepta $H_0$ a las 7 sem*: no detecta nada.

    La alternativa mal elegida (2.0) hace *perder* un riesgo real.
  ]
)

== El problema del SPRT clásico (II): síntomas neurológicos

#grid(columns: (58%, 42%), gutter: 1em,
  image("figures/sprt_neuro.png", width: 100%),
  [
    Mismos datos, otro evento. Acá el RR verdadero termina en *≈ 2.75*.

    - con $H_A: "RR"=2.0$ (cerca del real) → *señal a las 32 sem*.
    - con $H_A: "RR"=1.2$ (lejos del real) → *señal recién a las 65 sem*.

    Ambas detectan, pero elegir mal el RR *retrasa la señal* casi 8 meses: el tiempo hasta la alarma depende del RR fijado a priori, que casi nunca conocemos.
  ]
)

== MaxSPRT Poisson

#separador_seccion[MaxSPRT: vigilar sin adivinar el riesgo][La alternativa compuesta y el modelo Poisson]

== La idea de MaxSPRT

El SPRT clásico obliga a fijar el $#cRR$ de la alternativa antes de ver los datos. MaxSPRT usa una alternativa *compuesta*:

$ H_0: #cRR = 1 quad "vs." quad H_A: #cRR > 1 $

Clave: se *apoya en el mismo $"LLR"_t$ ya deducido*. En vez de un $#cRR$ fijo, se usa el $#cRR$ que maximiza la verosimilitud — estimado de los datos, no elegido a priori:

$ "LLR"_t = max_(#cRR >= 1) [ (1-#cRR)#cMt + c_t ln(#cRR) ] $

(razón de verosimilitud generalizada, Lorden 1973). No hay que adivinar el RR "interesante": el test busca el que mejor explica lo observado. Se paga algo de potencia, pero es mucho más robusto.

== Maximizando sobre RR

Derivando $"LLR"_t$ respecto de $#cRR$ e igualando a cero, el máximo se da en el $#cRR$ que iguala lo observado con lo esperado:

$ (partial "LLR"_t)/(partial #cRR) = -#cMt + c_t/#cRR = 0 quad => quad hat(#cRR)_t = c_t/#cMt $

Ese es el estimador sin restricción. Como la alternativa solo admite $#cRR >= 1$, el estimador restringido es

$ hat(#cRR)_(t,"restr") = max(1, c_t/#cMt) $

Por eso solo toma el valor $c_t/#cMt > 1$ cuando $c_t > #cMt$; si $c_t <= #cMt$, el máximo queda en el borde $#cRR=1$ y el LLR vale cero.

== MaxSPRT para datos Poisson

Sustituyendo $hat(#cRR)_t = c_t\/#cMt$:

$ "LLR"_t = (#cMt - c_t) + c_t ln(c_t / #cMt) $

y $"LLR"_t = 0$ si $c_t <= #cMt$ (no hay evidencia de exceso).

Regla, monitoreada continuamente:

- se sigue vigilando mientras $"LLR"_t$ no llegue a un valor crítico $V$
- se rechaza $H_0$ apenas $"LLR"_t >= V$ (señal)
- ya no hay límite inferior para "aceptar" $H_A$ - no interesa frenar si la vacuna resulta protectora
- sí hay un límite superior $T$ en la duración de la vigilancia (cantidad esperada de eventos bajo $H_0$), para no vigilar para siempre

== Variante binomial

#separador_seccion[Cuando no conocemos los casos esperados][La variante binomial: comparar eventos entre expuestos y controles]

== Y si no tenemos $mu_t$ confiable...

A veces no hay tasa basal esperada confiable. Se compara directamente tiempo expuesto vs. no expuesto (diseño pareado, en vez de Poisson):

- autocontrolado: mismo individuo, ventana expuesta vs. no expuesta
- con controles: expuestos vs. no expuestos, pareados por edad/sexo

Cada evento cae del lado expuesto o no. Sea $z$ la razón de oportunidades no expuestas sobre expuestas (tiempo, personas o controles).

$ H_0: P("expuesto") = 1/(z+1) quad quad H_A: P("expuesto") = "RR"/("RR"+z) $

Solo en un pareo 1:1, $z=1 arrow P("expuesto")=1/2$; con dos controles, $z=2$ y vale $1/3$.

El límite de vigilancia ya no se mide en eventos esperados ($T$) sino en eventos observados en total ($N$).

== MaxSPRT binomial: la verosimilitud

Sea $c_n$ la cantidad de esos $n$ eventos que cae del lado expuesto. Condicional en $n$, $C_n$ es binomial (prob. $p$ de caer expuesto):

$ P(C_n = c_n) = binom(n, c_n) p^(c_n) (1-p)^(n-c_n) $

Esa $p$ depende del RR: $p("RR")="RR"\/("RR"+z)$ y $1-p = z\/("RR"+z)$. Reemplazando:

$ P(C_n = c_n | "RR") = binom(n, c_n) ("RR"/("RR"+z))^(c_n) (z/("RR"+z))^(n-c_n) $

El cociente contra $H_0$ ($p=1\/(z+1)$) cancela $binom(n,c_n)$ y se maximiza sobre $"RR">1$:

$ "LR"_n("RR") = (("RR"/("RR"+z))^(c_n) (z/("RR"+z))^(n-c_n)) / ((1/(z+1))^(c_n) (z/(z+1))^(n-c_n)) $

== MaxSPRT binomial: maximización y estadístico

Maximizar la binomial sobre $p$ da el estimador de proporción $hat(p) = c_n\/n$, que corresponde a $hat("RR")_n = z c_n\/(n-c_n)$. Reemplazando, el numerador queda $ (c_n\/n)^(c_n) ((n-c_n)\/n)^(n-c_n) $ y el estadístico es:

$ "LLR"_n = c_n ln(c_n/n) + (n-c_n)ln((n-c_n)/n) - c_n ln(1/(z+1)) - (n-c_n) ln(z/(z+1)) $

(válido si $z c_n\/(n-c_n) > 1$; si no, $"LLR"_n = 0$: no hay exceso del lado expuesto).

Misma lógica que el caso Poisson —alternativa compuesta $"RR">1$, maximización de la verosimilitud, tablas de valores críticos— pero partiendo de la binomial en vez de la Poisson. La frontera binomial se calibra para $(alpha, N, z)$.

== Umbrales Poisson

#separador_seccion[¿Cómo se decide cuándo alertar?][Horizonte, umbral crítico, potencia y construcción exacta de las tablas del paper]

== El umbral crítico controla las falsas alarmas

$V$ controla el error de tipo I: el valor más chico tal que, bajo $H_0$, la probabilidad de cruzarlo alguna vez antes de $T$ sea $alpha$.

$ P_(H_0) (exists thin t <= T : "LLR"_t >= V) = alpha $

No hay fórmula cerrada: $alpha$ se calcula con una recursión numérica bajo $H_0$, y se ajusta $V$ hasta ese valor. Tabulado por $(alpha, T)$ una vez, no hace falta recalcularlo.

== El LLR baja entre eventos y salta cuando llega uno

#align(center)[
#cetz.canvas(length: 1.15cm, {
  import cetz.draw: *
  let azul = rgb("#2f6f9f")
  let oscuro = rgb("#1f2933")
  let naranja = rgb("#b3541e")
  let gris = rgb("#9aa5b1")

  // Ejes: la escala horizontal es media esperada acumulada.
  line((0, 0), (11.4, 0), stroke: 1.2pt + oscuro, mark: (end: ">"))
  line((0, 0), (0, 6.0), stroke: 1.2pt + oscuro, mark: (end: ">"))
  content((11.35, 0.18), anchor: "south-east", text(size: 13pt)[$mu_t$ (eventos esperados)])
  content((-0.2, 6.0), anchor: "south-east", text(size: 13pt)[$"LLR"_t$])

  // Frontera horizontal de rechazo y horizonte vertical.
  line((0, 4.55), (11.0, 4.55), stroke: (paint: naranja, thickness: 1.2pt, dash: "dashed"))
  content((0.15, 4.72), anchor: "south-west", text(size: 12pt, fill: naranja)[umbral $V$: rechazar $H_0$])
  line((9.7, 0), (9.7, 5.75), stroke: (paint: gris, thickness: 1.2pt, dash: "dashed"))
  content((9.7, -0.18), anchor: "north", text(size: 12pt, fill: oscuro)[$T=mu(t_"máx")$])
  content((9.7, 5.75), anchor: "south", text(size: 11pt, fill: oscuro)[fin si no hubo señal])

  // Trayectoria ilustrativa: descenso continuo y saltos en los eventos.
  line(
    (0, 0), (1.2, 0), (1.2, 1.25),
    (2.5, 0.65), (2.5, 2.15),
    (3.9, 1.35), (3.9, 3.15),
    (5.5, 2.15), (5.5, 4.05),
    (7.0, 3.05), (7.0, 5.05),
    stroke: 2.2pt + azul,
  )

  for (x, y0, y1) in ((1.2, 0, 1.25), (2.5, .65, 2.15), (3.9, 1.35, 3.15), (5.5, 2.15, 4.05)) {
    circle((x, y1), radius: .07, fill: azul, stroke: none)
  }
  circle((7.0, 5.05), radius: .11, fill: naranja, stroke: none)
  content((7.2, 5.15), anchor: "south-west", text(size: 12pt, fill: naranja)[evento $arrow.r$ señal])

  content((2.55, .35), anchor: "south-west", text(size: 11pt, fill: azul)[entre eventos: disminuye])
  content((4.05, 3.4), anchor: "west", text(size: 11pt, fill: azul)[evento: salta])
})
]

$mu_t$, la cantidad acumulada de eventos esperados, es una función estrictamente creciente. Por eso hay una correspondencia biyectiva entre el tiempo transcurrido y los eventos esperados. La vigilancia termina al cruzar $V$ o, si no hay señal, al llegar a $T$.

== El horizonte determina el umbral

Para $alpha = 0.05$, algunos ejemplos de valores críticos son:

#align(center)[
#table(
  columns: (1fr, 1fr),
  align: center,
  stroke: 0.5pt + luma(190),
  inset: 8pt,
  [*Horizonte $T$*], [*Umbral crítico $V$*],
  [$2$], [$approx 3.05$],
  [$10$], [$approx 3.47$],
  [$800$], [$approx 4.29$],
)
]

A mayor $T$, mayor $V$: se permiten más chances de cruzarlo por azar, así que hace falta una señal más fuerte para mantener el mismo $alpha$.

#text(weight: "bold", fill: rgb("#2f6f9f"))[
  Un horizonte mayor permite detectar incrementos de riesgo más moderados, pero exige acumular evidencia durante más tiempo.
]

== Por dentro: ¿cómo se calcula $V$?

¿Cómo se hace ese cálculo numérico, en la práctica?

Entre eventos, $c_t$ queda fijo y #cMt crece:

$ (partial "LLR"_t)/(partial #cMt) = 1 - c_t/#cMt < 0 quad "si " c_t > #cMt $

*El rechazo solo puede pasar en el instante de un evento nuevo* — alcanza con mirar ahí.

Dos herramientas para hacerlo exacto (no simulado):

- *Función $W$ de Lambert*: resuelve $w e^w = z$; da el instante donde $n$ eventos cruzarían $V$.
- *Incrementos independientes de Poisson*: lo que pasa en un tramo nuevo no depende de lo anterior.

== Paso 1a — normalizar la ecuación de cruce

Para cada $n=1,2,3,...$ buscamos $s_n$ tal que, con exactamente $n$ eventos ahí, el LLR da justo $V$:

$ (mu_(s_n) - n) + n ln(n/mu_(s_n)) = V $

Definimos $x=mu_(s_n)\/n$. Como el rechazo exige un exceso, buscamos $0<x<1$. Dividiendo por $n$:

$ x - 1 + ln(1/x) = V/n $

$ x - ln x = 1 + V/n $

== Paso 1b — llevarla a la forma de Lambert $W$

Primero aislamos el logaritmo y exponenciamos:

$ ln x = x - 1 - V/n $

$ x = e^(x-1-V/n) $

Luego multiplicamos por $e^(-x)$ y por $-1$:

$ x e^(-x) = e^(-1-V/n) $

$ (-x)e^(-x) = -e^(-1-V/n) $

Ahora el lado izquierdo tiene exactamente la forma $w e^w$, con $w=-x$.

== Paso 1c — invertir y elegir la rama correcta

La función de Lambert satisface $W(z)e^(W(z))=z$. Por lo tanto:

$ -x = W_0(-e^(-1-V/n)) $

y, volviendo a $mu_(s_n)=n x$:

#align(center)[*$ mu_(s_n) = -n W_0(-e^(-1-V\/n)) $*]

Se usa la rama principal $W_0$: produce $0<x<1$. La rama $W_(-1)$ da $x>=1$, fuera de la región de exceso del test unilateral.

Así obtenemos cada punto de cruce en la escala de media esperada, $mu_(s_n)$, sin simular trayectorias continuas.

== Lema: incrementos independientes de Poisson

Bajo $H_0$, sea $C_t$ un proceso de Poisson con media acumulada conocida $#cMt$. Para $s < t$:

$ C_t - C_s tilde "Poisson"(mu_t - mu_s), quad "independiente de " C_s $

El conteo en un tramo nuevo no depende de lo que ya pasó — el proceso "no tiene memoria" entre tramos disjuntos.

== Paso 2 — la recursión: ¿seguís vivo hasta $s_n$?

$pi_n (k) := P(C_(s_n)=k med, "no rechazado hasta la etapa " n)$, para $n>=1$ y $k=0,...,n-1$.

*Caso inicial*: en la primera etapa, seguir sin rechazo exige no haber observado eventos:

$ pi_1(0) = P(C_(s_1)=0) = e^(-mu_(s_1)) $

*Paso recursivo* — para $n>=2$, con $Delta_n = mu_(s_n)-mu_(s_(n-1))$, combinamos $pi_(n-1)$ con el incremento (Poisson($Delta_n$), independiente) y truncamos a $k<n$:

$ pi_n(k) &= sum_(j=0)^(min(k,n-2)) pi_(n-1)(j) med P("incremento"=k-j) \
&= sum_(j=0)^(min(k,n-2)) pi_(n-1)(j) med e^(-Delta_n) (Delta_n^(k-j))/((k-j)!), quad k=0,...,n-1 $

== Paso 3 — sumar todo lo rechazado

Sea $Pi_n := sum_(k=0)^(n-1) pi_n(k)$ la probabilidad total de no rechazar tras la etapa $n$. Definimos por separado $Pi_0:=1$: antes de la primera etapa, la vigilancia sigue activa con certeza.

Acá entra $T$: define cuándo parar,

$ n_max := min { n : mu_(s_n) >= T }. $

En la última etapa no avanzamos hasta $mu_(s_(n_max))$: usamos el incremento $T-mu_(s_(n_max-1))$ para terminar exactamente en el horizonte.

$Pi_n$ es la probabilidad de continuar sin rechazo después de la etapa $n$. Por lo tanto, $Pi_(n-1)-Pi_n$ es la probabilidad de rechazar precisamente durante la etapa $n$. Al sumar sobre todas las etapas obtenemos una suma telescópica:

$ alpha(V,T) = sum_(n=1)^(n_max) (Pi_(n-1)-Pi_n) = 1 - Pi_(n_max) $

== Ejemplo de la recursión

Ejemplo con $V=2$, $T=3$. Definimos $u_n:=min(mu_(s_n),T)$:

#grid(
  columns: (1fr, 1fr),
  gutter: 18pt,
  align(center)[
    #text(size: 14pt)[$mu_(s_7)=-7W_0(-e^(-1-2/7)) approx 2.949$]
  ],
  align(center)[
    #text(size: 14pt)[$mu_(s_8)=-8W_0(-e^(-1-2/8)) approx 3.590$]
  ],
)

#align(center)[
  #text(size: 14pt, weight: "bold")[$mu_(s_7)<T<mu_(s_8) quad arrow.r quad n_max=8$]
]

#align(center)[
#text(size: 14pt)[
#table(
  columns: 4,
  align: center,
  stroke: 0.5pt + luma(200),
  inset: 5pt,
  [*n*], [*$u_n$*], [*$Pi_n$*], [*$1-Pi_n$*],
  [1], [0.053], [0.9489], [0.0511],
  [2], [0.317], [0.9210], [0.0790],
  [3], [0.720], [0.9030], [0.0970],
  [$dots.v$], [$dots.v$], [$dots.v$], [$dots.v$],
  [8], [3.000 ($T$)], [0.8647], [*0.1353*],
)
]
]

#v(14pt)
#align(center)[
  #block(
    width: 54%,
    inset: 10pt,
    fill: rgb("#eef5fa"),
    stroke: 1pt + rgb("#2f6f9f"),
  )[
    #align(center)[
      #text(size: 20pt, weight: "bold", fill: rgb("#2f6f9f"))[
        $alpha(V=2,T=3)=1-Pi_8=0.1353$
      ]
    ]
  ]
]

== Paso 4 — invertir $alpha(V,T)$ para obtener $V$

La recursión anterior responde: dado un candidato $V$, ¿cuánto vale el error real $alpha(V,T)$ bajo $H_0$?

Para construir la Tabla 1 fijamos $T$ y el $alpha$ deseado, y resolvemos:

$ f(V) := alpha(V,T) - alpha_0 = 0 $

Subir $V$ hace más difícil rechazar, así que $alpha(V,T)$ es decreciente. Por lo tanto, es fácil buscar la raíz numéricamente.

Cada evaluación vuelve a ejecutar la cadena completa:

#align(center)[Lambert $W$ $arrow.r$ grilla $mu_(s_n)$ $arrow.r$ recursión Poisson $arrow.r$ masa absorbida $alpha(V,T)$.]

== Valores críticos para distintos horizontes

#align(center)[
#text(size: 13pt)[
#table(
  columns: (22%, 26%, 26%, 26%),
  align: center,
  stroke: 0.5pt + luma(190),
  inset: 4.5pt,
  [*Horizonte $T$*], [*$alpha=0.05$*], [*$alpha=0.01$*], [*$alpha=0.001$*],
  [$0.1$], [$2.04407$], [$4.11929$], [$6.57967$],
  [$0.5$], [$2.63793$], [$4.48374$], [$7.03447$],
  [$1$], [$2.85394$], [$4.67043$], [$7.17261$],
  [$2$], [$3.04698$], [$4.86222$], [$7.34145$],
  [$5$], [$3.29718$], [$5.09191$], [$7.56931$],
  [$10$], [$3.46795$], [$5.26051$], [$7.72486$],
  [$50$], [$3.81990$], [$5.60597$], [$8.06707$],
  [$100$], [$3.95232$], [$5.73897$], [$8.19940$],
  [$500$], [$4.22263$], [$6.01109$], [$8.47318$],
  [$1000$], [$4.32492$], [$6.11423$], [$8.57725$],
)
]
]

== Revisiones periódicas

#separador_seccion[De vigilancia continua a revisiones periódicas][Qué cambia al mirar cada semana o cada mes, en lugar de hacerlo ante cada evento]

== Discretización: el experimento

Hasta acá $V$ es exacto para vigilancia *continua*. En la práctica se mira cada tanto: $K$ miradas igualmente espaciadas en $#cMt$ acumulado, $#cMi = i T\/K$.

Se deja el *mismo* $V$ continuo fijo — no se recalibra por $K$, porque recalibrar escondería el efecto que se quiere mostrar. Simulando muchas trayectorias con ese $V$, se mide:

- $alpha$ empírico: fracción que cruza $V$ bajo $H_0$
- potencia empírica: fracción que cruza $V$ bajo un $"RR"$ verdadero, contra la potencia continua exacta

== Dos fuentes, tres series de reportes

Usamos datos públicos mensuales de 2016 a 2025:

#align(center)[
#table(
  columns: (28%, 42%, 30%),
  align: (left, left, left),
  stroke: 0.5pt + luma(190),
  inset: 8pt,
  [*Fuente*], [*Producto*], [*Evento buscado*],
  [FAERS / openFDA], [ELIQUIS], [Haemorrhage],
  [VAERS], [HPV9 (GARDASIL 9)], [Syncope],
  [VAERS], [MENB], [Pyrexia],
)
]

#v(16pt)

#align(center)[
  #text(size: 20pt, weight: "bold", fill: rgb("#2f6f9f"))[
    Cada fila representa reportes recibidos, no todas las personas expuestas ni todas las dosis administradas.
  ]
]

== ¿Contra qué comparamos cada serie?

#text(size: 15pt)[
*ELIQUIS:* usamos como referencia los reportes no-ELIQUIS del mismo mes.

$ E_t = "reportes ELIQUIS"_t times ("Haemorrhage en no-ELIQUIS"_t)/("reportes no-ELIQUIS"_t) $

*HPV9 y MENB:* usamos MNQ —vacuna meningocócica ACWY de la misma visita adolescente— como comparador contemporáneo.

$ E_t = "reportes de la vacuna"_t times ("evento en MNQ"_t)/("reportes MNQ"_t) $
]

#v(10pt)

En la pestaña *Binomial MENB vs MNQ* no construimos esperados Poisson: comparamos directamente las proporciones mensuales de reportes con Pyrexia, condicionando por los totales del mes.

== Qué vamos a mirar

La demo permite mostrar con datos reales:

- cómo evolucionan MaxSPRT y los SPRT clásicos;
- cuándo aparece una señal secuencial.

== Demostración

#separador_seccion[¡Veámoslo en acción!][Shiny R app]

== Cierre

#separador_seccion[Cerrando][Conclusiones principales]

== Para llevarse

- El SPRT clásico de Wald funciona, pero es demasiado sensible al RR elegido para la alternativa - un problema serio cuando no sabemos ese valor de antemano.
- MaxSPRT resuelve esto con una alternativa compuesta ($"RR">1$) y maximizando la verosimilitud en vez de fijar un valor.
- Existen versiones para datos Poisson y binomiales, cubriendo los diseños más usados en farmacovigilancia.
- Es un método que se sigue usando hoy en sistemas reales de vigilancia de vacunas y fármacos.

== Referencias

- Kulldorff, M., Davis, R. L., Kolczak, M., Lewis, E., Lieu, T. y Platt, R. (2011). _A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance_. Sequential Analysis, 30(1), 58-78.
- Centers for Disease Control and Prevention (2024). _Chapter 4: Vaccine Safety_. Epidemiology and Prevention of Vaccine-Preventable Diseases (Pink Book).
- Centers for Disease Control and Prevention (2024). _FY 2025 CDC Budget Overview_. Office of Financial Resources.
- Wald, A. (1945). _Sequential Tests of Statistical Hypotheses_. Annals of Mathematical Statistics.
- Lorden, G. (1973). _Open-Ended Tests for Koopman-Darmois Families_. Annals of Statistics.
- Lieu, T. A. et al. (2007). _Real-Time Vaccine Safety Surveillance for the Early Detection of Adverse Events_. Medical Care.
