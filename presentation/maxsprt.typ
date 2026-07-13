#import "deck_setup.typ": deck
#show: deck

// colores para distinguir el RR de la alternativa/estimado y los esperados mu
#let aRR = rgb("#c0392b")   // RR (alternativa / a estimar)
#let aMU = rgb("#1f8a70")   // mu (esperados bajo H0)
#let cRR = text(fill: aRR)[RR]
#let cMt = text(fill: aMU)[$mu_t$]
#let cMi = text(fill: aMU)[$mu_i$]

= MaxSPRT: vigilancia secuencial de seguridad de vacunas

Un test de razón de verosimilitud secuencial con alternativa compuesta

Santiago Dandois y Santiago Olszevicki

Basado en Kulldorff, Davis, Kolczak, Lewis, Lieu y Platt (2011), _A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance_, Sequential Analysis 30(1).

== El problema

Los ensayos clínicos no alcanzan para detectar eventos adversos raros: muestra chica y población más homogénea que la real. Una vez en el mercado, hay que seguir vigilando.

Importa *revisar la evidencia en forma continua/casi continua*, controlando los falsos positivos.


== ¿Qué es el riesgo relativo?

$ "RR" = ("casos observados")/("casos esperados sin efecto de la vacuna") $

- $"RR" = 1$: no hay exceso, todo como se esperaba.
- $"RR" = 2$: el doble de casos que los esperados.
- $"RR" = 1.2$: un 20% más de casos que lo esperado.

== SPRT clásico (Wald, 1945-47)

$C_t$: eventos adversos observados hasta $t$ (en $D$ días); $#cMt$ los esperados bajo $H_0$. El SPRT contrasta dos hipótesis *simples*:

$ H_0: C_t tilde "Poisson"(#cMt) quad (#cRR=1) quad "vs." quad H_A: C_t tilde "Poisson"(#cRR #cMt) quad (#cRR" fija, a priori") $

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

== La idea de MaxSPRT

El SPRT clásico obliga a fijar el $#cRR$ de la alternativa antes de ver los datos. MaxSPRT usa una alternativa *compuesta*:

$ H_0: #cRR = 1 quad "vs." quad H_A: #cRR > 1 $

Clave: se *monta sobre el mismo $"LLR"_t$ ya deducido*. En vez de un $#cRR$ fijo, se usa el $#cRR$ que maximiza la verosimilitud — estimado de los datos, no elegido a priori:

$ "LLR"_t = max_(#cRR > 1) [ (1-#cRR)#cMt + c_t ln(#cRR) ] $

(razón de verosimilitud generalizada, Lorden 1973). No hay que adivinar el RR "interesante": el test busca el que mejor explica lo observado. Se paga algo de potencia, pero es mucho más robusto.

== Maximizando sobre RR

Derivando $"LLR"_t$ respecto de $#cRR$ e igualando a cero, el máximo se da en el $#cRR$ que iguala lo observado con lo esperado:

$ (partial "LLR"_t)/(partial #cRR) = -#cMt + c_t/#cRR = 0 quad => quad hat(#cRR)_t = c_t/#cMt $

Es el estimador de máxima verosimilitud: *reemplaza al $#cRR$ fijo de la alternativa*. Sustituyéndolo en $"LLR"_t$ se obtiene el estadístico final.

== MaxSPRT para datos Poisson

Sustituyendo $hat(#cRR)_t = c_t\/#cMt$:

$ "LLR"_t = (#cMt - c_t) + c_t ln(c_t / #cMt) $

y $"LLR"_t = 0$ si $c_t < #cMt$ (no hay evidencia de exceso).

Regla, monitoreada continuamente:

- se sigue vigilando mientras $"LLR"_t$ no llegue a un valor crítico $V$
- se rechaza $H_0$ apenas $"LLR"_t >= V$ (señal)
- ya no hay límite inferior para "aceptar" $H_A$ - no interesa frenar si la vacuna resulta protectora
- sí hay un límite superior $T$ en la duración de la vigilancia (cantidad esperada de eventos bajo $H_0$), para no vigilar para siempre

== Y si no tenemos $mu_t$ confiable...

A veces no hay tasa basal esperada confiable. Se compara directamente tiempo expuesto vs. no expuesto (diseño pareado, en vez de Poisson):

- autocontrolado: mismo individuo, ventana expuesta vs. no expuesta
- con controles: expuestos vs. no expuestos, pareados por edad/sexo

Cada evento es como tirar una moneda: cae del lado expuesto o no. Sea $z$ la duración del período no expuesto sobre la del expuesto ($z=1$ para un pareo 1:1).

$ H_0: P("expuesto") = 1/(z+1) quad quad H_A: P("expuesto") = "RR"/("RR"+z) $

El límite de vigilancia ya no se mide en eventos esperados ($T$) sino en eventos observados en total ($N$).

== MaxSPRT binomial: la verosimilitud

Sea $c_n$ de los $n$ eventos el número del lado expuesto. Condicional en $n$, $C_n$ es binomial (prob. $p$ de caer expuesto):

$ P(C_n = c_n) = binom(n, c_n) p^(c_n) (1-p)^(n-c_n) $

Esa $p$ depende del RR: $p("RR")="RR"\/("RR"+z)$ y $1-p = z\/("RR"+z)$. Reemplazando:

$ P(C_n = c_n | "RR") = binom(n, c_n) ("RR"/("RR"+z))^(c_n) (z/("RR"+z))^(n-c_n) $

El cociente contra $H_0$ ($p=1\/(z+1)$) cancela $binom(n,c_n)$ y se maximiza sobre $"RR">1$:

$ "LR"_n("RR") = (("RR"/("RR"+z))^(c_n) (z/("RR"+z))^(n-c_n)) / ((1/(z+1))^(c_n) (z/(z+1))^(n-c_n)) $

== MaxSPRT binomial: maximización y estadístico

Maximizar la binomial sobre $p$ da el estimador de proporción $hat(p) = c_n\/n$, que corresponde a $hat("RR")_n = z c_n\/(n-c_n)$. Reemplazando, el numerador queda $ (c_n\/n)^(c_n) ((n-c_n)\/n)^(n-c_n) $ y el estadístico es:

$ "LLR"_n = c_n ln(c_n/n) + (n-c_n)ln((n-c_n)/n) - c_n ln(1/(z+1)) - (n-c_n) ln(z/(z+1)) $

(válido si $z c_n\/(n-c_n) > 1$; si no, $"LLR"_n = 0$: no hay exceso del lado expuesto).

Misma lógica que el caso Poisson —alternativa compuesta $"RR">1$, maximización de la verosimilitud, tablas de valores críticos— pero partiendo de la binomial en vez de la Poisson.

== Valores críticos y ejemplo numérico

$V$ controla el error de tipo I: el valor más chico tal que, bajo $H_0$, la probabilidad de cruzarlo alguna vez antes de $T$ sea $alpha$.

$ P_(H_0) (exists thin t <= T : "LLR"_t >= V) = alpha $

No hay fórmula cerrada: se enumeran numéricamente las trayectorias bajo $H_0$ y se ajusta $V$ hasta ese $alpha$. Tabulado por $(alpha, T)$ una vez, no hace falta recalcularlo.

== Valores críticos y ejemplo numérico

#align(center)[
  #image("figures/sprt_trayectoria.png", width: 65%)
]

== Un ejemplo concreto

Para $alpha = 0.05$, algunos valores críticos de la Tabla 1 del paper:

- $T=2$ → $V approx 3.05$
- $T=10$ → $V approx 3.47$
- $T=800$ (vigilancia larga, ≈2 años) → $V approx 4.29$

A mayor $T$, mayor $V$: se permiten más chances de cruzarlo por azar, así que hace falta una señal más fuerte para mantener el mismo $alpha$.

Con $T=2$ y $mu_t=1$: si $c_t=4$, $"LLR"_t = (1-4)+4ln(4) approx 2.55$, todavía no alcanza $V=3.05$. Con $c_t=5$: $"LLR"_t=(1-5)+5ln(5) approx 4.05$, cruza el límite y se rechaza $H_0$.

== Potencia y el trade-off de siempre

La potencia depende del RR verdadero y de $T$.

- $T$ chico: se termina antes, pero cuesta detectar riesgos moderados.
- $T$ grande: más potencia, pero hay que esperar más si la señal tarda.

El mismo trade-off de siempre entre muestra y potencia, definido en eventos esperados bajo $H_0$ en vez de un $n$ fijo. Como es vigilancia observacional, estirar $T$ casi solo cuesta cómputo.

== Simulación / Ejemplo (Santi D.)

#align(center + horizon)[
  #text(size: 24pt, fill: luma(150))[_Placeholder — simulación / ejemplo (Santi D.)_]
]

== Para llevarse

- El SPRT clásico de Wald funciona, pero es demasiado sensible al RR elegido para la alternativa - un problema serio cuando no sabemos ese valor de antemano.
- MaxSPRT resuelve esto con una alternativa compuesta ($"RR">1$) y maximizando la verosimilitud en vez de fijar un valor.
- Existen versiones para datos Poisson y binomiales, cubriendo los diseños más usados en farmacovigilancia.
- Es un método que se sigue usando hoy en sistemas reales de vigilancia de vacunas y fármacos.

== Referencias

- Kulldorff, M., Davis, R. L., Kolczak, M., Lewis, E., Lieu, T. y Platt, R. (2011). _A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance_. Sequential Analysis, 30(1), 58-78.
- Wald, A. (1945). _Sequential Tests of Statistical Hypotheses_. Annals of Mathematical Statistics.
- Lorden, G. (1973). _Open-Ended Tests for Koopman-Darmois Families_. Annals of Statistics.
- Lieu, T. A. et al. (2007). _Real-Time Vaccine Safety Surveillance for the Early Detection of Adverse Events_. Medical Care.
