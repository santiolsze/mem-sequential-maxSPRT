#import "deck_setup.typ": deck
#show: deck

= MaxSPRT: vigilancia secuencial de seguridad de vacunas

Un test de razón de verosimilitud secuencial con alternativa compuesta

Santiago Dandois y Santiago Olszevicki

Basado en Kulldorff, Davis, Kolczak, Lewis, Lieu y Platt (2011), _A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance_, Sequential Analysis 30(1).

== El problema

Los ensayos clínicos (fase 2/3) no alcanzan para detectar eventos adversos raros: la muestra es chica y la población suele ser más homogénea que la real.

Por eso, una vez que un fármaco o vacuna sale al mercado, hace falta seguir vigilando.

- Datos disponibles: reclamos de seguros de salud, historias clínicas electrónicas.
- Se quiere revisar la evidencia seguido (semanal, casi continuo), no una sola vez al final.
- Mirar los datos muchas veces sin corregir infla los falsos positivos. Hace falta un método pensado para eso.

Ahí entra el análisis secuencial.

== ¿Qué es el riesgo relativo?

El riesgo relativo (RR) compara cuántos casos de un evento se observan realmente contra cuántos se esperarían si la vacuna no tuviera ningún efecto:

$ "RR" = ("casos observados")/("casos esperados sin efecto de la vacuna") $

- $"RR" = 1$: no hay exceso, todo es como se esperaba.
- $"RR" = 2$: se observó el doble de casos que los esperados.
- $"RR" = 1.2$: un 20% más de casos que lo esperado.

Ejemplo: si históricamente se esperan 50 casos de fiebre en cierto período y se observan 100, $"RR"=2$.

Todo lo que sigue busca decidir, a medida que llegan los datos, si el RR verdadero es mayor a 1.

== SPRT clásico (Wald, 1945-47)

Sea $C_t$ el número de eventos adversos observados hasta el tiempo $t$, con media esperada $mu_t$ bajo $H_0$ (sin riesgo agregado).

$ H_0: "RR" = 1 quad "vs." quad H_A: "RR" = "RR fijo" $

El estadístico es la razón de verosimilitudes (Poisson):

$ "LLR"_t = (1-"RR") mu_t + c_t ln("RR") $

Regla de decisión, monitoreada continuamente:

- si $"LLR"_t >= ln[(1-beta)\/alpha]$: se rechaza $H_0$ (señal)
- si $"LLR"_t <= ln[beta\/(1-alpha)]$: se acepta $H_0$

Simple, pero hay que fijar un RR de antemano. Ese "detalle" resulta ser el problema.

== El problema del SPRT clásico

Dos vigilancias sobre los mismos datos, cambiando solo el RR de la alternativa: una llega antes al límite, la otra tarda mucho más (o ni llega).

== El problema del SPRT clásico

#align(center)[
  #image("figures/sprt_sensibilidad.png", width: 72%)
]

== Pediarix: un caso real de resultados contradictorios

Datos de Vaccine Safety Datalink, ≈650.000 chicos. Se mira si aumenta el riesgo de fiebre o de síntomas neurológicos en los 28 días post-vacunación.

*Fiebre:*
- con $H_A: "RR"=2.0$ → no se rechaza $H_0$ (7 semanas y se acepta)
- con $H_A: "RR"=1.2$ → se rechaza $H_0$ a las 13 semanas

*Síntomas neurológicos:*
- con $H_A: "RR"=1.2$ → señal a las 65 semanas
- con $H_A: "RR"=2.0$ → señal a las 32 semanas

Mismos datos, conclusión opuesta según qué RR se haya elegido a priori. Y en la práctica casi nunca se sabe de antemano cuál va a ser el riesgo relativo verdadero.

== La idea de MaxSPRT

En lugar de fijar un RR puntual, se usa una alternativa *compuesta*:

$ H_0: "RR" = 1 quad "vs." quad H_A: "RR" > 1 $

El estadístico ya no evalúa un único RR: maximiza la verosimilitud sobre todos los RR posibles (test de razón de verosimilitud generalizado, Lorden 1973).

$ "LR"_t = max_("RR">1) (P(C_t=c_t | H_A)) / (P(C_t=c_t | H_0)) $

Ventaja práctica: no hay que adivinar el riesgo relativo "interesante". El propio test busca el RR que mejor explica lo observado.

Se paga un precio en potencia comparado con adivinar el RR correcto de entrada, pero es mucho más robusto cuando no sabemos ese valor.

== Formalizando el riesgo relativo

Para trabajar con esto matemáticamente, dividimos el tiempo en intervalos $1, ..., t$. En cada uno se observan $n_i$ eventos, con $mu_i$ esperados bajo $H_0$, independientes entre sí:

$ n_i tilde "Poisson"(mu_i) "  bajo " H_0 quad quad n_i tilde "Poisson"("RR" mu_i) "  bajo " H_A $

La verosimilitud conjunta hasta $t$ es la productoria de esas densidades:

$ L_t("RR") = product_(i=1)^t (e^(-"RR" mu_i) ("RR" mu_i)^(n_i))/(n_i !) $

== Cociente de verosimilitud Poisson

Reagrupando la productoria (las exponenciales se suman, las potencias de RR se acumulan):

$ L_t("RR") = e^(-"RR" mu_t) "RR"^(c_t) product_(i=1)^t mu_i^(n_i)/n_i! $

con $mu_t = sum_i mu_i$ y $c_t = sum_i n_i$: lo esperado y lo observado, acumulados hasta $t$. El último factor no depende de RR, así que se cancela en cualquier cociente.

El cociente de verosimilitud contra $H_0$ ($"RR"=1$):

$ "LR"_t("RR") = L_t("RR")/L_t(1) = e^((1-"RR")mu_t) "RR"^(c_t) quad => quad "LLR"_t("RR") = (1-"RR")mu_t + c_t ln("RR") $

== Maximizando sobre RR

MaxSPRT no fija RR: lo maximiza sobre $"RR">1$. Derivando $"LLR"_t("RR")$ e igualando a cero, el máximo se da en el RR que coincide con lo observado sobre lo esperado:

$ (partial "LLR"_t)/(partial "RR") = -mu_t + c_t/"RR" = 0 quad => quad hat("RR")_t = c_t/mu_t $

Reemplazando $hat("RR")_t$ en $"LLR"_t("RR")$ se obtiene el estadístico final del test.

== MaxSPRT para datos Poisson

$ "LLR"_t = (mu_t - c_t) + c_t ln(c_t / mu_t) $

y $"LLR"_t = 0$ si $c_t < mu_t$ (no hay evidencia de exceso).

Regla, monitoreada continuamente:

- se sigue vigilando mientras $"LLR"_t$ no llegue a un valor crítico $V$
- se rechaza $H_0$ apenas $"LLR"_t >= V$ (señal)
- ya no hay límite inferior para "aceptar" $H_A$ - no interesa frenar si la vacuna resulta protectora
- sí hay un límite superior $T$ en la duración de la vigilancia (cantidad esperada de eventos bajo $H_0$), para no vigilar para siempre

== Valores críticos y ejemplo numérico

$V$ se elige para controlar el error de tipo I: es el valor más chico tal que, bajo $H_0$, la probabilidad de que $"LLR"_t$ llegue a cruzarlo alguna vez antes de $T$ sea exactamente $alpha$.

$ P_(H_0) (exists thin t <= T : "LLR"_t >= V) = alpha $

No hay fórmula cerrada: se calcula enumerando numéricamente las trayectorias posibles bajo $H_0$ y ajustando $V$ hasta lograr ese $alpha$ (el paper itera por interpolación). Una vez tabulado para cada par $(alpha, T)$, nadie más necesita recalcularlo.

== Valores críticos y ejemplo numérico

#align(center)[
  #image("figures/sprt_trayectoria.png", width: 65%)
]

== Un ejemplo concreto

Para $alpha = 0.05$, algunos valores críticos de la Tabla 1 del paper:

- $T=2$ → $V approx 3.05$
- $T=10$ → $V approx 3.47$
- $T=800$ (el caso de fiebre en Pediarix, más adelante) → $V approx 4.29$

A mayor $T$, mayor $V$: se permiten más chances de cruzarlo por azar, así que hace falta una señal más fuerte para mantener el mismo $alpha$.

Con $T=2$ y $mu_t=1$: si $c_t=4$, $"LLR"_t = (1-4)+4ln(4) approx 2.55$, todavía no alcanza $V=3.05$. Con $c_t=5$: $"LLR"_t=(1-5)+5ln(5) approx 4.05$, cruza el límite y se rechaza $H_0$.

== Potencia y el trade-off de siempre

La potencia depende del RR verdadero y de $T$ (cuánto se está dispuesto a vigilar).

- $T$ chico: se termina antes, pero cuesta detectar riesgos moderados.
- $T$ grande: más potencia, pero hay que estar dispuesto a esperar más si la señal tarda.

No es gratis: es el mismo trade-off de siempre entre tamaño de muestra y potencia, solo que acá se define en términos de eventos esperados bajo $H_0$ en vez de un $n$ fijo.

Como es vigilancia observacional (los datos se juntan igual, haya o no señal), estirar $T$ no cuesta mucho más que tiempo de cómputo.

== Y si no tenemos $mu_t$ confiable...

A veces no hay una tasa basal esperada bien establecida. Ahí conviene comparar directamente tiempo expuesto contra tiempo no expuesto: un diseño pareado, en vez de Poisson.

- autocontrolado: mismo individuo, ventana expuesta vs. no expuesta
- con controles: individuos expuestos vs. no expuestos, pareados por edad/sexo

Cada evento adverso es como tirar una moneda: cae del lado expuesto o del lado no expuesto. Sea $z$ la duración del período pareado no expuesto sobre la del expuesto (por ejemplo, $z=1$ para un pareo 1:1). Bajo $H_0$, $P("expuesto") = p = 1\/(z+1)$.

Como ya no hay una tasa Poisson de la que depender, el límite de vigilancia ya no se mide en eventos esperados ($T$) sino en eventos observados en total ($N$).

== MaxSPRT para datos binomiales

Sea $n$ el total de eventos observados hasta el momento y $c_n <= n$ los que cayeron del lado expuesto. Igual que en el caso Poisson, se maximiza la verosimilitud sobre $"RR">1$ (estimador $hat("RR")=z c_n\/(n-c_n)$), y se llega a:

$ "LLR"_n = c_n ln(c_n/n) + (n-c_n)ln((n-c_n)/n) - c_n ln(1/(z+1)) - (n-c_n) ln(z/(z+1)) $

(si $z c_n\/(n-c_n) > 1$; si no, $"LLR"_n = 0$).

Misma lógica que el caso Poisson —alternativa compuesta, maximización, tablas de valores críticos— pero con la binomial como distribución de base.

== Volviendo a Pediarix, ahora con MaxSPRT

Se aplica el MaxSPRT Poisson a los mismos datos, con $T$ = 800 eventos esperados para fiebre y 15 para síntomas neurológicos (≈2 años de vigilancia), $alpha=0.05$:

- *Fiebre*: señal a las *13 semanas* (97 casos observados vs. 69.7 esperados, $hat("RR")=1.39$, $"LLR"=4.78$)
- *Neurológicos*: señal a las *42 semanas* (15 casos vs. 5.5 esperados, $hat("RR")=2.7$, $"LLR"=5.51$)

En ambos casos el LLR observado queda por encima del $V approx 4.29$ (fiebre) y $V approx 3.56$ (neurológicos) que corresponden a esos $T$ - coherente con la señal.

Con $alpha=0.01$ o con $T$ más chico (≈3 meses) los resultados casi no cambian: 13 semanas para fiebre siempre, 32 a 42 semanas para neurológicos según el caso.

== Pediarix: qué significan estos resultados

Ya no hace falta elegir un RR de antemano ni justificar por qué se descartó tal o cual valor: el resultado no depende de una decisión arbitraria previa a mirar los datos.

- *Fiebre*: a las 82 semanas el riesgo relativo observado fue 1.16 (16% de exceso). Es un efecto secundario ya conocido de Pediarix, así que no sorprende.
- *Neurológicos*: el exceso de casos resultó estar ligado, al menos en parte, a un cambio en cómo se completaban los formularios de historia clínica - no necesariamente un efecto real de la vacuna.

Una señal estadística pide una investigación epidemiológica, no la reemplaza. (En la práctica se haría un solo análisis con parámetros fijados de antemano; acá se muestran varios nada más para comparar métodos.)

== Algunas cosas para tener en cuenta

- Una señal estadística no es lo mismo que una relación causal probada: hay que investigarla clínicamente.
- El sesgo depende mucho de qué grupo de comparación se use (histórico, pareado, autocontrol) - cada uno tiene sus propios problemas.
- El método asume homogeneidad de Poisson o proporción fija conocida bajo $H_0$; si esos supuestos fallan, las señales pueden ser espurias.
- MaxSPRT ya se usaba en producción en el momento de publicarse el paper (Vaccine Safety Datalink, CDC), no era solo una propuesta teórica.

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
