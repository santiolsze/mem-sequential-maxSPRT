#import "deck_setup.typ": deck
#show: deck

= MaxSPRT: vigilancia secuencial de seguridad de vacunas

Un test de razón de verosimilitud secuencial con alternativa compuesta

Santiago Olszevicki

Basado en Kulldorff, Davis, Kolczak, Lewis, Lieu y Platt (2011), _A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance_, Sequential Analysis 30(1).

== El problema

Los ensayos clínicos (fase 2/3) no alcanzan para detectar eventos adversos raros: la muestra es chica y la población suele ser más homogénea que la real.

Por eso, una vez que un fármaco o vacuna sale al mercado, hace falta seguir vigilando.

- Datos disponibles: reclamos de seguros de salud, historias clínicas electrónicas.
- Se quiere revisar la evidencia seguido (semanal, casi continuo), no una sola vez al final.
- Mirar los datos muchas veces sin corregir infla los falsos positivos. Hace falta un método pensado para eso.

Ahí entra el análisis secuencial.

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

Datos de Vaccine Safety Datalink, ~650.000 chicos. Se mira si aumenta el riesgo de fiebre o de síntomas neurológicos en los 28 días post-vacunación.

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

== Cociente de verosimilitud Poisson

Bajo $H_0$, $C_t tilde "Poisson"(mu_t)$; bajo $H_A$, $C_t tilde "Poisson"("RR" mu_t)$. El cociente de verosimilitudes ($c_t !$ y $mu_t^(c_t)$ se cancelan):

$ "LR"_t = (e^(-"RR" mu_t)("RR" mu_t)^(c_t) \/ c_t !) / (e^(-mu_t) mu_t^(c_t) \/ c_t !) = e^((1-"RR")mu_t) "RR"^(c_t) $

En log-verosimilitud, para un RR fijo:

$ "LLR"_t ("RR") = (1-"RR")mu_t + c_t ln("RR") $

MaxSPRT no fija RR: lo maximiza. Derivando e igualando a cero, el máximo se da en el RR que coincide con lo observado sobre lo esperado:

$ (partial "LLR"_t)/(partial "RR") = -mu_t + c_t/"RR" = 0 quad => quad hat("RR")_t = c_t/mu_t $

== MaxSPRT para datos Poisson

Reemplazando $hat("RR")_t = c_t\/mu_t$ en el $"LLR"_t$ anterior, se obtiene el estadístico del test (si $c_t >= mu_t$):

$ "LLR"_t = (mu_t - c_t) + c_t ln(c_t / mu_t) $

y $"LLR"_t = 0$ si $c_t < mu_t$ (no hay evidencia de exceso).

Regla:

- se sigue vigilando mientras $"LLR"_t$ no llegue a un valor crítico $V$
- se rechaza $H_0$ apenas $"LLR"_t >= V$ (señal)
- ya no hay límite inferior para "aceptar" $H_A$ - no interesa frenar si la vacuna resulta protectora
- sí hay un límite superior $T$ en la duración de la vigilancia (cantidad esperada de eventos bajo $H_0$), para no vigilar para siempre

== Valores críticos y ejemplo numérico

Los valores críticos $V$ se calculan numéricamente (no hay fórmula cerrada), tabulados para combinaciones de $alpha$ y $T$. El paper da estas tablas para que nadie tenga que recalcularlas.

== Valores críticos y ejemplo numérico

#align(center)[
  #image("figures/sprt_trayectoria.png", width: 65%)
]

== Potencia y el trade-off de siempre

La potencia depende del RR verdadero y de $T$ (cuánto se está dispuesto a vigilar).

- $T$ chico: se termina antes, pero cuesta detectar riesgos moderados.
- $T$ grande: más potencia, pero hay que estar dispuesto a esperar más si la señal tarda.

No es gratis: es el mismo trade-off de siempre entre tamaño de muestra y potencia, solo que acá se define en términos de eventos esperados bajo $H_0$ en vez de un $n$ fijo.

Como es vigilancia observacional (los datos se juntan igual, haya o no señal), estirar $T$ no cuesta mucho más que tiempo de cómputo.

== Y si no tenemos $mu_t$ confiable...

A veces no hay una tasa esperada poblacional bien establecida. Ahí se usa un diseño de comparación directa (binomial) en vez de Poisson:

- controles emparejados (mismo individuo, ventana expuesta vs. no expuesta) o
- individuos expuestos vs. no expuestos, emparejados por edad/sexo

Bajo $H_0$, la probabilidad de que un evento cayó del lado "expuesto" es $p$ conocido (0.5 si el emparejamiento es 1:1).

$ "LLR"_n = c_n ln(c_n/n) + (n-c_n)ln((n-c_n)/n) - c_n ln(1/(z+1)) - (n-c_n) ln(z/(z+1)) $

Misma lógica que el caso Poisson, mismo tipo de tablas de valores críticos, distinta distribución de base.

== Volviendo a Pediarix, ahora con MaxSPRT

Con $T approx$ 2 años de vigilancia esperada y $alpha=0.05$:

- *Fiebre*: señal a las *13 semanas* (97 casos observados vs. 69.7 esperados, $hat("RR")=1.39$)
- *Neurológicos*: señal a las *42 semanas* (15 casos vs. 5.5 esperados, $hat("RR")=2.7$)

Ya no hace falta elegir un RR de antemano ni justificar por qué se descartó tal o cual valor. El resultado no depende de una decisión arbitraria previa a mirar los datos.

(La fiebre leve es un efecto ya conocido de Pediarix, así que ese resultado no sorprende. Lo de síntomas neurológicos resultó estar ligado, al menos en parte, a un cambio en cómo se registraban los diagnósticos.)

== Algunas cosas para tener en cuenta

- Una señal estadística no es lo mismo que una relación causal probada: hay que investigarla clínicamente.
- El sesgo depende mucho de qué grupo de comparación se use (histórico, emparejado, autocontrol) - cada uno tiene sus propios problemas.
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
