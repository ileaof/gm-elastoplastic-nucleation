<h1 align="center">Formulação completa</h1>

<p align="center"><a href="FORMULATION.md">English</a> · <b>Português</b> · <a href="../README.pt-BR.md">voltar ao README</a></p>

Formulação de Ferreira (2024, *Physica B* 695, 416494; MRS Meeting 2026) do tensor de campo térmico na nucleação,
com a superfície de Gurtin–Murdoch estendida à elasto-plasticidade neste repositório, e a sua propagação ao
espaçamento dendrítico secundário (SDAS) pelo modelo de Rappaz–Boettinger modificado (MRB). As seções 1–4 seguem a
forma enviada pelo autor; as seções 5–8 transcrevem o que os scripts de `model/` calculam.

## Sumário

1. [Campo térmico a partir da 1ª lei](#1-campo-térmico-a-partir-da-1ª-lei)
2. [Tensor de campo térmico](#2-tensor-de-campo-térmico)
3. [Nucleação: energia livre, raios críticos e tensores](#3-nucleação-energia-livre-raios-críticos-e-tensores)
4. [Taxa de nucleação](#4-taxa-de-nucleação)
5. [Energia de superfície: Gurtin–Murdoch](#5-energia-de-superfície-gurtinmurdoch)
6. [Superfície elasto-plástica](#6-superfície-elasto-plástica)
7. [Solução numérica em função do gradiente térmico](#7-solução-numérica-em-função-do-gradiente-térmico)
8. [Espaçamento dendrítico secundário (MRB)](#8-espaçamento-dendrítico-secundário-mrb)
9. [Símbolos](#9-símbolos)

**Notação.** $`\mathbf{T}`$ (negrito) é o campo térmico, em K, e $`\nabla\mathbf{T}`$ o seu gradiente; $`T`$ é a
temperatura e $`\Delta T`$ o sub-resfriamento. $`\boldsymbol{\Gamma}`$ é o tensor de campo térmico (coeficiente de
Gibbs–Thomson generalizado) e $`\mathbf{A}`$ a área da nova superfície, normal ao campo $`\mathbf{T}`$.

---

## 1. Campo térmico a partir da 1ª lei

Pela 1ª lei, a energia $`E`$ é função de ponto, e $`Q`$ e $`W`$ dependem do caminho:

```math
\delta Q - \delta W = dE
```

Tomando $`dE \equiv \partial E_{,i}`$ em relação à variável primitiva $`T`$, definem-se dois gradientes de campo
térmico:

```math
\nabla\mathbf{T}_{\mathrm{pi}} = \nabla\left[\frac{\partial T}{\delta Q}\,(Q-W)\right] = \nabla\left[\frac{\partial T}{\delta Q}\,E\right]
\qquad\text{(independente do caminho)}
```

```math
\nabla\mathbf{T}_{\mathrm{pd}} = \nabla\left[\frac{\partial T}{\delta Q}\,Q\right]
= \nabla\mathbf{T}_{\mathrm{pi}} + \nabla\left[\frac{\partial T}{\delta Q}\,W\right]
\qquad\text{(dependente do caminho)}
```

Como $`E`$ e $`W`$ são independentes, os termos cruzados são desprezados:

```math
\nabla\mathbf{T} = \nabla\left[\frac{\partial T}{\delta Q}\,(E+W)\right]
\cong \nabla\left[\frac{\partial T}{\partial E}\,E\right] + \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

Com $`E = U + K + E_P + \Sigma_{other}`$ e a energia interna aberta nas suas formas,

```math
\nabla\mathbf{T} \cong
\nabla\left[\frac{\partial T}{\partial U}\left(U_{transl} + U_{vib} + U_{rot} + U_{elect} + U_{mag} + U_{other}\right)\right]
+ \nabla\left[\frac{\partial T}{\partial K}\,K\right]
+ \nabla\left[\frac{\partial T}{\partial E_P}\,E_P\right]
+ \nabla\left[\frac{\partial T}{\partial \Sigma_{other}}\,\Sigma_{other}\right]
+ \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

Só o gradiente do campo térmico desloca o equilíbrio sólido–líquido. No experimento de solidificação transiente,
$`\nabla\mathbf{T}`$ na isoterma liquidus é $`G_L = \dot T / V_L = T_R / V_L`$.

## 2. Tensor de campo térmico

Toda nova superfície criada ou deformada, de área $`\mathbf{A}`$, é normal ao campo térmico:

```math
\boldsymbol{\Gamma} = \mathbf{A}\cdot\nabla\mathbf{T} \cong \mathbf{A}\cdot\nabla\left[\frac{\partial T}{\partial E}\,E\right]
```

## 3. Nucleação: energia livre, raios críticos e tensores

Energia livre total do núcleo (calota de ângulo $`\theta`$), com $`\Delta G_V = \Delta S_V\,\Delta T`$ e
$`\Delta G_S = \gamma_{SL}`$, sendo $`\Delta S_V`$, $`\Delta T`$, $`\gamma_{SL}`$ e $`\theta`$ funções do raio:

```math
\Delta G = \left(\frac{1}{3}\pi r^3\,\Delta S_V\,\Delta T + \pi r^2\,\gamma_{SL}\right) f(\theta),
\qquad f(\theta) = 2 - 3\cos\theta + \cos^3\theta
```

O máximo em relação a $`r`$, $`\partial\Delta G/\partial r = 0`$, dá a equação parabólica

```math
\left[\left(\frac{\partial \Delta S_V}{\partial r}\Delta T + \Delta S_V\frac{\partial \Delta T}{\partial r}\right) f(\theta)
+ \Delta S_V\,\Delta T\,\frac{\partial f(\theta)}{\partial r}\right] r_C^2
+ 3\left[\left(\Delta S_V\,\Delta T + \frac{\partial \gamma_{SL}}{\partial r}\right) f(\theta)
+ \gamma_{SL}\frac{\partial f(\theta)}{\partial r}\right] r_C + 6\gamma_{SL}\, f(\theta) = 0
```

**Raio crítico de 2ª ordem**, em termos das energias livres de volume, de superfície e configuracional
($`\Delta G_C = \gamma_{SL}\,(\partial f/\partial r)/f`$):

```math
r_C^{2nd} = -\frac{3}{2}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
{\dfrac{\partial\left(\Delta S_V\,\Delta T\, f(\theta)\right)}{\partial r}}
= -\frac{3}{2}\,\frac{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}{\dfrac{\partial \Delta G_V}{\partial r}}
```

**Raio crítico de 1ª ordem** e a teoria clássica (CNT), recuperada quando $`\partial\Delta G_S/\partial r = \Delta G_C = 0`$:

```math
r_C^{1st} = -\frac{2\,\Delta G_S}{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C},
\qquad r_C^{CNT} = -\frac{2\,\Delta G_S}{\Delta G_V}
```

**Tensor de campo térmico de 2ª ordem:**

```math
\boldsymbol{\Gamma}^{2nd} = -\frac{3}{4}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
{\dfrac{1}{\Delta T}\dfrac{\partial\left(\Delta S_V\,\Delta T\, f(\theta)\right)}{\partial r}}
= -\frac{3}{4}\,\frac{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}{\dfrac{\partial \Delta S_V}{\partial r}}
= \mathbf{A}\cdot\nabla\mathbf{T}
```

**Tensor de campo térmico de 1ª ordem**, com entropias $`S \equiv G/\Delta T`$:

```math
\boldsymbol{\Gamma}^{1st} = -\frac{\Delta G_S\,\Delta T}{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}
= -\frac{\Delta G_S}{\Delta S_V + \dfrac{\partial \Delta S_S}{\partial r} + \Delta S_C}
```

Os dois têm a forma de Gibbs–Thomson:

```math
\boldsymbol{\Gamma}^{1st} = \tfrac{1}{2}\,\Delta T\, r_C^{1st}, \qquad \boldsymbol{\Gamma}^{2nd} = \tfrac{1}{2}\,\Delta T\, r_C^{2nd}
```

O $`\Gamma^{1st}`$ contém a entropia de superfície $`\Delta S_S = (1/\Delta T)\,\partial\gamma_{SL}/\partial r`$, isto é, a
derivada de $`\gamma_{SL}`$; o $`\Gamma^{2nd}`$ depende de $`\gamma_{SL}`$.

## 4. Taxa de nucleação

```math
I = \frac{D\,A}{\lambda^4}\,\frac{N}{V}\,\exp\left(\frac{\Delta G_C}{\Delta G_{C,Eq}}\right)
```

com $`D`$ o coeficiente de autodifusão, $`N/V`$ a densidade de modos da nova fase, $`\lambda`$ a distância de salto até a
nova superfície e $`A`$ a área da nova fase.

## 5. Energia de superfície: Gurtin–Murdoch

Na superfície $`\mathcal{S}`$ do núcleo, com gradiente de superfície $`\nabla_S`$ e projetor tangente
$`\mathbf{P} = \mathbf{I} - \mathbf{n}\otimes\mathbf{n}`$ (Gurtin e Murdoch, 1978; Ferreira, 2022):

```math
\nabla_S \cdot \boldsymbol{\Sigma} = \boldsymbol{\Theta}\,\mathbf{n} + \rho_0\,\ddot{\mathbf{u}}
```

```math
\boldsymbol{\Sigma} = \sigma\,\mathbf{I} + 2\,(\mu_0 - \sigma)\,\mathbf{E} + (\lambda_0 - \sigma)\,(\mathrm{tr}\,\mathbf{E})\,\mathbf{I} + \sigma\,\nabla_S\mathbf{u}
\qquad\text{em } \mathcal{S},
\qquad \mathbf{E} = \tfrac{1}{2}\left(D\mathbf{u} + D\mathbf{u}^{T}\right)
```

com $`\sigma`$ a tensão superficial residual, $`\rho_0`$ a massa específica de superfície e $`\lambda_0`$, $`\mu_0`$ os módulos
de Lamé da superfície. Tensão e energia de superfície se relacionam pela equação de Shuttleworth,

```math
\Sigma_{\alpha\beta} = \gamma\,\delta_{\alpha\beta} + \frac{\partial \gamma}{\partial \varepsilon_{\alpha\beta}},
```

e a energia livre de superfície na forma métrica (escalar de Ricci $`R = 2/r_0^2`$ na 2-esfera) é

```math
\Delta G_S = \gamma_0\,\Phi\!\left(r_0^2 R\right) - \frac{1}{2A}\int_{\mathcal{S}} \Sigma_{ij}\,g^{ij}\,dA .
```

**Solução esférica usada nos scripts.** Com as constantes de Lamé do volume $`\lambda`$, $`\mu`$:

```math
y(r) = (3\lambda + 2\mu)\,r + 2(\lambda_0 + \mu_0), \qquad y_0 = y(r_{eq})
```

```math
s(r) = -2\,\sigma_0 \ln\frac{y_0}{y}, \qquad \tau(r) = -\sigma_0 + s(r), \qquad
\gamma_{SL}(r) = \gamma_0\left(\frac{r_{eq}}{r}\right)^{2} - s(r)
```

```math
\frac{\partial \gamma_{SL}}{\partial r} = -\frac{2\gamma_0\,r_{eq}^2}{r^3} - \frac{2\sigma_0}{r\,(1+\alpha)},
\qquad \alpha = \frac{2\left(\sigma_0 + 2(\lambda_0+\mu_0)\right)}{r\,(3\lambda + 2\mu)}
```

$`\sigma_0 < 0`$ é uma força compressiva; nos scripts `sigma_0` guarda o seu módulo, de modo que $`\tau(r_{eq}) = -\sigma_0`$.
$`\gamma_0`$, $`\sigma_0`$ e $`r_{eq}`$ são valores de equilíbrio (sem deformação da área criada). Com
$`E_A = 2\ln(y/y_0)`$, a deformação de área da superfície, a solução elástica é $`s = C_A E_A`$ com $`C_A = \sigma_0`$.

## 6. Superfície elasto-plástica

**Cinemática** (plano tangente) e lei de Gurtin–Murdoch com a parte elástica:

```math
\mathbf{E} = \mathbf{E}^{e} + \mathbf{E}^{p}, \qquad \mathbf{E}^{p} = \mathbf{P}\,\mathbf{E}^{p}\,\mathbf{P},
\qquad
\boldsymbol{\Sigma} = \sigma\,\mathbf{I} + 2\,(\mu_0 - \sigma)\,\mathbf{E}^{e} + (\lambda_0 - \sigma)\,(\mathrm{tr}\,\mathbf{E}^{e})\,\mathbf{I} + \sigma\,\nabla_S\mathbf{u}
```

**Escoamento sobre a tensão em excesso** $`\tilde{\boldsymbol{\Sigma}} = \boldsymbol{\Sigma} - \sigma\,\mathbf{P}`$ (a tensão residual
sozinha não escoa):

```math
f_A = \left|\mathrm{tr}\,\tilde{\boldsymbol{\Sigma}}\right| - 2\,\sigma_y^{s}(\alpha_A) \le 0, \qquad \sigma_y^{s} = \sigma_{y0}^{s} + H_A\,\alpha_A
\qquad\text{(área: interface coerente → semicoerente)}
```

```math
f_D = \left\lVert \mathrm{dev}_S\, \tilde{\boldsymbol{\Sigma}} \right\rVert - \sigma_{yD}(\alpha_D) \le 0, \qquad \sigma_{yD} = \sigma_{yD0} + H_D\,\alpha_D
\qquad\text{(J}_2\text{ no plano tangente)}
```

**Lei de fluxo associativa e condições de Kuhn–Tucker:**

```math
\dot{\mathbf{E}}^{p} = \dot{\lambda}_A\,\mathrm{sign}\!\left(\mathrm{tr}\,\tilde{\boldsymbol{\Sigma}}\right)\tfrac{1}{2}\mathbf{P}
+ \dot{\lambda}_D\,\frac{\mathrm{dev}_S\,\tilde{\boldsymbol{\Sigma}}}{\lVert \mathrm{dev}_S\,\tilde{\boldsymbol{\Sigma}} \rVert},
\qquad \dot{\alpha}_k = \dot{\lambda}_k, \qquad \dot{\lambda}_k \ge 0,\quad f_k \le 0,\quad \dot{\lambda}_k f_k = 0
```

Na esfera, $`\mathbf{E} = \varepsilon\,\mathbf{P}`$ é equibiaxial, $`\mathrm{dev}_S\,\mathbf{E} = \mathbf{0}`$, e só o mecanismo de área atua.

**Return mapping** (passo $`n \rightarrow n+1`$, $`r`$ diminuindo a partir de $`r_{eq}`$):

```math
s^{tr} = C_A\left(E_A - E_{A,n}^{p}\right), \qquad f^{tr} = |s^{tr}| - \left(\sigma_{y0}^{s} + H_A\,\alpha_{A,n}\right)
```

```math
f^{tr} > 0:\quad \Delta\lambda = \frac{f^{tr}}{C_A + H_A},\quad s = s^{tr} - C_A\,\Delta\lambda\,\mathrm{sign}(s^{tr}),\quad
E_{A,n+1}^{p} = E_{A,n}^{p} + \Delta\lambda\,\mathrm{sign}(s^{tr}),\quad \alpha_{A,n+1} = \alpha_{A,n} + \Delta\lambda
```

com tangente consistente $`\partial s/\partial E_A = C_A H_A/(C_A + H_A)`$ e solução fechada para carga monotônica

```math
s = \mathrm{sign}(E_A)\left[\sigma_{y0}^{s} + \frac{H_A\,C_A}{C_A + H_A}\left(|E_A| - \frac{\sigma_{y0}^{s}}{C_A}\right)\right],
\qquad |E_A| > \sigma_{y0}^{s}/C_A .
```

**Área plástica e referência de equilíbrio.** A área plástica desloca a referência em que $`\gamma_0`$, $`\sigma_0`$ e
$`r_{eq}`$ foram definidos:

```math
\gamma_{SL} = \gamma_0\left(\frac{r_{eq,p}}{r}\right)^{2} - s, \qquad r_{eq,p} = r_{eq}\,\exp\!\left(E_A^{p}/2\right),
\qquad
\frac{\partial \gamma_{SL}}{\partial r} = -\frac{2\gamma_0\,r_{eq,p}^2}{r^3} + \frac{\gamma_0\,r_{eq,p}^2}{r^2}\frac{\partial E_A^{p}}{\partial r} - \frac{\partial s}{\partial E_A}\frac{\partial E_A}{\partial r}
```

Com plasticidade perfeita ($`H_A = 0`$), depois do escoamento:

```math
\gamma_{SL} = \gamma_0\,e^{\sigma_y^{s}/C_A} + \sigma_y^{s}, \qquad \Gamma^{2nd} = \frac{\gamma_{SL}}{|\Delta S_V|}
\quad\text{(independentes do raio e do gradiente térmico)}
```

**Energia da superfície plástica.** Energia livre $`\psi = \psi^{e}(\mathbf{E}^{e}) + \psi^{p}(\alpha_A)`$, com
$`\psi^{p} = \tfrac{1}{2}H_A\alpha_A^2`$ e força de encruamento $`K = \partial\psi^{p}/\partial\alpha_A = H_A\alpha_A`$. Pela
desigualdade de Clausius–Duhem, $`D = s\,\dot{E}_A^{p} - K\dot{\alpha}_A = \sigma_{y0}^{s}\,\dot{\alpha}_A \ge 0`$, e

```math
w_p = \int \left(\sigma_{y0}^{s} + H_A\,\alpha_A\right) d\alpha_A = \psi^{p} + q_{diss}, \qquad q_{diss} = \sigma_{y0}^{s}\,\alpha_A .
```

| Estado da superfície (`GM_DISSIPATION_MODE`) | $`\gamma_{SL}`$ |
|---|---|
| carregada, parte elástica (`stored`) | $`\gamma_0\,(r_{eq,p}/r)^2 - s`$ |
| carregada + energia armazenada (`thermo`) | $`\gamma_0\,(r_{eq,p}/r)^2 - s + \psi^{p}`$ |
| carregada, total (`barrier`) | $`\gamma_0\,(r_{eq,p}/r)^2 - s + w_p`$ |
| após a restauração elástica (`restored`) | $`\gamma_0`$ |
| após a restauração + trabalho plástico (`restored_barrier`) | $`\gamma_0 + w_p`$ |

## 7. Solução numérica em função do gradiente térmico

O script percorre o raio a partir do equilíbrio, $`r_i = r_{eq} - i\,\Delta r`$ com $`\Delta r = 0{,}005\,r_{eq}`$, e em cada
passo calcula $`\gamma_{SL}`$ e $`\partial\gamma_{SL}/\partial r`$ (seções 5–6) e, em seguida:

**Sub-resfriamento.** A raiz que o script obtém de `f_DT` coincide com a relação capilar

```math
\Delta T = \frac{2\,\gamma_{SL}}{r\,|\Delta S_V|}, \qquad \Delta S_V = -\frac{\Delta H\,\rho_s}{T_L}
```

**Gradiente térmico** (núcleo homogêneo, $`\theta = \pi`$, $`f = 4`$):

```math
\nabla T = \frac{\Delta T}{8\pi r}
```

**Ângulo de 2ª ordem** $`\theta_2`$: raiz de

```math
-\frac{3}{2}\,\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial\gamma_{SL}}{\partial r}\right)\dfrac{f(\theta)}{4} + \gamma_{SL}\,f'_r(\theta)}
{\dfrac{\partial(\Delta S_V\Delta T)}{\partial r}\dfrac{f(\theta)}{4} + \Delta S_V\,\Delta T\,f'_r(\theta)} - r_{het}(r,\theta) = 0,
\qquad
\frac{\partial(\Delta S_V\Delta T)}{\partial r} = -\frac{3}{2r}\left(\Delta S_V\,\Delta T + \frac{\partial\gamma_{SL}}{\partial r}\right)
```

com $`f'_r(\theta) = -3f(\theta) - (1-\cos\theta)(2-\cos\theta-\cos^2\theta)`$ (função `dfthetadr` do script) e o raio da calota

```math
r_{het}(r,\theta) = r\,\exp\!\left[\int_{\pi}^{\theta}\frac{\sin\theta'\,(1+\cos\theta')}{2-\cos\theta'-\cos^2\theta'}\,d\theta'\right]
```

**Tensor de 2ª ordem heterogêneo** (coluna `GT_het_2nd`):

```math
\nabla T_{het} = \frac{\Delta T}{4\pi\,r_{het}\,(1-\cos\theta_2)}, \qquad
\Gamma^{2nd}_{het} = 2\pi r^2\,(1-\cos\theta_2)\,\nabla T_{het} = \frac{r^2}{r_{het}}\,\frac{\Delta T}{2}
```

Como $`\theta_2 \approx \pi`$, $`r_{het} \approx r`$ e $`\Gamma^{2nd}_{het} \approx r\,\Delta T/2 = \gamma_{SL}/|\Delta S_V| = 4\pi r^2\,\nabla T`$.

**Estado num gradiente experimental.** As grandezas são interpoladas linearmente na curva $`\nabla T(r)`$ no gradiente
medido (médio ou local).

## 8. Espaçamento dendrítico secundário (MRB)

Para $`n`$ solutos (Al-3Cu-5Nb-0,1Fe: Cu, Nb, Fe), o modelo de Rappaz–Boettinger modificado é resolvido por iteração
de ponto fixo em SDAS:

```math
\mathrm{SDAS} = 5{,}5\,\left(M\,t_{SL}\right)^{1/3}, \qquad
M = \bar{\beta}\left(-\frac{\Gamma}{p_1}\right)\ln\frac{s_2}{t_3}
```

```math
p_1 = \sum_i \frac{m_{L,i}(1-k_{ef,i})(C_{f,i}-C_{0,i})}{D_{L,i}}, \qquad
s_2 = \sum_i \frac{m_{L,i}(1-k_{ef,i})\,C_{f,i}}{D_{L,i}}, \qquad
t_3 = \sum_i \frac{m_{L,i}(1-k_{ef,i})\,C_{0,i}}{D_{L,i}}
```

**Retrodifusão no sólido:**

```math
Fo_i = \frac{D_{S,i}\,t_{SL}}{(\mathrm{SDAS}/2)^2}, \qquad
\gamma_i = \frac{4\,Fo_i\,k_{0,i}}{4\,Fo_i\,k_{0,i} + 1}, \qquad
\beta_i = \frac{2\,\gamma_i\,Fo_i}{2\,Fo_i + \gamma_i}, \qquad
\bar{\beta} = \sum_i w_i\,\beta_i, \quad w_i = \frac{C_{0,i}}{\sum_j C_{0,j}}
```

(no script, $`\beta_3`$ usa o número de Fourier do soluto 2: $`\beta_3 = 2\gamma_3 Fo_2/(2Fo_3 + \gamma_3)`$).

**Partição efetiva** com o comprimento médio de difusão $`\bar\delta_i`$ (média de $`D_{L,i}/V`$ na faixa de velocidades):

```math
k_{ef,i} = \frac{k_{0,i}}{k_{0,i} + (1-k_{0,i})\exp\!\left(-V\,\bar\delta_i / D_{L,i}\right)}, \qquad
D_{L,i} = D_{0L,i}\,e^{-Q_{L,i}/RT_L}, \quad D_{S,i} = D_{0S,i}\,e^{-Q_{S,i}/RT_S}
```

O modelo RB original é recuperado com $`\bar\beta = 1`$ e $`k_{ef} = k_0`$.

**Leis do lingote** (Al-3Cu-5Nb-0,1Fe, horizontal; P em mm):

```math
V_L = 3{,}4\,P^{-0,56}\ \mathrm{mm\,s^{-1}}, \qquad t_{SL} = 4{,}93116\,P^{0,87606}\ \mathrm{s}, \qquad
G_L(P) = 10405{,}54 - 200{,}22\,P + 1{,}04692\,P^2\ \mathrm{K\,m^{-1}}
```

**SDAS explícito no gradiente térmico.** Invertendo $`G_L(P)`$ (monótono para $`P < 95{,}6`$ mm):

```math
P(\nabla T_L) = \frac{200{,}22 - \sqrt{40088{,}05 - 4{,}18768\,(10405{,}54 - \nabla T_L)}}{2{,}09384}, \qquad
t_{SL}(\nabla T_L) = 4{,}93116\,P(\nabla T_L)^{0,87606}
```

```math
\mathrm{SDAS}(\nabla T_L) = 5{,}5\left[M\,4{,}93116\,P(\nabla T_L)^{0,87606}\right]^{1/3}
```

e, com $`t_{SL} = \Delta T_{SL}/(\nabla T_L\,V_L)`$, a forma com expoente $`-1/3`$:

```math
\mathrm{SDAS} = 5{,}5\left[\frac{M\,\Delta T_{SL}}{\nabla T_L\,V_L}\right]^{1/3} = K\left(\nabla T_L\,V_L\right)^{-1/3}
```

**Tensor local.** Com $`\Gamma`$ calculado no gradiente medido em cada posição, $`\Gamma(P) = \Gamma^{2nd}(\nabla T_L(P))`$ entra
em $`M`$ ponto a ponto. Na superfície elástica, $`\Gamma^{2nd} \approx 4\pi r^2\,\nabla T`$ é proporcional ao gradiente;
na superfície elasto-plástica com referência atualizada e $`H_A = 0`$, $`\Gamma^{2nd}`$ é constante depois do escoamento.

## 9. Símbolos

| Símbolo | Significado | Unidade |
|---|---|---|
| $`\mathbf{T}`$, $`\nabla\mathbf{T}`$, $`G_L`$ | campo térmico, seu gradiente, gradiente na isoterma liquidus | K, K m⁻¹ |
| $`\Delta T`$ | sub-resfriamento | K |
| $`\Delta S_V`$ | entropia volumétrica de fusão | J m⁻³ K⁻¹ |
| $`\gamma_{SL}`$, $`\gamma_0`$ | energia de superfície sólido–líquido e o seu valor de equilíbrio | J m⁻² |
| $`\sigma_0`$, $`\tau`$, $`s`$ | tensão superficial residual, tensão superficial, tensão em excesso | N m⁻¹ |
| $`\lambda_0`$, $`\mu_0`$; $`\lambda`$, $`\mu`$ | módulos de Lamé da superfície; do volume | N m⁻¹; Pa |
| $`r`$, $`r_C`$, $`r_{eq}`$, $`r_{het}`$ | raio, raio crítico, raio de equilíbrio, raio da calota | m |
| $`\theta`$, $`f(\theta)`$ | ângulo de contato e fator geométrico | rad, – |
| $`\Gamma^{1st}`$, $`\Gamma^{2nd}`$ | tensores de campo térmico de 1ª e 2ª ordem | m K |
| $`E_A`$, $`E_A^p`$, $`\alpha_A`$ | deformação de área total e plástica, deformação plástica acumulada | – |
| $`\sigma_{y0}^s`$, $`H_A`$, $`C_A`$ | tensão de escoamento de área, endurecimento, módulo de área | N m⁻¹ |
| $`w_p`$, $`\psi^p`$, $`q_{diss}`$ | trabalho plástico, energia armazenada, calor dissipado | J m⁻² |
| $`t_{SL}`$, $`V_L`$, $`T_R`$ | tempo local de solidificação, velocidade da ponta, taxa de resfriamento | s, m s⁻¹, K s⁻¹ |
| $`m_{L,i}`$, $`k_{0,i}`$, $`k_{ef,i}`$ | inclinação liquidus, partição de equilíbrio e efetiva | K wt.%⁻¹, – |
| $`C_{0,i}`$, $`C_{f,i}`$ | concentração inicial e final do líquido | wt.% |
| $`D_{L,i}`$, $`D_{S,i}`$ | difusividade no líquido e no sólido | m² s⁻¹ |
| $`Fo_i`$, $`\beta_i`$, $`\bar\beta`$ | número de Fourier, parâmetro de retrodifusão e sua média ponderada | – |
