# Simulação de Sistemas de Controle — Controle e Servomecanismos

**João Pedro Pinheiro · Suelen F. Paulino · Priscila O. Valente**

Disciplina de Controle e Servomecanismos — Engenharia de Computação — FACOM/UFMS — 2026.2
Prof. Dr. Victor Leonardo Yoshimura

Repositório com os scripts das listas de exercícios de simulação da disciplina,
organizados em uma pasta por lista dentro de `scripts/`.

## Índice

- [Lista de Exercícios 4](#lista-de-exercícios-4)
  - [Sobre](#sobre)
  - [Scripts](#scripts)
- [Lista de Exercícios 2](#lista-de-exercícios-2)
  - [Sobre](#sobre-1)
  - [Scripts](#scripts-1)

---

## Lista de Exercícios 4

### Sobre

Scripts utilizados na quarta lista de exercícios, feita no **Scilab**. A lista tem três
exercícios:

1. Sistema massa-mola-amortecedor com a base submetida a oscilações senoidais: escolha
   de $k$ para $0{,}4 < \zeta < 0{,}8$, ganhos em $10$ e $20\,\mathrm{rad/s}$ com `horner`,
   diagramas de Bode e simulação das oscilações.
2. Sistema em realimentação unitária com planta $\hat g(s) = 1/(s(s+1)(s+2))$ e
   controlador $\hat c(s) = 6$: diagramas de Nyquist e de Black mostrando a estabilidade
   marginal e simulação da resposta à entrada $u = \sin(\sqrt{2}\,t)$, com a escolha do
   passo de tempo pelo teorema da amostragem.
3. Dois tanques comunicantes em espaço de estados: controlabilidade e observabilidade,
   realimentação de estados, realimentação dinâmica de saída com observador e funções de
   Lyapunov.

### Scripts

| Arquivo | Conteúdo |
| --- | --- |
| `scripts/lista_4/ex2_scripts.sce` | Margens de ganho e de fase, diagramas de Nyquist (`show_margins`) e de Black, simulação em malha fechada para $u = \sin(\sqrt{2}\,t)$ com a envoltória $\pm 0{,}64\,t$ e comparação de três períodos de amostragem |

---

## Lista de Exercícios 2

### Sobre

Scripts utilizados na segunda lista de exercícios de simulação. Todas as simulações
foram feitas no **Scilab 2026.1.0**, usando o **Xcos** nos exercícios em que o enunciado
pede a montagem por diagrama de blocos.

A lista tem quatro exercícios:

1. Sistema mecânico massa-mola-amortecedor, em espaço de estados (Xcos).
2. Dois tanques comunicantes, montados apenas com blocos de soma, ganho e
   integradores (Xcos), com a verificação do balanço de massa no Scilab.
3. Redução de ordem de uma planta de sexta ordem para quinta, quarta e terceira
   ordens pelo método dos polos dominantes, com comparação das respostas a
   entradas $u = t^i$.
4. Projeto de compensadores para uma planta de sexta ordem em realimentação
   unitária, pela equação diofantina e pelo lugar das raízes, incluindo o projeto
   sobre um modelo reduzido e a simulação dos três resultados na planta completa.

### Scripts

| Arquivo | Conteúdo |
| --- | --- |
| `scripts/lista_2/ex2_tanques.sce` | Gráficos dos níveis, das vazões e a verificação de $q - q_2 = C_1\dot h_1 + C_2\dot h_2$ |
| `scripts/lista_2/ex3_reducao_ordem.sce` | Reduções de 5ª, 4ª e 3ª ordens e comparação com a planta completa |
| `scripts/lista_2/ex4_regiao_omega.sce` | Região $\Omega$ de desempenho garantido |
| `scripts/lista_2/ex4_diofantina.sce` | Projeto do compensador pela equação diofantina |
| `scripts/lista_2/ex4_lugar_raizes.sce` | Traçado do lugar das raízes da planta completa |
| `scripts/lista_2/ex4_lugar_raizes_reduzida.sce` | Redução de ordem da planta e lugar das raízes do modelo reduzido |
| `scripts/lista_2/ex4_projeto_planta_e_planta_reduzida.sce` | Projeto dos compensadores por avanço de fase e simulação na planta completa |

Os scripts dos exercícios 1 e 2 dependem dos diagramas montados no Xcos: rode primeiro
a simulação do diagrama, que exporta os sinais para o Scilab, e depois o script.
