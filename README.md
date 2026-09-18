# Simulação de Sistemas de Controle — Lista de Exercícios 2

**João Pedro Pinheiro · Suelen F. Paulino · Priscila O. Valente**

Disciplina de Controle e Servomecanismos — Engenharia de Computação — FACOM/UFMS — 2026.2
Prof. Dr. Victor Leonardo Yoshimura

## Sobre

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

## Scripts

| Arquivo | Conteúdo |
| --- | --- |
| `scripts/ex2_tanques.sce` | Gráficos dos níveis, das vazões e a verificação de $q - q_2 = C_1\dot h_1 + C_2\dot h_2$ |
| `scripts/ex3_reducao_ordem.sce` | Reduções de 5ª, 4ª e 3ª ordens e comparação com a planta completa |
| `scripts/ex4_regiao_omega.sce` | Região $\Omega$ de desempenho garantido |
| `scripts/ex4_diofantina.sce` | Projeto do compensador pela equação diofantina |
| `scripts/ex4_lugar_raizes.sce` | Traçado do lugar das raízes da planta completa |
| `scripts/ex4_lugar_raizes_reduzida.sce` | Redução de ordem da planta e lugar das raízes do modelo reduzido |
| `scripts/ex4_projeto_planta_e_planta_reduzida.sce` | Projeto dos compensadores por avanço de fase e simulação na planta completa |

Os scripts dos exercícios 1 e 2 dependem dos diagramas montados no Xcos: rode primeiro
a simulação do diagrama, que exporta os sinais para o Scilab, e depois o script.
