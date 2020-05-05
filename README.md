# Modelagem Matemática do COVID-19

Grupo composto por alunos do Instituto de Matemática da UFRJ, com orientação de professores, com o objetivo de modelar a disseminação do Covid-19  na cidade do Rio de Janeiro.
Por meio do modelo epidemiológico SEIR-QAD ([1] Jia et. al.), utilizamos dados do Rio de Janeiro para encontrar a curva mais próxima dos dados por meio de técnicas de otimização ([2] DiffEqFlux.jl) e variamos os parâmetros do modelo para analisar outros cenários.
O modelo de equações diferenciais foi modelado com auxílio do pacote ([3] DifferentialEquations.jl)
Os resultados podem ser vistos no notebook Higher Order Bayes rules for ODEs.

# Grupo

- Rafael Benchimol Klausner
- Vitor Luiz Pinto de Piña Ferreira

# Próximos passos

- Utilizar quantificação de incertezas para estimar os parâmetros

# Referências

- ![[1] Jia et. al., 2020](https://arxiv.org/abs/2003.02985)
- ![[2] DiffEqFlux.jl](https://arxiv.org/abs/1902.02376)
- ![[3] DifferentialEquations.jl](https://zenodo.org/record/3740780#.XrGHlahKiUk)
