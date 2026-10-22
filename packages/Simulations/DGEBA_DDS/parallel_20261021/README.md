# Campanhas paralelas DGEBA/DDS — 21/10/2026

> **Public release note:** the private cluster-access controllers (`scripts/carbono.py`, `scripts/parallel_control.py`), job-submission logs and status reports were removed from this public repository. The commands below that reference them are kept as documentation of how the campaign was run; LAMMPS inputs, SLURM templates, analysis scripts and the formulation manifest are included.


23 formulações adicionais à referência 100 DGEBA + 50 DDS (job 175498). Cada campanha executa preparação, cura, equilíbrio e triagem de Tg sequencialmente na mesma alocação. A submissão fica registrada em `array_submission.json`; consultar `status.txt` para o último estado coletado.

## Recursos e protocolo

- Partição `grafeno`, conta `users`, QoS `normal`, 1 processo × 16 threads OpenMP, 768 MiB por campanha.
- Array de 23 tarefas, no máximo 23 simultâneas: 368 CPUs adicionais. Com a referência, teto de 384 CPUs. A fila pode disponibilizar menos recursos em cada instante.
- Limite de 36 h por tarefa, incluindo todas as etapas. O limite não é uma previsão de término.
- Preparação igual ao piloto de referência: minimização, aquecimento de 2 ps, NPT de 10 ps e teste reativo de 20 ps.
- Cura ampliada igual à referência: até dois ciclos reativos de 2 ns a 503 K/1 atm, separados por recozimento 503→800 K em 30 ps, NPT 800 K por 2 ns, retorno em 30 ps e NPT 503 K por 100 ps.
- Depois de aprovar a rede: NPT de 2 ns a 503 K/1 atm. Depois de aprovar estacionariedade: aquecimento 30 ps, NPT 800 K por 500 ps e 26 patamares de 800 a 300 K, com 50 ps de acomodação e 50 ps de amostragem cada.
- O protocolo gera estimativas preliminares dependentes da taxa de resfriamento. Uma tentativa pode terminar sem Tg resolvida ou com qualidade insuficiente.

## Como as receitas foram construídas

Fonte: réplica 1 publicada por Livraghi, `data.100dgeba0_50dds_1gcm3_1`. Identificamos moléculas por conectividade (o arquivo publicado usa IDs de molécula zero), selecionamos moléculas inteiras e reatribuímos IDs de átomos, moléculas e interações. Coordenadas, cargas, tipos, coeficientes e caixa dos átomos retidos permanecem iguais. Nenhuma interação ligada cruza a seleção. Não inserimos moléculas nem escalamos distâncias; a preparação NPT acomoda a densidade inicial menor.

As 23 caixas passaram por verificação de contagem/conectividade de moléculas, carga, conservação das propriedades físicas dos átomos retidos e comprimentos periódicos das ligações. A validação de execução cobre três extremos/representantes: D100_H040, D100_H049 e D087_H050. Os resultados desse teste curto são só verificação de execução.

## Conversão e critério de aprovação

Há 2 grupos epóxi por DGEBA e 4 hidrogênios amínicos reativos por DDS. A razão é `r = 4*N_DDS/(2*N_DGEBA)`. O máximo estequiométrico de conversão epóxi é `100*min(1,r)`.

O alvo é **90% da capacidade do reagente limitante**. Logo, com r=0,8, o alvo de conversão epóxi é 72%, pois 90% epóxi seriam impossíveis. Com r≥1, o alvo continua 90% epóxi. O número de eventos é inteiro e arredondado para cima. A conversão realmente atingida é registrada, inclusive eventuais ultrapassagens dentro de um bloco de 100 ps. Não tratar todas as receitas como se tivessem a mesma conversão epóxi.

Critérios estruturais: contagem esperada de átomos, |carga| <0,001 e, alvo atingido, número de reações não maior que a capacidade química e pelo menos 90% dos átomos no maior componente conectado. A conectividade é um critério operacional, não uma validação de percolação/tamanho finito.

Equilíbrio e ajuste de Tg usam os mesmos critérios da campanha original, com falhas explicitamente registradas. Nenhum resultado entra automaticamente no BO (`bo_eligible=false`). A composição, conversão real e protocolo devem acompanhar qualquer futura comparação. Uma realização por receita não quantifica a incerteza entre redes; candidatos promissores precisarão de confirmação com réplicas.

## Matriz submetida

| Índice | Formulação | DGEBA | DDS | Razão H/epóxi | Alvo epóxi (%) | Átomos |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | D100_H040 | 100 | 40 | 0.800000 | 72.00 | 6060 |
| 1 | D100_H041 | 100 | 41 | 0.820000 | 73.80 | 6089 |
| 2 | D100_H042 | 100 | 42 | 0.840000 | 75.60 | 6118 |
| 3 | D100_H043 | 100 | 43 | 0.860000 | 77.40 | 6147 |
| 4 | D100_H044 | 100 | 44 | 0.880000 | 79.20 | 6176 |
| 5 | D100_H045 | 100 | 45 | 0.900000 | 81.00 | 6205 |
| 6 | D100_H046 | 100 | 46 | 0.920000 | 82.80 | 6234 |
| 7 | D100_H047 | 100 | 47 | 0.940000 | 84.60 | 6263 |
| 8 | D100_H048 | 100 | 48 | 0.960000 | 86.40 | 6292 |
| 9 | D100_H049 | 100 | 49 | 0.980000 | 88.20 | 6321 |
| 10 | D099_H050 | 99 | 50 | 1.010101 | 90.00 | 6301 |
| 11 | D098_H050 | 98 | 50 | 1.020408 | 90.00 | 6252 |
| 12 | D097_H050 | 97 | 50 | 1.030928 | 90.00 | 6203 |
| 13 | D096_H050 | 96 | 50 | 1.041667 | 90.00 | 6154 |
| 14 | D095_H050 | 95 | 50 | 1.052632 | 90.00 | 6105 |
| 15 | D094_H050 | 94 | 50 | 1.063830 | 90.00 | 6056 |
| 16 | D093_H050 | 93 | 50 | 1.075269 | 90.00 | 6007 |
| 17 | D092_H050 | 92 | 50 | 1.086957 | 90.00 | 5958 |
| 18 | D091_H050 | 91 | 50 | 1.098901 | 90.00 | 5909 |
| 19 | D090_H050 | 90 | 50 | 1.111111 | 90.00 | 5860 |
| 20 | D089_H050 | 89 | 50 | 1.123596 | 90.00 | 5811 |
| 21 | D088_H050 | 88 | 50 | 1.136364 | 90.00 | 5762 |
| 22 | D087_H050 | 87 | 50 | 1.149425 | 90.00 | 5713 |

## Arquivos e acompanhamento

- `manifest.json`: receitas, parâmetros, auditoria e hashes dos arquivos comuns.
- `cases/ID/case.json`: composição e critérios de cada caixa.
- `cases/ID/data.initial`: estruturas geradas, disponíveis localmente e na Carbono; ignoradas no Git por volume.
- `common/`: cópia isolada dos inputs e scripts usados por esta campanha. A referência em execução não foi alterada.
- `results/`: cópia pontual dos logs e resumos coletados; não é sincronização contínua.

No diretório principal do repositório, com o ambiente Python/Paramiko já usado:

```powershell
& $py DGEBA_DDS_Livraghi/scripts/parallel_control.py status
& $py DGEBA_DDS_Livraghi/scripts/parallel_control.py collect
```

`stage`, `validate` e `submit` criam ou submetem trabalho; não repetir para consultar estado. O controlador impede a repetição das submissões registradas.

Para restaurar estruturas ignoradas em outra cópia, primeiro obter os arquivos fonte conforme o README principal e executar `prepare_parallel.py --restore-missing-data`; o hash da fonte e dos dados regenerados é verificado.

Diretório remoto: `$HOME/formulationlab-dds-livraghi-pilot-20261021/parallel_20261021`. Os jobs escrevem em `runs/ID/`; `phase.txt`, `campaign_summary.json`, `equilibrium_summary.json` e `tg_summary.json` registram as etapas e resultados.
