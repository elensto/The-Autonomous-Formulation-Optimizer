# Piloto DGEBA/4,4′-DDS na Carbono

> **Public release note:** the private cluster-access controllers (`scripts/carbono.py`, `scripts/parallel_control.py`), job-submission logs and status reports were removed from this public repository. The commands below that reference them are kept as documentation of how the campaign was run; LAMMPS inputs, SLURM templates, analysis scripts and the formulation manifest are included.


Esta pasta isola a rota de Livraghi et al., **Block Chemistry for Accurate Modeling of Epoxy Resins**, J. Phys. Chem. B 127 (2023), 7648–7662, DOI [10.1021/acs.jpcb.3c04724](https://doi.org/10.1021/acs.jpcb.3c04724). Dados e scripts originais: [Zenodo 7273800](https://zenodo.org/records/7273800), CC BY 4.0.

O objetivo imediato é verificar leitura do campo de forças, estabilidade inicial, atualização de cargas/topologia por `bond/react` e velocidade na Carbono. **O piloto curto não calcula Tg, não estabelece equilíbrio e não gera observações aprovadas para o BO.**

**Estado em 21/10/2026:** piloto 175482 concluído (8% de conversão). Validação de execução de toda a sequência concluída no job 175497. Campanha de cura **175498 em execução**, buscando 90%; propriedades **175499 aguardando seu sucesso**, com critérios adicionais de qualidade. Ver [campanha e próximos passos](reports/CAMPANHA.md) e [resultado do piloto](reports/RESULTADO_PILOTO.md). Ainda não há Tg ou densidade equilibrada aprovadas.

## Organização

- `sources/`: artigo PDF, metadados do Zenodo, cinco arquivos baixados e manifesto com checksums.
- `raw/`: extração dos arquivos originais, sem adaptações.
- `inputs/`: cópias pequenas e adaptadas para execução; não altera os originais.
- `scripts/`: auditoria local, build, SLURM, acesso direto e análise do piloto.
- `reports/`: auditoria estrutural, preflight e IDs de submissão.
- `results/job-ID/`: arquivos coletados do job; ficam fora do Git por padrão.

Os arquivos volumosos permanecem disponíveis localmente, mas são ignorados pelo Git. Nenhum commit ou push é feito automaticamente. Senhas não pertencem a esta pasta.

## Sistema escolhido

Original: `lammps/data/small/data.100dgeba0_50dds_1gcm3_1`.

| Quantidade | Valor |
| --- | --- |
| DGEBA n=0 | 100 moléculas |
| 4,4′-DDS | 50 moléculas |
| Átomos | 6.350 |
| Ligações iniciais | 6.700 |
| Grupos epóxi iniciais | 200 |
| H amínicos reativos | 200 |
| Razão H amínico/epóxi | 1 |
| Réplica | configuração inicial 1 |

O ID de molécula original é zero; a auditoria identifica 100 componentes de 49 átomos e 50 de 29 átomos pela conectividade. A carga total inicial é aproximadamente zero. A densidade inicial próxima de 0,994 g/cm³ é de preparação, não um resultado de equilíbrio. Esta primeira etapa usa uma receita, não uma série de formulações.

## Adaptações documentadas

1. Caminhos absolutos da máquina dos autores substituídos por arquivos locais.
2. Reservado espaço adicional para tipos e topologia que surgem na reação.
3. Coeficientes novos extraídos de `in.min`, sem alterar seus valores.
4. `BondingIDs` renomeado para `InitiatorIDs` nas duas cópias dos mapas.
5. `update_edges charges` substituído por `custom_charges no`, atualizando todas as cargas especificadas no template pós-reação. Conservação da carga é um critério de aprovação do diagnóstico.
6. Mantidos `lj/cut/coul/long 10`, mistura aritmética, `special_bonds amber`, PPPM 1e-4 e estilos harmonic dos autores. Não usar o bloco CHARMM de `lammps_Tg/prototype.in.lmp`.
7. Mantidos 503 K, 1 atm, raios reativos 2,8/3,1 Å e estabilização por 60 passos. O aquecimento curto do piloto parte de 50 K e não reproduz a preparação longa original.
8. Inserida a linha separadora ausente depois do cabeçalho `Dihedrals` no arquivo inicial publicado. Sem ela, o LAMMPS descarta o primeiro registro e conta 19.099 em vez de 19.100 diedros. Nenhum token de coordenada, carga, parâmetro ou topologia é alterado. Essa diferença de formatação foi identificada no primeiro job; ver `reports/input_adaptations.json`.

## O que o piloto executa

`inputs/in.pilot` contém:

1. Minimização CG, limitada a 2.000 iterações/10.000 avaliações.
2. Aquecimento NVT 50 → 503 K em 2 ps.
3. Relaxação NPT 503 K/1 atm por 10 ps.
4. Teste reativo NPT por 20 ps, com ambas as reações epóxi/amina.
5. Checkpoints, estrutura final e resumo automático.

O passo é 1 fs. Essas durações foram escolhidas para um diagnóstico limitado, **não como protocolo validado de cura/equilíbrio**. A conversão é contada pelo número de oxigênios do produto hidroxila (tipo 12), dividido pelos 200 epóxis iniciais. O resumo distingue execução bem-sucedida de reação efetivamente observada.

O arquivo `pilot_summary.json` informa conservação de carga, contagem de átomos, conversão, avisos e desempenho. `bo_eligible`, `density_equilibrated` e `tg_calculated` permanecem falsos. Se nenhuma reação ocorrer, a sintaxe pode ter passado, mas a aplicação das regras de reação ainda não está validada.

## Carbono

Diretório dedicado:

`$HOME/formulationlab-dds-livraghi-pilot-20261021`

Acesso direto por `<cluster-host>` na rede UFABC, usuário `<user>`. `scripts/carbono.py` lê `<credential file>` somente em memória e exige chave SSH previamente confiável. Não desativa a verificação da chave.

O módulo institucional LAMMPS 23Jun2022 Update 3 consultado não oferece `bond/react`. O build particular antigo também não inclui REACTION. Por isso `build.slurm` compila uma instalação separada com MOLECULE, KSPACE, REACTION e OPENMP, reutilizando as fontes já existentes em `formulationlab-upet-20261020/lammps-metatomic` e o ambiente de compilação existente. Esse script pressupõe esses recursos; não baixa novamente o código LAMMPS.

- Build: partição curta `grafite`, duas CPUs, 8 GiB, limite de 30 min.
- Piloto: `grafeno`, 16 CPUs OpenMP, 16 GiB, limite de 30 min.
- Conta/QoS explícitas: `users`/`normal`.
- Limites de partição consultados: `grafite` permite até duas CPUs por job; `grafeno`, até 16. A fila pode atrasar o início.

Os jobs rodam nos nós de computação; não há MD ou compilação no nó de login. Nenhuma GPU é necessária. Os demais ambientes são preservados.

## Comandos no Windows

Dependência do controlador local: `paramiko` (`scripts/requirements.txt`). A auditoria usa apenas a biblioteca padrão. Na preparação desta pasta foi usado o Python do runtime Codex, pois `python` no PATH abre o alias da Microsoft Store.

```powershell
$py = '<local-path>'
$env:PYTHONPATH = '<local-path>'
Set-Location '<local-path>'
& $py scripts/verify_sources.py
& $py scripts/prepare_inputs.py
& $py scripts/audit_source.py
& $py scripts/carbono.py stage
& $py scripts/carbono.py build
& $py scripts/carbono.py build_status
# Somente depois que a compilação tiver terminado com sucesso:
& $py scripts/carbono.py preflight
& $py scripts/carbono.py submit
& $py scripts/carbono.py status
& $py scripts/carbono.py collect
```

Para reconstruir os downloads em outra cópia do repositório: `& $py scripts/download_sources.py`. O script retoma transferências parciais do Zenodo e verifica seus checksums antes da extração.

Alternativamente, `& $py scripts/carbono.py submit --after-build` agenda o piloto com dependência de sucesso da compilação registrada, cancelando-o se a dependência falhar. Isso não significa que já executou. Consultar `status` e conferir os resultados antes de interpretar o teste como aprovado.

Cada `submit` cria um novo diretório e um novo job. Não repetir para consultar andamento; usar `status`. IDs ficam em `reports/`. Se o diretório temporário de dependências não existir em uma sessão futura, instalar Paramiko em um ambiente Python dedicado.

## Campanha de cura e propriedades

**Campanhas de composição em paralelo:** 23 receitas adicionais foram preparadas e submetidas no array **175540**, com até 368 CPUs além das 16 da referência. Ver a [matriz, protocolo e critérios de cura](parallel_20261021/README.md). O lote usa arquivos isolados e não altera o job de referência. Consultar `scripts/parallel_control.py status` e `collect` para acompanhar e atualizar os resultados locais.

A sequência já está implementada e submetida. Os jobs e critérios estão em [reports/CAMPANHA.md](reports/CAMPANHA.md). A execução curta de validação testa os comandos; sua duração não permite interpretar suas propriedades físicas.

Consultar e coletar os jobs existentes, no mesmo ambiente Python acima:

```powershell
& $py scripts/carbono.py campaign_status
& $py scripts/carbono.py properties_status
& $py scripts/carbono.py campaign_collect
& $py scripts/carbono.py properties_collect
```

Para uma NOVA campanha: `stage`, `validation_submit`, consultar `validation_status` até concluir, `campaign_submit`, `properties_submit`. Não repetir submissões para consultar andamento. A submissão da campanha exige o marcador de sucesso da validação completa.

O próximo resultado decisivo é a conversão/conectividade da rede. Só depois de obter propriedades estáveis desta referência devemos gerar outras configurações/proporções. O resfriamento implementado é uma triagem rápida dependente da taxa, não reprodução integral do protocolo publicado.

Alterar apenas o número de DDS muda tanto a proporção quanto a conversão máxima acessível. Para um BO de DGEBA/DDS, a feature de composição é `4*N_DDS/(2*N_DGEBA)`; o notebook atual de DGEBA/PGE/EDA precisa ser adaptado. Não transplantar dados de uma química para outra.

O throughput anterior de D230 não é um benchmark deste sistema. Não há promessa de tempo para uma campanha de Tg até medir a nova rota.
