# Ambiente no Mac

Desde 2026-11-09 o projeto é trabalhado só no Mac (Apple Silicon); o LAMMPS roda numa VM Multipass.

| Passo | Comando |
|---|---|
| Python (modelos, figuras, montagem do .docx) | `bash ambiente/setup_mac.sh` e depois `source .venv/bin/activate` |
| VM `lammps` + LAMMPS stable_30Sep2026 (MEAM, MANYBODY, OPENMP) | `bash ambiente/criar_vm_lammps.sh [CPUS] [MEMORIA]` — só com a máquina livre |
| Entrar na VM (o repositório fica montado em `~/gm-elastoplastic-nucleation`) | `multipass shell lammps` |
| Montar um manuscrito | `python paper/build_paper.py` → `paper/out/` (usa o `mathml2omml.xsl` do Word para Mac) |

Sem GPU no Apple Silicon: as rodadas antes feitas com `-sf gpu -pk gpu 1` usam `-sf omp -pk omp N`.

Pastas locais que não vão para o GitHub: `.venv/`, `paper/out/`, `privado/` (documentos editoriais
do artigo) e os estados `*.data`/`*.dump` da recristalização.
