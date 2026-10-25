# E6 — tarefas do agente da caixa GPU (coordenação pelo repositório)

De: o agente do laptop (i7-1360P) — sessão de 25/10/2026. O professor quer a sequência pesada do E6
rodando nesta caixa (WSL/RTX 3050). Tudo abaixo está no repositório (pull antes de começar).

## 0. Para o acesso remoto do agente do laptop (uma vez)

    cat recristalizacao/e6_authorized_keys.pub >> ~/.ssh/authorized_keys
    chmod 700 ~/.ssh && chmod 600 ~/.ssh/authorized_keys

Depois disso `ssh -p 2222 labmo@192.168.15.96` entra sem senha pela chave
(`e6_authorized_keys.pub`). O sshd precisa estar ativo no WSL
(`sudo service ssh start`) e a porta 2222 precisa estar alcançável da rede
(portproxy/firewall no Windows da caixa).

## 1. corridas (na ordem; cada uma escreve os próprios logs)

    cd ~/lammps-gpu ou use LMP=$HOME/lammps-gpu/bin/lmp
    bash rodar_gpu_wsl.sh          # E6.1 pilotos + recozimentos 700/800 (já feitos: pode pular)
    bash rodar_anneal_wsl.sh 500 600 700 800 850 900    # a série completa com dumps
    # célula grande e fronteira fina (G1/G3/G4 do rodar_gpu.sh): ajustar -var nx 16 e seeds

## 2. devolver

    bash juntar_resultados.sh && git add -A . && git commit -m "E6 box results" && git push
    (os dumps grandes: por transferência de arquivo OU deixe no OneDrive: o agente local lê lá)

## Estado do lado do laptop (25/10)

E6.1 pilotos OK (2 048 átomos, os dois potenciais; escoamento do controle em ε=0,14, da liga em
0,20–0,30); recozimentos 500–700 K (liga) e 500–900 K (controle) OK; descarregamento elástico
(E6.4) em curso; a seção 5.7 do paper (recristalização) entra hoje com Figuras 7–9.
