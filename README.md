<p align="center">
  <img src="https://i.pinimg.com/originals/a4/96/c2/a496c2b6bc5d7cfe0e0674f6598c38ad.gif" alt="Espaço animado" width="100%">
</p>

<h1 align="center"> 🪐 Integrais e Evolução Orbital </h1>

<p align="center">
  <em>Simulação gráfica e interativa baseada em mecânica celeste para o estudo secular da evolução de sistemas estelares binários sob transferência de massa.</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-140033?style=for-the-badge&logo=python&logoColor=00e5ff" alt="Python">
  <img src="https://img.shields.io/badge/SciPy-0a001a?style=for-the-badge&logo=scipy&logoColor=ff00aa" alt="SciPy">
  <img src="https://img.shields.io/badge/Matplotlib-140033?style=for-the-badge&logo=python&logoColor=00e5ff" alt="Matplotlib">
  <img src="https://img.shields.io/badge/Astrofísica-0a001a?style=for-the-badge&logo=spacetelescope&logoColor=ff00aa" alt="Astrofísica">
</p>

<hr>

### 🌌 Sobre o Projeto

Este repositório contém a simulação computacional desenvolvida para o estudo da dinâmica de sistemas binários de estrelas. Focada na resolução numérica das equações seculares de transferência de massa, a simulação demonstra como a perda de massa de uma estrela doadora altera o semieixo maior ($a$) e a excentricidade ($e$) da órbita ao longo de milhões de anos, utilizando a biblioteca `SciPy` e visualizações avançadas via `Matplotlib`.

### 🧮 O Modelo Matemático

A simulação integra numericamente as taxas de variação secular dos parâmetros orbitais:
- **Variação do Semieixo Maior ($a$):** Derivado da solução integral da Equação 25, modelando a expansão ou contração orbital em função da razão de massas ($q = M_1 / M_2$).
- **Variação da Excentricidade ($e$):** Derivado da Equação 26, rastreando como a transferência conservativa ou não-conservativa de momento angular molda o formato elíptico da órbita.

### 🛠️ Tecnologias Utilizadas

- **Python 3.x** — Linguagem principal de programação.
- **SciPy** — Utilizado através do `solve_ivp` para a integração numérica de equações diferenciais ordinárias (ODEs) ao longo do tempo.
- **NumPy** — Manipulação de arrays e modelagem matemática vetorial.
- **Matplotlib** — Renderização gráfica dos subplots temporais e da geometria analítica polar das órbitas.

### 🚀 Como executar a simulação localmente

1. Clone este repositório:
   ```bash
   git clone https://github.com/AruDaguer/integrais-evolucao-orbital.git

2. Acesse a pasta do projeto:
   ```bash
   cd integrais-evolucao-orbital

3. Instale as dependências necessárias:
   ```bash
   pip install numpy matplotlib scipy

4. Execute o script da simulação:
   ```bash
   python simulacao.py