import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# 1. PARÂMETROS FÍSICOS INICIAIS
massa_estrela_doadora = 2.0  # M1 em Massas Solares 
massa_estrela_acretora = 2.0 # M2 em Massas Solares
tamanho_orbita_inicial = 1.0  # Semieixo maior (a) em UA
formato_orbita_inicial = 0.2  # Excentricidade (e)

taxa_transferencia_massa = -1e-7  # M_dot_0 (negativa)
tempo_simulacao_anos = 5000000


# 2. MODELO MATEMÁTICO COMPLETO (EQUAÇÕES 39, 40 E A15 DE SEPINSKY ET AL. 2007)
def calcular_taxas_seculares(tempo, variaveis_orbitais):
    semieixo_a, excentricidade_e, massa_doadora, massa_acretora = (
        variaveis_orbitais
    )

    excentricidade_e = np.clip(excentricidade_e, 1e-4, 0.99)
    razao_massas_q = massa_doadora / massa_acretora

    # Rotação síncrona no periastro (f = Omega_1 / Omega_orb, P = 1.0)
    f_rotacao = 1.0

    # Posição adimensional de L1 no periastro (Equação A15)
    log_q = np.log10(razao_massas_q)
    X_L1 = (
        0.529
        + 0.231 * log_q
        - (f_rotacao**2) * (0.031 + 0.025 * excentricidade_e) * (1 + 0.4 * log_q)
    )

    # Razão entre a distância de L1 ao centro de massa do doador e o semieixo maior
    r_A1_sobre_a = X_L1 * (1 - excentricidade_e)

    factor_massa = taxa_transferencia_massa / massa_doadora
    raiz_e = np.sqrt(1 - excentricidade_e**2)

    # Equação 39 (Taxa secular do semieixo maior)
    termo_a = (
        excentricidade_e * r_A1_sobre_a
        + (razao_massas_q - 1) * (1 - excentricidade_e**2)
    )
    taxa_variacao_semieixo = (semieixo_a / np.pi) * factor_massa * (1 / raiz_e) * termo_a

    # Equação 40 (Taxa secular da excentricidade)
    termo_e = r_A1_sobre_a + 2 * (1 - excentricidade_e) * (razao_massas_q - 1)
    taxa_variacao_excentricidade = (
        (raiz_e / (2 * np.pi)) * factor_massa * termo_e
    )

    taxa_variacao_massa_doadora = taxa_transferencia_massa
    taxa_variacao_massa_acretora = -taxa_transferencia_massa

    return [
        taxa_variacao_semieixo,
        taxa_variacao_excentricidade,
        taxa_variacao_massa_doadora,
        taxa_variacao_massa_acretora,
    ]


# 3. RESOLUÇÃO NUMÉRICA
pontos_no_tempo = np.linspace(0, tempo_simulacao_anos, 1000)
estado_inicial = [
    tamanho_orbita_inicial,
    formato_orbita_inicial,
    massa_estrela_doadora,
    massa_estrela_acretora,
]

evolucao_orbital = solve_ivp(
    calcular_taxas_seculares,
    (0, tempo_simulacao_anos),
    estado_inicial,
    t_eval=pontos_no_tempo,
)

historico_semieixo = evolucao_orbital.y[0]
historico_excentricidade = evolucao_orbital.y[1]
historico_massa_doadora = evolucao_orbital.y[2]
historico_massa_acretora = evolucao_orbital.y[3]

# Cálculo da Razão de Massas (q = M_doadora / M_acretora) ao longo do tempo
historico_q = historico_massa_doadora / historico_massa_acretora


# 4. VISUALIZAÇÃO DOS GRÁFICOS (4 Subplots)
fig = plt.figure(figsize=(20, 5))

# Gráfico 1: Evolução do Semieixo Maior
ax1 = fig.add_subplot(1, 4, 1)
ax1.plot(evolucao_orbital.t, historico_semieixo, color='#1f77b4', lw=2)
ax1.set_title('Tamanho da Órbita (a)', fontweight='bold')
ax1.set_xlabel('Tempo (anos)')
ax1.set_ylabel('Semieixo Maior (UA)')
ax1.grid(True, linestyle='--', alpha=0.7)

# Gráfico 2: Evolução da Excentricidade
ax2 = fig.add_subplot(1, 4, 2)
ax2.plot(evolucao_orbital.t, historico_excentricidade, color='#d62728', lw=2)
ax2.set_title('Formato da Órbita (e)', fontweight='bold')
ax2.set_xlabel('Tempo (anos)')
ax2.set_ylabel('Excentricidade')
ax2.grid(True, linestyle='--', alpha=0.7)

# Gráfico 3: Razão de Massa x Excentricidade (NOVO)
ax3 = fig.add_subplot(1, 4, 3)
ax3.plot(historico_q, historico_excentricidade, color='#2ca02c', lw=2)
ax3.set_title('Razão de Massa (q) vs. Excentricidade (e)', fontweight='bold')
ax3.set_xlabel('Razão de Massas ($q = M_1 / M_2$)')
ax3.set_ylabel('Excentricidade ($e$)')
ax3.grid(True, linestyle='--', alpha=0.7)
# Inverte o eixo x opcionalmente se quiser mostrar a evolução temporal de q decrescente
ax3.invert_xaxis() 

# Gráfico 4: Geometria Analítica Polar (Órbita Inicial vs Final)
ax4 = fig.add_subplot(1, 4, 4, projection='polar')
angulo_theta = np.linspace(0, 2 * np.pi, 500)

raio_inicial = (
    tamanho_orbita_inicial * (1 - formato_orbita_inicial**2)
) / (1 + formato_orbita_inicial * np.cos(angulo_theta))
raio_final = (
    historico_semieixo[-1] * (1 - historico_excentricidade[-1] ** 2)
) / (1 + historico_excentricidade[-1] * np.cos(angulo_theta))

ax4.plot(
    angulo_theta,
    raio_inicial,
    label='Órbita Inicial',
    color='black',
    linestyle='--',
)
ax4.plot(
    angulo_theta,
    raio_final,
    label='Órbita Final (5 Milhões de Anos)',
    color='green',
    lw=2,
)
ax4.set_title('Comparativo Geométrico', fontweight='bold')
ax4.legend(loc='upper right', bbox_to_anchor=(1.4, 1.1))

plt.tight_layout()
plt.show()