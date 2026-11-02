import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. PARÂMETROS FÍSICOS INICIAIS
massa_estrela_doadora = 2.0     # M1 em Massas Solares
massa_estrela_acretora = 1.44   # M2 em Massas Solares
tamanho_orbita_inicial = 1.0    # Semieixo maior (a) em Unidades Astronômicas
formato_orbita_inicial = 0.2    # Excentricidade (e)

# Taxa de transferência de massa (negativa pois a doadora perde massa)
taxa_transferencia_massa = -1e-7 
tempo_simulacao_anos = 5000000

# 2. O MODELO MATEMÁTICO (EQUAÇÕES 25 E 26 INTEGRADAS)
def calcular_taxas_seculares(tempo, variaveis_orbitais):
    semieixo_a, excentricidade_e, massa_doadora, massa_acretora = variaveis_orbitais
    
    # Trava de segurança para a excentricidade não ultrapassar o limite físico (1.0)
    excentricidade_e = np.clip(excentricidade_e, 0.0, 0.99) 
    
    razao_massas_q = massa_doadora / massa_acretora
    
    # Derivada do semieixo maior (Resultado da integral da Equação 25)
    taxa_variacao_semieixo = (semieixo_a / np.pi) * (taxa_transferencia_massa / massa_doadora) * np.sqrt(1 - excentricidade_e**2) * (razao_massas_q - 1)
    
    # Derivada da excentricidade (Resultado da integral da Equação 26)
    taxa_variacao_excentricidade = (1 / np.pi) * (taxa_transferencia_massa / massa_doadora) * np.sqrt(1 - excentricidade_e**2) * (1 - excentricidade_e) * (razao_massas_q - 1)
    
    taxa_variacao_massa_doadora = taxa_transferencia_massa
    taxa_variacao_massa_acretora = -taxa_transferencia_massa 
    
    return [taxa_variacao_semieixo, taxa_variacao_excentricidade, taxa_variacao_massa_doadora, taxa_variacao_massa_acretora]

# 3. RESOLUÇÃO NUMÉRICA
pontos_no_tempo = np.linspace(0, tempo_simulacao_anos, 1000)
estado_inicial = [tamanho_orbita_inicial, formato_orbita_inicial, massa_estrela_doadora, massa_estrela_acretora]

# O SciPy calcula a evolução de todas as variáveis ao longo de 5 milhões de anos
evolucao_orbital = solve_ivp(calcular_taxas_seculares, (0, tempo_simulacao_anos), estado_inicial, t_eval=pontos_no_tempo)

historico_semieixo = evolucao_orbital.y[0]
historico_excentricidade = evolucao_orbital.y[1]

# 4. VISUALIZAÇÃO DOS GRÁFICOS
fig = plt.figure(figsize=(15, 5))

# Gráfico A: Evolução do Tamanho da Órbita
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(evolucao_orbital.t, historico_semieixo, color='#1f77b4', lw=2)
ax1.set_title('Evolução do Tamanho da Órbita (a)', fontweight='bold')
ax1.set_xlabel('Tempo (anos)')
ax1.set_ylabel('Semieixo Maior (UA)')
ax1.grid(True, linestyle='--', alpha=0.7)

# Gráfico B: Evolução do Formato da Órbita
ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(evolucao_orbital.t, historico_excentricidade, color='#d62728', lw=2)
ax2.set_title('Evolução do Formato da Órbita (e)', fontweight='bold')
ax2.set_xlabel('Tempo (anos)')
ax2.set_ylabel('Excentricidade')
ax2.grid(True, linestyle='--', alpha=0.7)

# Gráfico C: Geometria Analítica (Órbita Inicial vs Final)
ax3 = fig.add_subplot(1, 3, 3, projection='polar')
angulo_theta = np.linspace(0, 2*np.pi, 500)

# Equação polar da elipse
raio_inicial = (tamanho_orbita_inicial * (1 - formato_orbita_inicial**2)) / (1 + formato_orbita_inicial * np.cos(angulo_theta))
raio_final = (historico_semieixo[-1] * (1 - historico_excentricidade[-1]**2)) / (1 + historico_excentricidade[-1] * np.cos(angulo_theta))

ax3.plot(angulo_theta, raio_inicial, label='Órbita Inicial', color='black', linestyle='--')
ax3.plot(angulo_theta, raio_final, label='Órbita Final (5 Milhões de Anos)', color='green', lw=2)
ax3.set_title('Comparativo Geométrico da Órbita', fontweight='bold')
ax3.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1)) 

plt.tight_layout()
plt.show()