import matplotlib.pyplot as plt
import particle
import numpy as np

m_e = particle.Particle.from_pdgid(11).mass / 1000  # GeV
m_mu = particle.Particle.from_pdgid(13).mass / 1000  # GeV

alpha = np.linspace(0, np.pi, 10001)
for l, m_l, E_mu in [(r'e', m_e, (0.93, 11.1, 28.2)), (r'\mu', m_mu, (33.6, 50.2, 77.2))]:
    for E_mu in E_mu:
        E = E_mu + m_e
        p = np.sqrt(E_mu*E_mu - m_mu*m_mu)
        Ep = np.sqrt(2*E_mu*m_e + m_e*m_e + m_mu*m_mu)
        gamma = E / Ep
        beta = p / E
        pp = np.sqrt(Ep*Ep/4 - m_l*m_l)
        B = beta * Ep / (2 * pp)
        get_A = lambda alpha: np.sin(alpha) / (np.cos(alpha) + B)
        get_theta = lambda alpha: np.arctan(get_A(alpha) / gamma)
        theta = get_theta(alpha)
        theta += np.pi * (theta < 0)
        plt.plot(alpha, theta, label=r'$\mu^+e^- \to l^+l^-$, $E_\mu = $'.replace('l', l) + '%.2f' % E_mu + ' GeV')
    plt.xlabel(r'$\alpha$ [rad]')
    plt.ylabel(r'$\theta$ [rad]')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'theta_vs_alpha_{l.replace("\\", "")}.pdf')
    plt.close()

for E_mu in (0.11, 0.2, 1, 10, 100):
    E = E_mu + m_e
    p = np.sqrt(E_mu*E_mu - m_mu*m_mu)
    Ep = np.sqrt(2*E_mu*m_e + m_e*m_e + m_mu*m_mu)
    gamma = E / Ep
    beta = p / E
    pp = np.sqrt(Ep*Ep - (m_e + m_mu)**2) / 2
    B = beta * Ep / (2 * pp)
    get_A = lambda alpha: np.sin(alpha) / (np.cos(alpha) + B)
    get_theta = lambda alpha: np.arctan(get_A(alpha) / gamma)
    theta = get_theta(alpha)
    theta += np.pi * (theta < 0)
    plt.plot(alpha, theta, label=r'$\mu^+e^- \to \mu e$, $E_\mu = $' + '%.2f' % E_mu + ' GeV')
plt.xlabel(r'$\alpha$ [rad]')
plt.ylabel(r'$\theta$ [rad]')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('theta_vs_alpha_mue.pdf')
plt.close()
