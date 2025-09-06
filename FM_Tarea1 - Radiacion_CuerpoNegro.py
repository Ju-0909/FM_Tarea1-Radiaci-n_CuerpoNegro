import numpy as np
import matplotlib.pyplot as plt

#Constantes (SI)
h = 6.62607015e-34 # J s
c = 2.99792458e8   # m/s
k = 1.380649e-23   # J/K

def planck_lambda(lam, T):
    """Radiancia espectral B_lambda ( W sr ^-1 m ^-3)."""
    lam = np.asarray(lam)
    x = (h*c)/(lam*k*T)
    return (2*h*c**2)/(lam**5) * (1.0/np.expm1(x))

def wien_lambda(lam, T):
    lam = np.asarray(lam)
    x = (h*c)/(lam*k*T)
    return (2*h*c**2)/(lam**5) * np.exp(-x)

### Juan Durán ###
def rayleigh_jeans_lambda(lam, T):
    """Aproximacion de Rayleigh-Jeans: B_lambda (W sr ^-1 m ^-3)"""
    lam = np.asarray(lam)
    return (2*c*k*T) / (lam**4)

# Mallado de longitudes: 200 nm a 10 micras
lam = np.linspace(200e-9, 10e-6, 4000)
# Temperaturas a explorar
temperaturas = [1000, 2000, 4000] # K

#--- Gráfica 1: escala lineal ---
plt.figure(figsize=(7,5))
for T in temperaturas:
    plt.plot(lam*1e9, planck_lambda(lam,T), label=f"Planck {T} K")
    plt.plot(lam*1e9, wien_lambda(lam,T), linestyle='--',label=f"Wien {T} K")
    plt.plot(lam*1e9, rayleigh_jeans_lambda(lam,T), linestyle=':',label=f"Rayleigh-Jeans {T} K")
plt.xlabel("Longitud de onda (nm)")
plt.ylabel(r"$B_\lambda$ (W sr$^{-1}$ m$^{-3}$)")
plt.title("Planck vs. Wien vs. Rayleigh-Jeans (Escala lineal)")
plt.legend()
plt.tight_layout()
plt.show()

#--- Gráfica 2: eje Y logaritmico ---
plt.figure(figsize=(7,5))
for T in temperaturas:
    plt.semilogy(lam*1e9, planck_lambda(lam,T), label=f"Planck {T} K")
    plt.semilogy(lam*1e9, wien_lambda(lam,T), linestyle='--',label=f"Wien {T} K")
    plt.semilogy(lam*1e9, rayleigh_jeans_lambda(lam,T), linestyle=':',label=f"Rayleigh-Jeans {T} K")
plt.xlabel("Longitud de onda (nm)")
plt.ylabel(r"$B_\lambda$ (W sr$^{-1}$ m$^{-3}$)")
plt.title("Planck vs. Wien vs. Rayleigh-Jeans (Y log)")
plt.legend()
plt.tight_layout()
plt.show()

# (Opcional) Verificación desplazamiento de Wien
b = 2.897771955e-3 # m K
for T in temperaturas:
    lam_max = b/T
    print(f"T={T} K-> lambda_max{lam_max*1e9:.0f} nm")