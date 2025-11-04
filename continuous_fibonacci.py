import numpy as np
import matplotlib.pyplot as plt

phi = (1 + np.sqrt(5))/2
x = np.linspace(0, 10, 100)
f = phi**x

plt.plot(x, f, label=r'$f(x) = \varphi^x$')
plt.plot(x, phi**(x+2), '--', label=r'$f(x+2)$')
plt.plot(x, phi**(x+1) + phi**x, ':', label=r'$f(x+1) + f(x)$')
plt.legend()
plt.title("Continuous Fibonacci: Exact Match")
plt.savefig("continuous_fibonacci_proof.png")
plt.show()