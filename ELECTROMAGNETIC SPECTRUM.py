from sympy import symbols
from sympy.physics.optics import TWave
A,phi,f = symbols('A,phi,f')
w = TWave(0.2,50)
answer = w.angular_velocity
print(answer)
