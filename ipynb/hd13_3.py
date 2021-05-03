from sympy import symbols, Matrix, simplify, cos, sin, I, atan2, exp, re, collect, sign, sqrt, Abs, arg
from math import radians, degrees

# Horntíðnin er fasti
omega = 2  # rad/s

# Breytan `t` táknar tíma
t = symbols('t', real=True)

# Gildi á rásaeiningum er fasti
r1, r2, r3, r4, c1, l1 = 1, 1, 1, 2, 2, 2

yr1, yr2, yr3, yr4, yc1, yl1 = 1 / r1, 1 / r2, 1 / r3, 1 / r4, (c1 * omega * I), 1 / (l1 * omega * I)

# Aðferð 1: Gerum ráð fyrir að i = 1A og margföldum með 12<30-gráðu virkja alveg í lokinn
va, vb, vc = simplify(Matrix([
    [yr1 + yr2 + yc1, -yr2, -yc1],
    [-yr2, yr2 + yr3 + yl1, -yr3],
    [-yc1, -yr3, yc1 + yr3 + yr4]
]).LUsolve(Matrix([[1, 0, 0]]).T))  # Set inn i = 1 A hér

v = simplify(va - vb)

# Breytum á pólform
r, theta = Abs(v), degrees(arg(v))

# Hér er margfaldað með 12 og 30 gráður lagðar við
print(f"v(t) = {r * 12 * cos(omega * t + 30 + theta)}")

# Aðferð 2: Ímyndum okkur að straumurinn sé tvinntala
r = radians(30)  # `Sympy.exp` gerir ráð fyrir radíönum
i = 12 * exp(I*(omega*t + r))

va, vb, vc = simplify(Matrix([
    [yr1 + yr2 + yc1, -yr2, -yc1],
    [-yr2, yr2 + yr3 + yl1, -yr3],
    [-yc1, -yr3, yc1 + yr3 + yr4]
]).LUsolve(Matrix([[i, 0, 0]]).T))  # Hér setjum við inn strauminn á tvinntöluformi

# Tökum rauntölu-hluta útmerkisins. Athugði að það er ennþá á röngu formi!
v_out = re(simplify(va - vb))
print(f"v(t) = {v_out}")

# Hér þurfum að varpa úr a*cos(x) * b*sin(x) -> c*cos(x + theta)
# Finnum stuðla og notum aðferð frá https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
studlar = collect(v_out, [sin(2*t + r), cos(2*t + r)], evaluate=False)
b, a = studlar[sin(2*t + r)], studlar[cos(2*t + r)]

c, theta = sign(a) * sqrt(a**2 + b**2), atan2(-b, a)

print(f"v(t) = {c * cos(omega * t + degrees(r + theta))}")
