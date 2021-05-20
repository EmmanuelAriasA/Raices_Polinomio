import math
import sys
import cmath as cmath
from math import inf

# ----------------------------------------------------------------Algoritmo para sacar la derivada---------------------------------------------------------------------------------


def deriv(grado, coeficiente):
    k = 0
    derivado = []
    if((len(coeficiente)) == 1):
        derivado.append(0)
    else:
        while k < grado:
            derivado.append((grado - k)*(coeficiente[k]))
            k += 1
    return derivado

# ----------------------------------------------------------------Algoritmo de Horner-------------------------------------------------------------------------------------


def horner(grado, coeficientes, x):
    polinomio = coeficientes[grado]
    k = grado - 1
    while (k >= 0):
        polinomio = coeficientes[k] + (polinomio*x)
        k = k - 1
    return polinomio

# ----------------------------------------------------------------Algoritmo para sacar Bairstown---------------------------------------------------------------------------------


def Bairstow(coeficientes, r, s, grado, raiz, tolerancia):
    if(grado < 1):
        return None
    if((grado == 1) and (coeficientes[1] != 0)):
        raiz.append(float(-coeficientes[0])/float(coeficientes[1]))
        return None
    if(grado == 2):
        D = (coeficientes[1]**2.0)-(4.0)*(coeficientes[2])*(coeficientes[0])
        if(D < 0):
            X1 = (-coeficientes[1] - cmath.sqrt(D))/(2.0*coeficientes[2])
            X2 = (-coeficientes[1] + cmath.sqrt(D))/(2.0*coeficientes[2])
        else:
            X1 = (-coeficientes[1] - math.sqrt(D))/(2.0*coeficientes[2])
            X2 = (-coeficientes[1] + math.sqrt(D))/(2.0*coeficientes[2])
        raiz.append(X1)
        raiz.append(X2)
        return None
    n = len(coeficientes)
    b = [0]*len(coeficientes)
    c = [0]*len(coeficientes)
    b[n-1] = coeficientes[n-1]
    b[n-2] = coeficientes[n-2] + r*b[n-1]
    i = n - 3
    while(i >= 0):
        b[i] = coeficientes[i] + r*b[i+1] + s*b[i+2]
        i = i - 1
    c[n-1] = b[n-1]
    c[n-2] = b[n-2] + r*c[n-1]
    i = n - 3
    while(i >= 0):
        c[i] = b[i] + r*c[i+1] + s*c[i+2]
        i = i - 1
    Din = ((c[2]*c[2])-(c[3]*c[1]))**(-1.0)
    r = r + (Din)*((c[2])*(-b[1])+(-c[3])*(-b[0]))
    s = s + (Din)*((-c[1])*(-b[1])+(c[2])*(-b[0]))
    if(abs(b[0]) > tolerancia or abs(b[1]) > tolerancia):
        return Bairstow(coeficientes, r, s, grado, raiz, tolerancia)
    if (grado >= 3):
        Dis = ((r)**(2.0))+((4.0)*(1.0)*(s))
        X1 = (r - (cmath.sqrt(Dis)))/(2.0)
        X2 = (r + (cmath.sqrt(Dis)))/(2.0)
        raiz.append(X1)
        raiz.append(X2)
        return Bairstow(b[2:], r, s, grado-2, raiz, tolerancia)


# ------------------------------------------------------------------------------Ingresamos el grado del polinomio---------------------------------------------------------------------
grdo = int(input("Grado del polinomo: "))
coef = []

# ------------------------------------------------------------------------------Compara si el grado es correcto-----------------------------------------------------------------------
if grdo <= 0 or grdo == 0:
    print("Ingrese un grado valido: ")
    sys.exit()

# ------------------------------------------------------------------------------Ingresamos los coeficientes de nuestro polinomio-------------------------------------------------------
print("\n * * * * * Coeficientes empezando por el termino de mayor grado * * * * * ")
for i in range(grdo+1):
    coeficiente = float(input("Ingresa el coeficiente: "))
    coef.append(coeficiente)

# ------------------------------------------------------------------------------Ingresamos cifras significativas-----------------------------------------------------------------------
signi = int(input("Cifras significativas: "))

# ------------------------------------------------------------------------------Sacamos la derivada de nuestro polinomio---------------------------------------------------------------
# Sacamos la primera derivada mandando llamar nuestra funcion
primDev = deriv(grdo, coef)

# Sacamos la segunda derivada mandando llamar nuestra funcion
segDev = deriv(grdo-1, primDev)

# --------------------------------------------------------------------Sacamos los valores que necesita el Algoritmo de Bairstow--------------------------------------------------------
r = coef[-1]/coef[0]
s = r
r1 = primDev[-1]/primDev[0]
r2 = segDev[-1]/segDev[0]
s1 = r1
s2 = r2
# Sacamos la tolerancia (o sifras significativas de los valores)
tolerancia = 0.5*pow(10, 2-signi)

# --------------------------------------------------------------------Sacamos los valores que necesita el Algoritmo de Bairstow--------------------------------------------------------
# Reacomodamos nuestras listas para utilizar el metodo de Bairstow
coef.reverse()
primDev.reverse()
segDev.reverse()

# Listas donde vamos a guardar las raices
raizP = []
raizD1 = []
raizD2 = []

# --------------------------------------------------------------------Sacamos las raices imaginarios con el Algoritmo de Bairstow--------------------------------------------------------
Bairstow(coef, r, s, grdo, raizP, tolerancia)
Bairstow(primDev, r1, s1, grdo - 1, raizD1, tolerancia)
Bairstow(segDev, r2, s2, grdo - 2, raizD2, tolerancia)

# Separar raices complejas de las reales
raizD1 = [x for x in raizD1 if type(x) is not complex]
raizD2 = [x for x in raizD2 if type(x) is not complex]

# ----------------------------------------------------------------------------Sacamos los maximos y los minimos-------------------------------------------------------------------------
# Máximos y mínimos
maximos = []
minimos = []

c = 1/pow(10, signi)

if(len(raizD1) != 0):
    for x in raizD1:
        primer = horner(grdo-1, primDev, (x-c))
        segundo = horner(grdo-1, primDev, (x+c))

    if ((primer < 0) and (segundo > 0)):
        minimos.append(complex(x, horner(grdo, coef, x)))
        minimos.append((x, horner(grdo, coef, x)))

    elif ((primer > 0) and (segundo < 0)):
        maximos.append(complex(x, horner(grdo, coef, x)))
        maximos.append((x, horner(grdo, coef, x)))
else:
    maximos.append("La funcion no tiene maximos")
    minimos.append("La funcion no tiene minimos")

# ----------------------------------------------------------------------------Sacamos los puntos de inflexion----------------------------------------------------------------------------
inflexion = []

if(len(raizD2) != 0) and grdo > 2:
    for x in raizD2:
        inflexion.append((x, horner(grdo, coef, x)))
else:
    inflexion.append("La funcion no tiene puntos de inflexion")

# -----------------------------------------------------------------Sacamos puntos donde es creciente y decreciente------------------------------------------------------------------------
crece = []
decrece = []
raizD1.sort()

if grdo == 1:
    crece.append("La funcion no crece")
    decrece.append("La funcion no decrece")
else:
    if(len(raizD1) != 0):
        for i in range(1, len(raizD1)):
            if horner(grdo-1, primDev, raizD1[i]-c) < 0:
                decrece.append((raizD1[i-1], raizD1[i]))
            else:
                #crece.append(complex(raizD1[i-1], raizD1[i]))
                crece.append((raizD1[i-1], raizD1[i]))

        if horner(grdo-1, primDev, raizD1[-1]+c) < 0:
            #decrece.append(complex(raizD1[-1], inf))
            decrece.append((raizD1[-1], inf))
        else:
            #crece.append(complex(raizD1[-1], inf))
            crece.append((raizD1[-1], inf))

        if horner(grdo-1, primDev, raizD1[0] - c < 0):
            decrece.append((-inf, raizD1[0]))
        else:
            crece.append((-inf, raizD1[0]))
    else:
        crece.append("La funcion no crece")
        decrece.append("La funcion no decrece")

# ----------------------------------------------------------------------------Sacamos los intervalos de concavidad------------------------------------
arriba = []
abajo = []
raizD2.sort()

if(len(raizD2) != 0) and grdo > 2:
    for i in range(1, len(raizD2)):
        if horner(grdo-2, segDev, raizD2[i]-c) < 0:
            #abajo.append(complex(raizD2[i-1], raizD2[i]))
            abajo.append((raizD2[i-1], raizD2[i]))
        else:
            #arriba.append(complex(raizD2[i-1], raizD2[i]))
            arriba.append((raizD2[i-1], raizD2[i]))

    if horner(grdo-2, segDev, raizD2[-1] + c) < 0:
        #abajo.append(complex(raizD2[-1], inf))
        abajo.append((raizD2[-1], inf))
    else:
        #arriba.append(complex(raizD2[-1], inf))
        arriba.append((raizD2[-1], inf))

    if horner(grdo-2, segDev, raizD2[0] - c) < 0:
        #abajo.append(complex(-inf, raizD2[0]))
        abajo.append((-inf, raizD2[0]))
    else:
        #arriba.append(complex(-inf, raizD2[0]))
        arriba.append((-inf, raizD2[0]))
else:
    arriba.append("La funcion no es concava hacia arriba en ningun punto\n")
    abajo.append("La funcion no es concava hacia abajo en ningun punto\n")

# ---------------------------------------------------------Mostramos en pantalla todos los datos obtenidos---------------------------------------------
print("\nLas raices reales e imaginarias del polinomio son:\n")

for i in raizP:
    print(i)

# -----------------------------------------------------------------------------------------------------------------------------------------------------
print("\n")
print("\nLos maximos del polinomio son:")

for x in maximos:
    print(x, end=", ")

print("\nlos minimos del polinomio son:")

for x in minimos:
    print(x, end=", ")

# -----------------------------------------------------------------------------------------------------------------------------------------------------
print("\n")
print("\nLos puntos de inflexion del polinomio son:")

for x in inflexion:
    print(x, end=", ")

# -----------------------------------------------------------------------------------------------------------------------------------------------------
print("\n")
print("\nEs creciente en los intervalos: ")

for x in crece:
    print(x, end=", ")

print("\nEs decreciente en los intervalos:")

for x in decrece:
    print(x, end=", ")

# -----------------------------------------------------------------------------------------------------------------------------------------------------
print("\n")
print("\nEs concava hacia arriba en los intervalos: ")

for x in arriba:
    print(x, end=", ")

print("\nEs concava hacia abajo en los intervalos: ")

for x in abajo:
    print(x, end=", ")
