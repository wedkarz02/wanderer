import csv
import matplotlib.pyplot as plt

n = []
gp = []
mc = []
err = []

with open("dump/mc_verify.csv", "r") as f:
    reader = csv.reader(f, delimiter=";")
    next(reader)
    for row in reader:
        n.append(int(row[0]))
        gp.append(float(row[1]))
        mc.append(float(row[2]))
        err.append(float(row[3]))

plt.figure()
plt.yscale("log")
plt.plot(n, err, color="red",
         linestyle="-", marker="", label="Error")
plt.xlabel("Rozmiar macierzy (liczba skrzyżowań)")
plt.ylabel("Wartość błędu")
plt.title("Porównanie różnicy wyników metod numerycznych względem Monte Carlo",
          loc="center", wrap=True)
plt.grid()
plt.legend()
plt.savefig("img/mc_verify.png")
