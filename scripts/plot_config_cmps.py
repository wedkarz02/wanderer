import csv
import matplotlib.pyplot as plt

n_ns = []
jacobi_ns = []
seidel_ns = []
gauss_ns = []
gp_ns = []

with open("dump/config_time_no_sparse.csv", "r") as f:
    reader = csv.reader(f, delimiter=";")
    next(reader)
    for row in reader:
        n_ns.append(int(row[0]))
        jacobi_ns.append(float(row[1]))
        seidel_ns.append(float(row[2]))
        gauss_ns.append(float(row[3]))
        gp_ns.append(float(row[4]))

plt.figure()
plt.yscale("log")
plt.plot(n_ns, jacobi_ns, color="red",
         linestyle="-", marker="", label="Jacobi")
plt.plot(n_ns, seidel_ns, color="blue", linestyle="-",
         marker="", label="Gauss-Seidel")
plt.plot(n_ns, gauss_ns, color="green", linestyle="-",
         marker="", label="Gauss")
plt.plot(n_ns, gp_ns, color="orange", linestyle="-",
         marker="", label="Gauss (wybór częściowy)")
plt.xlabel("Rozmiar macierzy (liczba skrzyżowań)")
plt.ylabel("Czas wykonania (w ms)")
plt.title("Porównanie czasu wykonania metod rozwiązywania\nukładów równań liniowych",
          loc="center", wrap=True)
plt.grid()
plt.legend()
plt.savefig("img/config_time_no_sparse.png")

n_s = []
jacobi_s = []
seidel_s = []
gauss_s = []
gp_s = []

with open("dump/config_time_sparse.csv", "r") as f:
    reader = csv.reader(f, delimiter=";")
    next(reader)
    for row in reader:
        n_s.append(int(row[0]))
        jacobi_s.append(float(row[1]))
        seidel_s.append(float(row[2]))
        gauss_s.append(float(row[3]))
        gp_s.append(float(row[4]))

plt.figure()
plt.yscale("log")
plt.plot(n_s, jacobi_s, color="red",
         linestyle="-", marker="", label="Jacobi")
plt.plot(n_s, seidel_s, color="blue",
         linestyle="-", marker="", label="Gauss-Seidel")
plt.plot(n_s, gauss_s, color="green", linestyle="-",
         marker="", label="Gauss")
plt.plot(n_s, gp_s, color="orange", linestyle="-",
         marker="", label="Gauss (wybór częściowy)")
plt.xlabel("Rozmiar macierzy rzadkiej (liczba skrzyżowań)")
plt.ylabel("Czas wykonania (w ms)")
plt.title("Porównanie czasu wykonania metod rozwiązywania układów równań liniowych przy użyciu macierzy rzadkich",
          loc="center", wrap=True)
plt.grid()
plt.legend()
plt.savefig("img/config_time_sparse.png")
