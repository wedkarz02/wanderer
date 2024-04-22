import csv
import matplotlib.pyplot as plt

n_ns = []
jacobi_ns = []
seidel_ns = []
gauss_ns = []
gp_ns = []

with open("dump/default_no_sparse_time.csv", "r") as f:
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
plt.xlabel("Rozmiar macierzy")
plt.ylabel("Czas wykonania (w ms)")
plt.title("Porównanie czasu wykonania metod rozwiązywania\nukładów równań liniowych",
          loc="center", wrap=True)
plt.grid()
plt.legend()
plt.savefig("img/default_no_sparse_time.png")

n_s = []
jacobi_s = []
gauss_s = []
gp_s = []

with open("dump/default_sparse_time.csv", "r") as f:
    reader = csv.reader(f, delimiter=";")
    next(reader)
    for row in reader:
        n_s.append(int(row[0]))
        jacobi_s.append(float(row[1]))
        gauss_s.append(float(row[2]))
        gp_s.append(float(row[3]))

plt.figure()
plt.yscale("log")
plt.plot(n_ns, jacobi_s, color="red",
         linestyle="-", marker="", label="Jacobi")
plt.plot(n_ns, gauss_s, color="green", linestyle="-",
         marker="", label="Gauss")
plt.plot(n_ns, gp_s, color="orange", linestyle="-",
         marker="", label="Gauss (wybór częściowy)")
plt.xlabel("Rozmiar macierzy rzadkiej")
plt.ylabel("Czas wykonania (w ms)")
plt.title("Porównanie czasu wykonania metod rozwiązywania układów równań liniowych przy użyciu macierzy rzadkich",
          loc="center", wrap=True)
plt.grid()
plt.legend()
plt.savefig("img/default_sparse_time.png")
