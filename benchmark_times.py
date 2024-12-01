import subprocess
import time

# Funkcja do przechwytywania czasu wykonania
def benchmark_mpi(program, procs, repeats,m_size):
    results = []

    for i in range(repeats):
        print(f"  Wykonanie {i+1}/{repeats}...")
        args = str(m_size)
        print(program,args)
        # Uruchomienie programu za pomocą mpiexec
        result = subprocess.run(
            ["mpiexec","--oversubscribe", "-n", str(procs), program,"c",args],
            capture_output=True, 
            text=True
        )

        # Zbieranie ostatniej linii wyjścia, która zawiera czas wykonania
        time_output = result.stdout.strip().split('\n')[-1]
        # time_output = result
        print(f"    Zebrany czas: {time_output} sekund")
        results.append(time_output)

    return results

# Funkcja do przechwytywania czasu wykonania
def benchmark_seq(program, repeats, matrix_size):
    results = []
    args = " "+str(m_size)
    for i in range(repeats):
        print(f"  Wykonanie {i+1}/{repeats}...")

        # Uruchomienie programu za pomocą mpiexec
        result = subprocess.run(
            [program, args],
            capture_output=True, 
            text=True
        )

        # Zbieranie ostatniej linii wyjścia, która zawiera czas wykonania
        time_output = result.stdout.strip().split('\n')[-1]

        print(f"    Zebrany czas: {time_output} sekund")
        results.append(time_output)

    return results

# Parametry benchmarku
program_parrarel = "matmul" 
program_ortodox = "./eigen"  
program_panmicia = "./seq"

procs_list = [1,2,4,6,8,12,16,32]  # Liczba procesorów
repeats_list = [10]*50  # Liczba powtórzeń
repeats = 1000

matrix_size = list(range(10,100,10))
matrix_size += list(range(100,1100,100))

print(matrix_size)
# Zbieranie wyników
output_file = "times_comm2.csv"
i = 0
with open(output_file, "w") as f:
    f.write("processes,size,time\n")
    for repeat in repeats_list:
        i += 1
        print("powt: "+str(i)+"/"+str(len(repeats_list)))
        for procs in procs_list:
            for m_size in matrix_size:
                print(f"Uruchamianie benchmarku dla {procs} procesorów, rozmiaru macierzy {m_size}")

                # Zbieranie wyników dla określonej liczby procesorów
                times = benchmark_mpi(program_parrarel, procs, repeat, m_size)

                # Zapisanie wyników do pliku
                for timee in times:
                    f.write(f"{procs},{m_size},{timee}\n")
                time.sleep(1)

print(f"Benchmark zakończony. Wyniki zapisano w pliku {output_file}.")
