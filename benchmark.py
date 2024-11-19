import subprocess

# Funkcja do przechwytywania czasu wykonania
def benchmark_mpi(program, procs, repeats,m_size):
    results = []

    for i in range(repeats):
        print(f"  Wykonanie {i+1}/{repeats}...")
        args = " "+str(m_size)
        # Uruchomienie programu za pomocą mpiexec
        result = subprocess.run(
            ["mpiexec","--oversubscribe", "-n", str(procs), program,args],
            capture_output=True, 
            text=True
        )

        # Zbieranie ostatniej linii wyjścia, która zawiera czas wykonania
        time_output = result.stdout.strip().split('\n')[-1]

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

procs_list = [2,4,6,8,16,32]  # Liczba procesorów
repeats = 10000  # Liczba powtórzeń
matrix_size = range(100,1100,100)

# Zbieranie wyników
output_file = "times.csv"
with open(output_file, "w") as f:
    f.write("processes,size,time\n")

    for procs in procs_list:
        for m_size in matrix_size:
            print(f"Uruchamianie benchmarku dla {procs} procesorów, rozmiaru macierzy {m_size}")

            # Zbieranie wyników dla określonej liczby procesorów
            times = benchmark_mpi(program_parrarel, procs, repeats, m_size)

            # Zapisanie wyników do pliku
            for time in times:
                f.write(f"{procs},{m_size},{time}\n")
 
    for m_size in matrix_size:
        print(f"Uruchamianie benchmarku sekwencyjnie dla rozmiaru macierzy {m_size}")

        # Zbieranie wyników dla określonej liczby procesorów
        times = benchmark_seq(program_ortodox, repeats, m_size)

        # Zapisanie wyników do pliku
        for time in times:
            f.write(f"{0},{m_size},{time}\n")

    for m_size in matrix_size:
        print(f"Uruchamianie benchmarku sekwencyjnie dla rozmiaru macierzy {m_size}")

        # Zbieranie wyników dla określonej liczby procesorów
        times = benchmark_seq(program_panmicia, repeats, m_size)

        # Zapisanie wyników do pliku
        for time in times:
            f.write(f"{-1},{m_size},{time}\n")

print(f"Benchmark zakończony. Wyniki zapisano w pliku {output_file}.")

import subprocess

# Funkcja do przechwytywania czasu wykonania
def benchmark_mpi(program, procs, repeats,m_size):
    results = []

    for i in range(repeats):
        print(f"  Wykonanie {i+1}/{repeats}...")
        args = " "+str(m_size)
        # Uruchomienie programu za pomocą mpiexec
        result = subprocess.run(
            ["mpiexec", "-n", str(procs), program,args],
            capture_output=True, 
            text=True
        )

        # Zbieranie ostatniej linii wyjścia, która zawiera czas wykonania
        time_output = result.stdout.strip().split('\n')[-1]

        print(f"    Zebrany czas: {time_output}. Procesory = {procs}, m_size = {m_size}")
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

# procs_list = [2,4,6,8,16,32]  # Liczba procesorów
procs_list = [4,2,1]
repeats = 5000  # Liczba powtórzeń
matrix_size = range(100,1100,100)

# Zbieranie wyników
output_file = "times.csv"
with open(output_file, "w") as f:
    f.write("processes,size,time\n")

    for procs in procs_list:
        for m_size in matrix_size:
            print(f"Uruchamianie benchmarku dla {procs} procesorów, rozmiaru macierzy {m_size}")

            # Zbieranie wyników dla określonej liczby procesorów
            times = benchmark_mpi(program_parrarel, procs, repeats, m_size)

            # Zapisanie wyników do pliku
            for time in times:
                f.write(f"{procs},{m_size},{time}\n")
 
    for m_size in matrix_size:
        print(f"Uruchamianie benchmarku sekwencyjnie dla rozmiaru macierzy {m_size}")

        # Zbieranie wyników dla określonej liczby procesorów
        times = benchmark_seq(program_ortodox, repeats, m_size)

        # Zapisanie wyników do pliku
        for time in times:
            f.write(f"{0},{m_size},{time}\n")

    for m_size in matrix_size:
        print(f"Uruchamianie benchmarku sekwencyjnie dla rozmiaru macierzy {m_size}")

        # Zbieranie wyników dla określonej liczby procesorów
        times = benchmark_seq(program_panmicia, repeats, m_size)

        # Zapisanie wyników do pliku
        for time in times:
            f.write(f"{-1},{m_size},{time}\n")

print(f"Benchmark zakończony. Wyniki zapisano w pliku {output_file}.")
