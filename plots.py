import pandas as pd
import matplotlib.pyplot as plt


def plot_histogram(df):
    """
    Funkcja rysuje histogramy dla danych czasów wykonania.
    """
    plt.figure(figsize=(10, 6))

    # Rysowanie histogramu dla czasów wykonania
    plt.hist(df["time"], bins=100, color='skyblue', edgecolor='black', alpha=0.7)

    # Dodanie etykiet i tytułu
    plt.xlabel('Czas wykonania (s)')
    plt.ylabel('Liczba prób')
    plt.title('Histogram czasów wykonania dla różnych liczby procesów')

    # Wyświetlenie wykresu
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # Ścieżka do pliku CSV
    file_name = 'times.csv'
    
    # Odczytanie danych z pliku
    df = pd.read_csv(file_name)
    print(df)
    # Rysowanie histogramu wyników
    plot_histogram(df.loc[df["processes"] == 4])
