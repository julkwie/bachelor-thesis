from tkinter import *
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage 
from scipy.spatial.distance import pdist 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def get_words(sequence, word_size: int):
    """
    Generuje zbiór wszystkich słów o zadanym rozmiarze z danej sekwencji.

    sequence - sekwencja, z której generujemy listę słów
    word_size - rozmiar/długość słowa
    """
    words = [sequence[i:i+word_size] for i in range(len(sequence)-word_size+1)]
    return words

def create_full_word_set(first, second, word_size = 3):
    """
    Tworzy pełen zbiór unikalnych słów z obu sekwencji. 
    Wykorzystuje funkcję get_words do uzyskania listy słów.

    first - pierwsza sekwencja
    second - druga sekwencja
    word_size - rozmiar/długość słowa
    """
    x = get_words(first, word_size)
    y = get_words(second, word_size)

    full = set()
    full = list(set(x).union(set(y)))
    
    return full


def transform_to_vector(word_list, word_set):
    """
    Tworzy wektor reprezentujący sekwencję poprzez zliczanie wystąpień słów z zestawu.

    word_list - lista słów danej sekwencji
    word_set - zbiór unikalnych słów do uwzględnienia w wektorze
    """
    word_count_dict = {word: word_list.count(word) for word in word_set}
    vector = np.array(list(word_count_dict.values()))
    return vector


def euclidean_distance(vector1: np.array, vector2: np.array) -> float:
    """
    Oblicza odległość euklidesową między dwoma wektorami.

    vector1 - pierwszy wektor
    vector2 - drugi wektor
    """
    return np.linalg.norm(vector1 - vector2)

def program(first_sequence, second_sequence, word_size = 3) -> float:
    """
    Zwraca znormalizowaną wartość odległości euklidesowej między dwoma wektorami, 
    reprezentującymi sekwencję poprzez zliczanie wystąpień słów z zestawu.
    Wykorzystuje funkcję create_full_word_set do tworzenia pełen zbiór unikalnych słów z obu sekwencji. 
    Wykorzystuje funkcję get_words do uzyskania listy słów.

    first_sequence - pierwsza sekwencja
    second_sequence - druga sekwencja
    word_size - rozmiar/długość słowa, domyślnie 3
    """

    first_words = get_words(first_sequence, word_size)
    second_words = get_words(second_sequence, word_size)

    full_set = create_full_word_set(first_sequence, second_sequence, word_size)
    a = transform_to_vector(first_words, full_set)
    b = transform_to_vector(second_words, full_set)
    
    return euclidean_distance(a, b)/len(a)


def SimMatrix(seqs, word_size = 3):
    """
    Tworzy macierz odległości między sekwencjami.

    seqs - lista sekwencji
    word_size - rozmiar/długość słowa, domyślnie 3
    """

    n = len(seqs)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = program(seqs[i], seqs[j], word_size) 

    return A

def parse_fasta_text(fasta_text):
    """
    Przetwarza tekst FASTA na listę organizmów oraz listę sekwencji.

    fasta_text - tekst w formie FASTA
    """
    lines = fasta_text.strip().split('\n')
    organisms = []
    sequences = []
    sequence = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if sequence:
                sequences.append(sequence)
                sequence = ""
            organisms.append(line[1:])
        else:
            sequence += line
    
    if sequence:
        sequences.append(sequence)
    
    return organisms, sequences

def process_fasta_file():
    """
    Używając funkcji get_results tworzy dendrogram UPGMA i wyświetla go.
    Wgrywa i przetwarza plik wybrany z eksploratora plików w formacie FASTA używając funkcji parse_fasta_text. 
    """

    file_path = filedialog.askopenfilename(
        title="Select a FASTA file"
    )
    if file_path:
        try:
            with open(file_path, 'r') as file:
                fasta_text = file.read()
            organisms, sequences = parse_fasta_text(fasta_text)
            get_results(organisms, sequences)
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")

def process_fasta_text():
    """
    Używając funkcji get_results tworzy dendrogram UPGMA i wyświetla go.
    Wgrywa i przetwarza wpisany tekst w formacie FASTA używając funkcji parse_fasta_text.
    """
    fasta_text = text_area.get("1.0", END)
    if fasta_text.strip():
        try:
            organisms, sequences = parse_fasta_text(fasta_text)
            get_results(organisms, sequences)
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")
    else:
        messagebox.showwarning("Input Error", "The FASTA text input is empty.")


root = Tk()

root.title("Kmers app")
window_width = 1000
window_height = 600
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
center_x = int(screen_width/2 - window_width/2)
center_y = int(screen_height/2 - window_height/1.8)
root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')


label = Label(root, text = "Insert the word size:")
label.pack(pady=5)
entry_size = Entry(root, width= 60)
entry_size.pack(pady=10)

Label(root, text="Paste your FASTA formatted text below:").pack(pady=10)
text_area = Text(root, height=10, width=60)
text_area.pack(pady=10)

Button(root, text="Process FASTA Text", command=process_fasta_text).pack(pady=5)

Label(root, text="or").pack(pady=5)

Button(root, text="Open FASTA File", command=process_fasta_file).pack(pady=5)


def get_results(org, seq):
    """
    Używając funkcji SimMatrix tworzy dendrogram sekwencji.

    org - lista organizmów
    seq - lista sekwencji
    """
    result_window = Toplevel(root)
    result_window.title("Results")

    word_size = int(entry_size.get())
    input_list_names = org
    input_list_seq = seq

    
    X = pdist(SimMatrix(input_list_seq, word_size)) 
    dist = linkage(X, method="average")  #average = UPGMA method
    fig = plt.figure(figsize=(10, 5))
    dn = dendrogram(dist, labels = input_list_names, orientation = "left")

    plt.gca().margins(x=0.1, y=0.1)
    plt.tight_layout()
    plt.plot()

    canvas = FigureCanvasTkAgg(fig, master = result_window)
    canvas.draw()
    canvas.get_tk_widget().pack()

    root.plot_canvas = canvas

root.mainloop()
