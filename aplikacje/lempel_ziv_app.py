import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage 
from scipy.spatial.distance import pdist
from tkinter import *
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    
def get_words(sequence):
    """
    Generuje zbiór wszystkich słów z danej sekwencji.

    sequence - sekwencja
    """

    d = []
    l = len(sequence)
    i = 0
    k = 1
    while i<l:
        while sequence[i: i+k] in d and sequence[::-1][-(i+k):-(i)] in d and i+k<l:
            k+=1
        if sequence[i:i+k] not in d:
            d.append(sequence[i:i+k])
        i+=k
        k = 1
    return d

def program(seq1, seq2):
    """
    Zwraca znormalizowaną wartość compression distance using the Lempel–Ziv complexity estimation algorithm
    """

    if seq1 == seq2: return 0

    l1 = len(get_words(seq1))
    l2 = len(get_words(seq2))
    l3 = len(get_words(seq1+seq2))
    l4 = len(get_words(seq2+seq1))
    C = (((l3 - min(l1,l2))/max(l1,l2))+((l4 - min(l1,l2))/max(l1,l2)))/2
    return C

def SimMatrix(seqs):
    """
    Tworzy macierz odległości między sekwencjami.

    seqs - lista sekwencji
    """

    n = len(seqs)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = program(seqs[i], seqs[j])
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

root.title("Lempel-Ziv app")
window_width = 1000
window_height = 600
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
center_x = int(screen_width/2 - window_width/2)
center_y = int(screen_height/2 - window_height/1.8)
root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')


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

    input_list_names = org
    input_list_seq = seq

    
    X = pdist(SimMatrix(input_list_seq)) 
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