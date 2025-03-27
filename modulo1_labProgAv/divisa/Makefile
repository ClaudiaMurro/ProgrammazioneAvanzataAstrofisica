#COMANDI TERMINALE
#compilazione ed esecuzione: make run
#creo grafici: make grafici
#cancello file oggetto e a.out (lascio grafici e dati): make clean
#cancello grafici e dati: make clean_result
#cancello tutti i file prodotti: make clean_all
#-------------------------------------------------------------------------

CC = gcc
CFLAGS = -Wall -O2
LIBS = -lm
OBJ = decadimento.o main.o

# Regola di default
all: a.out

a.out: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o a.out $(LIBS)

# Compilazione di decadimento.o
decadimento.o: decadimento.c decadimento.h
	$(CC) $(CFLAGS) -c decadimento.c

# Compilazione di main.o
main.o: main.c decadimento.h
	$(CC) $(CFLAGS) -c main.c

# Creazione dell'eseguibile
a.out: $(OBJ)
	$(CC) $(CFLAGS) -o a.out $(OBJ)

# Esegui il programma
run: a.out
	./a.out

# Comando per generare il grafico
grafici: grafici.gp
	gnuplot grafici.gp

# Pulizia dei file compilati
clean:
	rm -f *.o a.out

# Pulizia dei file di risultato
clean_result: 
	rm -f *.dat *.png

# Pulizia di tutti i file prodotti
clean_all: 
	rm -f *.dat *.png *.o a.out

