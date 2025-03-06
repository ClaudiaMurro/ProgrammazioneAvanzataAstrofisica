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
OBJ = sph.o main.o

all: a.out

a.out: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o a.out $(LIBS)

sph.o: sph.c sph.h
	$(CC) $(CFLAGS) -c sph.c

main.o: main.c sph.h
	$(CC) $(CFLAGS) -c main.c

# Esegui il programma
run: a.out
	./a.out

# Comando per generare il grafico
grafici: grafici.gp animazione.gp
	gnuplot grafici.gp animazione.gp

# Pulizia dei file compilati
clean:
	rm -f *.o a.out

# Pulizia dei file di risultato
clean_result: 
	rm -f *.dat *.png *.gif

# Pulizia di tutti i file prodotti
clean_all: 
	rm -f *.dat *.png *.o a.out *.gif
