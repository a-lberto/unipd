# Metodo dei pannelli

## Costruzione profilo alare
Si definiscono gli assi unitari per le ascisse [0,1] e per le ordinate [-0.5,0.5], con un leggero incremento (ordine di 0.05) per evidenziare gli spigoli.

Viene inizializzato il vettore delle ascisse x, compreso tra 0 e c.

Vengono richieste le cifre del profilo NACA 4 cifre.

Vengono impostate le costanti Ai per lo spessore sulla linea media con correzione di Moran per A2.

Vengono definite le funzioni di linea di alzata media e linea dello spessore.

Si disegna la linea di alzata media tratteggiata.

Si disegnano gli spessori sommando e sottraendo dalla linea di alzata media i valori della funzione spessore.

## Discretizzazione con pannelli
Si sceglie il numero dei pannelli in cui suddividere il profilo dell'ala.

Si costruisce una circonferenza suddivisa in n angoli (n+1 punti).
I punti vengono proiettati verticalmente sugli spessori.
La scelta delle ascisse segue la distribuzione dei nodi Chebychev.
Le ordinate sono calcolate con la stessa modalità per il disegno del profilo d'ala.

Nel frattempo si crea un vettore che tiene nota dell'indice del nodo.

Si hanno ora gli n+1 punti che delimitano i n pannelli.

## Metodo dei pannelli con sorgenti (e vortice)
In ogni punto di controllo viene inserita una singolarità. Sorgente di portata q variabile da pannello a pannello, e un vortice di circolazione gamma, che è costante per ogni pannello.
Le incognite sono appunto le q(i) e gamma (n+1).


Si trova per ogni pannello l'angolo di inclinazione dello stesso, la sua lunghezza e il suo punto di mezzeria.